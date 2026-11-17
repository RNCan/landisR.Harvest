# Static geometry is cached by harvest_extension() and may also be prepared
# explicitly for repeated calls to harvest_step().
.harvest_prepare_cell_topology <- function(simulation, topology) {
  forest <- landisR::landis_csr(simulation)
  topology$partial_neighbors <- .Call("h05_cell_topology", topology$ids, topology$positions,
    topology$stand_offsets, as.integer(c(forest$nrow, forest$ncol)),
    PACKAGE="landisR.harvest")
  topology
}

.harvest_is_partial <- function(site_selection) {
  site_selection %in% c("PartialStandSpread", "PartialStandSpreading")
}

# Mixed Random/MaxCohortAge prescriptions for living cohorts.
harvest_plan_e4d <- function(simulation, topology, prescriptions,
                            history=NULL, year=simulation$year,
                            landis_seed, landis_rng_offset=0,
                            allocation_mode=c("corrected","landis_compatible"),
                            rules_by_prescription=NULL,
                            harvest_timestep=simulation$project$canonical$timestep, damage_history=NULL) {
  .h04_sim(simulation)
  forest <- landisR::landis_csr(simulation)
  allocation_mode <- match.arg(allocation_mode)
  valid_num <- function(x,lo,hi) is.numeric(x)&&length(x)==1L&&
    !is.na(x)&&is.finite(x)&&x>=lo&&x<=hi&&x==trunc(x)
  if(missing(landis_seed)||!valid_num(landis_seed,0,4294967295))
    stop("landis_seed must be a uint32",call.=FALSE)
  shared <- simulation$landis_rng
  if (!valid_num(landis_rng_offset, 0, if (is.null(shared)) 1e9 else 2^53-1))
    stop("Invalid landis_rng_offset", call. = FALSE)
  if (!is.null(shared)) {
    if (allocation_mode != "landis_compatible" || landis_seed != simulation$seed ||
        (landis_rng_offset != 0 && landis_rng_offset != shared$offset))
      stop("Shared LANDIS RNG needs landis_compatible and the simulation seed and offset.",
           call. = FALSE)
    landis_rng_offset <- shared$offset
  }
  if(!valid_num(year,0,1e9))stop("Invalid harvest year",call.=FALSE)
  if(is.null(topology$d2_neighbors)||!inherits(topology,"landisR_h04_topology")||
     !identical(topology$ids,as.integer(forest$cell_id)))
    stop("Call harvest_prepare_adjacency_d2(simulation,topology) first",call.=FALSE)
  required<-c("management_area","prescription","fraction","ranking","target_type",
    "minimum_age","maximum_age","begin_time","end_time","site_selection")
  if(!is.data.frame(prescriptions)||!nrow(prescriptions)||!all(required%in%names(prescriptions)))
    stop("prescriptions must provide E2 columns",call.=FALSE)
  p<-prescriptions
  integers<-c("management_area","minimum_age","maximum_age","begin_time","end_time")
  for(name in integers)p[[name]]<-.h04_integer(p[[name]],name,0L)
  for(name in c("prescription","ranking","target_type","site_selection"))
    p[[name]]<-as.character(p[[name]])
  if(anyNA(p$prescription)||any(!nzchar(p$prescription))||
     anyNA(p$ranking)||any(!p$ranking%in%c("MaxCohortAge","Random","Economic"))||
     anyNA(p$target_type)||any(p$target_type!="area")||
     anyNA(p$site_selection)||
     any(!p$site_selection %in% c("Complete", "CompleteStandSpread", "CompleteStandSpreading") &
         !.harvest_is_partial(p$site_selection) & p$site_selection != "PatchCutting"))
    stop("Expected MaxCohortAge, Random or Economic ranking and a supported site selection with area quotas",
         call.=FALSE)
  patch <- p$site_selection == "PatchCutting"
  spreading <- p$site_selection != "Complete" & !patch
  partial <- .harvest_is_partial(p$site_selection)
  if (any(patch)) {
    if (!is.numeric(p$patch_fraction) || !is.numeric(p$patch_size_ha))
      stop("PatchCutting requires numeric patch_fraction and patch_size_ha", call.=FALSE)
    patch_fraction <- p$patch_fraction[patch]
    patch_size_ha <- p$patch_size_ha[patch]
    if (any(!is.finite(patch_fraction) | patch_fraction <= 0 | patch_fraction > 1) ||
        any(!is.finite(patch_size_ha) | patch_size_ha <= 0))
      stop("Expected 0 < patch_fraction <= 1 and patch_size_ha > 0", call.=FALSE)
    patch_option <- if ("patch_option" %in% names(p)) as.character(p$patch_option[patch])
      else rep("", sum(patch))
    patch_option[is.na(patch_option)] <- ""
    if (any(!patch_option %in% c("", "AllowOverlap", "PatchSize", "PercentCut")))
      stop("Invalid patch_option", call.=FALSE)
  }
  block_sizes <- function(name) {
    values <- numeric(nrow(p))
    if (any(spreading)) {
      if (!name %in% names(p) || !is.numeric(p[[name]]))
        stop("Stand spreading requires numeric ", name, call.=FALSE)
      values[spreading] <- p[[name]][spreading]
      if (anyNA(values) || any(!is.finite(values) | values < 0))
        stop("Invalid ", name, call.=FALSE)
    }
    values
  }
  min_block_ha <- block_sizes("min_block_ha")
  max_block_ha <- block_sizes("max_block_ha")
  if (any(spreading & (max_block_ha <= 0 | max_block_ha < min_block_ha)))
    stop("Stand spreading requires max_block_ha > 0 and max_block_ha >= min_block_ha",
         call.=FALSE)
  if(!is.numeric(p$fraction)||anyNA(p$fraction)||any(!is.finite(p$fraction))||
     any(p$fraction<0|p$fraction>1)||
     any(p$maximum_age<p$minimum_age)||any(p$end_time<p$begin_time)||
     anyDuplicated(p[c("management_area","prescription")]))
    stop("Invalid prescription target/age/period/duplicate",call.=FALSE)
  if(any(!p$management_area%in%topology$management_areas))
    stop("Unknown management area in prescriptions",call.=FALSE)
  if(allocation_mode=="landis_compatible" && any(p$fraction==0))
    stop("LANDIS-II refuses zero percent implementation",call.=FALSE)
  optional_integer <- function(name,default) {
    val<-if(name%in%names(p))p[[name]] else rep(default,nrow(p))
    .h04_integer(val,name,0L)
  }
  minimum_biomass <- numeric(nrow(p))
  if ("minimum_biomass_mg_ha" %in% names(p)) {
    values <- p$minimum_biomass_mg_ha
    if ((!is.numeric(values) && !(is.logical(values) && all(is.na(values)))) ||
        any(is.nan(values)) || any(!is.na(values) & (!is.finite(values) | values < 0))) {
      stop("minimum_biomass_mg_ha must contain nonnegative finite numbers or NA", call.=FALSE)
    }
    minimum_biomass[!is.na(values)] <- values[!is.na(values)]
  }
  threshold<-optional_integer("stand_adjacency",0L)
  hist_gap<-optional_integer("minimum_time_since_last_harvest",0L)
  min_damage <- optional_integer("min_time_since_damage", 0L)
  track_damage <- any(min_damage > 0L) || !is.null(damage_history) ||
    !is.null(attr(history, "harvest_damage_history", exact = TRUE))
  if (allocation_mode == "landis_compatible")
    min_damage[p$site_selection == "Complete"] <- 0L
  set_aside <- optional_integer("adjacency_neighbor_set_aside", 0L)
  use_reserve <- allocation_mode == "corrected" &&
    (any(set_aside > 0L) || "set_aside_until" %in% names(history))
  adj_type<-if("adjacency_type"%in%names(p))as.character(p$adjacency_type) else
    rep("none",nrow(p))
  key<-tolower(gsub("[[:space:]_-]+","",trimws(adj_type)))
  key[!is.na(key)&key==""&threshold==0L]<-"none"
  codes<-c(none=0L,standage=1L,timesincelastharvested=2L)
  if(anyNA(key)||any(!key%in%names(codes))||any(threshold>0L & key=="none"))
    stop("Active stand adjacency requires StandAge or TimeSinceLastHarvested",call.=FALSE)
  kind<-as.integer(unname(codes[key]))
  species_names<-as.character(forest$species)
  forest_species<-if("forest_type_species"%in%names(p))
    as.character(p$forest_type_species) else rep(NA_character_,nrow(p))
  if(length(forest_species)!=nrow(p))stop("Invalid forest_type_species",call.=FALSE)
  forest_index<-match(forest_species,species_names)-1L
  absent<-is.na(forest_species)|!nzchar(forest_species)
  if(any(!absent & is.na(forest_index)))
    stop("Unknown forest_type_species",call.=FALSE)
  forest_index[absent]<- -1L
  forest_lo<-optional_integer("forest_type_age_min",1L)
  forest_hi<-optional_integer("forest_type_age_max",65535L)
  forest_pct<-if("forest_type_percent"%in%names(p))
    as.numeric(p$forest_type_percent) else rep(50,nrow(p))
  if(length(forest_pct)!=nrow(p)||anyNA(forest_pct)||any(!is.finite(forest_pct))||
     any(forest_lo>forest_hi)||any(!absent&(forest_pct<=0|forest_pct>100)))
    stop("Invalid Required ForestType rule",call.=FALSE)
  last<-rep.int(NA_integer_,nrow(topology$stands))
  reserve_until <- rep.int(NA_real_, length(last))
  if(!is.null(history)){
    if(!is.data.frame(history)||!all(c("stand","last_harvest_year")%in%names(history)))
      stop("history needs stand,last_harvest_year",call.=FALSE)
    hs<-.h04_integer(history$stand,"history stand",0L)
    if(anyDuplicated(hs)||length(hs)!=length(last)||
       !setequal(hs,topology$stands$stand))stop("history/topology mismatch",call.=FALSE)
    history_order <- match(topology$stands$stand, hs)
    values<-as.numeric(history$last_harvest_year[history_order])
    if(any(!is.na(values)&(!is.finite(values)|values<0|values>year|values!=trunc(values))))
      stop("Invalid last_harvest_year",call.=FALSE)
    last<-as.integer(values)
    if (use_reserve && "set_aside_until" %in% names(history)) {
      expiry <- history$set_aside_until[history_order]
      if (!is.numeric(expiry) ||
          any(!is.na(expiry) & (!is.finite(expiry) | expiry < 0 |
                               expiry > 1e9 + .Machine$integer.max | expiry != trunc(expiry))))
        stop("Invalid set_aside_until", call.=FALSE)
      reserve_until <- as.numeric(expiry)
    }
  }
  harvest_timestep <- .h04_integer(harvest_timestep, "harvest_timestep")
  if (length(harvest_timestep) != 1L || harvest_timestep %% simulation$project$canonical$timestep != 0L)
    stop("harvest_timestep must be a positive multiple of the succession timestep", call. = FALSE)
  repeats <- .harvest_repeat_options(p, rules_by_prescription, harvest_timestep)
  repeat_state <- .harvest_repeat_input(history, p, topology, repeats)
  if (!is.null(repeats) && is.null(rules_by_prescription))
    stop("Repeating prescriptions require treatment rules", call. = FALSE)
  # Without treatment rules this is an allocation-only preview.
  removal_rules <- NULL
  continuous_ids <- FALSE
  if (!is.null(rules_by_prescription)) {
    active <- p$fraction > 0 & year >= p$begin_time & year <= p$end_time
    names_needed <- unique(p$prescription[active])
    if (!is.null(repeats)) {
      due <- Filter(function(job) job$year <= year, repeat_state$jobs)
      repeat_rule_names <- ifelse(repeats$kind == 1L, paste0(p$prescription, "(1)"), p$prescription)
      due_prescriptions <- unique(vapply(due, function(job) job$prescription + 1L, integer(1L)))
      names_needed <- union(names_needed, repeat_rule_names[due_prescriptions])
    }
    if (!is.list(rules_by_prescription) || is.null(names(rules_by_prescription)) ||
        !all(names_needed %in% names(rules_by_prescription))) {
      stop("Need one named cohort-removal table per active prescription", call.=FALSE)
    }
    # Structured treatments opt into the continuous LANDIS-II journal.
    # Keep the historical per-step IDs for callers using only removal tables.
    continuous_ids <- allocation_mode == "landis_compatible" && is.null(repeats) &&
      any(vapply(rules_by_prescription, function(x) is.list(x) && !is.data.frame(x),
                 logical(1L)))
    removal_rules <- lapply(names_needed, function(name) {
      rules <- .harvest_treatment(rules_by_prescription[[name]], species_names)$rules
      list(species=as.integer(match(rules$species, species_names) - 1L),
           lo=as.numeric(rules$age_min), hi=as.numeric(rules$age_max),
           fraction=as.numeric(rules$fraction))
    })
    rule_index <- match(p$prescription, names_needed)
    distinct_rules <- unique(removal_rules)
    rule_lookup <- match(removal_rules, distinct_rules)
    rule_index <- rule_lookup[rule_index]
    rule_index[!active] <- 0L
    removal_rules <- distinct_rules
  }
  cfg<-list(ma=as.integer(match(p$management_area,topology$management_areas)-1L),
    fraction=as.numeric(p$fraction),
    ranking_code=as.integer(match(p$ranking, c("MaxCohortAge", "Random", "Economic")) - 1L),
    min_age=as.integer(p$minimum_age),
    max_age=as.integer(p$maximum_age),begin=as.integer(p$begin_time),
    end=as.integer(p$end_time),adjacent_limit=as.integer(threshold),
    adjacent_type=kind,min_harvest_gap=as.integer(hist_gap),
    species_count=as.integer(length(species_names)),
    forest_species=as.integer(forest_index),forest_min_age=as.integer(forest_lo),
    forest_max_age=as.integer(forest_hi),forest_percent=as.numeric(forest_pct))
  economic <- which(p$ranking == "Economic")
  if (length(economic)) {
    if (!is.list(p$economic) || length(p$economic) != nrow(p))
      stop("Economic requires an economic list column of species ranking tables", call.=FALSE)
    tables <- lapply(p$economic[economic], function(x) {
      fields <- c("species", "rank", "merchantable_age")
      if (!is.data.frame(x) || !all(fields %in% names(x)) ||
          anyNA(x$species) || anyDuplicated(x$species) || any(!x$species %in% species_names))
        stop("Invalid Economic species table", call.=FALSE)
      weights <- .h04_integer(x$rank, "Economic rank", 0L)
      ages <- .h04_integer(x$merchantable_age, "Economic merchantable age", 0L)
      if (any(weights > 100L) || any(ages > 65535L))
        stop("Invalid Economic rank or merchantable age", call.=FALSE)
      lookup <- match(species_names, x$species)
      weights <- weights[lookup]
      ages <- ages[lookup]
      ages[is.na(ages)] <- 0L
      weights <- ifelse(ages > 0, weights / pmax(1, ages), 0)
      list(weight = weights, age = ages)
    })
    distinct <- unique(tables)
    cfg$economic_index <- rep.int(-1L, nrow(p))
    cfg$economic_index[economic] <- match(tables, distinct) - 1L
    cfg$economic_weight <- do.call(rbind, lapply(distinct, `[[`, "weight"))
    cfg$economic_age <- do.call(rbind, lapply(distinct, `[[`, "age"))
  }
  use_biomass <- minimum_biomass > 0 & p$fraction > 0 &
    year >= p$begin_time & year <= p$end_time
  if (any(use_biomass)) cfg$minimum_biomass_mg_ha <- minimum_biomass
  if (!is.null(shared)) cfg$landis_rng_state <- shared$state
  excluded <- simulation$establishment$excluded_cells
  if (length(excluded))
    cfg$delisted_cells <- as.integer(match(excluded, topology$ids) - 1L)
  if (track_damage) {
    cfg$min_damage <- min_damage
    cfg$damage_compatible <- allocation_mode == "landis_compatible"
    internal <- .harvest_damage_positions(attr(history, "harvest_damage_history", exact = TRUE),
                                          topology$ids, year)
    fire_records <- landisR::landis_fire_history(simulation)
    if (!is.null(fire_records))
      fire_records <- fire_records[fire_records$year <= year, , drop = FALSE]
    fire <- .harvest_damage_positions(fire_records, topology$ids, year)
    external <- .harvest_damage_positions(damage_history, topology$ids, year)
    cfg$harvest_damage_positions <- internal$position
    cfg$harvest_damage_year <- internal$year
    cfg$external_damage_positions <- c(fire$position, external$position)
    cfg$external_damage_year <- c(fire$year, external$year)
  }
  if (use_reserve) {
    cfg$set_aside <- set_aside
    cfg$set_aside_until <- reserve_until
  }
  if (continuous_ids) {
    previous_id <- attr(history, "harvest_event_id", exact = TRUE)
    cfg$last_event_id <- if (is.null(previous_id)) 0L else previous_id
  }
  if (!is.null(repeats)) {
    cfg$repeat_kind <- repeats$kind
    cfg$repeat_interval <- repeats$interval
    cfg$repeat_count <- repeats$count
    cfg$repeat_jobs <- repeat_state$jobs
    cfg$last_event_id <- repeat_state$event_id
    cfg$current_rank <- repeat_state$rank
    cfg$repeat_end_year <- simulation$project$canonical$duration
    cfg$repeat_compatible <- allocation_mode == "landis_compatible"
    cfg$repeat_rule_index <- as.integer(rule_lookup[match(repeat_rule_names, names_needed)] - 1L)
  }
  if (any(spreading)) {
    cfg$spreading <- spreading
    cfg$min_block_ha <- min_block_ha
    cfg$max_block_ha <- max_block_ha
  }
  if (any(spreading | patch))
    cfg$cell_area_ha <- (as.numeric(forest$cell_length) / 100)^2
  if (any(patch)) {
    cfg$patch_fraction <- cfg$patch_size_ha <- numeric(nrow(p))
    cfg$patch_overlap <- rep(FALSE, nrow(p))
    cfg$patch_fraction[patch] <- patch_fraction
    cfg$patch_size_ha[patch] <- patch_size_ha
    cfg$patch_overlap[patch] <- patch_option == "AllowOverlap"
  }
  if (any(partial | patch)) {
    if (is.null(topology$partial_neighbors))
      topology <- .harvest_prepare_cell_topology(simulation, topology)
    cfg$partial <- partial
    cfg$cell_topology <- topology$partial_neighbors
    cfg$landis_compatible <- allocation_mode == "landis_compatible"
  }
  if (!is.null(removal_rules)) {
    cfg$removal_rules <- removal_rules
    cfg$rule_index <- as.integer(rule_index - 1L)
  }
  # One batched native pass over all existing shared living CSR states.
  nat<-.Call("h04e4d_plan",forest$states,forest$cell_state,
    topology$positions,topology$stand_offsets,topology$ma_offsets,
    topology$d2_neighbors$offsets,topology$d2_neighbors$indices,
    topology$d2_neighbors$first_position,cfg,last,
    as.numeric(year),as.numeric(landis_seed),as.numeric(landis_rng_offset),
    PACKAGE="landisR.harvest")
  s<-as.integer(nat$selected_stands);j<-as.integer(nat$prescription_indices)
  loc <- lapply(seq_along(s), function(row) {
    selected <- nat$partial_cells[[row]]
    if (is.null(selected)) {
      stand <- s[row]
      index <- seq.int(topology$stand_offsets[stand]+1L, topology$stand_offsets[stand+1L])
      selected <- topology$positions[index]
    }
    as.integer(topology$ids[selected+1L])
  })
  events<-data.frame(EventID=as.integer(nat$event_ids),Time=rep.int(as.integer(year),length(s)),
    ManagementArea=topology$stands$management_area[s],Prescription=p$prescription[j],
    Stand=topology$stands$stand[s],StandAge=as.integer(nat$event_age),
    StandRank=as.numeric(nat$event_ranks),NumberOfSites=as.integer(nat$event_sites),
    rng_draw_index=as.numeric(nat$draw_indices),
    rng_uniform=as.numeric(nat$draw_uniform),
    selection_probability=as.numeric(nat$selection_probability),
    stringsAsFactors=FALSE)
  events$prescription_row<-j;events$cells<-I(loc)
  if (!is.null(repeats)) {
    events$RepeatNumber <- as.integer(nat$repeat_numbers)
    events$treatment <- events$Prescription
    returning <- events$RepeatNumber > 0L
    events$Prescription[returning] <- paste0(events$Prescription[returning], "(",
                                            events$RepeatNumber[returning], ")")
    single_return <- returning & repeats$kind[j] == 1L
    events$treatment[single_return] <- events$Prescription[single_return]
  }
  history_out<-data.frame(management_area=topology$stands$management_area,
    stand=topology$stands$stand,last_harvest_year=as.integer(nat$last_harvest_year))
  if (use_reserve)
    history_out$set_aside_until <- as.numeric(nat$set_aside_until)
  if (!is.null(repeats))
    attr(history_out, "harvest_repeat_state") <- .harvest_repeat_output(nat, p, topology, repeats)
  if (continuous_ids)
    attr(history_out, "harvest_event_id") <- as.integer(nat$last_event_id)
  if (track_damage)
    attr(history_out, "harvest_damage_history") <- data.frame(
      cell_id = topology$ids[nat$harvest_damage_positions + 1L],
      year = as.integer(nat$harvest_damage_year))
  plan_type <- if (any(patch)) {
    "Random/MaxCohortAge with PatchCutting"
  } else if (any(partial)) {
    "Random/MaxCohortAge with PartialStandSpread"
  } else if (any(spreading)) {
    "Random/MaxCohortAge with CompleteStandSpread"
  } else {
    "H04-E4d pilot: Random/MaxCohortAge Complete, multiple prescriptions, per-rule adjacency/history/Required forest type"
  }
  plan <- list(type=plan_type,
    year=as.integer(year),events=events,cells=as.integer(unlist(loc,use.names=FALSE)),
    history=history_out,
    prescriptions=within(p,{allocated_sites<-as.integer(nat$allocated_sites);
      eligible_stands<-as.integer(nat$eligible_stands);
      attempted_sites<-as.integer(nat$attempted_sites)}),
    stand_ages=as.integer(nat$stand_age),
    landis_seed=as.numeric(landis_seed),
    landis_rng_offset_start=as.numeric(landis_rng_offset),
    landis_rng_offset_end=as.numeric(nat$rng_offset_end),
    draws_consumed=as.numeric(nat$draws_consumed))
  if (!is.null(shared))
    plan$landis_rng <- list(state = nat$landis_rng_state, offset = plan$landis_rng_offset_end)
  plan
}

harvest_step_e4d <- function(simulation,topology,prescriptions,rules_by_prescription,
                            history=NULL,year=simulation$year,landis_seed,
                            landis_rng_offset=0,
                            allocation_mode=c("corrected","landis_compatible"),
                            harvest_timestep=simulation$project$canonical$timestep, damage_history=NULL) {
  if (!is.list(rules_by_prescription) || is.null(names(rules_by_prescription)))
    stop("Need a named list of cohort-removal tables",call.=FALSE)
  allocation_mode <- match.arg(allocation_mode)
  plan<-harvest_plan_e4d(simulation,topology,prescriptions,history,year,
    landis_seed,landis_rng_offset,allocation_mode,rules_by_prescription,harvest_timestep,damage_history)
  .harvest_apply_plan_e4d(simulation, plan, rules_by_prescription, allocation_mode)
}

.harvest_apply_plan_e4d <- function(simulation, plan, rules_by_prescription, allocation_mode,
                                     batch = TRUE) {
  forest <- landisR::landis_csr(simulation)
  events <- plan$events
  chosen <- unique(as.character(events$Prescription))
  event_groups <- split(seq_len(nrow(events)), factor(events$Prescription, levels = chosen))
  preserve_order <- allocation_mode == "landis_compatible"
  # Only group independent edits and use the exact, order-preserving state key.
  # Zero-biomass cohorts are ignored by that key, so retain sequential
  # application when their stored representation might affect state reuse.
  selected <- if (batch && preserve_order && length(chosen) > 1L)
    as.integer(unlist(events$cells, use.names = FALSE)) else integer()
  batch <- length(selected) > 0L && anyDuplicated(selected) == 0L &&
    isTRUE(all(forest$states$biomass > 0))
  positions <- if (batch) match(selected, forest$cell_id) else integer()
  if (anyNA(positions)) stop("Unknown or inactive cell", call. = FALSE)
  event_sizes <- lengths(events$cells)
  event_starts <- c(0L, cumsum(event_sizes)) + 1L
  current <- simulation
  effects <- sites <- planting <- planted_species <- excluded <- list()
  replacement <- cleared <- journal <- vector("list", length(chosen))
  species <- as.character(forest$species)
  for (i in seq_along(chosen)) {
    pres <- chosen[i]
    k <- event_groups[[i]]
    target <- as.integer(unlist(events$cells[k], use.names = FALSE))
    if (!length(target)) next
    treatment_name <- if ("treatment" %in% names(events)) events$treatment[k[1L]] else pres
    treatment <- .harvest_treatment(rules_by_prescription[[treatment_name]], species)
    if (batch) {
      indices <- sequence(event_sizes[k], from = event_starts[k])
      res <- .harvest_remove_at_positions(forest, as.integer(forest$cell_id),
                                          positions[indices], species, treatment$rules)
      if (length(res$changed_cells)) {
        replacement[[i]] <- res$replacement
        cleared[[i]] <- res$clear_cells
        journal[[i]] <- .harvest_apply_metadata(res)
      }
    } else {
      res <- harvest_apply_cells(current, target, rules = treatment$rules,
                                 preserve_cohort_order = preserve_order)
      current <- res$simulation
    }
    if (length(treatment$plant)) {
      planted_species[[pres]] <- treatment$plant
      planting[[pres]] <- data.frame(
        cell_id = rep.int(target, length(treatment$plant)),
        species = rep(treatment$plant, each = length(target))
      )
    }
    if (treatment$prevent_establishment)
      excluded[[pres]] <- res$changed_cells
    effects[[pres]]<-res$effects
    z<-res$sites;z$prescription<-pres;sites[[pres]]<-z
  }
  if (batch && any(lengths(journal))) {
    current <- landisR::landis_apply_cohort_batch(
      simulation, cohorts = do.call(rbind, replacement),
      clear_cells = as.integer(unlist(cleared, use.names = FALSE)),
      event_metadata = Filter(Negate(is.null), journal),
      preserve_cohort_order = TRUE)
  }
  all_sites<-if(length(sites))do.call(rbind,sites)else data.frame()
  # Preserve selected-but-no-effect events as zero-site, zero-biomass log rows.
  # Sparse selected-cell mapping; never expand a landscape raster per event.
  ev<-plan$events
  ev$HarvestedSites<-integer(nrow(ev))
  ev$MgBiomassRemoved<-numeric(nrow(ev))
  if(nrow(ev)&&nrow(all_sites)){
    chosen_cells <- if (length(selected)) selected else as.integer(unlist(ev$cells, use.names = FALSE))
    event_of_cell<-rep.int(seq_len(nrow(ev)),lengths(ev$cells))
    at<-match(as.integer(all_sites$cell_id),chosen_cells)
    if(anyNA(at))stop("Harvest output includes unplanned cells",call.=FALSE)
    eidx<-event_of_cell[at]
    ev$HarvestedSites<-tabulate(eidx[all_sites$removed_g_m2>0],nbins=nrow(ev))
    total<-rowsum(as.numeric(all_sites$removed_g_m2),group=eidx,reorder=TRUE)
    if(nrow(total))ev$MgBiomassRemoved[as.integer(rownames(total))]<-
      total[,1L]*(as.numeric(forest$cell_length)^2/1e6)
  }
  credited <- integer(nrow(plan$prescriptions))
  if (nrow(ev)) {
    initial <- if ("RepeatNumber" %in% names(ev)) ev$RepeatNumber == 0L else rep(TRUE, nrow(ev))
    totals <- rowsum(ev$HarvestedSites[initial], ev$prescription_row[initial])
    credited[as.integer(rownames(totals))] <- totals[,1L]
  }
  if (!identical(credited, plan$prescriptions$allocated_sites))
    stop("Effective quota credit does not match applied harvest", call.=FALSE)
  plant <- do.call(rbind, planting)
  if (NROW(plant)) rownames(plant) <- NULL
  exclude <- as.integer(unlist(excluded, use.names = FALSE))
  if (length(planted_species) || length(excluded)) {
    ev$plant_species <- I(lapply(ev$Prescription, function(pres)
      as.character(planted_species[[pres]])))
    ev$establishment_excluded <- I(Map(function(pres, cells)
      intersect(cells, excluded[[pres]]), ev$Prescription, ev$cells))
  }
  if (NROW(plant) || length(exclude))
    current <- landisR::landis_schedule_establishment(current, plant, exclude)
  if (!is.null(plan$landis_rng)) current$landis_rng <- plan$landis_rng
  list(simulation=current,plan=plan,events=ev,history=plan$history,
    effects=if(length(effects))do.call(rbind,effects)else data.frame(),
    sites=all_sites,landis_rng_offset_end=plan$landis_rng_offset_end)
}
