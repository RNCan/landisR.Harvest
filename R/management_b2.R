# H04-B2: One prescription per management area, Complete, ClearCut.
# Ranking is MaxCohortAge or Economic; target is area or eligible stands.
# Never normalise cohort ages here: use CURRENT CSR ages.
harvest_plan <- function(simulation, topology, targets, economic=NULL,
                         allocation_mode=c("corrected","landis_compatible")) {
  .h04_sim(simulation)
  allocation_mode <- match.arg(allocation_mode)
  if (!inherits(topology,"landisR_h04_topology") ||
      !identical(topology$ids,as.integer(simulation$ctx$ids)))
    stop("Topology does not match active cell IDs/order",call.=FALSE)
  if (!is.data.frame(targets) ||
      !all(c("management_area","fraction","minimum_age","ranking","target_type") %in% names(targets)))
    stop("targets must have management_area, fraction, minimum_age, ranking, target_type",call.=FALSE)
  ma <- .h04_integer(targets$management_area,"target management_area")
  minimum <- .h04_integer(targets$minimum_age,"minimum_age",0)
  f <- targets$fraction
  if(!is.numeric(f)||anyNA(f)||any(!is.finite(f))||any(f<0|f>1)||
     anyDuplicated(ma)||!setequal(ma,topology$management_areas))
    stop("One valid target per management area, fraction in [0,1]",call.=FALSE)
  ranking <- as.character(targets$ranking)
  target_type <- as.character(targets$target_type)
  if(anyNA(ranking)||any(!ranking %in% c("MaxCohortAge","Economic")))
    stop("H04-B2 ranking supports only MaxCohortAge or Economic",call.=FALSE)
  if(anyNA(target_type)||any(!target_type %in% c("area","stands")))
    stop("H04-B2 target_type supports area or stands",call.=FALSE)
  if("site_selection" %in% names(targets) &&
     (anyNA(targets$site_selection)||any(targets$site_selection != "Complete")))
    stop("H04-B2 supports SiteSelection Complete only",call.=FALSE)
  if("cohorts_removed" %in% names(targets) &&
     (anyNA(targets$cohorts_removed)||any(targets$cohorts_removed != "ClearCut")))
    stop("H04-B2 supports ClearCut only",call.=FALSE)
  order_ma <- match(topology$management_areas,ma)
  species <- as.character(simulation$ctx$config$species$species)
  nma <- length(order_ma);nsp <- length(species)
  weights <- matrix(0,nrow=nma,ncol=nsp)
  merchant <- matrix(1L,nrow=nma,ncol=nsp) # unused species have zero economic rank
  if(any(ranking=="Economic")) {
    if(!is.data.frame(economic)||!all(c("species","rank","merchantable_age")%in%names(economic)))
      stop("Economic requires economic data.frame(species,rank,merchantable_age[,management_area])",call.=FALSE)
    if(!is.character(economic$species)||anyNA(economic$species)||any(!economic$species%in%species))
      stop("Unknown economic species",call.=FALSE)
    r <- economic$rank
    if(!is.numeric(r)||anyNA(r)||any(!is.finite(r))||any(r!=trunc(r)|r<0|r>100))
      stop("Economic ranks must be integer in [0,100]",call.=FALSE)
    eage <- .h04_integer(economic$merchantable_age,"merchantable_age",1)
    if("management_area" %in% names(economic)) {
      ema <- .h04_integer(economic$management_area,"economic management_area")
      if(any(!ema%in%topology$management_areas)||anyDuplicated(data.frame(ema,economic$species)))
        stop("Invalid/duplicate economic species by management area",call.=FALSE)
      target_m <- match(ema,topology$management_areas)
      target_s <- match(economic$species,species)
      weights[cbind(target_m,target_s)] <- r
      merchant[cbind(target_m,target_s)] <- eage
      if(any(!topology$management_areas[which(ranking[order_ma]=="Economic")] %in% ema))
        stop("Missing economic table for Economic management area",call.=FALSE)
    } else {
      if(anyDuplicated(economic$species)) stop("Duplicate economic species",call.=FALSE)
      for(m in seq_len(nma)) {
        weights[m,match(economic$species,species)] <- r
        merchant[m,match(economic$species,species)] <- eage
      }
    }
  }
  rn <- ranking[order_ma]; tn <- target_type[order_ma]
  native <- .Call("h04b2_plan",simulation$states,as.integer(simulation$cell_state),
    topology$positions,topology$stand_offsets,topology$ma_offsets,
    weights,merchant,as.integer(match(rn,c("MaxCohortAge","Economic"))-1L),
    as.numeric(f[order_ma]),minimum[order_ma],
    as.integer(match(tn,c("area","stands"))-1L),
    as.integer(rep(identical(allocation_mode,"landis_compatible"),nma)),
    PACKAGE="landisR.harvest")
  stands <- topology$stands
  stands$stand_age <- as.integer(native$stand_age)
  stands$rank <- as.numeric(native$stand_rank)
  stands$selected <- seq_len(nrow(stands)) %in% native$selected_stands
  areas <- data.frame(management_area=topology$management_areas,
    target_fraction=f[order_ma],target_type=tn,
    target_sites=ifelse(tn=="area",f[order_ma]*native$total,NA_real_),
    target_stands=ifelse(tn=="stands",native$target_stands,NA_integer_),
    allocated_sites=native$allocated,allocated_stands=native$selected_count,
    total_sites=native$total,eligible_stands=native$eligible)
  list(stands=stands,selected_stands=stands[native$selected_stands,,drop=FALSE],
    cells=topology$ids[native$selected_positions],areas=areas,
    states_evaluated=native$states_evaluated,
    allocation_mode=allocation_mode,
    type="H04-B2 deterministic Complete ClearCut (no harvest scheduler)")
}

harvest_complete_b2 <- function(simulation,topology,targets,economic=NULL,
                                allocation_mode=c("corrected","landis_compatible")) {
  plan <- harvest_plan(simulation,topology,targets,economic,allocation_mode)
  out <- harvest_apply_cells(simulation,plan$cells,clearcut=TRUE)
  if(!setequal(out$changed_cells,plan$cells))
    stop("Unexpected zero removal in selected stand",call.=FALSE)
  out$plan <- plan
  out
}
