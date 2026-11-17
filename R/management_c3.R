# H04-C3b: LANDIS-II v7 Random ranking with a seed + shared-stream offset.
# Exact on the validated small Complete/ClearCut LANDIS-II scenarios.
# Coupled exact RNG requires that ALL prior core/succession/extension draws
# are accounted for; landisR's simulation RNG must never be overwritten.
harvest_plan_c3 <- function(simulation, topology, targets,
                            economic=NULL, forest_type=NULL,
                            allocation_mode=c("corrected","landis_compatible"),
                            landis_seed, landis_rng_offset=0) {
  forest <- landisR::landis_csr(simulation)
  allocation_mode <- match.arg(allocation_mode)
  if(missing(landis_seed)) stop("landis_seed is required to reproduce LANDIS-II Random",call.=FALSE)
  if(!is.numeric(landis_seed) || length(landis_seed)!=1L || is.na(landis_seed) ||
      !is.finite(landis_seed) || landis_seed<0 || landis_seed>4294967295 ||
      landis_seed!=trunc(landis_seed)) stop("landis_seed must be a uint32",call.=FALSE)
  if(!is.numeric(landis_rng_offset) || length(landis_rng_offset)!=1L ||
      is.na(landis_rng_offset) || !is.finite(landis_rng_offset) ||
      landis_rng_offset<0 || landis_rng_offset>1e9 ||
      landis_rng_offset!=trunc(landis_rng_offset)) stop("Invalid landis_rng_offset",call.=FALSE)
  if(!is.data.frame(targets) || !all(c("management_area","ranking","fraction","target_type")%in%names(targets)))
    stop("targets must specify one ranking and target per MA",call.=FALSE)
  chosen <- as.character(targets$ranking)
  if(anyNA(chosen)||any(!chosen%in%c("MaxCohortAge","Economic","RegulateAges","Random")))
    stop("Unknown C3 ranking",call.=FALSE)
  # C2 performs validated CSR cohort qualification, ForestType, cohort ages,
  # economic weighting, topology checks, and age bounds. It is called only
  # with zero allocation, not for any harvest or random number generation.
  work <- targets
  work$ranking[work$ranking=="Random"] <- "MaxCohortAge"
  work$fraction[] <- 0
  base <- harvest_plan_c2(simulation,topology,work,economic,forest_type,allocation_mode)
  ord <- match(topology$management_areas,targets$management_area)
  ranking <- as.integer(match(targets$ranking[ord],
                    c("MaxCohortAge","Economic","RegulateAges","Random"))-1L)
  fraction <- as.numeric(targets$fraction[ord]); target_type <- as.character(targets$target_type[ord])
  if(anyNA(ranking)||anyNA(fraction)||any(!is.finite(fraction))||any(fraction<0|fraction>1))stop("Invalid target",call.=FALSE)
  if(anyNA(target_type)||any(!target_type%in%c("area","stands")))stop("Invalid target_type",call.=FALSE)
  nr <- as.integer(forest$nrow);nc <- as.integer(forest$ncol)
  first <- integer(nrow(topology$stands))
  for(s in seq_along(first)) {
    p <- seq.int(topology$stand_offsets[s]+1L,topology$stand_offsets[s+1L])
    cells <- topology$ids[topology$positions[p]+1L]
    first[s] <- as.integer(min(((cells-1L)%%nr)*nc + (cells-1L)%/%nr))
  }
  native <- .Call("h04c3_random_plan",first,topology$stand_offsets,topology$ma_offsets,
    base$stands$qualified,as.numeric(base$stands$rank),ranking,fraction,
    as.integer(match(target_type,c("area","stands"))-1L),
    rep.int(as.integer(identical(allocation_mode,"landis_compatible")),length(ranking)),
    as.numeric(landis_seed),as.numeric(landis_rng_offset),PACKAGE="landisR.harvest")
  base$stands$rank <- as.numeric(native$stand_rank)
  base$stands$selected <- seq_len(nrow(base$stands))%in%native$selected_stands
  base$selected_stands <- base$stands[native$selected_stands,,drop=FALSE]
  pos <- unlist(lapply(native$selected_stands,function(s){
    i <- seq.int(topology$stand_offsets[s]+1L,topology$stand_offsets[s+1L]);topology$positions[i]+1L
  }),use.names=FALSE)
  base$cells <- topology$ids[pos]
  base$areas <- data.frame(management_area=topology$management_areas,
    ranking=chosen[ord],target_fraction=fraction,target_type=target_type,
    target_sites=ifelse(target_type=="area",fraction*native$total,NA_real_),
    target_stands=ifelse(target_type=="stands",native$target_stands,NA_integer_),
    allocated_sites=native$allocated,allocated_stands=native$selected_count,
    total_sites=native$total,eligible_stands=native$eligible)
  base$landis_rng <- list(seed=as.numeric(landis_seed),
    start_offset=as.numeric(native$rng_start_offset),
    draws_used=as.numeric(native$rng_draws_used),
    next_offset=as.numeric(native$rng_next_offset),
    scope="LANDIS-II v7 seed+offset pilot; NOT landisR succession RNG")
  base$type <- "H04-C3b Complete/ClearCut pilot; LANDIS-II v7 Random MT19937 parity on qualified fixtures"
  base
}

harvest_complete_c3 <- function(simulation,topology,targets,economic=NULL,forest_type=NULL,
                                allocation_mode=c("corrected","landis_compatible"),
                                landis_seed,landis_rng_offset=0) {
  plan <- harvest_plan_c3(simulation,topology,targets,economic,forest_type,
                          allocation_mode,landis_seed,landis_rng_offset)
  out <- harvest_apply_cells(simulation,plan$cells,clearcut=TRUE)
  if(!setequal(out$changed_cells,plan$cells))stop("Unexpected zero removal",call.=FALSE)
  out$plan <- plan
  out
}
