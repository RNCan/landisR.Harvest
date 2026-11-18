# H04-D2: site-selection Complete adjacency and harvest-history pilot.
# Qualified against LANDIS-II v7 H04-D1/H04-D1b references.
# A single sparse static rook adjacency is cached once per topology; no
# simulation/R RNG mutation and no whole-landscape journal per harvest.
# AdjacencyNeighborSetAside has no effect for SiteSelection Complete in v7.
harvest_prepare_adjacency_d2 <- function(simulation, topology) {
  .harvest_prepare_adjacency(simulation, topology)
}

.harvest_prepare_adjacency <- function(simulation, topology, landis_order = FALSE) {
  .h04_sim(simulation)
  forest <- landisR::landis_csr(simulation)
  if(!inherits(topology,"landisR_h04_topology") ||
     !identical(topology$ids,as.integer(forest$cell_id)))
    stop("Topology does not match current landscape",call.=FALSE)
  ans <- .Call("h04d2_topology",as.integer(topology$ids),topology$positions,
    topology$stand_offsets,as.integer(c(forest$nrow, forest$ncol, landis_order)),
    PACKAGE="landisR.harvest")
  topology$d2_neighbors <- ans
  topology
}

harvest_plan_d2 <- function(simulation, topology, targets, history=NULL, year=simulation$year,
                            economic=NULL, forest_type=NULL,
                            allocation_mode=c("corrected","landis_compatible")) {
  .h04_sim(simulation)
  forest <- landisR::landis_csr(simulation)
  allocation_mode <- match.arg(allocation_mode)
  if(is.null(topology$d2_neighbors) ||
     !identical(topology$ids,as.integer(forest$cell_id)))
    stop("Prepare cached adjacency first: harvest_prepare_adjacency_d2",call.=FALSE)
  if(!is.data.frame(targets) ||
     !all(c("management_area","ranking","fraction","target_type") %in% names(targets)))
    stop("targets must specify one ranking, fraction and target_type per MA",call.=FALSE)
  if(anyNA(targets$ranking) || any(!targets$ranking %in%
     c("MaxCohortAge","Economic","RegulateAges")))
    stop("H04-D2 adjacency pilot supports deterministic ranking only; Random coupling deferred",call.=FALSE)
  if("site_selection"%in%names(targets) &&
     (anyNA(targets$site_selection)||any(targets$site_selection!="Complete")))
    stop("H04-D2 is qualified for SiteSelection Complete only",call.=FALSE)
  if(!is.numeric(year)||length(year)!=1L||is.na(year)||!is.finite(year)||
     year<0||year>1e9||year!=trunc(year))stop("year must be a nonnegative integer",call.=FALSE)
  nma <- length(topology$management_areas)
  i <- match(topology$management_areas,targets$management_area)
  if(anyNA(i) || anyDuplicated(targets$management_area) || nrow(targets)!=nma)
    stop("Exactly one target per management area",call.=FALSE)
  param <- function(key,default){
    if(key%in%names(targets))targets[[key]][i] else rep(default,nma)
  }
  threshold <- .h04_integer(param("stand_adjacency",0),"stand_adjacency",0)
  min_history <- .h04_integer(param("minimum_time_since_last_harvest",0),
                             "minimum_time_since_last_harvest",0)
  # Explicit canonicalization prevents a disabled adjacency constraint from
  # failing on blank/NA input, while rejecting missing types for active limits.
  # The lookup never silently changes a positive threshold's semantics.
  raw_type <- as.character(param("adjacency_type","StandAge"))
  raw_type <- trimws(raw_type)
  blank <- is.na(raw_type) | tolower(raw_type) %in% c("", "na", "nan", "null")
  if(any(blank & threshold > 0L))
    stop("Missing adjacency_type for active stand_adjacency in management area(s): ",
         paste(topology$management_areas[blank & threshold > 0L],collapse=","),call.=FALSE)
  raw_type[blank] <- "StandAge"
  key <- tolower(gsub("[[:space:]_-]+", "", raw_type))
  allowed <- c(standage="StandAge", timesincelastharvested="TimeSinceLastHarvested",
               none="none")
  bad <- is.na(key) | !(key %in% names(allowed))
  if(any(bad))
    stop("Invalid adjacency_type for management area(s) ",
         paste(topology$management_areas[bad],collapse=","),
         ": ",paste(shQuote(raw_type[bad]),collapse=", "),
         "; allowed: StandAge, TimeSinceLastHarvested, none",call.=FALSE)
  type <- unname(allowed[key])
  setaside <- .h04_integer(param("adjacency_neighbor_set_aside",0),
                          "adjacency_neighbor_set_aside",0)
  # Under LANDIS-II v7 CIL, set-aside writes are gated by a spreading
  # SiteSelection. In this Complete-only planner it is inert by design.
  if(!is.null(history)) {
    if(!is.data.frame(history)||!all(c("stand","last_harvest_year")%in%names(history)))
      stop("history requires stand,last_harvest_year",call.=FALSE)
    old <- .h04_integer(history$stand,"history stand")
    if(anyDuplicated(old)||length(old)!=nrow(topology$stands)||
       !setequal(old,topology$stands$stand))stop("history does not match topology",call.=FALSE)
    last <- as.numeric(history$last_harvest_year[match(topology$stands$stand,old)])
    if(any(!is.na(last) & (!is.finite(last)|last<0|last>year|last!=trunc(last))))
      stop("Invalid history year",call.=FALSE)
    last <- as.integer(last)
  } else last <- rep.int(NA_integer_,nrow(topology$stands))
  work <- targets
  work$fraction <- rep.int(0,nma)
  if("cohorts_removed"%in%names(work))work$cohorts_removed<-NULL
  # C2 does the cohort eligibility and ranking on a single CURRENT CSR state.
  # Set quota to zero so this baseline pass selects no stand and no cell.
  base <- harvest_plan_c2(simulation,topology,work,economic,forest_type,allocation_mode)
  target_type <- as.character(targets$target_type[i]);fraction <- as.numeric(targets$fraction[i])
  native <- .Call("h04d2_plan",topology$stand_offsets,topology$ma_offsets,
    topology$d2_neighbors$offsets,topology$d2_neighbors$indices,
    topology$d2_neighbors$first_position,
    as.integer(base$stands$stand_age),as.numeric(base$stands$rank),
    as.logical(base$stands$qualified),as.numeric(fraction),
    as.integer(match(target_type,c("area","stands"))-1L),
    rep.int(as.integer(identical(allocation_mode,"landis_compatible")),nma),
    threshold,as.integer(match(type,c("none","StandAge","TimeSinceLastHarvested"))-1L),
    min_history,last,as.numeric(year),PACKAGE="landisR.harvest")
  base$stands$selected <- seq_len(nrow(base$stands))%in%native$selected_stands
  base$stands$spatially_eligible <- as.logical(native$stand_eligible)
  base$selected_stands <- base$stands[native$selected_stands,,drop=FALSE]
  selected_positions <- unlist(lapply(native$selected_stands,function(s){
    j<-seq.int(topology$stand_offsets[s]+1L,topology$stand_offsets[s+1L]);topology$positions[j]+1L
  }),use.names=FALSE)
  base$cells<-topology$ids[selected_positions]
  base$areas<-data.frame(management_area=topology$management_areas,
    ranking=as.character(targets$ranking[i]),target_fraction=fraction,
    target_type=target_type,
    target_sites=ifelse(target_type=="area",fraction*native$total,NA_real_),
    target_stands=ifelse(target_type=="stands",native$target_stands,NA_integer_),
    allocated_sites=native$allocated,allocated_stands=native$selected_count,
    total_sites=native$total,eligible_stands=native$eligible)
  base$history<-data.frame(management_area=topology$stands$management_area,
    stand=topology$stands$stand,last_harvest_year=as.integer(native$last_harvest_year))
  base$year<-as.integer(year)
  base$adjacency<-data.frame(management_area=topology$management_areas,
    stand_adjacency=threshold,adjacency_type=type,
    minimum_time_since_last_harvest=min_history,
    adjacency_neighbor_set_aside=setaside,
    setaside_effective=FALSE)
  base$type<-"H04-D2 historical Complete adjacency/history; single prescription per MA; deterministic"
  base
}

harvest_step_d2 <- function(simulation,topology,targets,history=NULL,year=simulation$year,
                            rules=NULL,clearcut=FALSE,economic=NULL,forest_type=NULL,
                            allocation_mode=c("corrected","landis_compatible")) {
  if(!isTRUE(clearcut)&&is.null(rules))
    stop("Specify clearcut=TRUE or cohort-removal rules",call.=FALSE)
  plan <- harvest_plan_d2(simulation,topology,targets,history,year,economic,
                          forest_type,allocation_mode)
  ans <- harvest_apply_cells(simulation,plan$cells,rules=rules,clearcut=clearcut)
  if(length(plan$cells)&&!setequal(plan$cells,ans$changed_cells))
    stop("The selected treatment did not affect every planned cell",call.=FALSE)
  ans$plan<-plan
  ans$history<-plan$history
  ans
}
