# H04-A2: one prescription per MA, MaxCohortAge, MinimumAge, Complete, ClearCut.
# Input maps have already been aligned to the core canonical cell IDs.
.h04_integer <- function(x, label, min=1) {
  if(!is.numeric(x)||anyNA(x)||any(!is.finite(x))||any(x!=trunc(x))||any(x<min|x>.Machine$integer.max)) stop("Invalid ",label,call.=FALSE)
  as.integer(x)
}
.h04_sim <- function(s) {
  if (!inherits(s, "landisR_simulation") || !identical(s$format, "landisR/simulation-v1"))
    stop("Harvest requires a landisR simulation with format landisR/simulation-v1.", call. = FALSE)
}
.harvest_check_management_columns <- function(management) {
  if (!is.data.frame(management) ||
      !all(c("cell_id", "management_area", "stand") %in% names(management)))
    stop("management requires cell_id, management_area, stand", call. = FALSE)
}
harvest_prepare_management <- function(simulation, management) {
  .h04_sim(simulation)
  forest <- landisR::landis_csr(simulation)
  .harvest_check_management_columns(management)
  ids <- as.integer(forest$cell_id)
  if(!length(ids)) stop("Empty landscape")
  cell <- .h04_integer(management$cell_id,"cell_id"); ma <- .h04_integer(management$management_area,"management_area"); st <- .h04_integer(management$stand,"stand")
  if(length(cell)!=length(ids)||anyDuplicated(cell)||!setequal(cell,ids)) stop("Exactly one management row per active cell is required")
  # Each stand belongs to exactly one management area.
  pairs <- unique(data.frame(stand=st,management_area=ma)); if(anyDuplicated(pairs$stand)) stop("Stand spans multiple management areas")
  p <- match(cell,ids); o <- order(ma,st,p); ma <- ma[o]; st <- st[o]; p <- p[o]
  start <- c(TRUE,diff(st)!=0L | diff(ma)!=0L); starts <- which(start)
  counts <- diff(c(starts,length(st)+1L))
  stands <- data.frame(management_area=ma[starts],stand=st[starts],sites=counts)
  ma_start <- which(c(TRUE,diff(stands$management_area)!=0L))
  structure(list(ids=ids,positions=as.integer(p-1L),stand_offsets=as.integer(c(0,cumsum(counts))),
    ma_offsets=as.integer(c(ma_start-1L,nrow(stands))),stands=stands,management_areas=stands$management_area[ma_start]),class="landisR_h04_topology")
}
harvest_plan_complete <- function(simulation, topology, targets) {
  .h04_sim(simulation)
  forest <- landisR::landis_csr(simulation)
  if(!inherits(topology,"landisR_h04_topology")||!identical(topology$ids,as.integer(forest$cell_id))) stop("Topology does not match active cell IDs/order")
  if(!is.data.frame(targets)||!all(c("management_area","fraction","minimum_age")%in%names(targets))) stop("targets requires management_area, fraction, minimum_age")
  ma <- .h04_integer(targets$management_area,"target management_area"); minimum <- .h04_integer(targets$minimum_age,"minimum_age",0)
  f <- targets$fraction
  if(!is.numeric(f)||anyNA(f)||any(!is.finite(f))||any(f<0|f>1)||anyDuplicated(ma)||!setequal(ma,topology$management_areas)) stop("One target per MA; fraction must be in [0,1]")
  k <- match(topology$management_areas,ma)
  n <- .Call("h04_plan",forest$states,as.integer(forest$cell_state),topology$positions,topology$stand_offsets,topology$ma_offsets,as.numeric(f[k]),minimum[k],PACKAGE="landisR.harvest")
  stands <- topology$stands; stands$stand_age <- n$stand_age; stands$rank <- n$stand_age
  stands$selected <- seq_len(nrow(stands)) %in% n$selected_stands
  list(stands=stands,selected_stands=stands[n$selected_stands,,drop=FALSE],cells=topology$ids[n$selected_positions],
    areas=data.frame(management_area=topology$management_areas,target_fraction=f[k],target_sites=f[k]*n$total,allocated_sites=n$allocated,total_sites=n$total,eligible_stands=n$eligible),states_evaluated=n$states_evaluated)
}
harvest_complete <- function(simulation, topology, targets) {
  plan <- harvest_plan_complete(simulation,topology,targets)
  result <- harvest_apply_cells(simulation,plan$cells,clearcut=TRUE)
  # Selected sites must all yield positive removal in this restricted pilot.
  if(!setequal(result$changed_cells,plan$cells)) stop("Unexpected zero removal in selected stand")
  result$plan <- plan
  result
}
