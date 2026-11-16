# H04-C2 candidate: deterministic Complete / ClearCut, one prescription per MA.
# Adds RegulateAges, inclusive MaximumAge, and ForestTypeTable.
# Random is deliberately not supported until .NET RNG equivalence is qualified.
# Use CURRENT CSR cohort ages. No LANDIS-II initialization-age correction here.
harvest_plan_c2 <- function(simulation, topology, targets,
                            economic=NULL, forest_type=NULL,
                            allocation_mode=c("corrected","landis_compatible")) {
  .h04_sim(simulation)
  allocation_mode <- match.arg(allocation_mode)
  if(!inherits(topology,"landisR_h04_topology") ||
     !identical(topology$ids,as.integer(simulation$ctx$ids)))
    stop("Topology does not match simulation",call.=FALSE)
  if(!is.data.frame(targets) ||
     !all(c("management_area","fraction","ranking","target_type")%in%names(targets)))
    stop("targets requires management_area, fraction, ranking, target_type",call.=FALSE)
  ma <- .h04_integer(targets$management_area,"management_area")
  nma <- length(topology$management_areas)
  if(anyDuplicated(ma)||!setequal(ma,topology$management_areas))
    stop("Exactly one target per management area required",call.=FALSE)
  j <- match(topology$management_areas,ma)
  f <- as.numeric(targets$fraction[j])
  if(!is.numeric(targets$fraction)||anyNA(f)||any(!is.finite(f))||any(f<0|f>1))
    stop("fraction must be in [0,1]",call.=FALSE)
  minimum <- if("minimum_age"%in%names(targets))
    .h04_integer(targets$minimum_age[j],"minimum_age",0) else rep.int(0L,nma)
  maximum <- if("maximum_age"%in%names(targets))
    .h04_integer(targets$maximum_age[j],"maximum_age",0) else rep.int(65535L,nma)
  if(any(maximum<minimum))stop("maximum_age must be >= minimum_age",call.=FALSE)
  ranking <- as.character(targets$ranking[j]); target_type <- as.character(targets$target_type[j])
  if(anyNA(ranking)||any(!ranking%in%c("MaxCohortAge","Economic","RegulateAges")))
    stop("H04-C2 supports only MaxCohortAge, Economic, RegulateAges. Random RNG remains unqualified.",call.=FALSE)
  if(anyNA(target_type)||any(!target_type%in%c("area","stands")))
    stop("target_type must be area or stands",call.=FALSE)
  if("site_selection"%in%names(targets) &&
     (anyNA(targets$site_selection)||any(targets$site_selection!="Complete")))
    stop("H04-C2 is qualified for SiteSelection Complete only",call.=FALSE)
  if("cohorts_removed"%in%names(targets) &&
     (anyNA(targets$cohorts_removed)||any(targets$cohorts_removed!="ClearCut")))
    stop("H04-C2 is qualified for ClearCut only",call.=FALSE)

  species <- as.character(simulation$ctx$config$species$species)
  nsp <- length(species)
  weights <- matrix(0,nrow=nma,ncol=nsp)
  merchant <- matrix(1L,nrow=nma,ncol=nsp)
  if(any(ranking=="Economic")) {
    if(!is.data.frame(economic)||
       !all(c("species","rank","merchantable_age")%in%names(economic)))
      stop("Economic requires species,rank,merchantable_age",call.=FALSE)
    if(anyNA(economic$species)||any(!economic$species%in%species))
      stop("Unknown Economic species",call.=FALSE)
    er <- as.numeric(economic$rank)
    if(anyNA(er)||any(!is.finite(er))||any(er!=trunc(er)|er<0|er>100))
      stop("Invalid Economic rank",call.=FALSE)
    em <- .h04_integer(economic$merchantable_age,"merchantable_age",1)
    if("management_area"%in%names(economic)) {
      mea <- .h04_integer(economic$management_area,"economic management_area")
      if(any(!mea%in%topology$management_areas)||
         anyDuplicated(data.frame(mea,economic$species)))
        stop("Invalid Economic species/management area",call.=FALSE)
      ii <- match(mea,topology$management_areas)
      jj <- match(economic$species,species)
      weights[cbind(ii,jj)] <- er
      merchant[cbind(ii,jj)] <- em
      if(any(!topology$management_areas[ranking=="Economic"]%in%mea))
        stop("Economic species table missing for a management area",call.=FALSE)
    } else {
      if(anyDuplicated(economic$species))stop("Duplicate Economic species",call.=FALSE)
      for(k in seq_len(nma)) {
        weights[k,match(economic$species,species)] <- er
        merchant[k,match(economic$species,species)] <- em
      }
    }
  }

  rule_ma <- rule_kind <- rule_lo <- rule_hi <- rule_highest <- pair_rule <- pair_sp <- integer()
  rule_pct <- numeric()
  if(!is.null(forest_type)) {
    if(!is.data.frame(forest_type)||
       !all(c("management_area","rule_id","inclusion","species","age_min","age_max","percent","highest")%in%names(forest_type)))
      stop("forest_type requires management_area,rule_id,inclusion,species,age_min,age_max,percent,highest",call.=FALSE)
    fm <- .h04_integer(forest_type$management_area,"forest_type management_area")
    rid <- .h04_integer(forest_type$rule_id,"forest_type rule_id")
    lo <- .h04_integer(forest_type$age_min,"forest_type age_min")
    hi <- .h04_integer(forest_type$age_max,"forest_type age_max")
    inc <- as.character(forest_type$inclusion)
    fs <- as.character(forest_type$species)
    highest <- as.logical(forest_type$highest)
    percent <- as.numeric(forest_type$percent)
    if(any(!fm%in%topology$management_areas)||any(hi<lo)|any(hi>65535L)||
       anyNA(inc)||any(!inc%in%c("Required","Optional","Forbidden"))||
       anyNA(fs)||any(!fs%in%species)||anyNA(highest)||
       anyNA(percent)||any(!is.finite(percent))||any(percent<=0|percent>100))
      stop("Invalid forest_type row",call.=FALSE)
    key <- paste(fm,rid,sep="::")
    distinct <- unique(key)
    rule_idx <- match(key,distinct)
    first <- match(distinct,key)
    for(col in list(inc,lo,hi,percent,highest)) if(any(vapply(
        split(col,rule_idx),function(x)length(unique(x))!=1L,logical(1))))
      stop("Each forest_type rule_id must have coherent attributes",call.=FALSE)
    if(anyDuplicated(data.frame(key,fs)))
      stop("Duplicate species in forest_type rule",call.=FALSE)
    rule_ma <- as.integer(match(fm[first],topology$management_areas)-1L)
    rule_kind <- as.integer(match(inc[first],c("Required","Optional","Forbidden"))-1L)
    rule_lo <- lo[first];rule_hi <- hi[first];rule_pct <- percent[first]
    rule_highest <- as.integer(highest[first])
    pair_rule <- as.integer(rule_idx-1L)
    pair_sp <- as.integer(match(fs,species)-1L)
  }

  native <- .Call("h04c2_plan",simulation$states,as.integer(simulation$cell_state),
    topology$positions,topology$stand_offsets,topology$ma_offsets,
    weights,merchant,as.integer(match(ranking,c("MaxCohortAge","Economic","RegulateAges"))-1L),
    f,minimum,maximum,as.integer(match(target_type,c("area","stands"))-1L),
    rep.int(as.integer(identical(allocation_mode,"landis_compatible")),nma),
    rule_ma,rule_kind,rule_lo,rule_hi,rule_pct,rule_highest,pair_rule,pair_sp,
    as.integer(c(simulation$ctx$nrow,simulation$ctx$ncol)),
    as.integer(simulation$ctx$ids),
    PACKAGE="landisR.harvest")
  stands <- topology$stands
  stands$stand_age <- as.integer(native$stand_age)
  stands$rank <- as.numeric(native$stand_rank)
  stands$qualified <- as.logical(native$stand_qualified)
  stands$selected <- seq_len(nrow(stands))%in%native$selected_stands
  areas <- data.frame(management_area=topology$management_areas,
    ranking=ranking,target_fraction=f,target_type=target_type,
    target_sites=ifelse(target_type=="area",f*native$total,NA_real_),
    target_stands=ifelse(target_type=="stands",native$target_stands,NA_integer_),
    allocated_sites=native$allocated,allocated_stands=native$selected_count,
    total_sites=native$total,eligible_stands=native$eligible)
  list(stands=stands,selected_stands=stands[native$selected_stands,,drop=FALSE],
       cells=topology$ids[native$selected_positions],areas=areas,
       states_evaluated=native$states_evaluated,
       allocation_mode=allocation_mode,
       type="H04-C2 deterministic Complete/ClearCut pilot; Random not supported")
}

harvest_complete_c2 <- function(simulation,topology,targets,economic=NULL,
                                forest_type=NULL,
                                allocation_mode=c("corrected","landis_compatible")) {
  plan <- harvest_plan_c2(simulation,topology,targets,economic,forest_type,allocation_mode)
  out <- harvest_apply_cells(simulation,plan$cells,clearcut=TRUE)
  if(!setequal(out$changed_cells,plan$cells))
    stop("Unexpected zero removal in selected stand",call.=FALSE)
  out$plan <- plan
  out
}
