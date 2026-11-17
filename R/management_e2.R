# H04-E2: multiple prescriptions, weighted remaining-target random selection.
# Pilot is deliberately scoped to MaxCohortAge + SiteSelection Complete.
# A separate MT19937 stream matches isolated LANDIS-II v7 at the supplied
# seed and draw offset. The landisR succession RNG is NEVER modified.
harvest_plan_e2 <- function(simulation, topology, prescriptions,
                            year=simulation$year,landis_seed,
                            landis_rng_offset=0,
                            allocation_mode=c("corrected","landis_compatible")) {
  .h04_sim(simulation)
  forest <- landisR::landis_csr(simulation)
  allocation_mode<-match.arg(allocation_mode)
  if(missing(landis_seed))stop("landis_seed is required",call.=FALSE)
  valid_scalar<-function(x,lo,hi){is.numeric(x)&&length(x)==1L&&!is.na(x)&&is.finite(x)&&x>=lo&&x<=hi&&x==trunc(x)}
  if(!valid_scalar(landis_seed,0,4294967295))stop("landis_seed must be a uint32",call.=FALSE)
  if(!valid_scalar(landis_rng_offset,0,1e9))stop("Invalid landis_rng_offset",call.=FALSE)
  if(!valid_scalar(year,0,1e9))stop("Invalid harvest year",call.=FALSE)
  if(!inherits(topology,"landisR_h04_topology") ||
     !identical(topology$ids,as.integer(forest$cell_id)))
     stop("Topology does not match current landscape",call.=FALSE)
  required<-c("management_area","prescription","fraction","ranking","target_type",
               "minimum_age","maximum_age","begin_time","end_time","site_selection")
  if(!is.data.frame(prescriptions)||!all(required%in%names(prescriptions))||!nrow(prescriptions))
     stop("prescriptions must provide management_area, prescription, fraction, ranking, target_type, minimum_age, maximum_age, begin_time, end_time, site_selection",call.=FALSE)
  p<-prescriptions
  p$management_area<-.h04_integer(p$management_area,"management_area")
  p$minimum_age<-.h04_integer(p$minimum_age,"minimum_age",0L)
  p$maximum_age<-.h04_integer(p$maximum_age,"maximum_age",0L)
  p$begin_time<-.h04_integer(p$begin_time,"begin_time",0L)
  p$end_time<-.h04_integer(p$end_time,"end_time",0L)
  p$prescription<-as.character(p$prescription)
  p$ranking<-as.character(p$ranking)
  p$target_type<-as.character(p$target_type)
  p$site_selection<-as.character(p$site_selection)
  if(anyNA(p$prescription)||any(!nzchar(p$prescription))||
      anyNA(p$ranking)||any(p$ranking!="MaxCohortAge")||
      anyNA(p$site_selection)||any(p$site_selection!="Complete")||
      anyNA(p$target_type)||any(p$target_type!="area"))
     stop("E2 pilot supports named MaxCohortAge, SiteSelection Complete, area quotas only",call.=FALSE)
  if(any(p$maximum_age<p$minimum_age)||any(p$end_time<p$begin_time)||
     !is.numeric(p$fraction)||anyNA(p$fraction)||any(!is.finite(p$fraction))||
     any(p$fraction<0|p$fraction>1)||
     anyDuplicated(p[c("management_area","prescription")]))
     stop("Invalid prescription ages/periods/fractions or repeated MA+prescription",call.=FALSE)
  if(any(!p$management_area%in%topology$management_areas) ||
     any(!topology$management_areas%in%p$management_area))
     stop("Every MA needs at least one known prescription",call.=FALSE)
  if(identical(allocation_mode,"landis_compatible")&&any(p$fraction==0))
     stop("LANDIS-II v7 rejects explicitly specified 0% HarvestImplementations rows",call.=FALSE)
  # ONE CSR cohort pass to obtain the age/ranking of every stand (no harvest).
  baseline<-data.frame(management_area=topology$management_areas,
    ranking="MaxCohortAge",fraction=0,target_type="area",
    minimum_age=0L,maximum_age=65535L,site_selection="Complete")
  base<-harvest_plan_c2(simulation,topology,baseline,allocation_mode=allocation_mode)
  nr<-as.integer(forest$nrow);nc<-as.integer(forest$ncol)
  first<-integer(nrow(topology$stands))
  for(s in seq_along(first)){
    i<-seq.int(topology$stand_offsets[s]+1L,topology$stand_offsets[s+1L])
    ids<-topology$ids[topology$positions[i]+1L]
    first[s]<-as.integer(min(((ids-1L)%%nr)*nc+(ids-1L)%/%nr))
  }
  nat<-.Call("h04e2_plan",topology$stand_offsets,topology$ma_offsets,
    first,as.integer(base$stands$stand_age),as.numeric(base$stands$rank),
    as.integer(match(p$management_area,topology$management_areas)-1L),
    as.numeric(p$fraction),rep.int(0L,nrow(p)),
    as.integer(p$minimum_age),as.integer(p$maximum_age),
    as.integer(p$begin_time),as.integer(p$end_time),
    as.numeric(year),as.numeric(landis_seed),as.numeric(landis_rng_offset),
    identical(allocation_mode,"landis_compatible"),PACKAGE="landisR.harvest")
  s<-as.integer(nat$selected_stands);j<-as.integer(nat$prescription_indices)
  events<-data.frame(EventID=seq_along(s),Time=rep.int(as.integer(year),length(s)),
    ManagementArea=topology$stands$management_area[s],
    Prescription=p$prescription[j],Stand=topology$stands$stand[s],
    StandAge=as.integer(nat$ages),StandRank=as.numeric(base$stands$rank[s]),
    NumberOfSites=topology$stands$sites[s],
    rng_draw_index=as.numeric(nat$draw_indices),
    rng_uniform=as.numeric(nat$draw_uniform),
    selection_probability=as.numeric(nat$draw_probability),stringsAsFactors=FALSE)
  events$prescription_row<-j
  locations<-vector("list",length(s))
  for(i in seq_along(s)){
    k<-seq.int(topology$stand_offsets[s[i]]+1L,topology$stand_offsets[s[i]+1L])
    locations[[i]]<-topology$ids[topology$positions[k]+1L]
  }
  events$cells<-I(locations)
  selected_cells<-as.integer(unlist(locations,use.names=FALSE))
  prescription_summary<-p
  prescription_summary$allocated_sites<-as.integer(nat$allocated_sites)
  prescription_summary$allocated_stands<-as.integer(nat$allocated_stands)
  prescription_summary$eligible_stands<-as.integer(nat$eligible_stands)
  prescription_summary$target_sites<-as.numeric(p$fraction)*as.numeric(nat$total_sites[match(p$management_area,topology$management_areas)])
  base$stands$selected<-seq_len(nrow(base$stands))%in%s
  base$events<-events
  base$selected_stands<-base$stands[s,,drop=FALSE]
  base$cells<-selected_cells
  base$prescriptions<-prescription_summary
  base$draws_consumed<-as.numeric(nat$draws_consumed)
  base$landis_rng_offset_start<-as.numeric(landis_rng_offset)
  base$landis_rng_offset_end<-as.numeric(nat$rng_offset_end)
  base$landis_seed<-as.numeric(landis_seed)
  base$year<-as.integer(year)
  base$type<-"H04-E2 Multi-prescription MaxCohortAge Complete area pilot; stochastic coupled RNG not yet qualified"
  base
}

# Applies disjoint stands by prescription using the already validated H03
# living-cohort engine. This does not duplicate the landscape per event.
harvest_step_e2 <- function(simulation,topology,prescriptions,rules_by_prescription,
                            year=simulation$year,landis_seed,landis_rng_offset=0,
                            allocation_mode=c("corrected","landis_compatible")) {
  plan<-harvest_plan_e2(simulation,topology,prescriptions,year,
             landis_seed,landis_rng_offset,allocation_mode)
  names_needed<-unique(as.character(plan$events$Prescription))
  if(!is.list(rules_by_prescription)||is.null(names(rules_by_prescription))||
     any(!names_needed%in%names(rules_by_prescription)))
     stop("rules_by_prescription requires one named cohort-removal data frame per selected prescription",call.=FALSE)
  updated<-simulation;effects<-list();sites<-list()
  for(pres in names_needed){
    rows<-which(plan$events$Prescription==pres)
    cells<-as.integer(unlist(plan$events$cells[rows],use.names=FALSE))
    r<-rules_by_prescription[[pres]]
    result<-harvest_apply_cells(updated,cells,rules=r,clearcut=FALSE)
    updated<-result$simulation
    effects[[pres]]<-result$effects
    z<-result$sites;z$prescription<-pres;sites[[pres]]<-z
  }
  list(simulation=updated,plan=plan,effects=if(length(effects))do.call(rbind,effects)else data.frame(),
       sites=if(length(sites))do.call(rbind,sites)else data.frame(),
       landis_rng_offset_end=plan$landis_rng_offset_end)
}
