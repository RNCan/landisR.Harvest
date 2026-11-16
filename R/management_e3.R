# H04-E3 initial native qualification: multiple stand prescriptions PLUS
# per-prescription rook adjacency/history and one Required ForestType rule.
# Exact LANDIS-II first-year sequences are the scientific target.
# SiteSelection Complete + MaxCohortAge only; no Random/Spreading/global RNG.
harvest_plan_e3 <- function(simulation, topology, prescriptions,
                            history=NULL, year=simulation$year,
                            landis_seed, landis_rng_offset=0,
                            allocation_mode=c("corrected","landis_compatible")) {
  .h04_sim(simulation)
  allocation_mode <- match.arg(allocation_mode)
  valid_num <- function(x,lo,hi) is.numeric(x)&&length(x)==1L&&
    !is.na(x)&&is.finite(x)&&x>=lo&&x<=hi&&x==trunc(x)
  if(missing(landis_seed)||!valid_num(landis_seed,0,4294967295))
    stop("landis_seed must be a uint32",call.=FALSE)
  if(!valid_num(landis_rng_offset,0,1e9))stop("Invalid landis_rng_offset",call.=FALSE)
  if(!valid_num(year,0,1e9))stop("Invalid harvest year",call.=FALSE)
  if(is.null(topology$d2_neighbors)||!inherits(topology,"landisR_h04_topology")||
     !identical(topology$ids,as.integer(simulation$ctx$ids)))
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
     anyNA(p$ranking)||any(p$ranking!="MaxCohortAge")||
     anyNA(p$target_type)||any(p$target_type!="area")||
     anyNA(p$site_selection)||any(p$site_selection!="Complete"))
    stop("E3 pilot supports MaxCohortAge, Complete, area quotas only",call.=FALSE)
  if(!is.numeric(p$fraction)||anyNA(p$fraction)||any(!is.finite(p$fraction))||
     any(p$fraction<0|p$fraction>1)||
     any(p$maximum_age<p$minimum_age)||any(p$end_time<p$begin_time)||
     anyDuplicated(p[c("management_area","prescription")]))
    stop("Invalid prescription target/age/period/duplicate",call.=FALSE)
  if(any(!p$management_area%in%topology$management_areas)||
     any(!topology$management_areas%in%p$management_area))
    stop("Missing or unknown management area",call.=FALSE)
  if(allocation_mode=="landis_compatible" && any(p$fraction==0))
    stop("LANDIS-II refuses zero percent implementation",call.=FALSE)
  optional_integer <- function(name,default) {
    val<-if(name%in%names(p))p[[name]] else rep(default,nrow(p))
    .h04_integer(val,name,0L)
  }
  threshold<-optional_integer("stand_adjacency",0L)
  hist_gap<-optional_integer("minimum_time_since_last_harvest",0L)
  adj_type<-if("adjacency_type"%in%names(p))as.character(p$adjacency_type) else
    rep("none",nrow(p))
  key<-tolower(gsub("[[:space:]_-]+","",trimws(adj_type)))
  key[!is.na(key)&key==""&threshold==0L]<-"none"
  codes<-c(none=0L,standage=1L,timesincelastharvested=2L)
  if(anyNA(key)||any(!key%in%names(codes))||any(threshold>0L & key=="none"))
    stop("Active stand adjacency requires StandAge or TimeSinceLastHarvested",call.=FALSE)
  kind<-as.integer(unname(codes[key]))
  species_names<-as.character(simulation$ctx$config$species$species)
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
  if(!is.null(history)){
    if(!is.data.frame(history)||!all(c("stand","last_harvest_year")%in%names(history)))
      stop("history needs stand,last_harvest_year",call.=FALSE)
    hs<-.h04_integer(history$stand,"history stand",0L)
    if(anyDuplicated(hs)||length(hs)!=length(last)||
       !setequal(hs,topology$stands$stand))stop("history/topology mismatch",call.=FALSE)
    values<-as.numeric(history$last_harvest_year[match(topology$stands$stand,hs)])
    if(any(!is.na(values)&(!is.finite(values)|values<0|values>year|values!=trunc(values))))
      stop("Invalid last_harvest_year",call.=FALSE)
    last<-as.integer(values)
  }
  cfg<-list(ma=as.integer(match(p$management_area,topology$management_areas)-1L),
    fraction=as.numeric(p$fraction),min_age=as.integer(p$minimum_age),
    max_age=as.integer(p$maximum_age),begin=as.integer(p$begin_time),
    end=as.integer(p$end_time),adjacent_limit=as.integer(threshold),
    adjacent_type=kind,min_harvest_gap=as.integer(hist_gap),
    species_count=as.integer(length(species_names)),
    forest_species=as.integer(forest_index),forest_min_age=as.integer(forest_lo),
    forest_max_age=as.integer(forest_hi),forest_percent=as.numeric(forest_pct))
  # One batched native pass over all existing shared living CSR states.
  nat<-.Call("h04e3_plan",simulation$states,simulation$cell_state,
    topology$positions,topology$stand_offsets,topology$ma_offsets,
    topology$d2_neighbors$offsets,topology$d2_neighbors$indices,
    topology$d2_neighbors$first_position,cfg,last,
    as.numeric(year),as.numeric(landis_seed),as.numeric(landis_rng_offset),
    PACKAGE="landisR.harvest")
  s<-as.integer(nat$selected_stands);j<-as.integer(nat$prescription_indices)
  loc<-lapply(s,function(k){i<-seq.int(topology$stand_offsets[k]+1L,topology$stand_offsets[k+1L]);
    as.integer(topology$ids[topology$positions[i]+1L])})
  events<-data.frame(EventID=seq_along(s),Time=rep.int(as.integer(year),length(s)),
    ManagementArea=topology$stands$management_area[s],Prescription=p$prescription[j],
    Stand=topology$stands$stand[s],StandAge=as.integer(nat$event_age),
    StandRank=as.numeric(nat$event_age),NumberOfSites=as.integer(nat$event_sites),
    rng_draw_index=as.numeric(nat$draw_indices),
    rng_uniform=as.numeric(nat$draw_uniform),
    selection_probability=as.numeric(nat$selection_probability),
    stringsAsFactors=FALSE)
  events$prescription_row<-j;events$cells<-I(loc)
  history_out<-data.frame(management_area=topology$stands$management_area,
    stand=topology$stands$stand,last_harvest_year=as.integer(nat$last_harvest_year))
  list(type="H04-E3 pilot: MaxCohortAge Complete, multiple prescriptions, per-rule adjacency/history/Required forest type",
    year=as.integer(year),events=events,cells=as.integer(unlist(loc,use.names=FALSE)),
    history=history_out,
    prescriptions=within(p,{allocated_sites<-as.integer(nat$allocated_sites);
      eligible_stands<-as.integer(nat$eligible_stands)}),
    stand_ages=as.integer(nat$stand_age),
    landis_seed=as.numeric(landis_seed),
    landis_rng_offset_start=as.numeric(landis_rng_offset),
    landis_rng_offset_end=as.numeric(nat$rng_offset_end),
    draws_consumed=as.numeric(nat$draws_consumed))
}

harvest_step_e3 <- function(simulation,topology,prescriptions,rules_by_prescription,
                            history=NULL,year=simulation$year,landis_seed,
                            landis_rng_offset=0,
                            allocation_mode=c("corrected","landis_compatible")) {
  plan<-harvest_plan_e3(simulation,topology,prescriptions,history,year,
    landis_seed,landis_rng_offset,allocation_mode)
  chosen<-unique(as.character(plan$events$Prescription))
  if(!is.list(rules_by_prescription)||is.null(names(rules_by_prescription))||
     any(!chosen%in%names(rules_by_prescription)))
    stop("Need one named cohort-removal table per selected prescription",call.=FALSE)
  current<-simulation;effects<-list();sites<-list()
  for(pres in chosen){
    k<-which(plan$events$Prescription==pres)
    target<-as.integer(unlist(plan$events$cells[k],use.names=FALSE))
    res<-harvest_apply_cells(current,target,rules=rules_by_prescription[[pres]],clearcut=FALSE)
    current<-res$simulation;effects[[pres]]<-res$effects
    z<-res$sites;z$prescription<-pres;sites[[pres]]<-z
  }
  list(simulation=current,plan=plan,history=plan$history,
    effects=if(length(effects))do.call(rbind,effects)else data.frame(),
    sites=if(length(sites))do.call(rbind,sites)else data.frame(),
    landis_rng_offset_end=plan$landis_rng_offset_end)
}
