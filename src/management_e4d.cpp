#include <Rcpp.h>
#include <landisR/progress.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "harvest_math.hpp"
#include "management_e4d_math.hpp"
using namespace Rcpp;

// A single fused pass over living CSR cohort states calculates stand age and
// the stand-level presence of the specific Required ForestType rules in E3.
// The forest stands/cells are never duplicated. This is scoped to one Required
// species/age/percentage rule per prescription, no general ForestType DSL yet.
RcppExport SEXP h04e4d_plan(SEXP states_, SEXP cell_state_, SEXP pos_, SEXP so_,
 SEXP mo_, SEXP no_, SEXP ni_, SEXP first_, SEXP config_, SEXP history_,
 SEXP year_, SEXP seed_, SEXP offset_) {
 BEGIN_RCPP
 List states(states_),cfg(config_);
 using Clock = std::chrono::steady_clock;
 using Seconds = std::chrono::duration<double>;
 const bool timing = cfg.containsElementNamed("timing") && as<bool>(cfg["timing"]);
 const auto start = timing ? Clock::now() : Clock::time_point{};
 IntegerVector cs(cell_state_),pos(pos_),so(so_),mo(mo_),no(no_),ni(ni_),first(first_),history(history_);
 IntegerVector ma=as<IntegerVector>(cfg["ma"]),lo=as<IntegerVector>(cfg["min_age"]),
   hi=as<IntegerVector>(cfg["max_age"]),beg=as<IntegerVector>(cfg["begin"]),
   end=as<IntegerVector>(cfg["end"]),adj=as<IntegerVector>(cfg["adjacent_limit"]),
   kind=as<IntegerVector>(cfg["adjacent_type"]),gap=as<IntegerVector>(cfg["min_harvest_gap"]),
   ft_species=as<IntegerVector>(cfg["forest_species"]),
   ft_lo=as<IntegerVector>(cfg["forest_min_age"]),
   ft_hi=as<IntegerVector>(cfg["forest_max_age"]),
   ranking_code=as<IntegerVector>(cfg["ranking_code"]);
 NumericVector fraction=as<NumericVector>(cfg["fraction"]),ft_pct=as<NumericVector>(cfg["forest_percent"]),
   yearv(year_),seedv(seed_),offsetv(offset_);
 const int ns=so.size()-1,nm=mo.size()-1,np=ma.size();
 if(ns<1||nm<1||np<1||so[0]!=0||mo[0]!=0||mo[nm]!=ns||
    no.size()!=ns+1||no[0]!=0||no[ns]!=ni.size()||so[ns]!=pos.size()||
    first.size()!=ns||history.size()!=ns||cs.size()<1||
    lo.size()!=np||hi.size()!=np||beg.size()!=np||end.size()!=np||adj.size()!=np||ranking_code.size()!=np||
    kind.size()!=np||gap.size()!=np||ft_species.size()!=np||ft_lo.size()!=np||ft_hi.size()!=np||
    fraction.size()!=np||ft_pct.size()!=np||yearv.size()!=1||seedv.size()!=1||offsetv.size()!=1||
    !std::isfinite(yearv[0])||yearv[0]!=std::floor(yearv[0])||yearv[0]<0||yearv[0]>1e9||
    !std::isfinite(seedv[0])||seedv[0]!=std::floor(seedv[0])||seedv[0]<0||seedv[0]>4294967295.0||
    !std::isfinite(offsetv[0])||offsetv[0]!=std::floor(offsetv[0])||offsetv[0]<0||offsetv[0]>9007199254740991.0)
     stop("Invalid H04-E3 native inputs");
 IntegerVector coff=as<IntegerVector>(states["offset"]),sp=as<IntegerVector>(states["species_id"]);
 NumericVector cohort_age=as<NumericVector>(states["age"]),biomass=as<NumericVector>(states["biomass"]);
 const int nstates=coff.size()-1;
 if(nstates<1||coff[0]!=0||coff[nstates]!=sp.size()||
    sp.size()!=cohort_age.size()||sp.size()!=biomass.size())stop("Invalid living CSR");
 for (int st=0; st<nstates; ++st) {
   if (coff[st]<0 || coff[st+1]<coff[st] || coff[st+1]>sp.size())
     stop("Invalid living CSR offsets");
 }
 // Validate species against the fixed catalog, NOT against species still
 // represented by living cohorts. An extinction or a full clearcut may remove
 // the highest-ID species from the CSR without invalidating ForestType rules.
 const int species_count=as<int>(cfg["species_count"]);
 if(species_count<1)stop("Invalid E3 species catalog size");
 for(int q=0;q<sp.size();++q) {
   if(sp[q]<0||sp[q]>=species_count||!std::isfinite(cohort_age[q])||cohort_age[q]<0||
      !std::isfinite(biomass[q])||biomass[q]<0)stop("Invalid living cohort");
 }
 const bool has_spread=cfg.containsElementNamed("spreading");
 LogicalVector spreading(np,false),partial(np,false),patch_overlap(np,false);
 NumericVector min_block(np,0.0),max_block(np,0.0),patch_fraction(np,0.0),patch_size(np,0.0);
 const double cell_area_ha=cfg.containsElementNamed("cell_area_ha") ?
   as<double>(cfg["cell_area_ha"]) : 1.0;
 if(cfg.containsElementNamed("patch_fraction")){
   patch_fraction=as<NumericVector>(cfg["patch_fraction"]);
   patch_size=as<NumericVector>(cfg["patch_size_ha"]);
   patch_overlap=as<LogicalVector>(cfg["patch_overlap"]);
   if(patch_fraction.size()!=np||patch_size.size()!=np||patch_overlap.size()!=np)
     stop("Invalid PatchCutting parameters");
 }
 if(has_spread){
   spreading=as<LogicalVector>(cfg["spreading"]);
   min_block=as<NumericVector>(cfg["min_block_ha"]);
   max_block=as<NumericVector>(cfg["max_block_ha"]);
   if(cfg.containsElementNamed("partial")) partial=as<LogicalVector>(cfg["partial"]);
   if(spreading.size()!=np||min_block.size()!=np||max_block.size()!=np||partial.size()!=np)
     stop("Invalid stand-spread parameters");
 }
 IntegerVector set_aside(np,0);
 std::vector<double> reserve_until;
 if(cfg.containsElementNamed("set_aside")){
   set_aside=as<IntegerVector>(cfg["set_aside"]);
   reserve_until=as<std::vector<double>>(cfg["set_aside_until"]);
   if(set_aside.size()!=np||reserve_until.size()!=static_cast<std::size_t>(ns))
     stop("Invalid adjacency reserve");
   for(double &value:reserve_until) if(NumericVector::is_na(value)) value=-1;
 }
 h04e4d::RepeatPlan repeat;
 if (cfg.containsElementNamed("last_event_id")) {
   repeat.last_event_id = as<int>(cfg["last_event_id"]);
   if (repeat.last_event_id < 0) stop("Invalid previous harvest event ID");
 }
 IntegerVector repeat_kind(np), repeat_interval(np), repeat_count(np), repeat_rule_index(np,-1);
 std::vector<int> due_rule;
 if (cfg.containsElementNamed("repeat_kind")) {
   repeat.enabled = true;
   repeat.compatible = as<bool>(cfg["repeat_compatible"]);
   repeat.end_year = as<int>(cfg["repeat_end_year"]);
   repeat.current_rank = as<double>(cfg["current_rank"]);
   repeat_kind = as<IntegerVector>(cfg["repeat_kind"]);
   repeat_interval = as<IntegerVector>(cfg["repeat_interval"]);
   repeat_count = as<IntegerVector>(cfg["repeat_count"]);
   repeat_rule_index = as<IntegerVector>(cfg["repeat_rule_index"]);
   if (repeat_kind.size()!=np || repeat_interval.size()!=np ||
       repeat_count.size()!=np || repeat_rule_index.size()!=np ||
       repeat.last_event_id<0 || !std::isfinite(repeat.current_rank))
     stop("Invalid repeat settings");
   repeat.rule_index = repeat_rule_index.begin();
   due_rule.assign(ns,-1);
   std::vector<unsigned char> scheduled(ns,0);
   List jobs=cfg["repeat_jobs"];
   for (int i=0;i<jobs.size();++i) {
     List row=jobs[i];
     const int j=as<int>(row["prescription"]);
     const int next_year=as<int>(row["year"]);
     const int iteration=as<int>(row["iteration"]);
     const int id=as<int>(row["event_id"]);
     const double rank=as<double>(row["rank"]);
     if(j<0 || j>=np || ma[j]<0 || ma[j]>=nm || repeat_kind[j]<1 ||
        next_year<0 || iteration<1 || iteration>repeat_count[j] ||
        id<1 || id>repeat.last_event_id || !std::isfinite(rank))
       stop("Invalid pending harvest");
     h04e4d::RepeatJob job{j,next_year,iteration,id,rank,
       as<std::vector<int>>(row["stands"]),{}};
     List cells=row["cells"];
     if (job.stands.empty() || static_cast<std::size_t>(cells.size())!=job.stands.size())
       stop("Invalid repeat harvest footprint");
     for (int k=0;k<cells.size();++k) {
       const int stand=job.stands[k];
       if(stand<mo[ma[j]] || stand>=mo[ma[j]+1] || scheduled[stand])
         stop("Invalid or duplicate pending stand");
       scheduled[stand]=1;
       if (so[stand]<0 || so[stand]>=so[stand+1] || so[stand+1]>pos.size())
         stop("Invalid scheduled stand offsets");
       auto indices=as<std::vector<int>>(cells[k]);
       for (int cell:indices)
         if (cell<0 || cell>=cs.size() ||
             !std::binary_search(pos.begin()+so[stand],pos.begin()+so[stand+1],cell))
           stop("Repeat cell does not belong to the scheduled stand");
       if (next_year<=yearv[0]) {
         if (repeat_rule_index[j]<0) stop("Missing rules for repeat harvest");
         due_rule[stand]=repeat_rule_index[j];
       }
       job.cells.push_back(std::move(indices));
     }
     repeat.jobs.push_back(std::move(job));
   }
 }
 h04e4d::DamageHistory damage;
 std::vector<int> listed_sites;
 if (cfg.containsElementNamed("delisted_cells")) {
   IntegerVector cells = cfg["delisted_cells"];
   damage.delisted.reserve(cells.size());
   for (int cell : cells) {
     if (cell < 0 || cell >= cs.size()) stop("Invalid delisted harvest cell");
     damage.delisted.insert(cell);
   }
   listed_sites.resize(ns, 0);
 }
 IntegerVector min_damage(np, 0);
 if (cfg.containsElementNamed("min_damage")) {
   damage.enabled = true;
   damage.compatible = as<bool>(cfg["damage_compatible"]);
   min_damage = as<IntegerVector>(cfg["min_damage"]);
   if (min_damage.size() != np) stop("Invalid minimum time since damage");
   const int initial = damage.compatible ? -100 : std::numeric_limits<int>::min();
   damage.harvest.assign(cs.size(), initial);
   damage.external.assign(cs.size(), initial);
   const auto load_damage = [&](std::vector<int>& target, const char* position_key, const char* year_key) {
     IntegerVector positions = cfg[position_key];
     IntegerVector years = cfg[year_key];
     if (positions.size() != years.size()) stop("Invalid damage history vectors");
     for (int i = 0; i < positions.size(); ++i) {
       const int cell = positions[i], year = years[i];
       if (cell < 0 || cell >= cs.size() || year < 0 || year > yearv[0])
         stop("Invalid damage history record");
       target[cell] = std::max(target[cell], year);
     }
   };
   load_damage(damage.harvest, "harvest_damage_positions", "harvest_damage_year");
   load_damage(damage.external, "external_damage_positions", "external_damage_year");
 }
 std::vector<h04e4d::Prescription> params;params.reserve(np);
 struct Filter {int species,lo,hi;};
 std::vector<Filter> filters;
 std::vector<int> filter_of(np,-1);
 for(int j=0;j<np;++j){
   params.push_back({ma[j],lo[j],hi[j],beg[j],end[j],adj[j],kind[j],gap[j],ranking_code[j],fraction[j],
     spreading[j]==TRUE,min_block[j],max_block[j],partial[j]==TRUE,patch_fraction[j],patch_size[j],patch_overlap[j]==TRUE,set_aside[j],
     repeat_kind[j],repeat_interval[j],repeat_count[j],min_damage[j]});
   if(ft_species[j]>=0){
     if(ft_species[j]>=species_count||ft_lo[j]<1||ft_hi[j]<ft_lo[j]||
         !std::isfinite(ft_pct[j])||ft_pct[j]<=0||ft_pct[j]>100)
       stop("Invalid E3 Required ForestType rule");
     int k=-1;
     for(std::size_t f=0;f<filters.size();++f)
       if(filters[f].species==ft_species[j]&&filters[f].lo==ft_lo[j]&&filters[f].hi==ft_hi[j]){k=f;break;}
     if(k<0){k=filters.size();filters.push_back({ft_species[j],ft_lo[j],ft_hi[j]});}
     filter_of[j]=k;
   }else if(ft_species[j]!=-1) stop("Invalid E3 forest species sentinel");
 }
 const auto prescription_offsets = h04e4d::prescription_offsets(as<std::vector<int>>(mo), params);
 std::vector<std::vector<int>> prescriptions_by_area(nm);
 for (int j = 0; j < np; ++j) prescriptions_by_area[ma[j]].push_back(j);

 IntegerVector economic_index(np, -1);
 NumericMatrix economic_weight;
 IntegerMatrix economic_age;
 int economic_count = 0;
 if (cfg.containsElementNamed("economic_index")) {
   economic_index = as<IntegerVector>(cfg["economic_index"]);
   economic_weight = as<NumericMatrix>(cfg["economic_weight"]);
   economic_age = as<IntegerMatrix>(cfg["economic_age"]);
   economic_count = economic_weight.nrow();
   if (economic_index.size() != np || economic_weight.ncol() != species_count ||
       economic_age.nrow() != economic_count || economic_age.ncol() != species_count)
     stop("Invalid Economic ranking tables");
   for (int j = 0; j < np; ++j)
     if (ranking_code[j] == 2 && (economic_index[j] < 0 || economic_index[j] >= economic_count))
       stop("Missing Economic ranking table");
 }
 std::vector<double> economic_cache(static_cast<std::size_t>(economic_count) * nstates,
                                    std::numeric_limits<double>::quiet_NaN());
 const auto economic_rank = [&](int table, int state) {
   double& rank = economic_cache[static_cast<std::size_t>(table) * nstates + state];
   if (std::isnan(rank)) {
     rank = 0;
     for (int q = coff[state]; q < coff[state + 1]; ++q) {
       const int minimum = economic_age(table, sp[q]);
       if (minimum > 0 && cohort_age[q] >= minimum)
         rank += economic_weight(table, sp[q]) * cohort_age[q];
     }
   }
   return rank;
 };
 // Cache whether a compressed cell state yields a positive removal under
 // each distinct prescription. Disjoint stands keep these predictions valid
 // until application, with no trial edits to the landscape.
 const bool credit_effective = cfg.containsElementNamed("removal_rules");
 using SpeciesRules = std::vector<std::vector<harvest_h03::Rule>>;
 std::vector<SpeciesRules> removal_rules;
 std::vector<std::vector<int>> active_prescriptions(nm);
 IntegerVector rule_index(np, -1);
 if (credit_effective) {
   List tables = cfg["removal_rules"];
   rule_index = as<IntegerVector>(cfg["rule_index"]);
   if (rule_index.size()!=np) stop("Invalid quota rule index");
   for (int j=0; j<tables.size(); ++j) {
     List table = tables[j];
     IntegerVector species = table["species"];
     NumericVector lower = table["lo"], upper = table["hi"], proportion = table["fraction"];
     if (species.size()!=lower.size() || lower.size()!=upper.size() ||
         upper.size()!=proportion.size()) stop("Invalid quota removal rules");
     SpeciesRules rules(species_count);
     for (int k=0; k<species.size(); ++k) {
       if (species[k]<0 || species[k]>=species_count) stop("Invalid quota rule species");
       rules[species[k]].push_back({lower[k],upper[k],proportion[k]});
     }
     removal_rules.push_back(std::move(rules));
   }
   for (int j=0; j<np; ++j) {
     const bool active = fraction[j]>0 && yearv[0]>=beg[j] && yearv[0]<=end[j];
     if ((active && rule_index[j]<0) || rule_index[j]>=tables.size())
       stop("Missing quota removal rules for active prescription");
     if (ma[j]<0 || ma[j]>=nm) stop("Invalid quota management area");
     if (active) active_prescriptions[ma[j]].push_back(j);
   }
 }
 for (int i:due_rule)
   if(i>=static_cast<int>(removal_rules.size())) stop("Missing rules for repeat harvest");
 std::vector<signed char> cache_removal(static_cast<std::size_t>(nstates)*removal_rules.size(),-1);
 std::vector<int> effective_sites(credit_effective ? prescription_offsets.back() : 0,0);
 const bool use_biomass = cfg.containsElementNamed("minimum_biomass_mg_ha");
 NumericVector minimum_biomass;
 std::vector<double> cache_biomass;
 if (use_biomass) {
   minimum_biomass = as<NumericVector>(cfg["minimum_biomass_mg_ha"]);
   if (minimum_biomass.size() != np) stop("Invalid minimum stand biomass vector");
   for (double value : minimum_biomass)
     if (!std::isfinite(value) || value < 0) stop("Invalid minimum stand biomass");
   cache_biomass.resize(nstates, 0.0);
 }
 std::vector<int> cache_age(nstates,-1);
 std::vector<unsigned char> cache_flag(static_cast<std::size_t>(nstates)*filters.size(),0);
 std::vector<int> age(ns,0);
 std::vector<double> ranks(prescription_offsets.back(),0);
 std::vector<unsigned char> qualifies(prescription_offsets.back(),1);
 std::vector<int> counts(filters.size(),0);
 int management_area=0;
 const auto analyse_start = timing ? Clock::now() : Clock::time_point{};
 landis_progress::Reporter progress;
 progress.start(0, 0.55, "Evaluating & ranking stands", ns);
 for(int s=0;s<ns;++s){
   progress.advance(s);
   if(s%4096==0)checkUserInterrupt();
   if(so[s]>=so[s+1]||so[s]<0||so[s+1]>pos.size())stop("Invalid E3 stand offsets");
   std::fill(counts.begin(),counts.end(),0);
   while (management_area+1<nm && s>=mo[management_area+1]) ++management_area;
   double sum_age=0, sum_biomass=0;
   for(int k=so[s];k<so[s+1];++k){
     int position=pos[k];if(position<0||position>=cs.size())stop("Invalid E3 active cell position");
     const bool listed = damage.listed(position);
     if (!listed && !use_biomass) continue;
     int st=cs[position]-1;if(st<0||st>=nstates)stop("Invalid E3 shared state reference");
     if(cache_age[st]<0){
       int ma=0;
       for(int q=coff[st];q<coff[st+1];++q) {
         ma=std::max(ma,static_cast<int>(cohort_age[q]));
         if (use_biomass) cache_biomass[st] += biomass[q];
         for(std::size_t f=0;f<filters.size();++f){
           const auto &x=filters[f];
           if(sp[q]==x.species&&cohort_age[q]>=x.lo&&cohort_age[q]<=x.hi)
             cache_flag[static_cast<std::size_t>(st)*filters.size()+f]=1;
         }
       }
       cache_age[st]=ma;
     }
     if (use_biomass) sum_biomass += cache_biomass[st];
     if (!listed) continue;
     if (!listed_sites.empty()) ++listed_sites[s];
     sum_age+=cache_age[st];
     for(std::size_t f=0;f<filters.size();++f)
       counts[f]+=cache_flag[static_cast<std::size_t>(st)*filters.size()+f];
     const auto evaluate_removal = [&](int index){
       auto &has_removal=cache_removal[static_cast<std::size_t>(index)*nstates+st];
       if (has_removal<0) {
         has_removal=0;
         const auto &rules=removal_rules[index];
         for (int q=coff[st]; q<coff[st+1]; ++q) {
           const double f=harvest_h03::removal_fraction(cohort_age[q],rules[sp[q]]);
           if (harvest_h03::removed(biomass[q],f)>0) {
             has_removal=1;
             break;
           }
         }
       }
       return has_removal;
     };
     for (int j : active_prescriptions[management_area])
       effective_sites[prescription_offsets[j] + (s - mo[management_area])] += evaluate_removal(rule_index[j]);
     for (int j : prescriptions_by_area[management_area])
       if (economic_index[j] >= 0 && yearv[0] >= beg[j] && yearv[0] <= end[j] && fraction[j] > 0)
         ranks[prescription_offsets[j] + (s - mo[management_area])] += economic_rank(economic_index[j], st);
     if (!due_rule.empty() && due_rule[s]>=0) evaluate_removal(due_rule[s]);
   }
   const int total = listed_sites.empty() ? so[s+1]-so[s] : listed_sites[s];
   if (total > 0) age[s]=static_cast<int>(sum_age/total);
   // The denominator includes all active managed cells, including empty ones.
   const double biomass_mg_ha = sum_biomass / (100.0 * (so[s+1]-so[s]));
   for (int j : prescriptions_by_area[management_area]) {
     const auto ix = prescription_offsets[j] + (s - mo[management_area]);
     if (total == 0) qualifies[ix] = 0;
     if(filter_of[j]>=0 && total>0){
       int matches=counts[filter_of[j]];
       qualifies[ix]=static_cast<unsigned char>(100.0*matches+1e-12>=ft_pct[j]*total);
     }
     if (use_biomass && biomass_mg_ha < minimum_biomass[j]) qualifies[ix] = 0;
     if (qualifies[ix])
       ranks[ix] = economic_index[j] >= 0 ? ranks[ix] / total : age[s];
     else ranks[ix] = 0;
   }
 }
 const auto analyse_end = timing ? Clock::now() : Clock::time_point{};
 std::vector<int> normalized_history=as<std::vector<int>>(history);
 for(int &v:normalized_history)if(v==NA_INTEGER)v=-1;
 h04e4d::PartialTopology cell_topology;
 IntegerVector raster_positions, neighbors, stand_of_cell;
 if(cfg.containsElementNamed("cell_topology")){
   List cells=cfg["cell_topology"];
   raster_positions=as<IntegerVector>(cells["positions"]);
   neighbors=as<IntegerVector>(cells["neighbors"]);
   stand_of_cell=as<IntegerVector>(cells["stand"]);
   if(raster_positions.size()!=cs.size()||stand_of_cell.size()!=cs.size()||
      neighbors.size()!=cs.size()*h04e4d::cell_directions)
     stop("Invalid partial-spread cell topology");
   cell_topology={raster_positions.begin(),neighbors.begin(),stand_of_cell.begin(),
                  cs.begin(),rule_index.begin(),
                  credit_effective ? cache_removal.data() : nullptr,
                  static_cast<std::size_t>(nstates),static_cast<std::size_t>(cs.size()),
                  as<bool>(cfg["landis_compatible"])};
 } else if (repeat.enabled || damage.enabled || !damage.delisted.empty()) {
   cell_topology = {pos.begin(),nullptr,nullptr,cs.begin(),rule_index.begin(),
                    credit_effective ? cache_removal.data() : nullptr,
                    static_cast<std::size_t>(nstates),0,repeat.compatible};
 }
 const bool shared = cfg.containsElementNamed("landis_rng_state");
 const auto offset = landis_rng::Stream::offset(offsetv[0]);
 h04c3::Stream rng(static_cast<std::uint32_t>(seedv[0]), shared ? 0 : offset);
 if (shared) {
   rng.restore(as<std::string>(cfg["landis_rng_state"]));
   rng.draws = offset;
 }
 const auto select_start = timing ? Clock::now() : Clock::time_point{};
 h04e4d::Result out=h04e4d::select(
   as<std::vector<int>>(so),as<std::vector<int>>(mo),
   as<std::vector<int>>(no),as<std::vector<int>>(ni),
   as<std::vector<int>>(first),age,std::move(ranks),qualifies,params,normalized_history,
   static_cast<int>(yearv[0]),rng,effective_sites,cell_area_ha,cell_topology,reserve_until,repeat,std::move(damage),listed_sites,
   cfg.containsElementNamed("landis_economic_order") && as<bool>(cfg["landis_economic_order"]),
   progress ? &progress : nullptr);
 progress.phase(0.97, "Collecting selected stands");
 const auto select_end = timing ? Clock::now() : Clock::time_point{};
 IntegerVector selected(out.events.size()),pres(out.events.size()),maevents(out.events.size()),
   event_age(out.events.size()),sites(out.events.size()),event_ids(out.events.size());
 NumericVector draw(out.events.size()),uniform(out.events.size()),probability(out.events.size()),event_rank(out.events.size());
 List partial_cells(out.events.size());
 IntegerVector repeat_numbers(out.events.size());
 for(std::size_t k=0;k<out.events.size();++k){
   const auto &e=out.events[k];selected[k]=e.stand+1;pres[k]=e.pres+1;
   event_ids[k]=e.id;maevents[k]=e.ma+1;event_age[k]=e.age;sites[k]=e.number_of_sites;event_rank[k]=e.rank;
   repeat_numbers[k]=e.repeat_number;
   if(e.cell_end>e.cell_begin || e.repeat_number>0)
     partial_cells[k]=IntegerVector(out.selected_cells.begin()+e.cell_begin,
                                    out.selected_cells.begin()+e.cell_end);
   draw[k]=static_cast<double>(e.draw);uniform[k]=e.uniform;probability[k]=e.probability;
 }
 IntegerVector updated=wrap(out.last_harvest);
 for(int k=0;k<updated.size();++k)if(updated[k]<0)updated[k]=NA_INTEGER;
 NumericVector expiry=wrap(out.set_aside_until);
 for(double &value:expiry) if(value<0) value=NA_REAL;
 List result=List::create(_["selected_stands"]=selected,_["prescription_indices"]=pres,
   _["management_indices"]=maevents,_["event_ids"]=event_ids,_["event_age"]=event_age,
   _["event_sites"]=sites,_["partial_cells"]=partial_cells,_["event_ranks"]=event_rank,_["draw_indices"]=draw,
   _["draw_uniform"]=uniform,_["selection_probability"]=probability,
   _["last_harvest_year"]=updated,_["set_aside_until"]=expiry,_["allocated_sites"]=wrap(out.allocated),
   _["attempted_sites"]=wrap(out.attempted),
   _["eligible_stands"]=wrap(out.eligible),_["stand_age"]=wrap(age),
   _["last_event_id"]=out.last_event_id,
   _["rng_offset_end"]=static_cast<double>(out.offset_end),
   _["draws_consumed"]=static_cast<double>(out.offset_end-offset));
 if (shared) result["landis_rng_state"] = rng.encode();
 if (repeat.enabled) {
   List jobs(out.repeat_jobs.size());
   for (std::size_t i=0;i<out.repeat_jobs.size();++i) {
     const auto& job=out.repeat_jobs[i];
     jobs[i]=List::create(_["prescription"]=job.prescription,_["year"]=job.year,
       _["iteration"]=job.iteration,_["event_id"]=job.event_id,_["rank"]=job.rank,
       _["stands"]=wrap(job.stands),_["cells"]=wrap(job.cells));
   }
   result["repeat_numbers"]=repeat_numbers;
   result["repeat_jobs"]=jobs;
   result["current_rank"]=out.current_rank;
 }
 if (!out.damage_year.empty()) {
   std::vector<int> visited, years;
   for (std::size_t cell = 0; cell < out.damage_year.size(); ++cell) {
     if (out.damage_year[cell] < 0) continue;
     visited.push_back(static_cast<int>(cell));
     years.push_back(out.damage_year[cell]);
   }
   result["harvest_damage_positions"] = visited;
   result["harvest_damage_year"] = years;
 }
 if (timing) {
   result["timing"] = NumericVector::create(
     _["validation_et_caches"] = Seconds(analyse_start - start).count(),
     _["analyse_peuplements"] = Seconds(analyse_end - analyse_start).count(),
     _["preparation_selection"] = Seconds(select_start - analyse_end).count(),
     _["classement_et_selection"] = Seconds(select_end - select_start).count(),
     _["construction_sorties"] = Seconds(Clock::now() - select_end).count()
   );
 }
 return result;
 END_RCPP
}
