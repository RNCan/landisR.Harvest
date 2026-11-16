#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "management_e3_math.hpp"
using namespace Rcpp;

// A single fused pass over living CSR cohort states calculates stand age and
// the stand-level presence of the specific Required ForestType rules in E3.
// The forest stands/cells are never duplicated. This is scoped to one Required
// species/age/percentage rule per prescription, no general ForestType DSL yet.
RcppExport SEXP h04e3_plan(SEXP states_, SEXP cell_state_, SEXP pos_, SEXP so_,
 SEXP mo_, SEXP no_, SEXP ni_, SEXP first_, SEXP config_, SEXP history_,
 SEXP year_, SEXP seed_, SEXP offset_) {
 BEGIN_RCPP
 List states(states_),cfg(config_);
 IntegerVector cs(cell_state_),pos(pos_),so(so_),mo(mo_),no(no_),ni(ni_),first(first_),history(history_);
 IntegerVector ma=as<IntegerVector>(cfg["ma"]),lo=as<IntegerVector>(cfg["min_age"]),
   hi=as<IntegerVector>(cfg["max_age"]),beg=as<IntegerVector>(cfg["begin"]),
   end=as<IntegerVector>(cfg["end"]),adj=as<IntegerVector>(cfg["adjacent_limit"]),
   kind=as<IntegerVector>(cfg["adjacent_type"]),gap=as<IntegerVector>(cfg["min_harvest_gap"]),
   ft_species=as<IntegerVector>(cfg["forest_species"]),
   ft_lo=as<IntegerVector>(cfg["forest_min_age"]),
   ft_hi=as<IntegerVector>(cfg["forest_max_age"]);
 NumericVector fraction=as<NumericVector>(cfg["fraction"]),ft_pct=as<NumericVector>(cfg["forest_percent"]),
   yearv(year_),seedv(seed_),offsetv(offset_);
 const int ns=so.size()-1,nm=mo.size()-1,np=ma.size();
 if(ns<1||nm<1||np<1||so[0]!=0||mo[0]!=0||mo[nm]!=ns||
    no.size()!=ns+1||no[0]!=0||no[ns]!=ni.size()||so[ns]!=pos.size()||
    first.size()!=ns||history.size()!=ns||cs.size()<1||
    lo.size()!=np||hi.size()!=np||beg.size()!=np||end.size()!=np||adj.size()!=np||
    kind.size()!=np||gap.size()!=np||ft_species.size()!=np||ft_lo.size()!=np||ft_hi.size()!=np||
    fraction.size()!=np||ft_pct.size()!=np||yearv.size()!=1||seedv.size()!=1||offsetv.size()!=1||
    !std::isfinite(yearv[0])||yearv[0]!=std::floor(yearv[0])||yearv[0]<0||yearv[0]>1e9||
    !std::isfinite(seedv[0])||seedv[0]!=std::floor(seedv[0])||seedv[0]<0||seedv[0]>4294967295.0||
    !std::isfinite(offsetv[0])||offsetv[0]!=std::floor(offsetv[0])||offsetv[0]<0||offsetv[0]>1e9)
     stop("Invalid H04-E3 native inputs");
 IntegerVector coff=as<IntegerVector>(states["offset"]),sp=as<IntegerVector>(states["species_id"]);
 NumericVector cohort_age=as<NumericVector>(states["age"]),biomass=as<NumericVector>(states["biomass"]);
 const int nstates=coff.size()-1;
 if(nstates<1||coff[0]!=0||coff[nstates]!=sp.size()||
    sp.size()!=cohort_age.size()||sp.size()!=biomass.size())stop("Invalid living CSR");
 // Validate species against the fixed catalog, NOT against species still
 // represented by living cohorts. An extinction or a full clearcut may remove
 // the highest-ID species from the CSR without invalidating ForestType rules.
 const int species_count=as<int>(cfg["species_count"]);
 if(species_count<1)stop("Invalid E3 species catalog size");
 for(int q=0;q<sp.size();++q) {
   if(sp[q]<0||sp[q]>=species_count||!std::isfinite(cohort_age[q])||cohort_age[q]<0||
      !std::isfinite(biomass[q])||biomass[q]<0)stop("Invalid living cohort");
 }
 std::vector<h04e3::Prescription> params;params.reserve(np);
 struct Filter {int species,lo,hi;};
 std::vector<Filter> filters;
 std::vector<int> filter_of(np,-1);
 for(int j=0;j<np;++j){
   params.push_back({ma[j],lo[j],hi[j],beg[j],end[j],adj[j],kind[j],gap[j],fraction[j]});
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
 std::vector<int> cache_age(nstates,-1);
 std::vector<unsigned char> cache_flag(static_cast<std::size_t>(nstates)*filters.size(),0);
 std::vector<int> age(ns,0);
 std::vector<double> ranks(static_cast<std::size_t>(np)*ns,0);
 std::vector<unsigned char> qualifies(static_cast<std::size_t>(np)*ns,1);
 std::vector<int> counts(filters.size(),0);
 for(int s=0;s<ns;++s){
   if(s%4096==0)checkUserInterrupt();
   if(so[s]>=so[s+1]||so[s]<0||so[s+1]>pos.size())stop("Invalid E3 stand offsets");
   std::fill(counts.begin(),counts.end(),0);
   double sum_age=0;
   for(int k=so[s];k<so[s+1];++k){
     int position=pos[k];if(position<0||position>=cs.size())stop("Invalid E3 active cell position");
     int st=cs[position]-1;if(st<0||st>=nstates)stop("Invalid E3 shared state reference");
     if(cache_age[st]<0){
       int ma=0;
       for(int q=coff[st];q<coff[st+1];++q) {
         ma=std::max(ma,static_cast<int>(cohort_age[q]));
         for(std::size_t f=0;f<filters.size();++f){
           const auto &x=filters[f];
           if(sp[q]==x.species&&cohort_age[q]>=x.lo&&cohort_age[q]<=x.hi)
             cache_flag[static_cast<std::size_t>(st)*filters.size()+f]=1;
         }
       }
       cache_age[st]=ma;
     }
     sum_age+=cache_age[st];
     for(std::size_t f=0;f<filters.size();++f)
       counts[f]+=cache_flag[static_cast<std::size_t>(st)*filters.size()+f];
   }
   age[s]=static_cast<int>(sum_age/(so[s+1]-so[s]));
   for(int j=0;j<np;++j){
     int ix=j*ns+s;
     if(filter_of[j]>=0){
       int matches=counts[filter_of[j]],total=so[s+1]-so[s];
       qualifies[ix]=static_cast<unsigned char>(100.0*matches+1e-12>=ft_pct[j]*total);
     }
     ranks[ix]=qualifies[ix]?age[s]:0.0;
   }
 }
 std::vector<int> normalized_history=as<std::vector<int>>(history);
 for(int &v:normalized_history)if(v==NA_INTEGER)v=-1;
 h04e3::Result out=h04e3::select(
   as<std::vector<int>>(so),as<std::vector<int>>(mo),
   as<std::vector<int>>(no),as<std::vector<int>>(ni),
   as<std::vector<int>>(first),age,ranks,qualifies,params,normalized_history,
   static_cast<int>(yearv[0]),static_cast<std::uint32_t>(seedv[0]),
   static_cast<std::uint64_t>(offsetv[0]));
 IntegerVector selected(out.events.size()),pres(out.events.size()),maevents(out.events.size()),
   event_age(out.events.size()),sites(out.events.size());
 NumericVector draw(out.events.size()),uniform(out.events.size()),probability(out.events.size());
 for(std::size_t k=0;k<out.events.size();++k){
   const auto &e=out.events[k];selected[k]=e.stand+1;pres[k]=e.pres+1;
   maevents[k]=e.ma+1;event_age[k]=e.age;sites[k]=e.number_of_sites;
   draw[k]=static_cast<double>(e.draw);uniform[k]=e.uniform;probability[k]=e.probability;
 }
 IntegerVector updated=wrap(out.last_harvest);
 for(int k=0;k<updated.size();++k)if(updated[k]<0)updated[k]=NA_INTEGER;
 return List::create(_["selected_stands"]=selected,_["prescription_indices"]=pres,
   _["management_indices"]=maevents,_["event_age"]=event_age,
   _["event_sites"]=sites,_["draw_indices"]=draw,
   _["draw_uniform"]=uniform,_["selection_probability"]=probability,
   _["last_harvest_year"]=updated,_["allocated_sites"]=wrap(out.allocated),
   _["eligible_stands"]=wrap(out.eligible),_["stand_age"]=wrap(age),
   _["rng_offset_end"]=static_cast<double>(out.offset_end),
   _["draws_consumed"]=static_cast<double>(out.offset_end-static_cast<std::uint64_t>(offsetv[0])));
 END_RCPP
}
