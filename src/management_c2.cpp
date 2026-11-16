#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include "management_b2_math.hpp"
using namespace Rcpp;
namespace {
struct Signature {
  std::vector<unsigned char> species, rule;
};
inline int safe_int(double x) {
  if(!R_finite(x)||x<0||x>65535||x!=std::trunc(x))stop("Invalid cohort age");
  return static_cast<int>(x);
}
inline int banker(double x) {
  if(!std::isfinite(x)||x>2e9||x<0)stop("RegulateAges score out of range");
  return h04b2::round_even(x);
}
}
// H04-C2 deterministic qualified pilot: MaxCohortAge, Economic, RegulateAges,
// MinimumAge/MaximumAge and ForestTypeTable. Single Complete/ClearCut prescription
// per MA. Legacy FTT replicates OBSERVED small-stand presence behaviour only;
// it does not claim a complete reverse-engineering of the original DLL.
// IMPORTANT: rule indices, species and site positions are zero-based.
RcppExport SEXP h04c2_plan(SEXP states_,SEXP cs_,SEXP positions_,SEXP so_,SEXP mo_,
 SEXP weights_,SEXP merchant_,SEXP ranking_,SEXP fraction_,SEXP minage_,SEXP maxage_,
 SEXP target_type_,SEXP mode_,SEXP rule_ma_,SEXP rule_kind_,SEXP rule_lo_,SEXP rule_hi_,
 SEXP rule_pct_,SEXP rule_highest_,SEXP pairs_rule_,SEXP pairs_species_,SEXP dims_,SEXP ids_) {
 BEGIN_RCPP
 List states(states_);
 IntegerVector cs(cs_),pos(positions_),so(so_),mo(mo_),ranking(ranking_),minage(minage_),
  maxage(maxage_),target_type(target_type_),mode(mode_),rule_ma(rule_ma_),rule_kind(rule_kind_),
  rule_lo(rule_lo_),rule_hi(rule_hi_),rule_highest(rule_highest_),pairs_rule(pairs_rule_),pairs_sp(pairs_species_);
 NumericVector fractions(fraction_),rule_pct(rule_pct_);
 NumericMatrix weights(weights_);IntegerMatrix merchant(merchant_);
 IntegerVector dims(dims_),ids(ids_);
 if(dims.size()!=2||dims[0]<1||dims[1]<1)stop("Invalid landscape dimensions");
 IntegerVector off=as<IntegerVector>(states["offset"]),sp=as<IntegerVector>(states["species_id"]);
 NumericVector ages=as<NumericVector>(states["age"]),biomass=as<NumericVector>(states["biomass"]);
 const int ns=off.size()-1,nst=so.size()-1,nma=mo.size()-1,nsp=weights.ncol(),nr=rule_ma.size();
 if(ids.size()!=cs.size())stop("Active cell IDs and state IDs do not match");
 if(ns<1||nst<1||nma<1||nsp<1||off[0]!=0||off[ns]!=sp.size()||sp.size()!=ages.size()||ages.size()!=biomass.size()||so[0]!=0||so[nst]!=pos.size()||mo[0]!=0||mo[nma]!=nst||weights.nrow()!=nma||merchant.nrow()!=nma||merchant.ncol()!=nsp||ranking.size()!=nma||minage.size()!=nma||maxage.size()!=nma||target_type.size()!=nma||mode.size()!=nma||fractions.size()!=nma||rule_kind.size()!=nr||rule_lo.size()!=nr||rule_hi.size()!=nr||rule_pct.size()!=nr||rule_highest.size()!=nr||pairs_rule.size()!=pairs_sp.size()) stop("Invalid H04-C2 CSR/topology/target/rule dimensions");
 for(int s=0;s<ns;++s)if(off[s]<0||off[s+1]<off[s]||off[s+1]>sp.size())stop("Invalid CSR offset");
 for(int j=0;j<sp.size();++j)if(sp[j]<0||sp[j]>=nsp||safe_int(ages[j])==0||!R_finite(biomass[j])||biomass[j]<=0)stop("Invalid cohort state");
 for(int m=0;m<nma;++m)if(mo[m]<0||mo[m+1]<=mo[m]||mo[m+1]>nst||!R_finite(fractions[m])||fractions[m]<0||fractions[m]>1||minage[m]<0||maxage[m]<minage[m]||ranking[m]<0||ranking[m]>2||target_type[m]<0||target_type[m]>1||mode[m]<0||mode[m]>1) stop("Invalid MA settings");
 std::vector<std::vector<int>> rule_sp(nr),ma_rules(nma);
 for(int r=0;r<nr;++r) {
  if(rule_ma[r]<0||rule_ma[r]>=nma||rule_kind[r]<0||rule_kind[r]>2||rule_lo[r]<1||rule_hi[r]<rule_lo[r]||rule_hi[r]>65535||rule_highest[r]<0||rule_highest[r]>1||!R_finite(rule_pct[r])||rule_pct[r]<=0||rule_pct[r]>100)stop("Invalid forest type rule");
  ma_rules[rule_ma[r]].push_back(r);
 }
 for(int j=0;j<pairs_rule.size();++j){int r=pairs_rule[j],s=pairs_sp[j];if(r<0||r>=nr||s<0||s>=nsp)stop("Invalid forest type species pair");rule_sp[r].push_back(s);}
 for(int r=0;r<nr;++r)if(rule_sp[r].empty())stop("Rule without species");
 std::vector<int> max_state(ns,-1),size(nst),nonempty(nst),stand_age(nst),first_pos(nst);
 NumericVector stand_rank(nst); LogicalVector stand_qualified(nst);
 IntegerVector allocated(nma),total(nma),eligible(nma),target_stands(nma),selected_count(nma);
 std::vector<int> selected_stands,selected_positions;int evaluated=0;
 for(int m=0;m<nma;++m){
  std::vector<double> econ(ns,std::numeric_limits<double>::quiet_NaN());
  std::unordered_map<int,Signature> signatures;
  std::unordered_map<int,int> counts;
  // Reused, preallocated scratch buffers for all stands in this management area.
  std::vector<int> rule_count(nr,0),species_count(nsp,0);
  for(int s=mo[m];s<mo[m+1];++s){
   if(s%2048==0)checkUserInterrupt();
   if(so[s]<0||so[s+1]<=so[s]||so[s+1]>pos.size())stop("Invalid stand offsets");
   size[s]=so[s+1]-so[s];total[m]+=size[s];first_pos[s]=std::numeric_limits<int>::max();
   double sum_age=0,sum_econ=0;
   std::fill(rule_count.begin(),rule_count.end(),0);
   std::fill(species_count.begin(),species_count.end(),0);
   for(int j=so[s];j<so[s+1];++j){
    int p=pos[j];if(p<0||p>=cs.size())stop("Invalid cell position");
    // Original LANDIS-II stand order follows raster row-major discovery,
    // not numeric stand ID or core CSR column-major order.
    if(ids[p]<1 || static_cast<long long>(ids[p])>static_cast<long long>(dims[0])*dims[1])stop("Invalid canonical cell ID");
    int zero=ids[p]-1;
    int scan=(zero % dims[0])*dims[1]+(zero/dims[0]);
    first_pos[s]=std::min(first_pos[s],scan);
    int st=cs[p]-1;if(st<0||st>=ns)stop("Invalid CSR state reference");
    if(max_state[st]<0){int age=0;for(int q=off[st];q<off[st+1];++q)age=std::max(age,static_cast<int>(ages[q]));max_state[st]=age;++evaluated;}
    sum_age+=max_state[st];if(max_state[st]>0)++nonempty[s];
    if(ranking[m]==1){
     if(std::isnan(econ[st])){double v=0;for(int q=off[st];q<off[st+1];++q){int k=sp[q],age=merchant(m,k);if(weights(m,k)>0&&age>0&&ages[q]>=age)v+=weights(m,k)*ages[q]/age;}econ[st]=v;++evaluated;}
     sum_econ+=econ[st];
    }
    if(!ma_rules[m].empty()){
     auto it=signatures.find(st);
     if(it==signatures.end()){
      Signature sig;sig.species.assign(nsp,0);sig.rule.assign(nr,0);
      for(int q=off[st];q<off[st+1];++q)sig.species[sp[q]]=1;
      for(int r:ma_rules[m]) {
       for(int q=off[st];q<off[st+1]&&!sig.rule[r];++q){
        if(ages[q]<rule_lo[r]||ages[q]>rule_hi[r])continue;
        for(int species:rule_sp[r])if(species==sp[q]){sig.rule[r]=1;break;}
       }
      }
      it=signatures.emplace(st,std::move(sig)).first;++evaluated;
     }
     for(int k=0;k<nsp;++k)species_count[k]+=it->second.species[k];
     for(int r:ma_rules[m])rule_count[r]+=it->second.rule[r];
    }
   }
   stand_age[s]=static_cast<int>(sum_age/size[s]);
   if(ranking[m]==0)stand_rank[s]=stand_age[s];
   else if(ranking[m]==1)stand_rank[s]=sum_econ/size[s];
   else counts[stand_age[s]]++;
   bool qualified=(stand_age[s]>=minage[m] && stand_age[s]<=maxage[m]);
   bool has_optional=false,optional_ok=false;
   for(int r:ma_rules[m]){
    bool matched=false;
    if(rule_highest[r]){
     int maxcount=0;for(int x:species_count)maxcount=std::max(maxcount,x);
     matched=rule_count[r]>0&&rule_count[r]>=maxcount;
    }else if(mode[m])matched=rule_count[r]>0;
    else matched=rule_count[r]>0 && 100.0*rule_count[r]>=rule_pct[r]*size[s]-1e-12;
    if(rule_kind[r]==0 && !matched)qualified=false;
    else if(rule_kind[r]==1){has_optional=true;optional_ok|=matched;}
    else if(rule_kind[r]==2 && matched)qualified=false;
   }
   if(has_optional&&!optional_ok)qualified=false;
   stand_qualified[s]=qualified;
  }
  if(ranking[m]==2){
   for(int s=mo[m];s<mo[m+1];++s){double x=static_cast<double>(counts[stand_age[s]])/static_cast<double>(mo[m+1]-mo[m])*std::exp(std::floor(stand_age[s]/10.0));stand_rank[s]=banker(x);}
  }
  std::vector<int> order;
  for(int s=mo[m];s<mo[m+1];++s){if(!stand_qualified[s])stand_rank[s]=0; if(stand_rank[s]>0)order.push_back(s);}
  std::sort(order.begin(),order.end(),[&](int a,int b){if(stand_rank[a]!=stand_rank[b])return stand_rank[a]>stand_rank[b];return first_pos[a]<first_pos[b];});
  eligible[m]=order.size();
  int n=0;int target_st=0;int have=0;
  if(target_type[m]==1){
   target_st=h04b2::round_even(fractions[m]*order.size());
   if(mode[m] && fractions[m]>0) n=target_st==0 ? 1 : order.size();
   else n=target_st;
   n=std::min(n,static_cast<int>(order.size()));
  }
  target_stands[m]=target_st;
  for(int s:order){
   if(fractions[m]==0)break;
   if(target_type[m]==1 && selected_count[m]>=n)break;
   if(target_type[m]==0 && have>=fractions[m]*total[m])break;
   if(nonempty[s]!=size[s])stop("Selected stand includes empty site: H04-C2 unsupported");
   selected_stands.push_back(s+1);++selected_count[m];have+=size[s];
   for(int j=so[s];j<so[s+1];++j)selected_positions.push_back(pos[j]+1);
  }
  allocated[m]=have;
 }
 return List::create(_["stand_age"]=stand_age,_["stand_rank"]=stand_rank,
   _["stand_qualified"]=stand_qualified,_["selected_stands"]=selected_stands,
   _["selected_positions"]=selected_positions,_["allocated"]=allocated,
   _["total"]=total,_["eligible"]=eligible,_["target_stands"]=target_stands,
   _["selected_count"]=selected_count,_["states_evaluated"]=evaluated);
 END_RCPP
}
