#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "management_b2_math.hpp"
using namespace Rcpp;
// One compact CSR state shared across cells. Cohort species_id is zero based.
// Multiple MAs may use different economic species-rank tables.
RcppExport SEXP h04b2_plan(SEXP states_,SEXP cell_state_,SEXP positions_,SEXP stand_offsets_,SEXP ma_offsets_,
 SEXP ranks_,SEXP merchant_,SEXP ranking_,SEXP fractions_,SEXP minimum_age_,SEXP target_type_,SEXP mode_) {
 BEGIN_RCPP
 List states(states_);
 IntegerVector cs(cell_state_),pos(positions_),so(stand_offsets_),mo(ma_offsets_),
   ranking(ranking_),minimum_age(minimum_age_),target_type(target_type_),mode(mode_);
 NumericMatrix weights(ranks_); IntegerMatrix merchant(merchant_);
 NumericVector fractions(fractions_);
 IntegerVector off=as<IntegerVector>(states["offset"]),sp=as<IntegerVector>(states["species_id"]);
 NumericVector ages=as<NumericVector>(states["age"]),biomass=as<NumericVector>(states["biomass"]);
 int ns=off.size()-1,nst=so.size()-1,nma=mo.size()-1,nsp=weights.ncol();
 if(ns<1||nst<1||nma<1||nsp<1||off[0]!=0||off[ns]!=ages.size()||sp.size()!=ages.size()||biomass.size()!=ages.size()||so[0]!=0||so[nst]!=pos.size()||mo[0]!=0||mo[nma]!=nst||weights.nrow()!=nma||merchant.nrow()!=nma||merchant.ncol()!=nsp||ranking.size()!=nma||minimum_age.size()!=nma||target_type.size()!=nma||mode.size()!=nma||fractions.size()!=nma) stop("Invalid H04-B2 CSR/topology/target dimensions");
 for(int q=0;q<ns;++q)if(off[q]<0||off[q+1]<off[q]||off[q+1]>ages.size())stop("Invalid CSR offsets");
 for(int q=0;q<ages.size();++q)if(sp[q]<0||sp[q]>=nsp||!R_finite(ages[q])||ages[q]<1||ages[q]>65535||ages[q]!=std::trunc(ages[q])||!R_finite(biomass[q])||biomass[q]<=0)stop("Invalid living cohort");
 std::vector<int> max_age(ns,-1),max_age_used(ns,0),size(nst),nonempty(nst,0),stand_age(nst);
 NumericVector stand_rank(nst); std::vector<int> selected_stands,selected_positions;
 IntegerVector allocated(nma),total(nma),eligible(nma),target_stands(nma),selected_count(nma);
 int evaluated=0;
 for(int m=0;m<nma;++m) {
   if(mo[m]<0||mo[m+1]<=mo[m]||mo[m+1]>nst||minimum_age[m]<0||!R_finite(fractions[m])||fractions[m]<0||fractions[m]>1||ranking[m]<0||ranking[m]>1||target_type[m]<0||target_type[m]>1||mode[m]<0||mode[m]>1)stop("Invalid MA allocation");
   // Economic ranks must be calculated separately for each MA's prescription.
   std::vector<double> econ(ns,std::numeric_limits<double>::quiet_NaN());
   for(int s=mo[m];s<mo[m+1];++s) {
     if(s%4096==0)checkUserInterrupt();
     if(so[s]<0||so[s+1]<=so[s]||so[s+1]>pos.size())stop("Invalid stand offsets");
     size[s]=so[s+1]-so[s];total[m]+=size[s];
     double age_sum=0,rank_sum=0;
     for(int j=so[s];j<so[s+1];++j) {
       int p=pos[j];if(p<0||p>=cs.size())stop("Invalid active cell position");
       int state=cs[p]-1;if(state<0||state>=ns)stop("Invalid state reference");
       if(max_age[state]<0) {
         int a=0;for(int q=off[state];q<off[state+1];++q)a=std::max(a,static_cast<int>(ages[q]));
         max_age[state]=a; ++evaluated;
       }
       age_sum+=max_age[state];if(max_age[state]>0)++nonempty[s];
       if(ranking[m]==1) {
         if(std::isnan(econ[state])) {
           double r=0;
           for(int q=off[state];q<off[state+1];++q) {
             int species=sp[q],ma=merchant(m,species); double w=weights(m,species);
             if(w>0 && ma>0 && ages[q]>=ma)r+=w*ages[q]/ma;
           }
           econ[state]=r;++evaluated;
         }
         rank_sum+=econ[state];
       }
     }
     stand_age[s]=static_cast<int>(age_sum/size[s]);
     stand_rank[s]=(ranking[m]==0) ? stand_age[s] : rank_sum/size[s];
     // Unlike some legacy C# logs, the R plan keeps the PRE-HARVEST age.
     if(stand_age[s]<minimum_age[m])stand_rank[s]=0;
   }
   std::vector<double> ranks(nst);for(int s=mo[m];s<mo[m+1];++s)ranks[s]=stand_rank[s];
   h04b2::Allocation allocation=h04b2::allocate(ranks,size,mo[m],mo[m+1],fractions[m],0,target_type[m]==1,mode[m]==1);
   eligible[m]=allocation.eligible;allocated[m]=allocation.sites;selected_count[m]=allocation.selected.size();target_stands[m]=allocation.target_stands;
   for(int s:allocation.selected) {
     if(nonempty[s]!=size[s])stop("Selected stand contains empty cells: H04-B2 not qualified");
     selected_stands.push_back(s+1);
     for(int j=so[s];j<so[s+1];++j)selected_positions.push_back(pos[j]+1);
   }
 }
 return List::create(_["stand_age"]=stand_age,_["stand_rank"]=stand_rank,
   _["selected_stands"]=selected_stands,_["selected_positions"]=selected_positions,
   _["total"]=total,_["allocated"]=allocated,_["eligible"]=eligible,
   _["target_stands"]=target_stands,_["selected_count"]=selected_count,
   _["states_evaluated"]=evaluated);
 END_RCPP
}
