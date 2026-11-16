#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "management_math.hpp"
using namespace Rcpp;
// All indices are zero-based except cell_state, which is the core's one-based index.
RcppExport SEXP h04_plan(SEXP states_,SEXP cell_state_,SEXP positions_,SEXP stand_offsets_,SEXP ma_offsets_,SEXP targets_,SEXP minima_) {
 BEGIN_RCPP
 List states(states_); IntegerVector cs(cell_state_),pos(positions_),so(stand_offsets_),mo(ma_offsets_),minimum(minima_); NumericVector targets(targets_);
 IntegerVector off=as<IntegerVector>(states["offset"]); NumericVector ages=as<NumericVector>(states["age"]),bio=as<NumericVector>(states["biomass"]);
 int ns=off.size()-1,nst=so.size()-1,nma=mo.size()-1;
 if(ns<1||nst<1||nma<1||off[0]!=0||off[ns]!=ages.size()||bio.size()!=ages.size()||so[0]!=0||so[nst]!=pos.size()||mo[0]!=0||mo[nma]!=nst||targets.size()!=nma||minimum.size()!=nma) stop("Invalid H04 layout");
 std::vector<int> maxage(ns,-1),standage(nst),sizes(nst),nonempty(nst);
 int evaluated=0;
 for(int s=0;s<nst;++s){
   if(s%4096==0) checkUserInterrupt();
   if(so[s]<0||so[s+1]<=so[s]||so[s+1]>pos.size()) stop("Invalid stand offsets");
   sizes[s]=so[s+1]-so[s]; double sum=0;
   for(int j=so[s];j<so[s+1];++j){
     int p=pos[j]; if(p<0||p>=cs.size()) stop("Invalid active cell position");
     int k=cs[p]-1; if(k<0||k>=ns) stop("Invalid state index");
     if(maxage[k]<0){
       if(off[k]<0||off[k+1]<off[k]||off[k+1]>ages.size()) stop("Invalid CSR offset");
       int a=0; for(int q=off[k];q<off[k+1];++q){
         if(!R_finite(ages[q])||ages[q]<1||ages[q]>65535||ages[q]!=std::trunc(ages[q])||!R_finite(bio[q])||bio[q]<=0) stop("Invalid living cohort");
         a=std::max(a,static_cast<int>(ages[q]));
       } maxage[k]=a; ++evaluated;
     }
     sum+=maxage[k]; if(maxage[k]>0) ++nonempty[s];
   }
   standage[s]=static_cast<int>(sum/sizes[s]);
 }
 std::vector<int> selected_stands,selected_positions,allocated(nma),total(nma),eligible(nma);
 for(int m=0;m<nma;++m){
   if(mo[m]<0||mo[m+1]<=mo[m]||mo[m+1]>nst||!R_finite(targets[m])||targets[m]<0||targets[m]>1||minimum[m]<0) stop("Invalid MA allocation");
   for(int s=mo[m];s<mo[m+1];++s){total[m]+=sizes[s];if(standage[s]>0&&standage[s]>=minimum[m]) ++eligible[m];}
   h04::Plan p=h04::allocate(standage,sizes,mo[m],mo[m+1],targets[m],minimum[m]);
   allocated[m]=p.count;
   for(int s:p.selected){
     // Empty-cell coverage and zero-removal allocation are not qualified yet.
     if(nonempty[s]!=sizes[s]) stop("Selected stand contains empty cells: H04-A2 coverage not qualified");
     selected_stands.push_back(s+1);
     for(int j=so[s];j<so[s+1];++j) selected_positions.push_back(pos[j]+1);
   }
 }
 return List::create(_["stand_age"]=standage,_["selected_stands"]=selected_stands,_["selected_positions"]=selected_positions,_["allocated"]=allocated,_["total"]=total,_["eligible"]=eligible,_["states_evaluated"]=evaluated);
 END_RCPP
}
