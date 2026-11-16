#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <cstdint>
#include "management_b2_math.hpp"
#include "management_c3_math.hpp"
using namespace Rcpp;
// Separate LANDIS-II RNG stream for H04-C3 reference parity. It is NOT the
// landisR succession RNG. The caller must supply the number of draws already
// consumed by OTHER modules for fully coupled simulations; the future shared
// core RNG will replace this explicit seed+offset boundary.
RcppExport SEXP h04c3_random_plan(SEXP first_, SEXP so_, SEXP mo_, SEXP qualifying_,
                                  SEXP base_rank_, SEXP ranking_, SEXP fraction_,
                                  SEXP type_, SEXP legacy_, SEXP seed_, SEXP skip_) {
 BEGIN_RCPP
 IntegerVector first(first_),so(so_),mo(mo_),ranking(ranking_),type(type_),legacy(legacy_);
 LogicalVector qualifying(qualifying_);
 NumericVector base_rank(base_rank_),fractions(fraction_);
 const int nst=so.size()-1,nma=mo.size()-1;
 if(nst<1||nma<1||first.size()!=nst||qualifying.size()!=nst||base_rank.size()!=nst||
   ranking.size()!=nma||fractions.size()!=nma||type.size()!=nma||legacy.size()!=nma||
   so[0]!=0||mo[0]!=0||mo[nma]!=nst)stop("Invalid C3 input dimensions");
 NumericVector seeds(seed_),skips(skip_);
 if(seeds.size()!=1||skips.size()!=1||!std::isfinite(seeds[0])||!std::isfinite(skips[0])||
    seeds[0]<0||seeds[0]>4294967295.0||seeds[0]!=std::floor(seeds[0])||
    skips[0]<0||skips[0]>1000000000.0||skips[0]!=std::floor(skips[0]))
    stop("Invalid MT19937 seed or existing RNG draw offset");
 h04c3::Stream stream(static_cast<std::uint32_t>(seeds[0]),static_cast<std::uint64_t>(skips[0]));
 NumericVector ranks(nst);
 IntegerVector target_stands(nma),selected_count(nma),eligible(nma),total(nma),allocated(nma);
 std::vector<int> selected;
 for(int m=0;m<nma;++m){
  if(mo[m]<0||mo[m+1]>nst||mo[m]>=mo[m+1]||ranking[m]<0||ranking[m]>3||
     type[m]<0||type[m]>1||legacy[m]<0||legacy[m]>1||
     !std::isfinite(fractions[m])||fractions[m]<0||fractions[m]>1)stop("Invalid C3 management settings");
  std::vector<int> raster_order,ordered;
  for(int s=mo[m];s<mo[m+1];++s){
    if(so[s]<0||so[s+1]<=so[s]||so[s+1]>so[nst]||first[s]<0)stop("Invalid stand topology");
    total[m]+=(so[s+1]-so[s]);raster_order.push_back(s);
  }
  std::sort(raster_order.begin(),raster_order.end(),[&](int a,int b){return first[a]<first[b];});
  for(size_t j=1;j<raster_order.size();++j)
    if(first[raster_order[j]]==first[raster_order[j-1]])stop("Duplicate first raster cell");
  if(ranking[m]==3){
    // Rank all stands, including unqualified ones (the reference does this).
    auto p=h04c3::landis_shuffle_ranks(static_cast<int>(raster_order.size()),stream);
    for(size_t j=0;j<p.size();++j)ranks[raster_order[j]]=p[j];
  } else for(int s:raster_order)ranks[s]=base_rank[s];
  for(int s:raster_order)if(qualifying[s] && ranks[s]>0)ordered.push_back(s);
  eligible[m]=ordered.size();
  std::sort(ordered.begin(),ordered.end(),[&](int a,int b){
    if(ranks[a]!=ranks[b])return ranks[a]>ranks[b];
    return first[a]<first[b];
  });
  int goal=static_cast<int>(ordered.size());
  if(type[m]==1){
    target_stands[m]=h04b2::round_even(fractions[m]*ordered.size());
    if(legacy[m]&&fractions[m]>0)goal=(target_stands[m]==0?1:ordered.size());
    else goal=target_stands[m];
    goal=std::min(goal,static_cast<int>(ordered.size()));
  }
  for(int s:ordered){
    if(fractions[m]==0)break;
    if(type[m]==1 && selected_count[m]>=goal)break;
    if(type[m]==0 && allocated[m]>=fractions[m]*total[m])break;
    selected.push_back(s+1);
    ++selected_count[m]; allocated[m]+=so[s+1]-so[s];
    // Site starting-position selection in LANDIS-II consumes one draw even
    // for a stand consisting of exactly one cell.
    (void)stream.next_double();
  }
 }
 return List::create(_["stand_rank"]=ranks,_["selected_stands"]=selected,
   _["allocated"]=allocated,_["eligible"]=eligible,_["total"]=total,
   _["target_stands"]=target_stands,_["selected_count"]=selected_count,
   _["rng_start_offset"]=skips[0],_["rng_draws_used"]=static_cast<double>(stream.draws)-skips[0],
   _["rng_next_offset"]=static_cast<double>(stream.draws));
 END_RCPP
}
