#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <vector>
#include "management_b2_math.hpp"
#include "management_partial.hpp"
using namespace Rcpp;

// Read-only, once-per-topology rook adjacency cache. Sparse pairs only.
// Active IDs are LANDIS-R core canonical column-major row+(col-1)*nr.
RcppExport SEXP h04d2_topology(SEXP ids_, SEXP pos_, SEXP so_, SEXP dims_) {
 BEGIN_RCPP
 IntegerVector ids(ids_), pos(pos_), so(so_), dims(dims_);
 if((dims.size()!=2 && dims.size()!=3) || dims[0]<1 || dims[1]<1 || so.size()<2 || so[0]!=0 ||
    so[so.size()-1]!=pos.size())stop("Invalid adjacency topology");
 const int nr=dims[0],nc=dims[1],nst=so.size()-1;
 std::unordered_map<int,int> index;index.reserve(ids.size()*2+1);
 IntegerVector first(nst);
 for(int s=0;s<nst;++s){
   int mn=std::numeric_limits<int>::max();
   if(so[s]>=so[s+1]||so[s]<0||so[s+1]>pos.size())stop("Empty/invalid stand");
   for(int j=so[s];j<so[s+1];++j){
      int p=pos[j];if(p<0||p>=ids.size())stop("Invalid active position");
      int id=ids[p];if(id<1||static_cast<double>(id)>static_cast<double>(nr)*nc)stop("Invalid cell ID");
      if(!index.emplace(id,s).second)stop("Duplicate cell in stand topology");
      int row=(id-1)%nr,col=(id-1)/nr;
      mn=std::min(mn,row*nc+col);
   }
   first[s]=mn;
 }
 if(index.size()!=static_cast<std::size_t>(ids.size()))stop("Some active cells are not assigned to stands");
 std::vector<std::pair<int,int>> edges;edges.reserve(4*ids.size());
 for(const auto &it:index){
   int id=it.first,s=it.second,row=(id-1)%nr,col=(id-1)/nr;
   int neighbor[4]={row>0?id-1:0,row+1<nr?id+1:0,col>0?id-nr:0,col+1<nc?id+nr:0};
   for(int n:neighbor){
     auto jt=index.find(n);if(n!=0&&jt!=index.end()&&jt->second!=s)edges.emplace_back(s,jt->second);
   }
 }
 std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
 IntegerVector offsets(nst+1),values(edges.size());
 for(const auto &e:edges)++offsets[e.first+1];
 for(int s=0;s<nst;++s)offsets[s+1]+=offsets[s];
 for(std::size_t j=0;j<edges.size();++j)values[j]=edges[j].second;
 List ans = List::create(_["offsets"]=offsets,_["indices"]=values,_["first_position"]=first);
 if (dims.size() == 3 && dims[2] != 0) {
   // Store the first time each side of a stand boundary was seen in the
   // LANDIS-II raster scan. Each cell checks the upper then left side.
   using Edge = std::tuple<int, int, std::uint64_t>;
   std::vector<Edge> discovery;
   discovery.reserve(edges.size());
   for (const auto& cell : index) {
     const int id = cell.first;
     const int row = (id - 1) % nr, col = (id - 1) / nr, stand = cell.second;
     const int earlier[] = {row > 0 ? id - 1 : 0, col > 0 ? id - nr : 0};
     for (int direction = 0; direction < 2; ++direction) {
       const auto neighbor = index.find(earlier[direction]);
       if (neighbor == index.end() || neighbor->second == stand) continue;
       const auto order = 2 * (static_cast<std::uint64_t>(row) * nc + col) + direction;
       discovery.emplace_back(stand, neighbor->second, order);
       discovery.emplace_back(neighbor->second, stand, order);
     }
   }
   std::sort(discovery.begin(), discovery.end());
   discovery.erase(std::unique(discovery.begin(), discovery.end(), [](const Edge& a, const Edge& b) {
     return std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b);
   }), discovery.end());
   if (discovery.size() != edges.size()) stop("Inconsistent stand adjacency order");
   IntegerVector ordered(values.size());
   for (int s = 0; s < nst; ++s) {
     const auto begin = discovery.begin() + offsets[s], end = discovery.begin() + offsets[s + 1];
     std::sort(begin, end, [](const Edge& a, const Edge& b) { return std::get<2>(a) < std::get<2>(b); });
     for (int i = offsets[s]; i < offsets[s + 1]; ++i) ordered[i] = std::get<1>(discovery[i]);
   }
   ans["landis_order"] = ordered;
 }
 return ans;
 END_RCPP
}

// Build only for partial selectors. Use an active-cell lookup rather
// than an array sized by the bounding raster of a sparse landscape.
RcppExport SEXP h05_cell_topology(SEXP ids_, SEXP positions_, SEXP offsets_, SEXP dimensions_) {
  BEGIN_RCPP
  IntegerVector ids(ids_), positions(positions_), offsets(offsets_), dims(dimensions_);
  const int nr=dims[0], nc=dims[1], count=ids.size();
  IntegerVector raster_positions=clone(positions), stand_of_cell(count);
  std::unordered_map<int,int> active;
  active.reserve(count);
  for(int i=0;i<count;++i) active.emplace(ids[i],i);
  for(int stand=0;stand+1<offsets.size();++stand){
    const auto begin=raster_positions.begin()+offsets[stand];
    const auto end=raster_positions.begin()+offsets[stand+1];
    for(auto cell=begin;cell!=end;++cell) stand_of_cell[*cell]=stand;
    std::sort(begin,end,[&](int a,int b){
      const int left=ids[a]-1,right=ids[b]-1;
      return left%nr<right%nr || (left%nr==right%nr && left<right);
    });
  }
  constexpr int directions=h04e4d::cell_directions;
  const int dr[directions]={-1,1,0,0,-1,-1,1,1};
  const int dc[directions]={0,0,-1,1,-1,1,-1,1};
  IntegerMatrix neighbors(directions,count);
  std::fill(neighbors.begin(),neighbors.end(),-1);
  for(int i=0;i<count;++i){
    if(i%65536==0) checkUserInterrupt();
    const int row=(ids[i]-1)%nr, col=(ids[i]-1)/nr;
    for(int d=0;d<directions;++d){
      const int r=row+dr[d], c=col+dc[d];
      if(r<0||r>=nr||c<0||c>=nc) continue;
      const long long id=1LL+r+static_cast<long long>(c)*nr;
      if(id>std::numeric_limits<int>::max()) continue;
      const auto found=active.find(static_cast<int>(id));
      if(found!=active.end()) neighbors(d,i)=found->second;
    }
  }
  return List::create(_["positions"]=raster_positions,
                      _["neighbors"]=neighbors,_["stand"]=stand_of_cell);
  END_RCPP
}

// Historical Complete-site-selection adjacency/harvest-history planner.
// Eligibility is snapshotted within each MA; recent-harvest history is updated
// after that MA so the next MA observes it, like reference LANDIS-II v7.
// AdjacencyNeighborSetAside is NOT triggered by SiteSelection Complete in v7.
RcppExport SEXP h04d2_plan(SEXP so_,SEXP mo_,SEXP no_,SEXP ni_,SEXP first_,
 SEXP age_,SEXP rank_,SEXP qualifies_,SEXP fractions_,SEXP type_,SEXP legacy_,
 SEXP adj_thr_,SEXP adj_kind_,SEXP min_hist_,SEXP last_,SEXP year_) {
 BEGIN_RCPP
 IntegerVector so(so_),mo(mo_),no(no_),ni(ni_),first(first_),age(age_),type(type_),
  legacy(legacy_),adj_thr(adj_thr_),adj_kind(adj_kind_),min_hist(min_hist_),last(last_);
 NumericVector rank(rank_),frac(fractions_),yearv(year_);
 LogicalVector qualifies(qualifies_);
 const int nst=so.size()-1,nma=mo.size()-1;
 if(nst<1||nma<1||so[0]!=0||mo[0]!=0||mo[nma]!=nst||no.size()!=nst+1||
    no[0]!=0||no[nst]!=ni.size()||first.size()!=nst||age.size()!=nst||
    rank.size()!=nst||qualifies.size()!=nst||last.size()!=nst||
    frac.size()!=nma||type.size()!=nma||legacy.size()!=nma||adj_thr.size()!=nma||
    adj_kind.size()!=nma||min_hist.size()!=nma||yearv.size()!=1||!R_finite(yearv[0])||
    yearv[0]<0||yearv[0]>1e9||yearv[0]!=std::floor(yearv[0]))stop("Invalid D2 inputs");
 const int year=static_cast<int>(yearv[0]);
 for(int s=0;s<nst;++s){
  if(so[s]<0||so[s+1]<=so[s]||no[s]<0||no[s+1]<no[s]||
     no[s+1]>ni.size()||first[s]<0||age[s]<0||!R_finite(rank[s])||
     (last[s]!=NA_INTEGER&&(last[s]<0||last[s]>year)))stop("Invalid D2 stand/topology/history");
 }
 for(int j=0;j<ni.size();++j)if(ni[j]<0||ni[j]>=nst)stop("Invalid D2 neighbor");
 for(int m=0;m<nma;++m){
   if(mo[m]<0||mo[m+1]<=mo[m]||!R_finite(frac[m])||frac[m]<0||frac[m]>1||
      (type[m]!=0&&type[m]!=1)||(legacy[m]!=0&&legacy[m]!=1)||adj_thr[m]<0||
      adj_kind[m]<0||adj_kind[m]>2||min_hist[m]<0)stop("Invalid D2 MA settings");
 }
 IntegerVector updated=clone(last),allocated(nma),eligible(nma),total(nma),
   target_stands(nma),selected_count(nma);
 LogicalVector spatially_eligible(nst);
 std::vector<int> selected;
 selected.reserve(nst);
 for(int m=0;m<nma;++m){
   if(m%256==0)checkUserInterrupt();
   std::vector<int> ranked;ranked.reserve(mo[m+1]-mo[m]);
   for(int s=mo[m];s<mo[m+1];++s){
      total[m]+=so[s+1]-so[s];
      if(!qualifies[s]||rank[s]<=0)continue;
      if(last[s]!=NA_INTEGER && year-last[s]<min_hist[m])continue;
      bool allowed=true;
      if(adj_thr[m]>0 && adj_kind[m]!=0){
       for(int k=no[s];k<no[s+1];++k){
         int other=ni[k];
         if(adj_kind[m]==1 && age[other]<adj_thr[m]){allowed=false;break;}
         if(adj_kind[m]==2 && updated[other]!=NA_INTEGER && year-updated[other]<adj_thr[m]){
            allowed=false;break;
         }
       }
      }
      if(allowed){spatially_eligible[s]=true;ranked.push_back(s);}
   }
   std::sort(ranked.begin(),ranked.end(),[&](int a,int b){
      if(rank[a]!=rank[b])return rank[a]>rank[b];
      return first[a]<first[b];
   });
   eligible[m]=ranked.size();
   int goal=0;
   if(type[m]==1){
      // H04-B2 default and original %Stands are preserved as in C2.
      target_stands[m]=h04b2::round_even(frac[m]*ranked.size());
      if(legacy[m] && frac[m]>0)goal=target_stands[m]==0 ? 1 : ranked.size();
      else goal=target_stands[m];
      goal=std::min(goal,static_cast<int>(ranked.size()));
   }
   for(int s:ranked){
     if(frac[m]==0)break;
     if(type[m]==1 && selected_count[m]>=goal)break;
     if(type[m]==0 && allocated[m]>=frac[m]*total[m])break;
     selected.push_back(s+1);
     ++selected_count[m];allocated[m]+=so[s+1]-so[s];
   }
   // Historical v7: eligibility for this MA was already computed before
   // cutting; only the NEXT MA observes the new last-harvest timestamps.
   for(int j=0;j<selected_count[m];++j){
      int s=selected[selected.size()-selected_count[m]+j]-1;
      updated[s]=year;
   }
 }
 return List::create(_["selected_stands"]=selected,_["last_harvest_year"]=updated,
  _["stand_eligible"]=spatially_eligible,_["eligible"]=eligible,
  _["allocated"]=allocated,_["total"]=total,
  _["target_stands"]=target_stands,_["selected_count"]=selected_count);
 END_RCPP
}
