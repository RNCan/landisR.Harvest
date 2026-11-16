#ifndef LANDISR_H04E3_MATH_HPP
#define LANDISR_H04E3_MATH_HPP
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include "management_e2_math.hpp"
namespace h04e3 {
struct Prescription {
  int ma, min_age, max_age, first, last, adjacent_limit, adjacent_type, min_harvest_gap;
  double target_fraction;
};
struct Event {
  int stand, pres, ma, age, number_of_sites;
  std::uint64_t draw;
  double uniform, probability;
};
struct Result {
  std::vector<Event> events;
  std::vector<int> last_harvest, allocated, eligible, total_sites;
  std::uint64_t offset_end;
};
inline Result select(
  const std::vector<int>& so,const std::vector<int>& mo,
  const std::vector<int>& no,const std::vector<int>& ni,
  const std::vector<int>& first,const std::vector<int>& age,
  const std::vector<double>& ranks, const std::vector<unsigned char>& qualifies,
  const std::vector<Prescription>& ps, const std::vector<int>& last,
  int year,std::uint32_t seed,std::uint64_t offset) {
  // matrices 'ranks' and 'qualifies': prescription-major, i.e. j*ns+s.
  const int ns=static_cast<int>(so.size())-1,nm=static_cast<int>(mo.size())-1,np=static_cast<int>(ps.size());
  if(ns<1||nm<1||np<1||so[0]!=0||mo[0]!=0||mo[nm]!=ns||
      no.size()!=so.size()||no[0]!=0||no[ns]!=static_cast<int>(ni.size())||
      first.size()!=static_cast<std::size_t>(ns)||age.size()!=static_cast<std::size_t>(ns)||
      ranks.size()!=static_cast<std::size_t>(ns*np)||qualifies.size()!=static_cast<std::size_t>(ns*np)||
      last.size()!=static_cast<std::size_t>(ns)||year<0)
    throw std::runtime_error("invalid E3 landscape/plan dimensions");
  Result out;out.last_harvest=last;out.allocated.assign(np,0);out.eligible.assign(np,0);out.total_sites.assign(nm,0);
  for(int s=0;s<ns;++s){
    if(so[s]>=so[s+1]||no[s]>no[s+1]||first[s]<0||age[s]<0||
      (last[s]>=0&&last[s]>year))throw std::runtime_error("invalid E3 stand/history");
  }
  for(int n:ni)if(n<0||n>=ns)throw std::runtime_error("invalid E3 neighbor index");
  for(const auto &p:ps){
    if(p.ma<0||p.ma>=nm||p.min_age<0||p.max_age<p.min_age||p.first<0||p.last<p.first||
      p.adjacent_limit<0||p.adjacent_type<0||p.adjacent_type>2||p.min_harvest_gap<0||
      p.target_fraction<0||p.target_fraction>1||!std::isfinite(p.target_fraction))
      throw std::runtime_error("invalid E3 prescription");
  }
  for(int m=0;m<nm;++m) {
    if(mo[m]>=mo[m+1])throw std::runtime_error("empty management area");
    for(int s=mo[m];s<mo[m+1];++s)out.total_sites[m]+=so[s+1]-so[s];
  }
  h04c3::Stream rng(seed,offset);
  std::vector<unsigned char> used(ns,0);
  for(int m=0;m<nm;++m){
    std::vector<int> active;
    std::vector<std::vector<int>> ordered;
    std::vector<double> goals;
    for(int j=0;j<np;++j){
      const auto &p=ps[j];
      if(p.ma!=m||p.target_fraction==0||year<p.first||year>p.last)continue;
      std::vector<int> candidates;
      for(int s=mo[m];s<mo[m+1];++s){
        int ix=j*ns+s;
        if(!qualifies[ix]||ranks[ix]<=0||!std::isfinite(ranks[ix])||
           age[s]<p.min_age||age[s]>p.max_age)continue;
        if(last[s]>=0&&year-last[s]<p.min_harvest_gap)continue;
        bool ok=true;
        if(p.adjacent_type&&p.adjacent_limit){
          for(int k=no[s];k<no[s+1];++k){
            int neighbor=ni[k];
            if(p.adjacent_type==1&&age[neighbor]<p.adjacent_limit){ok=false;break;}
            if(p.adjacent_type==2&&out.last_harvest[neighbor]>=0&&
               year-out.last_harvest[neighbor]<p.adjacent_limit){ok=false;break;}
          }
        }
        if(ok)candidates.push_back(s);
      }
      std::sort(candidates.begin(),candidates.end(),[&](int a,int b){
        if(ranks[j*ns+a]!=ranks[j*ns+b])return ranks[j*ns+a]>ranks[j*ns+b];
        return first[a]<first[b];
      });
      out.eligible[j]=candidates.size();
      double goal=out.total_sites[m]*p.target_fraction;
      if(!candidates.empty()&&goal>0){active.push_back(j);goals.push_back(goal);ordered.push_back(std::move(candidates));}
    }
    std::vector<int> newly_cut;
    while(!active.empty()){
      std::vector<double> weights(active.size(),0.0);
      double total=0.0;
      for(std::size_t k=0;k<active.size();++k){
        bool available=false;
        for(int s:ordered[k])if(!used[s]){available=true;break;}
        if(available)weights[k]=std::max(0.0,(goals[k]-out.allocated[active[k]])/goals[k]);
        total+=weights[k];
      }
      if(total<=0)break;
      const std::uint64_t draw=rng.draws;
      const double u=rng.next_double();
      double cumulative=0;int kchosen=-1;
      for(std::size_t k=0;k<active.size();++k){
        cumulative+=weights[k]/total;if(u<cumulative){kchosen=static_cast<int>(k);break;}
      }
      if(kchosen<0)for(int k=static_cast<int>(active.size())-1;k>=0;--k)
        if(weights[k]>0){kchosen=k;break;}
      if(kchosen<0)throw std::runtime_error("E3 random quota selection failed");
      int j=active[kchosen],s=-1;
      for(int x:ordered[kchosen])if(!used[x]){s=x;break;}
      if(s<0)throw std::runtime_error("E3 selected exhausted prescription");
      used[s]=1;newly_cut.push_back(s);
      const int cells=so[s+1]-so[s];out.allocated[j]+=cells;
      out.events.push_back({s,j,m,age[s],cells,draw,u,weights[kchosen]/total});
    }
    // H04-D2 observed semantics: within-MA candidate snapshot, then
    // only subsequent management areas observe current-year history.
    for(int s:newly_cut)out.last_harvest[s]=year;
  }
  out.offset_end=rng.draws;
  return out;
}
} // namespace
#endif
