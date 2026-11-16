#ifndef LANDISR_H04E2_MATH_HPP
#define LANDISR_H04E2_MATH_HPP
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include "management_c3_math.hpp"

namespace h04e2 {
struct Prescription {
  int ma, minimum_age, maximum_age, first_year, last_year, target_type;
  double target_fraction;
};
struct Event {
  int stand, prescription, management_area, age;
  double selection_u, probability;
  std::uint64_t rng_draw_index;
};
struct Result {
  std::vector<Event> events;
  std::vector<int> sites, stands, eligible, total_sites;
  std::uint64_t next_offset;
};
inline int nearest_even(double x) {
  // R and .NET System.Math.Round both use ties-to-even in the validated fixtures.
  if(x<0 || !std::isfinite(x))throw std::runtime_error("invalid rounded quota");
  double floorx=std::floor(x), rem=x-floorx;
  return static_cast<int>(rem>.5 || (rem==.5 && std::fmod(floorx,2.0)==1.0) ? floorx+1 : floorx);
}
inline Result select(const std::vector<int>& offsets,
                     const std::vector<int>& ma_offsets,
                     const std::vector<int>& raster_first,
                     const std::vector<int>& ages,
                     const std::vector<double>& ranks,
                     const std::vector<Prescription>& prescriptions,
                     int year,std::uint32_t seed,std::uint64_t offset,
                     bool landis_compatible) {
  const int ns=static_cast<int>(offsets.size())-1, nm=static_cast<int>(ma_offsets.size())-1;
  if(ns<1 || nm<1 || ma_offsets[0]!=0 || ma_offsets[nm]!=ns || offsets[0]!=0 ||
     ages.size()!=static_cast<std::size_t>(ns) || ranks.size()!=static_cast<std::size_t>(ns) ||
     raster_first.size()!=static_cast<std::size_t>(ns) || year<0 || prescriptions.empty())
       throw std::runtime_error("invalid E2 topology/age/rank configuration");
  Result out;const int np=static_cast<int>(prescriptions.size());
  out.sites.assign(np,0);out.stands.assign(np,0);out.eligible.assign(np,0);out.total_sites.assign(nm,0);
  for(int s=0;s<ns;++s) {
    if(offsets[s]>=offsets[s+1] || ages[s]<0 || raster_first[s]<0 || !std::isfinite(ranks[s]))
      throw std::runtime_error("invalid stand or rank");
  }
  for(int m=0;m<nm;++m){
    if(ma_offsets[m]>=ma_offsets[m+1])throw std::runtime_error("empty management area");
    for(int s=ma_offsets[m];s<ma_offsets[m+1];++s)out.total_sites[m]+=offsets[s+1]-offsets[s];
  }
  for(const auto& p:prescriptions){
    if(p.ma<0||p.ma>=nm || !std::isfinite(p.target_fraction) ||p.target_fraction<0||p.target_fraction>1||
       p.target_type<0||p.target_type>1||p.minimum_age<0||p.maximum_age<p.minimum_age||
       p.first_year<0||p.last_year<p.first_year)
       throw std::runtime_error("invalid E2 prescription");
    if(landis_compatible && p.target_fraction==0)
       throw std::runtime_error("LANDIS-II rejects an explicit 0% HarvestImplementations row");
  }
  h04c3::Stream rng(seed,offset);
  std::vector<bool> used(static_cast<std::size_t>(ns),false);
  for(int m=0;m<nm;++m){
    std::vector<int> active;
    std::vector<double> goals;
    std::vector<int> available_count;
    std::vector<std::vector<int>> ranked;
    for(int j=0;j<np;++j){
      const auto& p=prescriptions[j];if(p.ma!=m || p.target_fraction<=0||year<p.first_year||year>p.last_year)continue;
      std::vector<int> order;
      for(int s=ma_offsets[m];s<ma_offsets[m+1];++s)
        if(ranks[s]>0 && ages[s]>=p.minimum_age && ages[s]<=p.maximum_age) order.push_back(s);
      std::sort(order.begin(),order.end(),[&](int a,int b){
        if(ranks[a]!=ranks[b])return ranks[a]>ranks[b];
        return raster_first[a]<raster_first[b];
      });
      if(order.empty())continue;
      out.eligible[j]=static_cast<int>(order.size());
      const double target=p.target_type==0?out.total_sites[m]*p.target_fraction:
                  nearest_even(order.size()*p.target_fraction);
      if(target<=0)continue;
      active.push_back(j);goals.push_back(target);
      available_count.push_back(0);ranked.emplace_back(std::move(order));
    }
    while(!active.empty()){
      std::vector<double> weight(active.size(),0.0);
      double total=0;
      for(std::size_t j=0;j<active.size();++j){
        int pj=active[j];const auto &p=prescriptions[pj];
        bool next_found=false;
        for(int s:ranked[j])if(!used[s]){next_found=true;break;}
        if(next_found){
          const double allocated=p.target_type==0?out.sites[pj]:out.stands[pj];
          weight[j]=std::max(0.0,(goals[j]-allocated)/goals[j]);
        }
        total+=weight[j];
      }
      if(total<=0)break;
      const std::uint64_t draw=rng.draws;
      const double u=rng.next_double();
      double cumulative=0;std::size_t chosen=active.size();
      for(std::size_t j=0;j<active.size();++j){
        cumulative+=weight[j]/total;
        if(u<cumulative){chosen=j;break;}
      }
      if(chosen==active.size()){
        // Equivalent to assigning the final positive probability the rest of [0,1).
        for(int j=static_cast<int>(active.size())-1;j>=0;--j)if(weight[j]>0){chosen=j;break;}
      }
      if(chosen==active.size())throw std::runtime_error("E2 invalid selector probability");
      int ps=active[chosen], s=-1;
      for(int cand:ranked[chosen])if(!used[cand]){s=cand;break;}
      if(s<0)throw std::runtime_error("selected prescription exhausted");
      used[s]=true;
      out.stands[ps]+=1;out.sites[ps]+=offsets[s+1]-offsets[s];
      out.events.push_back(Event{s,ps,m,ages[s],u,weight[chosen]/total,draw});
    }
  }
  out.next_offset=rng.draws;
  return out;
}
}
#endif
