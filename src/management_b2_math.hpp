#ifndef LANDISR_H04B2_MATH_HPP
#define LANDISR_H04B2_MATH_HPP
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
namespace h04b2 {
// Math.Round(x) in C# defaults to ties-to-even for non-negative x.
inline int round_even(double x) {
  if(!std::isfinite(x)||x<0||x>2147480000.) throw std::runtime_error("invalid stands target");
  double lo=std::floor(x),delta=x-lo;
  if(std::fabs(delta-.5)<1e-12) return static_cast<int>(lo)+(static_cast<long long>(lo)%2LL==1LL ? 1:0);
  return static_cast<int>(lo)+(delta>.5 ? 1:0);
}
struct Allocation {std::vector<int> selected; int sites=0; int target_stands=0; int eligible=0;};
// method 0 = target area, 1 = target stands; legacy only applies to site selection Complete.
inline Allocation allocate(const std::vector<double>& rank,const std::vector<int>& size,
 int first,int last,double fraction,int minimum_rank,bool stands,bool legacy) {
  Allocation out;
  if(!std::isfinite(fraction)||fraction<0||fraction>1 || first<0 || last<first || last>(int)rank.size()) throw std::runtime_error("invalid allocation arguments");
  if(fraction==0) return out;
  std::vector<int> eligible;
  int total=0;
  for(int s=first;s<last;++s) {
    if(size[s]<=0||rank[s]<0||!std::isfinite(rank[s])) throw std::runtime_error("invalid stand rank or size");
    total+=size[s];
    if(rank[s]>minimum_rank) eligible.push_back(s);
  }
  std::sort(eligible.begin(),eligible.end(),[&](int a,int b){return rank[a]>rank[b];});
  for(size_t i=1;i<eligible.size();++i)
    if(eligible[i-1]!=eligible[i] && rank[eligible[i-1]]==rank[eligible[i]])
      throw std::runtime_error("Eligible stand rank tie: not qualified in H04-B2");
  out.eligible=static_cast<int>(eligible.size());
  if(stands) {
    out.target_stands=round_even(fraction*out.eligible);
    // Exact tested legacy fixtures: target rounded to zero still selects one;
    // positive target keeps choosing until all qualified stands are exhausted.
    int n=legacy ? (out.target_stands==0 ? 1 : out.eligible) : out.target_stands;
    n=std::min(n,out.eligible);
    for(int i=0;i<n;++i){out.selected.push_back(eligible[i]);out.sites+=size[eligible[i]];}
  } else {
    const double target=fraction*total;
    for(int s:eligible){if(out.sites>=target)break;out.selected.push_back(s);out.sites+=size[s];}
  }
  return out;
}
}
#endif
