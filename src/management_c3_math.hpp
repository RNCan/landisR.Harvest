#ifndef LANDISR_H04C3_MATH_HPP
#define LANDISR_H04C3_MATH_HPP
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>
#include <landisR/mt19937.hpp>

// LANDIS-II v7 core / Troschuetz.Random MT19937 generator parity.
// The underlying MT19937 is seeded with a 32-bit integer. NextDouble uses
// the upper 31 bits of ONE uint32, unlike std::uniform_real_distribution.
namespace h04c3 {
using Stream = landis_rng::Stream;
// Fenwick tree removes the selected rank without changing the order of
// remaining values, equivalent to Landis.Model.shuffle List.RemoveAt.
class Fenwick {
  std::vector<int> tree;
  int n;
public:
  explicit Fenwick(int count):tree(count+1,0),n(count){
    for(int i=1;i<=n;++i){tree[i]++;int j=i+(i&-i);if(j<=n)tree[j]+=tree[i];}
  }
  int kth(int k) const {
    if(k<1)throw std::runtime_error("invalid kth");
    int idx=0,step=1;
    while(step<=n/2)step*=2;
    for(;step>0;step/=2){int next=idx+step;if(next<=n && tree[next]<k){idx=next;k-=tree[next];}}
    return idx+1;
  }
  void remove(int i){for(;i<=n;i+=(i&-i))tree[i]--;}
};
inline std::vector<int> landis_shuffle_ranks(int n,Stream &rng){
  Fenwick f(n);std::vector<int> p;p.reserve(n);
  for(int r=n;r>=1;--r){int k=f.kth(rng.next(r)+1);p.push_back(k);f.remove(k);}
  return p;
}
}
#endif
