#ifndef LANDISR_H04_MANAGEMENT_MATH_HPP
#define LANDISR_H04_MANAGEMENT_MATH_HPP
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
namespace h04 {
struct Plan {std::vector<int> selected; int count=0;};
// Equal-size cells: area target in cell units; full stands may overshoot.
inline Plan allocate(const std::vector<int>& age,const std::vector<int>& size,
                     int first,int last,double fraction,int minimum) {
  Plan p; if(fraction==0) return p;
  int total=0; std::vector<int> candidates;
  for(int s=first;s<last;++s){total+=size[s]; if(age[s]>0 && age[s]>=minimum) candidates.push_back(s);}
  std::sort(candidates.begin(),candidates.end(),[&](int a,int b){return age[a]>age[b];});
  // Deliberately do not invent LANDIS tie-breaking. Reject eligible ties.
  for(size_t j=1;j<candidates.size();++j) if(age[candidates[j]]==age[candidates[j-1]]) throw std::runtime_error("Eligible stand-age tie: H04-A2 does not yet support ties");
  const double target=fraction*total;
  for(int s:candidates){if(p.count>=target) break; p.selected.push_back(s);p.count+=size[s];}
  return p;
}
}
#endif
