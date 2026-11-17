#ifndef LANDISR_H04E4D_MATH_HPP
#define LANDISR_H04E4D_MATH_HPP
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "management_e2_math.hpp"
#include "management_patch.hpp"
namespace h04e4d {
struct Prescription {
  int ma, min_age, max_age, first, last, adjacent_limit, adjacent_type, min_harvest_gap, ranking;
  double target_fraction;
  bool spread = false;
  double min_block_ha = 0.0, max_block_ha = 0.0;
  bool partial = false;
  double patch_fraction = 0.0;
  double patch_area = 0.0;
  bool patch_overlap = false;
  int set_aside = 0;
  int repeat_kind = 0, repeat_interval = 0, repeat_count = 0;
  int min_damage = 0;
};
struct Event {
  int stand, pres, ma, age, number_of_sites, id;
  double rank;
  std::uint64_t draw;
  double uniform, probability;
  std::size_t cell_begin, cell_end;
  int repeat_number = 0;
};
// Pending work contains identifiers and the footprint of a harvested block.
struct RepeatJob {
  int prescription, year, iteration, event_id;
  double rank;
  std::vector<int> stands;
  std::vector<std::vector<int>> cells;
};
struct RepeatPlan {
  std::vector<RepeatJob> jobs;
  int last_event_id = 0;
  double current_rank = 0;
  bool enabled = false;
  bool compatible = false;
  int end_year = std::numeric_limits<int>::max();
  const int* rule_index = nullptr;
};
struct Result {
  std::vector<Event> events;
  std::vector<int> selected_cells;
  std::vector<int> last_harvest, allocated, attempted, eligible, total_sites;
  std::vector<double> set_aside_until;
  std::uint64_t offset_end;
  std::vector<RepeatJob> repeat_jobs;
  int last_event_id = 0;
  double current_rank = 0;
  std::vector<int> damage_year;
};
// Each prescription stores entries only for stands in its own area.
inline std::vector<std::size_t> prescription_offsets(
    const std::vector<int>& management_offsets, const std::vector<Prescription>& prescriptions) {
  std::vector<std::size_t> offsets(prescriptions.size() + 1, 0);
  for (std::size_t j = 0; j < prescriptions.size(); ++j) {
    const int m = prescriptions[j].ma;
    if (m < 0 || static_cast<std::size_t>(m) >= management_offsets.size() - 1 ||
        management_offsets[m] > management_offsets[m + 1])
      throw std::runtime_error("Invalid prescription management area");
    offsets[j + 1] = offsets[j] + management_offsets[m + 1] - management_offsets[m];
  }
  return offsets;
}
inline Result select(
  const std::vector<int>& so,const std::vector<int>& mo,
  const std::vector<int>& no,const std::vector<int>& ni,
  const std::vector<int>& first,const std::vector<int>& age,
  std::vector<double> ranks, const std::vector<unsigned char>& qualifies,
  const std::vector<Prescription>& ps, const std::vector<int>& last,
  int year, h04c3::Stream& rng,
  const std::vector<int>& effective_sites = {}, double cell_area_ha = 1.0,
  const PartialTopology& cell_topology = {},
  const std::vector<double>& reserve_until = {},
  const RepeatPlan& repeat = {}, DamageHistory damage = {},
  const std::vector<int>& listed_sites = {}) {
  const auto offsets = prescription_offsets(mo, ps);
  const int ns=static_cast<int>(so.size())-1,nm=static_cast<int>(mo.size())-1,np=static_cast<int>(ps.size());
  if(ns<1||nm<1||np<1||so[0]!=0||mo[0]!=0||mo[nm]!=ns||
      no.size()!=so.size()||no[0]!=0||no[ns]!=static_cast<int>(ni.size())||
      first.size()!=static_cast<std::size_t>(ns)||age.size()!=static_cast<std::size_t>(ns)||
      ranks.size()!=offsets.back()||qualifies.size()!=offsets.back()||
      last.size()!=static_cast<std::size_t>(ns)||year<0||
      !std::isfinite(cell_area_ha)||cell_area_ha<=0||
    (!effective_sites.empty() && effective_sites.size()!=ranks.size()) ||
    (!reserve_until.empty() && reserve_until.size()!=static_cast<std::size_t>(ns)))
    throw std::runtime_error("invalid E3 landscape/plan dimensions");
  Result out;out.last_harvest=last;out.allocated.assign(np,0);out.attempted.assign(np,0);out.eligible.assign(np,0);out.total_sites.assign(nm,0);
  out.set_aside_until=reserve_until;
  for(double& expiry:out.set_aside_until) if(expiry<=year) expiry=-1;
  std::vector<unsigned char> pending;
  std::vector<std::vector<const RepeatJob*>> due;
  if (repeat.enabled) {
    pending.assign(ns, 0);
    due.resize(nm);
    for (const auto& job : repeat.jobs) {
      for (int stand : job.stands) pending[stand] = 1;
      if (job.year <= year) due[ps[job.prescription].ma].push_back(&job);
      else out.repeat_jobs.push_back(job);
    }
    for (auto& jobs : due)
      std::stable_sort(jobs.begin(), jobs.end(), [](const RepeatJob* a, const RepeatJob* b){
        if (a->prescription != b->prescription) return a->prescription < b->prescription;
        return a->year < b->year;
      });
  }
  const auto reserved=[&](int stand){
    if (!pending.empty() && pending[stand]) return true;
    return !out.set_aside_until.empty() && out.set_aside_until[stand]>year;
  };
  for(int s=0;s<ns;++s){
    if(so[s]>=so[s+1]||no[s]>no[s+1]||first[s]<0||age[s]<0||
      (last[s]>=0&&last[s]>year))throw std::runtime_error("invalid E3 stand/history");
  }
  for(int n:ni)if(n<0||n>=ns)throw std::runtime_error("invalid E3 neighbor index");
  for(const auto &p:ps){
    if(p.ma<0||p.ma>=nm||p.min_age<0||p.max_age<p.min_age||p.first<0||p.last<p.first||
      p.adjacent_limit<0||p.adjacent_type<0||p.adjacent_type>2||p.min_harvest_gap<0||p.ranking<0||p.ranking>2||
      p.target_fraction<0||p.target_fraction>1||!std::isfinite(p.target_fraction)||
      (p.spread && (!std::isfinite(p.min_block_ha)||!std::isfinite(p.max_block_ha)||
                    p.min_block_ha<0||p.max_block_ha<=0||p.max_block_ha<p.min_block_ha)) ||
      !std::isfinite(p.patch_fraction) || p.patch_fraction<0 || p.patch_fraction>1 ||
      (p.patch_fraction>0 && (!std::isfinite(p.patch_area) || p.patch_area<=0 ||
                              cell_topology.cell_count==0)))
      throw std::runtime_error("invalid E3 prescription");
  }
  const auto stand_size = [&](int stand) {
    return listed_sites.empty() ? so[stand+1]-so[stand] : listed_sites[stand];
  };
  for(int m=0;m<nm;++m) {
    if(mo[m]>=mo[m+1])throw std::runtime_error("empty management area");
    for(int s=mo[m];s<mo[m+1];++s)out.total_sites[m]+=so[s+1]-so[s];
  }
  std::vector<std::vector<int>> prescriptions_by_area(nm);
  for (int j = 0; j < np; ++j)
    if (ps[j].target_fraction > 0 && year >= ps[j].first && year <= ps[j].last)
      prescriptions_by_area[ps[j].ma].push_back(j);
  std::vector<unsigned char> used(ns,0);
  std::vector<std::size_t> seen;
  std::vector<int> frontier, block;
  std::size_t attempt_id=0;
  int event_id = repeat.last_event_id;
  double current_rank = repeat.current_rank;
  const auto update_reserve = [&](int stand, int ma, int duration){
    if (duration <= 0 || out.set_aside_until.empty()) return;
    const double expiry = static_cast<double>(year) + duration;
    for (int k = no[stand]; k < no[stand+1]; ++k) {
      const int neighbor = ni[k];
      if (neighbor >= mo[ma] && neighbor < mo[ma+1] && !used[neighbor])
        out.set_aside_until[neighbor] = std::max(out.set_aside_until[neighbor], expiry);
    }
  };
  const auto cutoff_for = [year](const Prescription& p) {
    return p.min_damage > 0 ? year - p.min_damage : std::numeric_limits<int>::max();
  };
  const auto enqueue = [&](RepeatJob job){
    const std::int64_t next = static_cast<std::int64_t>(year) + ps[job.prescription].repeat_interval;
    if (next > std::numeric_limits<int>::max() ||
        (!repeat.compatible && next > repeat.end_year)) return;
    job.year = static_cast<int>(next);
    out.repeat_jobs.push_back(std::move(job));
  };
  PartialBlock partial_block(cell_topology);
  PatchCut patch_cut(cell_topology, std::any_of(ps.begin(), ps.end(),
      [](const Prescription& p) { return p.patch_fraction > 0; }));
  std::vector<std::size_t> cell_ends;
  std::vector<int> filtered_cells;
  if (std::any_of(ps.begin(), ps.end(), [](const Prescription& p) { return p.spread; }))
    seen.resize(ns,0);
  for(int m=0;m<nm;++m){
    // Scientific Random mode: each active Random prescription receives a
    // uniform random permutation of all stands in this MA, in raster order.
    // This is NOT a claim of event-by-event .NET RNG stream equivalence.
    std::vector<int> raster;
    for (int j : prescriptions_by_area[m]) {
      if (ps[j].ranking != 1) continue;
      if (raster.empty()) {
        raster.resize(mo[m + 1] - mo[m]);
        std::iota(raster.begin(), raster.end(), mo[m]);
        std::sort(raster.begin(), raster.end(), [&](int a, int b) { return first[a] < first[b]; });
      }
      auto permutation=h04c3::landis_shuffle_ranks(static_cast<int>(raster.size()),rng);
      for(std::size_t k=0;k<raster.size();++k)ranks[offsets[j] + (raster[k] - mo[m])]=permutation[k];
    }
    std::vector<int> active;
    std::vector<std::vector<int>> ordered;
    std::vector<std::vector<unsigned char>> spread_eligibility;
    std::vector<double> goals;
    for (int j : prescriptions_by_area[m]) {
      const auto &p=ps[j];
      std::vector<int> candidates;
      for(int s=mo[m];s<mo[m+1];++s){
        const auto ix = offsets[j] + (s - mo[m]);
        if(reserved(s)||!qualifies[ix]||ranks[ix]<=0||!std::isfinite(ranks[ix])||
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
        if(ranks[offsets[j] + (a - mo[m])]!=ranks[offsets[j] + (b - mo[m])])return ranks[offsets[j] + (a - mo[m])]>ranks[offsets[j] + (b - mo[m])];
        return first[a]<first[b];
      });
      out.eligible[j]=candidates.size();
      double goal=out.total_sites[m]*p.target_fraction;
      if(!candidates.empty()&&goal>0){
        active.push_back(j);
        goals.push_back(goal);
        spread_eligibility.emplace_back();
        if(p.spread || p.patch_fraction>0){
          auto &eligible=spread_eligibility.back();
          eligible.resize(mo[m+1]-mo[m],0);
          for(int stand:candidates) eligible[stand-mo[m]]=1;
        }
        ordered.push_back(std::move(candidates));
      }
    }
    if (repeat.enabled) {
      for (const RepeatJob* job : due[m]) {
        const int j = job->prescription;
        const auto& prescription = ps[j];
        const int cutoff = cutoff_for(prescription);
        const bool single = prescription.repeat_kind == 1;
        const double rank = repeat.compatible ? current_rank : job->rank;
        const auto draw = rng.draws;
        int id = job->event_id;
        for (int stand : job->stands) used[stand] = 1;
        for (std::size_t i = 0; i < job->stands.size(); ++i) {
          const int stand = job->stands[i];
          const int size = stand_size(stand);
          if (size == 0) continue;
          const auto start = out.selected_cells.size();
          // A compatible periodic spread calls the original selector while
          // the stand is still set aside. That selector returns an empty block.
          if (single || (prescription.spread && !repeat.compatible)) {
            const auto& cells = job->cells[i];
            for (int cell : cells)
              if (damage.listed(cell) && (damage.compatible || damage.eligible(cell, cutoff)))
                out.selected_cells.push_back(cell);
          } else if (prescription.patch_fraction > 0) {
            patch_cut.select(stand, so[stand], so[stand+1], cell_area_ha,
                             prescription.patch_fraction, prescription.patch_area,
                             prescription.patch_overlap, rng, damage, cutoff);
            out.selected_cells.insert(out.selected_cells.end(), patch_cut.cells.begin(), patch_cut.cells.end());
          } else if (!prescription.spread) {
            // The complete-stand selector is replayed at each periodic return.
            for (int k = so[stand]; k < so[stand+1]; ++k) {
              const int cell = cell_topology.raster_positions[k];
              if (damage.eligible(cell, cutoff)) out.selected_cells.push_back(cell);
            }
          }
          const auto finish = out.selected_cells.size();
          if (!single && i == 0 && finish > start) id = ++event_id;
          bool positive = false;
          for (auto k = start; k < finish; ++k) {
            const int cell = out.selected_cells[k];
            const bool removed = cell_topology.rule_credit(repeat.rule_index[j], cell) > 0;
            positive = positive || removed;
            damage.visit(cell, year, removed);
          }
          // A compatible SingleRepeat bypasses the original selector and
          // therefore leaves the last-harvest date unchanged.
          if (repeat.compatible ? !single && finish > start : positive)
            out.last_harvest[stand] = year;
          if (positive) update_reserve(stand, m, prescription.set_aside);
          out.events.push_back({stand,j,m,age[stand],size,id,rank,draw,
                                std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN(),start,finish,job->iteration});
        }
        if (!single && job->iteration < prescription.repeat_count) {
          RepeatJob next = *job;
          ++next.iteration;
          enqueue(std::move(next));
        }
      }
    }
    std::vector<int> newly_cut;
    std::vector<std::size_t> cursor(active.size(),0);
    std::vector<double> weights(active.size(),0.0);
    while(!active.empty()){
      std::fill(weights.begin(),weights.end(),0.0);
      double total=0.0;
      for(std::size_t k=0;k<active.size();++k){
        const auto &candidates=ordered[k];
        const auto &eligible=spread_eligibility[k];
        while(cursor[k]<candidates.size()){
          const int stand=candidates[cursor[k]];
          if(!used[stand] && !reserved(stand) && (eligible.empty() || eligible[stand-mo[m]])) break;
          ++cursor[k];
        }
        if(cursor[k]<candidates.size())
          weights[k]=std::max(0.0,(goals[k]-out.allocated[active[k]])/goals[k]);
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
      const int j=active[kchosen], initial=ordered[kchosen][cursor[kchosen]];
      const auto &prescription=ps[j];
      const int cutoff = cutoff_for(prescription);
      const bool filtered = !damage.delisted.empty() ||
        (damage.enabled && prescription.min_damage > 0);
      block.clear();
      block.push_back(initial);
      const std::vector<int>* subset=nullptr;
      if(prescription.patch_fraction>0){
        patch_cut.select(initial,so[initial],so[initial+1],cell_area_ha,
                         prescription.patch_fraction,prescription.patch_area,
                         prescription.patch_overlap,rng,damage,cutoff);
        if(patch_cut.cells.empty()){
          spread_eligibility[kchosen][initial-mo[m]]=0;
          continue;
        }
        cell_ends.assign(1,patch_cut.cells.size());
        subset=&patch_cut.cells;
      }else if(prescription.spread){
        if (filtered && !prescription.partial) {
          filtered_cells.clear();
          subset = &filtered_cells;
        }
        auto &eligible=spread_eligibility[kchosen];
        const auto lower_priority=[&](int a,int b){
          if(ranks[offsets[j] + (a - mo[m])]!=ranks[offsets[j] + (b - mo[m])]) return ranks[offsets[j] + (a - mo[m])]<ranks[offsets[j] + (b - mo[m])];
          return first[a]>first[b];
        };
        // A block grows through the highest ranked neighbor of any
        // member. Trial stands are reserved only after the minimum is met.
        ++attempt_id;
        seen[initial]=attempt_id;
        frontier.clear();
        int block_sites=0;
        cell_ends.clear();
        if(prescription.partial){
          partial_block.reset();
          subset=&partial_block.cells;
        }
        for(std::size_t i=0;i<block.size();++i){
          const int stand=block[i];
          if(prescription.partial){
            const int entry=i==0 ? partial_block.start(so[stand],so[stand+1],rng,damage,cutoff)
                                 : partial_block.entry(so[stand],so[stand+1],damage,cutoff);
            if(entry<0){
              eligible[stand-mo[m]] = 0;
              block.pop_back();
              if (filtered && !damage.compatible && !frontier.empty()) {
                std::pop_heap(frontier.begin(),frontier.end(),lower_priority);
                block.push_back(frontier.back());
                frontier.pop_back();
                --i; // Retry this slot with the next accessible branch.
                continue;
              }
              break;
            }
            partial_block.append(entry,stand,cell_area_ha,prescription.max_block_ha,damage,cutoff);
            cell_ends.push_back(partial_block.cells.size());
            block_sites=static_cast<int>(partial_block.cells.size());
          } else if (filtered) {
            for (int k = so[stand]; k < so[stand+1]; ++k) {
              const int cell = cell_topology.raster_positions[k];
              if (damage.eligible(cell, cutoff)) filtered_cells.push_back(cell);
            }
            cell_ends.push_back(filtered_cells.size());
            block_sites = static_cast<int>(filtered_cells.size());
          } else {
            block_sites+=so[stand+1]-so[stand];
          }
          const bool empty_member = filtered && cell_ends.back() ==
            (i == 0 ? 0 : cell_ends[i-1]);
          if (empty_member && !damage.compatible) {
            eligible[stand-mo[m]] = 0;
            if (i == 0) break;
          }
          if(block_sites*cell_area_ha>=prescription.max_block_ha) break;
          for(int k=no[stand]; (!empty_member || damage.compatible) && k<no[stand+1]; ++k){
            const int neighbor=ni[k];
            if(neighbor<mo[m]||neighbor>=mo[m+1]||used[neighbor]||reserved(neighbor)||
               !eligible[neighbor-mo[m]]||seen[neighbor]==attempt_id) continue;
            seen[neighbor]=attempt_id;
            frontier.push_back(neighbor);
            std::push_heap(frontier.begin(),frontier.end(),lower_priority);
          }
          if(frontier.empty()) break;
          std::pop_heap(frontier.begin(),frontier.end(),lower_priority);
          block.push_back(frontier.back());
          frontier.pop_back();
        }
        if(block_sites == 0 || block_sites*cell_area_ha<prescription.min_block_ha){
          // A failed block stays available to other prescriptions in the MA.
          for(int stand:block) eligible[stand-mo[m]]=0;
          continue;
        }
      } else if (filtered) {
        filtered_cells.clear();
        for (int k = so[initial]; k < so[initial+1]; ++k) {
          const int cell = cell_topology.raster_positions[k];
          if (damage.eligible(cell, cutoff)) filtered_cells.push_back(cell);
        }
        if (filtered_cells.empty()) {
          ++cursor[kchosen];
          continue;
        }
        cell_ends.assign(1, filtered_cells.size());
        subset = &filtered_cells;
      }
      ++event_id;
      current_rank = ranks[offsets[j] + (initial - mo[m])];
      RepeatJob scheduled{j, year, 1, event_id, current_rank, {}, {}};
      bool has_cut = false;
      // Every member in a committed block is already selected when the
      // buffer is written, so no member is reserved by its own block.
      for (std::size_t i = 0; i < block.size(); ++i)
        if (!subset || cell_ends[i] != (i == 0 ? 0 : cell_ends[i-1])) used[block[i]] = 1;
      for(std::size_t i=0;i<block.size();++i){
        const int stand=block[i];
        const int stand_sites=stand_size(stand);
        const auto cell_begin=out.selected_cells.size();
        int attempted=stand_sites;
        int credit=effective_sites.empty() ? stand_sites : effective_sites[offsets[j] + (stand - mo[m])];
        if(subset){
          const auto start=i==0 ? 0 : cell_ends[i-1];
          const auto end=cell_ends[i];
          attempted=static_cast<int>(end-start);
          credit=0;
          for(auto k=start;k<end;++k){
            const int cell=(*subset)[k];
            out.selected_cells.push_back(cell);
            credit+=cell_topology.credit(j,cell);
          }
        }
        if (attempted == 0) {
          if (!spread_eligibility[kchosen].empty()) spread_eligibility[kchosen][stand-mo[m]] = 0;
          continue;
        }
        newly_cut.push_back(stand);
        if (damage.enabled) {
          for (int k = 0; k < attempted; ++k) {
            const int cell = subset ? out.selected_cells[cell_begin+k]
                                    : cell_topology.raster_positions[so[stand]+k];
            damage.visit(cell, year, cell_topology.credit(j, cell) > 0);
          }
        }
        if (credit > 0) update_reserve(stand, m, prescription.set_aside);
        if (prescription.repeat_kind > 0 &&
            (credit > 0 || (repeat.compatible && prescription.spread))) {
          has_cut = has_cut || credit > 0;
          scheduled.stands.push_back(stand);
          scheduled.cells.emplace_back();
          // Preserve the winning cells of each member of the block. The
          // compatible single return associates its sites with the entry stand.
          const bool record_cells = prescription.repeat_kind == 1
            ? (!repeat.compatible || i == 0)
            : (prescription.spread && !repeat.compatible);
          if (record_cells) {
            auto& saved = scheduled.cells.back();
            for (int k = 0; k < attempted; ++k) {
              const int cell = subset ? out.selected_cells[cell_begin+k]
                                      : cell_topology.raster_positions[so[stand]+k];
              if (prescription.repeat_kind == 2 || cell_topology.credit(j, cell)) saved.push_back(cell);
            }
          }
        }
        out.attempted[j]+=attempted;
        out.allocated[j]+=credit;
        out.events.push_back({stand,j,m,age[stand],stand_sites,event_id,ranks[offsets[j] + (stand - mo[m])],
                              draw,u,weights[kchosen]/total,cell_begin,out.selected_cells.size()});
      }
      if (has_cut) enqueue(std::move(scheduled));
    }
    // H04-D2 observed semantics: within-MA candidate snapshot, then
    // only subsequent management areas observe current-year history.
    for(int s:newly_cut)out.last_harvest[s]=year;
  }
  out.offset_end=rng.draws;
  out.damage_year = std::move(damage.harvest);
  out.last_event_id = event_id;
  out.current_rank = current_rank;
  return out;
}
} // namespace
#endif
