#ifndef LANDISR_MANAGEMENT_PARTIAL_HPP
#define LANDISR_MANAGEMENT_PARTIAL_HPP

#include <cstddef>
#include <limits>
#include <vector>
#include <unordered_set>
#include "management_c3_math.hpp"

namespace h04e4d {
constexpr int cell_directions = 8;

// Cached active-cell indices. Offsets for raster_positions are the
// same as for the original stand topology.
struct PartialTopology {
  const int* raster_positions = nullptr;
  const int* neighbors = nullptr;
  const int* stand_of_cell = nullptr;
  const int* cell_states = nullptr;
  const int* rule_index = nullptr;
  const signed char* removal = nullptr;
  std::size_t state_count = 0;
  std::size_t cell_count = 0;
  bool compatible = false;

  int rule_credit(int rule, int cell) const {
    if (!removal) return 1;
    return removal[static_cast<std::size_t>(rule) * state_count + cell_states[cell] - 1];
  }

  int credit(int prescription, int cell) const {
    return removal ? rule_credit(rule_index[prescription], cell) : 1;
  }
};

struct DamageHistory {
  std::vector<int> harvest, external;
  std::unordered_set<int> delisted;

  bool listed(int cell) const {
    return delisted.empty() || delisted.count(cell) == 0;
  }
  bool enabled = false, compatible = false;

  bool eligible(int cell, int cutoff) const {
    return listed(cell) &&
      (!enabled || cutoff == std::numeric_limits<int>::max() ||
       (harvest[cell] <= cutoff && external[cell] <= cutoff));
  }

  void visit(int cell, int year, bool positive) {
    if (enabled && (compatible || positive)) harvest[cell] = year;
  }
};

class PartialBlock {
  const PartialTopology& topology;
  std::vector<std::size_t> visited;
  std::vector<int> queue;
  std::size_t attempt = 0;

public:
  std::vector<int> cells;

  explicit PartialBlock(const PartialTopology& topology)
    : topology(topology), visited(topology.cell_count, 0) {}

  void reset() {
    ++attempt;
    cells.clear();
  }

  int entry(int begin, int end, const DamageHistory& damage, int cutoff) const {
    for (int i = begin; i < end; ++i) {
      const int cell = topology.raster_positions[i];
      if (!damage.listed(cell) || (!damage.compatible && !damage.eligible(cell, cutoff))) continue;
      const int* neighbors = topology.neighbors + static_cast<std::size_t>(cell) * cell_directions;
      for (int d = 0; d < cell_directions; ++d) {
        if (neighbors[d] >= 0 && visited[neighbors[d]] == attempt) return cell;
      }
    }
    return -1;
  }

  int start(int begin, int end, h04c3::Stream& rng, const DamageHistory& damage, int cutoff) {
    const bool filter_damage = damage.enabled && !damage.compatible &&
      cutoff != std::numeric_limits<int>::max();
    if (filter_damage || !damage.delisted.empty()) {
      queue.clear();
      for (int i = begin; i < end; ++i) {
        const int cell = topology.raster_positions[i];
        if (damage.listed(cell) && (!filter_damage || damage.eligible(cell, cutoff)))
          queue.push_back(cell);
      }
      const int count = static_cast<int>(queue.size());
      const int choices = topology.compatible && count > 1 ? count - 1 : count;
      return count == 0 ? -1 : queue[rng.next(choices)];
    }
    const int count = end - begin;
    // LANDIS-II v7 samples floor(u * (count - 1)).
    const int choices = topology.compatible && count > 1 ? count - 1 : count;
    return topology.raster_positions[begin + rng.next(choices)];
  }

  void append(int first_cell, int stand, double cell_area, double max_area,
              const DamageHistory& damage, int cutoff) {
    queue.clear();
    queue.push_back(first_cell);
    visited[first_cell] = attempt;
    for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
      const int cell = queue[cursor];
      if (damage.eligible(cell, cutoff)) cells.push_back(cell);
      // Any queued cells left at this point will not be used for entry
      // into another stand, because the block has reached its target.
      if (cells.size() * cell_area >= max_area) break;
      const int* neighbors = topology.neighbors + static_cast<std::size_t>(cell) * cell_directions;
      for (int d = 0; d < cell_directions; ++d) {
        const int next = neighbors[d];
        if (next >= 0 && topology.stand_of_cell[next] == stand && visited[next] != attempt &&
            damage.listed(next) && (damage.compatible || damage.eligible(next, cutoff))) {
          visited[next] = attempt;
          queue.push_back(next);
        }
      }
    }
  }
};
}
#endif
