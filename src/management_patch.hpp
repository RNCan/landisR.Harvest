#ifndef LANDISR_MANAGEMENT_PATCH_HPP
#define LANDISR_MANAGEMENT_PATCH_HPP

#include "management_partial.hpp"

namespace h04e4d {

// The Fenwick tree preserves LANDIS-II's candidate-list order with
// indexed selection and deletion in logarithmic time. Pending cells
// are appended in queue order when the current patch ends.
class PatchCandidates {
  std::vector<int> slots;
  std::vector<int> tree;
  std::vector<std::size_t> position;
  int remaining = 0;

  static std::size_t low_bit(std::size_t i) { return i & (-i); }

  int prefix(std::size_t i) const {
    int count = 0;
    for (; i > 0; i -= low_bit(i)) count += tree[i];
    return count;
  }

public:
  explicit PatchCandidates(std::size_t cell_count) : position(cell_count, 0) {}

  int size() const { return remaining; }
  void clear_position(int cell) { position[cell] = 0; }

  void reset(const int* cells, int count) {
    slots.resize(count + 1);
    tree.resize(count + 1);
    remaining = count;
    for (int j = 0; j < count; ++j) {
      const std::size_t i = static_cast<std::size_t>(j) + 1;
      slots[i] = cells[j];
      position[cells[j]] = i;
      tree[i] = static_cast<int>(low_bit(i));
    }
  }

  void append(int cell) {
    const auto i = slots.size();
    slots.push_back(cell);
    position[cell] = i;
    tree.push_back(1 + prefix(i - 1) - prefix(i - low_bit(i)));
    ++remaining;
  }

  bool remove(int cell) {
    const auto index = position[cell];
    if (index == 0) return false;
    position[cell] = 0;
    --remaining;
    for (auto i = index; i < tree.size(); i += low_bit(i)) --tree[i];
    return true;
  }

  int draw(h04c3::Stream& rng, bool compatible) {
    const int choices = compatible && size() > 1 ? size() - 1 : size();
    int target = rng.next(choices) + 1;
    std::size_t index = 0;
    std::size_t bit = 1;
    while (bit <= (tree.size() - 1) / 2) bit *= 2;
    for (; bit > 0; bit /= 2) {
      const auto next = index + bit;
      if (next < tree.size() && tree[next] < target) {
        index = next;
        target -= tree[index];
      }
    }
    const int cell = slots[index + 1];
    remove(cell);
    return cell;
  }
};

class PatchCut {
  const PartialTopology& topology;
  PatchCandidates available;
  std::vector<int> pending, eligible;

public:
  std::vector<int> cells;

  PatchCut(const PartialTopology& topology, bool enabled)
    : topology(topology), available(enabled ? topology.cell_count : 0) {}

  void select(int stand, int begin, int end, double cell_area,
              double fraction, double patch_area, bool overlap, h04c3::Stream& rng,
              const DamageHistory& damage, int cutoff) {
    cells.clear();
    int listed_count = end - begin;
    const bool filter_damage = damage.enabled && !damage.compatible &&
      cutoff != std::numeric_limits<int>::max();
    if (filter_damage || !damage.delisted.empty()) {
      eligible.clear();
      for (int i = begin; i < end; ++i) {
        const int cell = topology.raster_positions[i];
        if (!damage.listed(cell)) {
          --listed_count;
          available.clear_position(cell);
        } else if (filter_damage && !damage.eligible(cell, cutoff)) {
          available.clear_position(cell);
        } else {
          eligible.push_back(cell);
        }
      }
      if (eligible.empty()) return;
      available.reset(eligible.data(), static_cast<int>(eligible.size()));
    } else {
      available.reset(topology.raster_positions + begin, end - begin);
    }
    const double target = listed_count * cell_area * fraction;
    int start = available.draw(rng, topology.compatible);
    while ((!topology.compatible || available.size() > 0) && cells.size() * cell_area < target) {
      const auto patch_begin = cells.size();
      pending.clear();
      pending.push_back(start);
      std::size_t cursor = 0;
      for (; cursor < pending.size();) {
        const int cell = pending[cursor++];
        const auto* neighbors = topology.neighbors + static_cast<std::size_t>(cell) * cell_directions;
        for (int d = 0; d < cell_directions; ++d) {
          const int next = neighbors[d];
          if (next >= 0 && topology.stand_of_cell[next] == stand && available.remove(next))
            pending.push_back(next);
        }
        if (damage.eligible(cell, cutoff)) cells.push_back(cell);
        if ((cells.size() - patch_begin) * cell_area >= patch_area ||
            cells.size() * cell_area >= target ||
            (topology.compatible && available.size() == 0)) break;
      }

      for (; cursor < pending.size(); ++cursor) available.append(pending[cursor]);
      if (!overlap) {
        for (auto i = patch_begin; i < cells.size(); ++i) {
          const auto* neighbors = topology.neighbors
              + static_cast<std::size_t>(cells[i]) * cell_directions;
          for (int d = 0; d < cell_directions; ++d) {
            const int neighbor = neighbors[d];
            if (neighbor >= 0 && topology.stand_of_cell[neighbor] == stand)
              available.remove(neighbor);
          }
        }
      }
      // The v7 selector draws the next start after closing the patch,
      // even when the stand target has just been reached.
      if (available.size() == 0 || (!topology.compatible && cells.size() * cell_area >= target)) break;
      start = available.draw(rng, topology.compatible);
    }
  }
};
}
#endif
