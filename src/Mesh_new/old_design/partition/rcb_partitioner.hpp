#pragma once

#include <array>
#include <vector>
#include "Mesh_new/structured/index_map.hpp"

namespace fem {

/// Represents an axis-aligned box in integer index space
template<int Dim>
struct Box {
    std::array<int, Dim> lower;
    std::array<int, Dim> upper;
};

template<int Dim> using Index = std::array<int, Dim>;

/// Recursive Coordinate Bisection partitioner
template<int Dim>
class RCBPartitioner {
public:
    /// Construct with number of ghost layers (default = 1)
    explicit constexpr RCBPartitioner(int ghost_layers = 1) noexcept
        : ghost_layers{ghost_layers} {}

    /// Partition a box into `num_parts` regions
    /// Optionally control which dimensions are splitable (default: all true)
    std::vector<StructuredIndexMap<Dim>>
    partition(const Index<Dim>& global_lower, 
              const Index<Dim>& global_upper,
              int num_partitions, 
              std::span<const bool, Dim> splitable) const
    {
        static_assert(Dim > 0, "Dim must be positive");
        Index<Dim> global_shape = diff(global_upper, global_lower);

        std::vector<Box<Dim>> boxes;
        split_recursive({global_lower, global_upper}, num_partitions, splitable, boxes);

        std::vector<StructuredIndexMap<Dim>> maps;
        maps.reserve(num_partitions);
        for (const auto& b : boxes)
            maps.emplace_back(global_shape, b.lower,
                              diff(b.upper, b.lower),
                              ghost_layers);
        return maps;
    }

    std::vector<StructuredIndexMap<Dim>> partition(
        const Index<Dim>& gl,
        const Index<Dim>& gu,
        int parts) const
    {
        static constexpr auto all = make_true_array();
        return partition(gl, gu, parts, std::span<const bool, Dim>{all});
    }

private:
    int ghost_layers;

    static constexpr Index<Dim> diff(const Index<Dim>& a,
                                   const Index<Dim>& b) {
      Index<Dim> r{};
      for (int d = 0; d < Dim; ++d) r[d] = a[d] - b[d];
      return r;
    }

    static void split_recursive(Box<Dim> box, int parts,
                                std::span<const bool, Dim> mask,
                                std::vector<Box<Dim>>& out) {
      if (parts == 1) { out.push_back(box); return; }

      int axis = longest_axis(box, mask);
      if (axis < 0) throw std::runtime_error("No splittable axis");

      int mid = (box.lower[axis] + box.upper[axis]) / 2;
      Box<Dim> left{box}, right{box};
      left.upper[axis] = mid;
      right.lower[axis] = mid;

      int l = parts / 2;
      split_recursive(left,  l, mask, out);
      split_recursive(right, parts - l, mask, out);
    }

    static int longest_axis(const Box<Dim>& b,
                            std::span<const bool, Dim> mask) {
        int best = -1, len = -1;
        for (int d = 0; d < Dim; ++d)
            if (mask[d] && (b.upper[d] - b.lower[d] > len))
                best = d, len = b.upper[d] - b.lower[d];
        return best;
    }

    static constexpr std::array<bool, Dim> make_true_array() {
        std::array<bool, Dim> out{};
        for (int d = 0; d < Dim; ++d)
            out[d] = true;
        return out;
    }
};

} // namespace fem
