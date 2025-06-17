#pragma once

#include <Kokkos_Core.hpp>
#include <unordered_map>
#include <tuple>
#include "Mesh_new/Blades.hpp"
#include "Mesh_new/IndexView.hpp"

namespace fem {

template <int N>
class GrassmanIndex {
public:
  using coord_type = Kokkos::Array<int, N>;

  explicit GrassmanIndex(const coord_type& logical_extent, int halo_width = 0)
        : logical_vertex_extent_(logical_extent), halo_width_(halo_width)
    {}


  // --- API for the LOGICAL VIEW (used by Mesh/Connectivity) ---
  template <typename Blade>
  auto get_logical_extent() const -> coord_type {
    coord_type extent;
    for (int i = 0; i < N; ++i) {
        extent[i] = logical_vertex_extent_[i] - ((Blade::Mask >> i) & 1);
    }
    return extent;
  }

  // --- API for the STORAGE VIEW (used by Layout) ---
  template <typename Blade>
  auto get_allocated_extent() const -> coord_type {
      auto extent = get_logical_extent<Blade>();
      for (int i = 0; i < N; ++i) {
          extent[i] += 2 * halo_width_;
      }
      return extent;
  }


  KOKKOS_INLINE_FUNCTION auto get_storage_offset() const -> coord_type {
      coord_type offset;
      for (int i = 0; i < N; ++i) offset[i] = halo_width_;
      return offset;
  }

private:
  coord_type logical_vertex_extent_;
  int halo_width_;
};

} // namespace fem
