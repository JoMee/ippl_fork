#pragma once
#include "Mesh_new/structured/GrassmanIndex.hpp"

namespace fem {
namespace Detail {
template<typename From, typename To, typename Indexer>
struct StructuredIncidenceProvider;

// Specialization for finding the 2 Vertices of an x-aligned Edge in 2D
template<typename Indexer>
struct StructuredIncidenceProvider<Blade<0>, Blade<>, Indexer> {
  // Returns the logical coordinates of the 2 vertices for a given edge's logical coordinate.
  KOKKOS_INLINE_FUNCTION auto get(const Kokkos::Array<int, Indexer::coord_type::size()>& edge_coord) const
  -> Kokkos::Array<Kokkos::Array<int, Indexer::coord_type::size()>, 2>
  {
    return {edge_coord, Kokkos::Array<int, 2>{edge_coord[0] + 1, edge_coord[1]}};
  }
};
// Specialization for finding the 2 Vertices of a y-aligned Edge in 2D
template<typename Indexer>
struct StructuredIncidenceProvider<Blade<1>, Blade<>, Indexer> {
  // Returns the logical coordinates of the 2 vertices for a given edge's logical coordinate.
  KOKKOS_INLINE_FUNCTION auto get(const Kokkos::Array<int, Indexer::coord_type::size()>& edge_coord) const
  -> Kokkos::Array<Kokkos::Array<int, Indexer::coord_type::size()>, 2>
  {
    return {edge_coord, Kokkos::Array<int, 2>{edge_coord[0], edge_coord[1] + 1}};
  }
};
}
template <typename From, typename To, int Dim>
struct Connectivity {
  static auto get_provider(const GrassmanIndex<Dim>& /*indexer*/) {
    return Detail::StructuredIncidenceProvider<From, To, GrassmanIndex<Dim>>{};
  }
};
}
