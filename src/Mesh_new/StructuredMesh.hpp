#pragma once

#include "Mesh_new/GrassmanIndex.hpp"
#include "Mesh_new/StructuredConnectivity.hpp"

namespace fem {

template <int N>
class StructuredMesh {
public:
  using coord_type = Kokkos::Array<int, N>;

  StructuredMesh(const coord_type& global_shape,
                 const coord_type& offset,
                 const coord_type& local_extent)
    : index(global_shape, offset, local_extent)
  {}

  template <typename Blade>
  auto get_index_space() const {
    return index.template get_subspace<Blade>();
  }

  template <typename FromBlade, typename ToBlade>
  auto get_incidence() const {
    return connectivity.template get<FromBlade, ToBlade>();
  }

  template <typename Blade>
  std::size_t entity_count() const {
    auto view = get_index_space<Blade>();
    std::size_t count = 1;
    for (int i = 0; i < N; ++i)
      count *= view.extent[i];
    return count;
  }

private:
  GrassmanIndex<N> index;
  StructuredConnectivity<N> connectivity; 
};

} // namespace fem

