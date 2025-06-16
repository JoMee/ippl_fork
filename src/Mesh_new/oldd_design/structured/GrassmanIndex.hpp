#pragma once

#include "Mesh_new/index_space.h"
#include <Kokkos_Core.hpp>
#include "Mesh_new/structured/Blades.hpp"
#include "Mesh_new/structured/BladeDispatch.hpp"

namespace fem {

template<int N>
class GrassmanIndex {

public:
  using index_type = Kokkos::Array<int, N + 1>;
  using coord_type = Kokkos::Array<int, N>;

  using local_id_type = index_type;
  using global_id_type = index_type;

  static constexpr int mask(const index_type& idx) {
    return idx[N];
  }

  GrassmanIndex(const coord_type& global_shape,
                const coord_type& offset,
                const coord_type& local_extent)
    : vertex_shape(global_shape), 
      vertex_offset(offset),
      vertex_extent(local_extent)
    {}

    

private:
  coord_type vertex_shape;
  coord_type vertex_offset;
  coord_type vertex_extent;

};

}
