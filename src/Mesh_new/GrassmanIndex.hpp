#pragma once

#include <Kokkos_Core.hpp>
#include <unordered_map>
#include <tuple>
#include "Mesh_new/Blades.hpp" // contains Blade<D...> and GenerateBlades<N>
#include "Mesh_new/IndexView.hpp" // helper to wrap coordinate bounds per Blade

namespace fem {

template <int N>
class GrassmanIndex {
public:
  using coord_type = Kokkos::Array<int, N>;

  GrassmanIndex() = default;

  GrassmanIndex(const coord_type& global_shape,
                const coord_type& offset,
                const coord_type& local_extent)
    : global_shape(global_shape),
      vertex_offset(offset),
      vertex_extent(local_extent)
  {
    // Build index ranges for all possible blades
    build_all_index_spaces();
  }

  template <typename Blade>
  auto get_subspace() const -> IndexView<N> {
    constexpr uint32_t mask = Blade::Mask;
    auto it = index_spaces.find(mask);
    if (it == index_spaces.end())
      throw std::runtime_error("No index space for given Blade.");
    return it->second;
  }

private:
  coord_type global_shape, vertex_offset, vertex_extent;

  std::unordered_map<uint32_t, IndexView<N>> index_spaces;

  void build_all_index_spaces() {
    using AllBlades = typename Detail::GenerateBlades<N>::type;
    [&]<typename... Blades>(std::tuple<Blades...>) {
      (insert_subspace<Blades>(), ...);
    }(AllBlades{});
  }

  template <typename Blade>
  void insert_subspace() {
    constexpr uint32_t mask = Blade::Mask;

    // Offset logic: this is oversimplified for now
    IndexView<N> view;
    for (int i = 0; i < N; ++i) {
      view.offset[i] = vertex_offset[i];
      view.extent[i] = vertex_extent[i] - ((mask >> i) & 1); 
    }
    index_spaces[mask] = view;
  }

};

} // namespace fem
