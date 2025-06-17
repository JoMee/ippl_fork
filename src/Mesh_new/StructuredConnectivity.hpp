#include <Kokkos_Core.hpp>

namespace fem {

template <typename FromBlade, typename ToBlade, int N>
std::vector<std::array<int, N>> generate_offsets() {
  // Dummy version: just returns a single offset [0, 0, ..., 0]
  // Replace later with oriented subcube offsets derived from mask differences.
  return { std::array<int, N>{} };
}

template<int N>
class StructuredConnectivity {
public:
  using index_type = Kokkos::Array<int, N>;

  // Compact object for a given (FromBlade, ToBlade) mapping
  template<typename FromBlade, typename ToBlade>
  struct IncidenceView {
    using from_type = FromBlade;
    using to_type = ToBlade;

    // Neighbor offsets from an entity of FromBlade type to ToBlade neighbors
    static constexpr auto neighbor_offsets = generate_offsets<FromBlade, ToBlade, N>();

    KOKKOS_INLINE_FUNCTION
    int num_neighbors() const {
      return neighbor_offsets.size();
    }

    // Query the neighbor at `neighbor_idx` for a given local_id (entity coord)
    KOKKOS_INLINE_FUNCTION
    index_type neighbor(const index_type& id, int neighbor_idx) const {
      index_type out = id;
      for (int d = 0; d < N; ++d)
        out[d] += neighbor_offsets[neighbor_idx][d];
      return out;
    }
  };

  // Retrieve a view for incidence between FromBlade → ToBlade
  template<typename FromBlade, typename ToBlade>
  KOKKOS_INLINE_FUNCTION
  IncidenceView<FromBlade, ToBlade> get_incidence() const {
    return {};
  }

};


}
