#include <concepts>
#include "Mesh_new/index_space.h"

namespace fem {

template<typename T>
concept Mesh = requires(const T& mesh, int entity_dim, int from_dim, int to_dim) {
    // The main index space for a given entity dimension
    { mesh.get_index_space(entity_dim) } -> LocalIndexSpace;

    // Get the incidence between entity sets (incidence)
    // The returned object must have a `get_neighbors(local_id)` method
    // For this, we require from_dim > to_dim. 
    { mesh.get_incidence(from_dim, to_dim) -> Connectivity};

    // Returns the number of entities owned by this mesh (without ghost entities)
    { mesh.entity_count(entity_dim) } -> std::integral;
};

template<typename T>
concept Connectivity = requires(const T& conn, int local_id) {
    { conn.get_neighbors(local_id) } -> std::ranges::input_range;
    { conn.neighbor_count(local_id) } -> std::integral;
};

}

/*

class StructuredConnectivity {
public:
    // Store grid dimensions, etc...

    struct NeighborView {
        Kokkos::View<const local_id_type*, DeviceType> entries;
        int begin;
        int end;

        KOKKOS_INLINE_FUNCTION int size() const { return end - begin; }
        KOKKOS_INLINE_FUNCTION local_id_type operator[](int i) const { return entries(begin + i); }
    };

    KOKKOS_INLINE_FUNCTION
    NeighborView get_neighbors(local_id_type lid) const {
        // No memory lookup!
        // Convert lid to logical index.
        // Calculate neighbor of logical index, etc., and their lids.
        // Return the result in the NeighborView struct.
    }
};
*/
