#include <concepts>
#include "Mesh_new/index_space.h"
#include "Mesh_new/structured/Blades.hpp"

namespace fem {

template<typename T>
concept Mesh = requires(const T& mesh,
                        typename T::global_index g,
                        typename T::local_index l) {
  // Metadata
  { T::topological_dimension } -> std::convertible_to<int>;

  // The index space associated with a given Blade or entity tag
  { mesh.template get_index_space<Blade<0>>() } -> LocalIndexSpace;

  // Entity count for that entity type (without ghosts)
  { mesh.template entity_count<Blade<0>>() } -> std::integral;

  // Incidence between two types
  { mesh.template get_incidence<Blade<0,1>, Blade<0>>() }; 
    -> IncidenceView<typename M::local_index>;
};

template<typename I, typename Index>
concept IncidenceView = requires(const I& incidence,
                                 const Index& from_idx,
                                 int i) {

    // number of neighbors for a given entity
    { incidence.degree(from_idx) } -> std::integral;

    // Returns (sign, neighbor index)
    { incidence.neighbor(from_idx, i) } 
        -> std::same_as<std::pair<int, Index>>;

};

}
