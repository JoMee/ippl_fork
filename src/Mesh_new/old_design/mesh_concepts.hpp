#pragma once
#include <concepts>
#include "incidence.hpp"

namespace fem {

template<class M>
concept BasicMeshModel = requires(const M& mesh,
                                  typename M::global_index g,
                                  int source_dim,
                                  int target_dim) {
    typename M::global_index;

    { M::topological_dimension } -> std::convertible_to<int>;
    { mesh.incidence(g, source_dim, target_dim) }
        -> std::same_as<IncidenceView<typename M::global_index>>;
};

} // namespace fem 
