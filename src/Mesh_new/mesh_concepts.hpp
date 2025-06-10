#pragma once
#include <concepts>
#include "incidence.hpp"

namespace fem {

// ---- Tier-0: topology only ---------------------------------------------------------------
template<class M>
concept BasicMeshModel = requires(const M& mesh,
                                  typename M::global_index g,
                                  int src_dim, int tgt_dim)
{
    typename M::global_index;
    { M::topological_dimension } -> std::convertible_to<int>;

    { mesh.entity_dimension(g) } -> std::same_as<int>;
    { mesh.incidence(g, src_dim, tgt_dim) }
        -> std::same_as<IncidenceView<typename M::global_index>>;
};

// ---- Tier-1: geometry available ----------------------------------------------------------
template<class M>
concept MetricMeshModel = BasicMeshModel<M> &&
    requires(const M& mesh, typename M::global_index g)
{
    typename M::vec_type;
    typename M::scalar_type;
    { mesh.centroid(g) } noexcept -> std::same_as<typename M::vec_type>;
    { mesh.measure(g) }  noexcept -> std::same_as<typename M::scalar_type>;
};

} // namespace fem 
