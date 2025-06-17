#pragma once

#include "Mesh_new/GrassmanIndex.hpp"
#include "Mesh_new/Connectivity.hpp"

namespace fem {

template <int Dim>
class StructuredCartesianMesh {
public:
    explicit StructuredCartesianMesh(const Kokkos::Array<int, Dim>& extents)
        : logical_extents_(extents) {}

    auto get_extents() const { return logical_extents_; }

    // This method now requires the caller (the Layout) to provide the indexer.
    // This removes the runtime branch and makes the Mesh stateless.
    template <typename FromBlade, typename ToBlade>
    auto get_incidence(const GrassmanIndex<Dim>& indexer) const {
        return Connectivity<FromBlade, ToBlade, Dim>::get_functor(indexer);
    }

private:
    Kokkos::Array<int, Dim> logical_extents_;
};

} // namespace fem

