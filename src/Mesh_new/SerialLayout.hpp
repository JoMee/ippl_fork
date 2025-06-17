#pragma once

#include "Mesh_new/StructuredMesh.hpp"

namespace fem {

template <typename MeshType>
class Layout;

// Specialization for a single-rank layout on a structured mesh
template <int Dim>
class Layout<StructuredCartesianMesh<Dim>> {
public:
    static constexpr int DIM = Dim;
    using MeshType = StructuredCartesianMesh<Dim>;

    explicit Layout(std::shared_ptr<MeshType> mesh, int halo_width = 0)
        : mesh_(mesh),
          indexer_(std::make_shared<GrassmanIndex<Dim>>(mesh->get_extents(), halo_width))
    {
        std::cout << "Layout created. It owns the single indexer instance." << std::endl;
    }

    template <typename Blade>
    auto get_alloc_extent() const {
        return indexer_->template get_allocated_extent<Blade>();
    }

    // A helper for operators to get a connectivity functor.
    // It calls the mesh's method, passing its own indexer.
    template <typename FromBlade, typename ToBlade>
    auto get_incidence() const {
        return mesh_->template get_incidence<FromBlade, ToBlade>(*indexer_);
    }

    // Trivial halo exchange for the playground
    void fill_halo() const {
        std::cout << "  -> (Layout) Performing trivial halo exchange." << std::endl;
    }

    const auto& get_mesh() const { return mesh_; }
    const auto& get_indexer() const { return *indexer_; }

private:
    std::shared_ptr<MeshType> mesh_;
    std::shared_ptr<const GrassmanIndex<Dim>> indexer_;
};

} // namespace fem
