#pragma once

#include "Mesh_new/MeshPolicies.hpp"

namespace fem {

// Specialization for a single-rank layout on a structured mesh
template <int Dim, typename MeshPolicy>
class Layout {
public:

    static constexpr int DIM = Dim;
    using Policy = MeshPolicy;
    using MeshType = typename Policy::template MeshType<Dim>;
    using IndexerType = typename Policy::template IndexerType<Dim>;

    explicit Layout(std::shared_ptr<MeshType> mesh, int halo_width = 0)
        : mesh_(mesh),
          indexer_(std::make_shared<IndexerType>(mesh->get_extents(), halo_width))
    {
        std::cout << "Generic Layout created using policy: "
                  << typeid(MeshPolicy).name() << std::endl;
    }

    template <typename GroupTag>
    auto get_alloc_extent() const {
        return indexer_->template get_allocated_extent<GroupTag>();
    }

    /**
     * @brief Gets a connectivity provider between two element groups.
     *
     * @tparam FromGroupTag The compile-time tag for the source element group.
     * @tparam ToGroupTag The compile-time tag for the target element group.
     * @return An instance of the correct connectivity provider for this policy.
     */
    template <typename FromGroupTag, typename ToGroupTag>
    auto get_incidence() const {
        // The policy provides the specific provider type based on the generic tags.
        using ProviderType = typename Policy::template ConnectivityProvider<
            FromGroupTag,
            ToGroupTag,
            IndexerType
        >;

        return ProviderType{};
    }

    template <typename GroupTag>
    auto get_geometry() const {
        return Policy::template create_geometry_provider<GroupTag, Dim>(*mesh_);
    }

    void fillHalo() const {
        Policy::LayoutStrategy::fillHalo(*this);
    }

    const auto& get_mesh() const { return *mesh_; }
    const auto& get_indexer() const { return *indexer_; }

private:
    std::shared_ptr<MeshType> mesh_;
    std::shared_ptr<const IndexerType> indexer_;
};

} // namespace fem
