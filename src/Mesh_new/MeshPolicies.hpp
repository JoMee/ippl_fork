#pragma once

// Include the existing components that this policy will use.
#include "Mesh_new/structured/StructuredMesh.hpp"
#include "Mesh_new/structured/GrassmanIndex.hpp"
#include "Mesh_new/structured/Connectivity.hpp"
#include "Mesh_new/structured/Geometry.hpp"
#include "Mesh_new/Blades.hpp"
#include "Mesh_new/LayoutStrategy.hpp"

namespace fem {

/**
 * @brief A policy that bundles all types and services for a structured Cartesian mesh.
 *
 * This struct serves as a compile-time configuration that tells generic classes
 * like Layout and Form how to behave for a structured grid. It contains no data itself,
 * only type definitions and static factories.
 */
struct StructuredCartesianPolicy {

    template <int Dim>
    using MeshType = StructuredCartesianMesh<Dim>;

    template <int Dim>
    using IndexerType = GrassmanIndex<Dim>;

    using LayoutStrategy = SerialLayoutStrategy;

    template <int Dim, int k>
    struct StorageModel {
        // The GroupTags for this policy are the Blade types themselves.
        using GroupTagTuple = typename BladesForGrade<Dim, k>::type;
    };

    template <typename BladeType, int Dim, typename MeshType>
    static auto create_geometry_provider(const MeshType& mesh) {
        using ProviderType = GeometryProvider<BladeType, Dim>;
        return ProviderType{mesh.get_spacing()};
    }

    template <typename From, typename To, typename Indexer>
    using ConnectivityProvider = Detail::StructuredIncidenceProvider<From, To, Indexer>;

    template <typename BladeType, int Dim>
    using GeometryProvider = Detail::StructuredGeometryProvider<BladeType, Dim>;

};

} // namespace fem
