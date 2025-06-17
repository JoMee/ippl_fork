// =========================================================================
// CMakeLists.txt
// =========================================================================
/*
cmake_minimum_required(VERSION 3.16)
project(FemFrameworkPlayground CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Kokkos REQUIRED)

add_executable(playground main.cpp)

target_link_libraries(playground PRIVATE Kokkos::kokkos)
*/


// =========================================================================
// Blades.hpp
// =========================================================================
#pragma once

#include <tuple>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace fem {
namespace detail {

template<int... Dims>
struct Blade {
    static constexpr int Dim = sizeof...(Dims);
    static constexpr uint32_t Mask = ((1u << Dims) | ... | 0u);
};

// Helper to generate all 2^N blades for N dimensions
template <typename T, int... N>
auto concat_tuples(std::integer_sequence<int, N...>, T t) {
    return std::tuple_cat(std::get<N>(t)...);
}

template <int N, typename T = std::tuple<std::tuple<detail::Blade<>>>>
struct blade_generator {
    using Next = typename blade_generator<N - 1>::type;
    using type = decltype(concat_tuples(
        std::make_integer_sequence<int, std::tuple_size<Next>::value>(),
        std::make_tuple(std::make_tuple(
            std::get<std::make_integer_sequence<int, std::tuple_size<Next>::value>>(Next){},
            detail::Blade<std::get<std::make_integer_sequence<int, std::tuple_size<Next>::value>>(Next)::Dims..., N - 1>{}
        )...)
    ));
};
template <> struct blade_generator<0> { using type = std::tuple<std::tuple<detail::Blade<>>>; };


// Helper to filter a tuple of blades by their dimension (grade)
template<int k, typename TupleOfBlades, typename ResultTuple = std::tuple<>>
struct FilterBladesByGrade;

template<int k, typename... Result>
struct FilterBladesByGrade<k, std::tuple<>, std::tuple<Result...>> {
    using type = std::tuple<Result...>;
};

template<int k, typename Head, typename... Tail, typename... Result>
struct FilterBladesByGrade<k, std::tuple<Head, Tail...>, std::tuple<Result...>> {
    using next_tuple = std::conditional_t<
        (Head::Dim == k),
        std::tuple<Result..., Head>,
        std::tuple<Result...>
    >;
    using type = typename FilterBladesByGrade<k, std::tuple<Tail...>, next_tuple>::type;
};

} // namespace detail

// Public-facing alias
template<int... Dims>
using Blade = detail::Blade<Dims...>;


// Public-facing helper to get all blades of a specific grade k in N dimensions
template<int N, int k>
struct BladesForGrade {
    using AllBlades = typename detail::blade_generator<N>::type;
    using type = typename detail::FilterBladesByGrade<k, AllBlades>::type;
};


// =========================================================================
// GrassmanIndex.hpp
// =========================================================================
#pragma once

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <unordered_map>

namespace fem {

// A simple container for subspace info.
template <int N>
struct IndexView {
    Kokkos::Array<long, N> offset;
    Kokkos::Array<long, N> extent;
};

// The evolved, unified indexer for structured grids.
template <int N>
class GrassmanIndex {
public:
    using coord_type = Kokkos::Array<long, N>;

    explicit GrassmanIndex(const coord_type& logical_extent, int halo_width = 0)
        : logical_vertex_extent_(logical_extent), halo_width_(halo_width)
    {}

    // --- API for the LOGICAL VIEW (used by Mesh/Connectivity) ---

    template <typename B>
    auto get_logical_extent() const -> coord_type {
        coord_type extent;
        for (int i = 0; i < N; ++i) {
            extent[i] = logical_vertex_extent_[i] - ((B::Mask >> i) & 1);
        }
        return extent;
    }

    template <typename B>
    KOKKOS_INLINE_FUNCTION long logical_coord_to_flat(const coord_type& ijk) const {
        auto logical_extents = get_logical_extent<B>();
        long flat_index = 0;
        long stride = 1;
        for (int i = 0; i < N; ++i) {
            flat_index += ijk[i] * stride;
            stride *= logical_extents[i];
        }
        return flat_index;
    }

    // --- API for the STORAGE VIEW (used by Layout) ---

    template <typename B>
    auto get_allocated_extent() const -> coord_type {
        auto extent = get_logical_extent<B>();
        for (int i = 0; i < N; ++i) {
            extent[i] += 2 * halo_width_;
        }
        return extent;
    }

    KOKKOS_INLINE_FUNCTION auto get_storage_offset() const -> coord_type {
        coord_type offset;
        for (int i = 0; i < N; ++i) offset[i] = halo_width_;
        return offset;
    }

private:
    coord_type logical_vertex_extent_;
    int halo_width_;
};


// =========================================================================
// Connectivity.hpp
// =========================================================================
#pragma once

#include "GrassmanIndex.hpp"

namespace fem {
namespace detail {
    // A placeholder functor. In a real system, this would contain the
    // arithmetic to resolve incidences using the provided indexer.
    template<typename IndexerType>
    struct DummyIncidenceFunctor {
        IndexerType indexer;
        KOKKOS_INLINE_FUNCTION void operator()() const {}
    };
} // namespace detail

// Compile-time connectivity provider
template <typename FromBlade, typename ToBlade, int Dim>
struct Connectivity {
    // Generic fallback
    static auto get_functor(const GrassmanIndex<Dim>& indexer) {
        return detail::DummyIncidenceFunctor{indexer};
    }
};


// =========================================================================
// Mesh.hpp
// =========================================================================
#pragma once

#include "Connectivity.hpp"
#include <memory>

namespace fem {

// Represents the pure topology of the local grid.
template <int Dim>
class StructuredCartesianMesh {
public:
    explicit StructuredCartesianMesh(const Kokkos::Array<long, Dim>& extents)
        : logical_extents_(extents), indexer_(nullptr) {}

    auto get_extents() const { return logical_extents_; }

    template <typename FromBlade, typename ToBlade>
    auto get_incidence() const {
        if (!indexer_) {
            throw std::runtime_error("Mesh's indexer not set by a Layout.");
        }
        return Connectivity<FromBlade, ToBlade, Dim>::get_functor(*indexer_);
    }

    void set_indexer(std::shared_ptr<const GrassmanIndex<Dim>> indexer) {
        indexer_ = indexer;
    }

private:
    Kokkos::Array<long, Dim> logical_extents_;
    std::shared_ptr<const GrassmanIndex<Dim>> indexer_;
};


// =========================================================================
// SerialLayout.hpp
// =========================================================================
#pragma once

#include "Mesh.hpp"

namespace fem {

// Forward declaration
template <typename MeshType>
class Layout;

// Specialization for a single-rank layout on a structured mesh
template <int Dim>
class Layout<StructuredCartesianMesh<Dim>> {
public:
    static constexpr int DIM = Dim;
    using MeshType = StructuredCartesianMesh<Dim>;

    explicit Layout(std::shared_ptr<MeshType> mesh, int halo_width = 0)
        : mesh_(mesh)
    {
        // Create the single indexer instance and inject it into the mesh.
        auto indexer = std::make_shared<GrassmanIndex<Dim>>(mesh->get_extents(), halo_width);
        indexer_ = indexer;
        mesh_->set_indexer(indexer);
        std::cout << "Layout created. Mesh is now configured with a shared indexer." << std::endl;
    }

    template <typename B>
    auto get_alloc_extent() const {
        return indexer_->template get_allocated_extent<B>();
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


// =========================================================================
// Fields.hpp
// =========================================================================
#pragma once

#include "SerialLayout.hpp"
#include <iostream>

namespace fem {

// --- Private Implementation Detail: ComponentField ---
namespace detail {
template <typename T, typename B, typename LT>
class ComponentField {
public:
    using view_type = Kokkos::View<T**>; // Simplified for 2D playground

    explicit ComponentField(std::shared_ptr<const LT> layout) : layout_(layout) {
        auto extent = layout->template get_alloc_extent<B>();
        view_ = view_type("component", extent[0], extent[1]);
        std::cout << "    - ComponentField for Blade<" << typeid(B).name() << "> allocated with size "
                  << extent[0] << "x" << extent[1] << std::endl;
    }

    ComponentField& operator=(T scalar) {
        Kokkos::deep_copy(view_, scalar);
        return *this;
    }

    void fillHalo() { layout_->fill_halo(); }
    view_type& view() { return view_; }
    const view_type& view() const { return view_; }
private:
    std::shared_ptr<const LT> layout_;
    view_type view_;
};
} // namespace detail


// --- Public-Facing Abstraction: Form ---
template <int k, typename T, typename LayoutType>
class Form {
private:
    template <typename> struct ComponentTupleFromBlades;
    template <typename... BT>
    struct ComponentTupleFromBlades<std::tuple<BT...>> {
        using type = std::tuple<detail::ComponentField<T, BT, LayoutType>...>;
    };

    using BladeTuple = typename BladesForGrade<LayoutType::DIM, k>::type;
    using ComponentTuple = typename ComponentTupleFromBlades<BladeTuple>::type;

    std::shared_ptr<const LayoutType> layout_;
    ComponentTuple components_;

    static ComponentTuple create_components(std::shared_ptr<const LayoutType> layout) {
        return std::apply(
            [&](auto... blade_args) {
                return std::make_tuple(detail::ComponentField<T, decltype(blade_args), LayoutType>(layout)...);
            }, BladeTuple{});
    }

public:
    explicit Form(std::shared_ptr<const LayoutType> layout)
        : layout_(layout), components_(create_components(layout)) {}

    Form& operator=(T scalar) {
        std::apply([&](auto&... component) { (component = scalar, ...); }, components_);
        return *this;
    }

    void fillHalo() {
        std::cout << "Form<" << k << ">: Orchestrating halo exchange..." << std::endl;
        std::apply([](auto&... component) { (component.fillHalo(), ...); }, components_);
    }

    template <typename B>
    const auto& get_component() const {
        return std::get<detail::ComponentField<T, B, LayoutType>>(components_);
    }
};

} // namespace fem


// =========================================================================
// main.cpp
// =========================================================================

#include "Fields.hpp"
#include <Kokkos_Core.hpp>
#include <iostream>

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        std::cout << "--- FEM Framework Playground ---" << std::endl;

        constexpr int Dim = 2;
        using T = double;
        using MeshType = fem::StructuredCartesianMesh<Dim>;
        using LayoutType = fem::Layout<MeshType>;

        // 1. Create the Mesh (pure topology)
        std::cout << "\n1. Creating Mesh..." << std::endl;
        auto mesh = std::make_shared<MeshType>(fem::Kokkos::Array<long, Dim>{10, 10});

        // 2. Create the Layout (owns indexer, injects it into mesh)
        //    Using 1 halo layer for demonstration.
        std::cout << "\n2. Creating Layout..." << std::endl;
        auto layout = std::make_shared<LayoutType>(mesh, 1);

        // 3. Create the user-facing Form object.
        //    This triggers allocation of its internal ComponentFields.
        std::cout << "\n3. Creating Form<1> (Vector Field)..." << std::endl;
        fem::Form<1, T, LayoutType> E(layout);

        // 4. Perform a simple operation.
        std::cout << "\n4. Assigning scalar value 12.3 to Form..." << std::endl;
        E = 12.3;

        // 5. Call a layout-dependent operation.
        std::cout << "\n5. Calling fillHalo()..." << std::endl;
        E.fillHalo();

        // 6. Inspect a value to verify everything is working.
        std::cout << "\n6. Verifying component data..." << std::endl;
        const auto& Ex_comp = E.get_component<fem::Blade<0>>();
        auto Ex_host = Kokkos::create_mirror_view(Ex_comp.view());
        Kokkos::deep_copy(Ex_host, Ex_comp.view());

        // We access (1,1) which corresponds to the logical (0,0) point's data
        // inside the haloed view.
        std::cout << "   Value of Ex component at a point: " << Ex_host(1, 1) << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
