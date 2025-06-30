#pragma once

#include "Mesh_new/Fields.hpp"
#include <tuple>

namespace fem {

struct P_r_Lambda {};
struct Q_r_Lambda {};

namespace Detail {
    template <typename GroupTagTuple, typename T, typename LayoutType>
    struct ComponentsFromGroups;

    template <typename... GroupTags, typename T, typename LayoutType>
    struct ComponentsFromGroups<std::tuple<GroupTags...>, T, LayoutType> {
        using type = std::tuple<
            Detail::ComponentField<T, GroupTags, LayoutType>...
        >;
    };
} // namespace Detail


// --- Generic SpaceTraits Template ---
// This is the primary template. It is specialized for each
// finite element space we want to support.
template <typename Family, int r, typename T, typename LayoutType>
struct SpaceTraits;


// --- Specialization for P_1 Lagrange Elements (hard-coded example) ---
template <typename T, typename LayoutType>
struct SpaceTraits<P_r_Lambda, 1, T, LayoutType> {
private:
    // Get the tuple of GroupTags for 0-cells (vertices) directly from your existing utility.
    // For any policy, this should correspond to the single vertex group.
    using VertexGroups = typename LayoutType::Policy::template StorageModel<LayoutType::DIM, 0>::GroupTagTuple;

public:
    using ComponentTuple = typename Detail::ComponentsFromGroups<VertexGroups, T, LayoutType>::type;

    template <typename... Cs> struct FormFromTuple;
    template <typename... Cs> struct FormFromTuple<std::tuple<Cs...>> { using type = Form<Cs...>; };
    using FormType = typename FormFromTuple<ComponentTuple>::type;

    template <typename GroupTag>
    static constexpr int get_num_dofs_per_entity() {
        return (GroupTag::Dim == 0) ? 1 : 0;
    }
};


template <typename T, typename LayoutType>
struct SpaceTraits<P_r_Lambda, 2, T, LayoutType> {
private:
    using VertexGroups = typename LayoutType::Policy::template StorageModel<LayoutType::DIM, 0>::GroupTagTuple;
    using EdgeGroups   = typename LayoutType::Policy::template StorageModel<LayoutType::DIM, 1>::GroupTagTuple;

    using AllDofGroups = typename Detail::Concat<VertexGroups, EdgeGroups>::type;

public:
    using ComponentTuple = typename Detail::ComponentsFromGroups<AllDofGroups, T, LayoutType>::type;

    template <typename... Cs> struct FormFromTuple;
    template <typename... Cs> struct FormFromTuple<std::tuple<Cs...>> { using type = Form<Cs...>; };
    using FormType = typename FormFromTuple<ComponentTuple>::type;

    template <typename GroupTag>
    static constexpr int get_num_dofs_per_entity() {
        if constexpr (GroupTag::Dim == 0) { 
            return 1; 
        }
        if constexpr (GroupTag::Dim == 1) { 
            return 1;
        }
        return 0;
    }
};

} // namespace fem

