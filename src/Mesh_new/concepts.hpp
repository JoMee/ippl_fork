#pragma once
#include <array>
#include <optional>
#include <concepts>

namespace mesh {

template<class M>
concept IndexMap = requires(const M& m,
                            typename M::GlobalIndex g,
                            typename M::LocalIndex  l) {
    { m.global_to_local(g) } -> std::same_as<std::optional<typename M::LocalIndex>>;
    { m.local_to_global(l) } -> std::same_as<typename M::GlobalIndex>;
    { m.owns(g) }           -> std::same_as<bool>;
    { m.is_ghost(g) }       -> std::same_as<bool>;
};

template<class M>
concept GeometryMap = requires(const M& m, typename M::GlobalIndex g) {
    { m.centroid(g) } -> std::same_as<typename M::Vec>;
    { m.measure (g) } -> std::same_as<typename M::Scalar>;
    { M::Dim }        -> std::convertible_to<int>;
};

template<class M>
concept MeshModel = IndexMap<M> && GeometryMap<M> && requires(
    const M& m,
    typename M::GlobalIndex src,
    int  tgt_dim,
    typename M::GlobalIndex* out,
    std::size_t max_n,
    std::size_t out_n) {
        { m.incident(src, tgt_dim, out, max_n, out_n) } -> std::same_as<void>;
};

