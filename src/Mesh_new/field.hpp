#pragma once
#include "mesh/concepts.hpp"
#include <Kokkos_Core.hpp>

namespace mesh {

/// Simple 1‑component field container

template<MeshModel M, typename T>
class Field {
public:
    using LocalIndex = typename M::LocalIndex;
    using View = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>;

    Field(const M& mesh, std::size_t n_local) : data_("f", n_local), mesh_(mesh) {}

    KOKKOS_INLINE_FUNCTION       T& operator[](LocalIndex i)       { return data_(i); }
    KOKKOS_INLINE_FUNCTION const T& operator[](LocalIndex i) const { return data_(i); }

    View view() const { return data_; }

    void update_ghosts() { /* TODO: pass to HaloExchangePlan */ }

private:
    View data_;
    const M& mesh_;
};
