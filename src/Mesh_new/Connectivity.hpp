#pragma once

#include "Mesh_new/GrassmanIndex.hpp"

namespace fem {
namespace Detail {
    // A placeholder functor. In a real system, this would contain the
    // arithmetic to resolve incidences using the provided indexer.
    template<typename IndexerType>
    struct DummyIncidenceFunctor {
        IndexerType indexer;
        KOKKOS_INLINE_FUNCTION void operator()() const {}
    };
} // namespace Detail

// Compile-time connectivity provider
template <typename FromBlade, typename ToBlade, int Dim>
struct Connectivity {
    // Generic fallback
    static auto get_functor(const GrassmanIndex<Dim>& indexer) {
        return Detail::DummyIncidenceFunctor{indexer};
    }
};

} // namespace fem
