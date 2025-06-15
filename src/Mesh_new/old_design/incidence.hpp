#pragma once
#include <Kokkos_Core.hpp>
#include "common.hpp"

namespace fem {

// device-copyable read-only view of an incidence list  (CSR storage lives elsewhere)
template<class Index>
struct IncidenceView {
    const Index*       ids   = nullptr;          // global indices of incident entities
    const Orientation* sign  = nullptr;          // nullptr ⇒ all +1
    int                count = 0;

    KOKKOS_INLINE_FUNCTION const Index operator[](int i) const { return ids[i]; }
    KOKKOS_INLINE_FUNCTION Orientation orientation(int i) const { return sign ? sign[i] : 1; }
};

} // namespace fem 

