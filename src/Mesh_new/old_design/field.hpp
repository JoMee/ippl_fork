#pragma once
#include <Kokkos_Core.hpp>
#include "field_layout.hpp"

namespace fem {

// lightweight data container concept
template<class F>
concept Field = requires(F& f, typename F::local_index idx) {
    typename F::value_type;
    typename F::device_type;
    typename F::memory_space;
    typename F::local_index;
    typename F::layout_type;

    { f.layout() } noexcept -> std::same_as<const typename F::layout_type&>;
    { f.view() }   noexcept;                        // a Kokkos::View
    { f(idx) }           -> std::same_as<typename F::value_type&>;  // host convenience
    { typename F::HostMirror{} };
};

} // namespace fem 

