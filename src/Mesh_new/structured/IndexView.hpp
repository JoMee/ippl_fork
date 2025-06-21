#pragma once

#include <Kokkos_Core.hpp>

template <int N>
struct IndexView {
  Kokkos::Array<int, N> offset;
  Kokkos::Array<int, N> extent;

  KOKKOS_INLINE_FUNCTION
  Kokkos::Array<int, N> upper() const {
    Kokkos::Array<int, N> result;
    for (int i = 0; i < N; ++i)
      result[i] = offset[i] + extent[i];
    return result;
  }
};

