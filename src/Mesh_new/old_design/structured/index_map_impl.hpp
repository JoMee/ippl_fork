#pragma once
#include <array>
#include <cstdint>

namespace fem {

// Multi-dimensional flattening (row-major)
template<int Dim>
inline std::uint64_t flatten(const std::array<int, Dim>& idx,
                             const std::array<int, Dim>& shape)
{
    std::uint64_t linear = 0;
    std::uint64_t stride = 1;
    for (int d = 0; d < Dim; ++d) {
        linear += idx[d] * stride;
        stride *= shape[d];
    }
    return linear;
}

template<int Dim>
inline std::array<int, Dim> unflatten(std::uint64_t gid,
                                      const std::array<int, Dim>& shape)
{
    std::array<int, Dim> idx = {};
    for (int d = 0; d < Dim; ++d) {
        idx[d] = gid % shape[d];
        gid /= shape[d];
    }
    return idx;
}

template<int Dim>
inline bool within_bounds(const std::array<int, Dim>& idx,
                          const std::array<int, Dim>& lower,
                          const std::array<int, Dim>& upper)
{
    for (int d = 0; d < Dim; ++d)
        if (idx[d] < lower[d] || idx[d] >= upper[d])
            return false;
    return true;
}

} // namespace fem

