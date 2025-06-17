#pragma once
#include "Mesh_new/Blades.hpp"
namespace fem {
namespace Detail {
    // Functor providing geometric info for a specific family of elements (Blade).
    template<typename BladeType, int Dim>
    struct GeometryProvider {
        Kokkos::Array<double, Dim> spacing;
        KOKKOS_INLINE_FUNCTION double get_primal_measure(const Kokkos::Array<long, Dim>& /*logical_coord*/) const {
            if constexpr (BladeType::Dim == 1) { // An edge
                if constexpr (BladeType::Mask == (1u << 0)) return spacing[0]; // x-edge
                if constexpr (BladeType::Mask == (1u << 1)) return spacing[1]; // y-edge
                if constexpr (BladeType::Mask == (1u << 2)) return spacing[2]; // z-edge
            }
            return 1.0; // Default for vertices, etc.
        }
    };
}
template <typename BladeType, int Dim>
struct Geometry {
    static auto get_provider(const Kokkos::Array<double, Dim>& spacing) {
        return Detail::GeometryProvider<BladeType, Dim>{spacing};
    }
};

}
