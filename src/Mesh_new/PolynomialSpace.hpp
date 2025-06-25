#pragma once

#include <cstdint>
#include <type_traits> // For std::is_same_v

namespace fem {

// NOTE: This is just a placeholder for very specific elements in specific dimensions.
// We should replace this with a compile-time generated struct templating mechanism.

struct P_r_Lambda {}; // For simplex elements (e.g., triangles, tetrahedra).
struct Q_r_Family {}; // For tensor-product elements (e.g., quadrilaterals, hexahedra).

// At runtime, each element in the mesh will have a type ID.
// This allows a single function space to operate on a hybrid mesh.
// NOTE: For a robust, extensible design, this enum would eventually be replaced
// by a system of struct tags to allow for compile-time dispatch. For this
// prototype, the enum is fine for now.
enum class ElementShape {
    Vertex,
    Line,
    Triangle,
    Quadrilateral,
    Tetrahedron,
    Hexahedron
};

template <typename Family, int r, int ElementDim>
class PolynomialSpace;


template <int r, int ElementDim>
class PolynomialSpace<P_r_Lambda, r, ElementDim> {
public:
    PolynomialSpace() = default;

    int get_num_coefficients(ElementShape shape) const {
        if constexpr (r == 1) {
            // For P_1 (linear) polynomials on an n-simplex, there are n+1 coeffs.
            switch (shape) {
                case ElementShape::Vertex:      return 1; // 0D simplex
                case ElementShape::Line:        return 2; // 1D simplex
                case ElementShape::Triangle:    return 3; // 2D simplex
                case ElementShape::Tetrahedron: return 4; // 3D simplex
                default:                        return -1; // Not a simplex
            }
        }
        return -1; 
    }
};

template <int r, int ElementDim>
class PolynomialSpace<Q_r_Family, r, ElementDim> {
public:
    PolynomialSpace() = default;

    int get_num_coefficients(ElementShape shape) const {
        if constexpr (r == 1) {
            // For Q_1 (bi/tri-linear) elements on an n-cube, there are 2^n coeffs.
            switch (shape) {
                case ElementShape::Vertex:        return 1; // 0D cube
                case ElementShape::Line:          return 2; // 1D cube
                case ElementShape::Quadrilateral: return 4; // 2D cube
                case ElementShape::Hexahedron:    return 8; // 3D cube
                default:                          return -1; // Not a cube
            }
        }
        return -1;
    }
};

}
