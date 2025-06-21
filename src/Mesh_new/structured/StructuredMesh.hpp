#pragma once

namespace fem {

template <int Dim>
class StructuredCartesianMesh {
public:
    explicit StructuredCartesianMesh(const Kokkos::Array<int, Dim>& extents,
                                     const Kokkos::Array<double, Dim>& spacing)
        : logical_extents_(extents), grid_spacing_(spacing) {}

    auto get_extents() const { return logical_extents_; }

    auto get_spacing() const {return grid_spacing_; }

private:
    Kokkos::Array<int, Dim> logical_extents_;
    Kokkos::Array<double, Dim> grid_spacing_;
};

} // namespace fem

