#include "Ippl.h"

#include <array>
#include <iostream>
#include <random>
#include <typeinfo>

#include "Mesh_new/FunctionSpace.hpp"

#include <iostream>
#include <cassert>

#include "Utility/ParameterList.h"

using namespace fem;

void test() {

    constexpr int Dim = 2;
    constexpr int k = 1; // 1-form
    constexpr int r = 1; // Using linear polynomials
    using T = double;
    using MeshPolicy = fem::StructuredCartesianPolicy;
    using Family = fem::Q_r_Family; // The correct family for structured grids

    Kokkos::Array<int, Dim> extents = {10, 20}; // A 10x20 grid
    Kokkos::Array<double, Dim> spacing = {0.1, 0.1};
    auto mesh = std::make_shared<fem::StructuredCartesianMesh<Dim>>(extents, spacing);

    using LayoutType = fem::Layout<Dim, MeshPolicy>;
    auto layout = std::make_shared<LayoutType>(mesh, 1); // Use a halo of width 1

    using FunctionSpaceType = fem::FunctionSpace<Family, k, r, T, LayoutType>;
    FunctionSpaceType Vh_1(layout);

    auto example_1_form = Vh_1.create_form();

    example_1_form.fillHalo();

}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
      test();
    }
    ippl::finalize();

    return 0;
}
