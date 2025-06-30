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
    constexpr int r = 2; 
    using T = double;
    using MeshPolicy = fem::StructuredCartesianPolicy;
    using Family = fem::P_r_Lambda;

    std::cout << "--- Setting up Mesh and Layout ---\n";
    Kokkos::Array<int, Dim> extents = {10, 20}; // A 10x20 vertex grid
    Kokkos::Array<double, Dim> spacing = {0.1, 0.1};
    auto mesh = std::make_shared<fem::StructuredCartesianMesh<Dim>>(extents, spacing);

    using LayoutType = fem::Layout<Dim, MeshPolicy>;
    auto layout = std::make_shared<LayoutType>(mesh, 1); // Use a halo of width 1

    std::cout << "\n--- Setting up P2 FunctionSpace ---\n";
    using P2_Lagrange_Space = fem::FunctionSpace<Family, r, T, LayoutType>;
    P2_Lagrange_Space Vh_P2(layout);

    std::cout << "\n--- Creating a P2 Lagrange Field ---\n";
    auto my_p2_field = Vh_P2.create_form();
    std::cout << "Successfully created a P2 Lagrange field.\n";

    my_p2_field.fillHalo();

  
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
      test();
    }
    ippl::finalize();

    return 0;
}
