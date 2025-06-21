#include "Ippl.h"

#include <array>
#include <iostream>
#include <random>
#include <typeinfo>

#include "Mesh_new/structured/StructuredMesh.hpp"
#include "Mesh_new/Layout.hpp"
#include "Mesh_new/Fields.hpp"
#include "Mesh_new/MeshPolicies.hpp"

#include <iostream>
#include <cassert>

#include "Utility/ParameterList.h"

using namespace fem;


void test() {

    constexpr int Dim = 3;
    using T = double;
    const int halo_width = 1;

    using MeshPolicy = StructuredCartesianPolicy;

    using LayoutType = Layout<Dim, MeshPolicy>;

    using MeshType = typename LayoutType::MeshType;

    auto mesh = std::make_shared<MeshType>(
        Kokkos::Array<int, Dim>{10, 10, 10},
        Kokkos::Array<double, Dim>{1.0, 1.0, 1.0}
    );

    auto layout = std::make_shared<LayoutType>(mesh, halo_width);

    Form<2, T, LayoutType> E(layout);

    E.fillHalo();

}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
      test();
    }
    ippl::finalize();

    return 0;
}
