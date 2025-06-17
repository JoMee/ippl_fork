#include "Ippl.h"

#include <array>
#include <iostream>
#include <random>
#include <typeinfo>

#include "Mesh_new/GrassmanIndex.hpp"
#include "Mesh_new/StructuredMesh.hpp"
#include "Mesh_new/SerialLayout.hpp"
#include "Mesh_new/Fields.hpp"

#include <iostream>
#include <cassert>

#include "Utility/ParameterList.h"

using namespace fem;


void test() {
  constexpr int Dim = 3;
  using T = double;
  using MeshType = StructuredCartesianMesh<Dim>;
  using LayoutType = Layout<MeshType>;

  auto mesh = std::make_shared<MeshType>(Kokkos::Array<int, Dim>{10, 10, 10},
                                         Kokkos::Array<double, Dim>{1.0,1.0,1.0});
  auto layout = std::make_shared<LayoutType>(mesh, 1);
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
