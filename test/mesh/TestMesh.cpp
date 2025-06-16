#include "Ippl.h"

#include <array>
#include <iostream>
#include <random>
#include <typeinfo>

#include "Mesh_new/GrassmanIndex.hpp"
#include "Mesh_new/StructuredMesh.hpp"

#include <iostream>
#include <cassert>

#include "Utility/ParameterList.h"

using namespace fem;


void test() {
  constexpr int N = 2;

  // Create GrassmanIndex for a 4x4 grid
  StructuredMesh<N> mesh(
      {4, 4},  // global_shape
      {0, 0},  // offset
      {4, 4}   // local_extent
  );

  std::cout << mesh.entity_count<Blade<>>() << std::endl;
  std::cout << std::endl;

}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
      test();
    }
    ippl::finalize();

    return 0;
}
