#include "Ippl.h"

#include <array>
#include <iostream>
#include <random>
#include <typeinfo>

#include "Mesh_new/mesh.hpp"

#include <iostream>
#include <cassert>

#include "Utility/ParameterList.h"

using namespace fem;

void test_2d_decomposition() {
  constexpr int Dim = 2;
  using IndexMap = StructuredIndexMap<Dim>;

  std::array<int, Dim> global_lower = {0, 0};
  std::array<int, Dim> global_upper = {64, 64};
  int num_parts = 5;

  RCBPartitioner<Dim> partitioner(1); // ghost layers = 1

  auto index_maps = partitioner.partition(global_lower, global_upper, num_parts);
  assert(index_maps.size() == static_cast<std::size_t>(num_parts));

  for (IndexMap map : index_maps) {
    std::array<int, Dim> local_shape_ = map.local_shape();
    std::cout << "[" << local_shape_[0] << ", "<< local_shape_[1] << "]" << std::endl;

  }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
      test_2d_decomposition();
    }
    ippl::finalize();

    return 0;
}
