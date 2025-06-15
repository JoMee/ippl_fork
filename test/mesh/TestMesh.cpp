#include "Ippl.h"

#include <array>
#include <iostream>
#include <random>
#include <typeinfo>

#include "Mesh_new/structured/GrassmanIndex.hpp"

#include <iostream>
#include <cassert>

#include "Utility/ParameterList.h"

using namespace fem;

template<typename Blade>
struct PrintDim {
    static int call() {
        return Blade::Dim;
    }
};


void test() {
  constexpr auto table = BladeDispatchTable<3, int>::create<PrintDim>();
  uint32_t m = Blade<0,2,1>::Mask;

  std::cout << table(m) << std::endl;



  GrassmanIndex<3>::coord_type global = {10, 10, 10};
  GrassmanIndex<3>::coord_type offset = {2, 2, 2};
  GrassmanIndex<3>::coord_type local  = {5, 5, 5};

  GrassmanIndex<3> index(global, offset, local);

  GrassmanIndex<3>::local_id_type point {10,10,10, 1};
  std::cout << index.mask(point) << std::endl;


}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
      test();
    }
    ippl::finalize();

    return 0;
}
