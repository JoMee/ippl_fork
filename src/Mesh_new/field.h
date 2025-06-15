#include <concepts>
#include <kokkos/Kokkos_Core.hpp>
#include "Mesh_new/index_space.h"

namespace fem {

template<typename T>
concept Field = requires(const T& field,
                         typename T::local_id_type lid,
                         int component) {
    typename T::value_type;
    typename T::device_type;

    // The data is just a Kokkos::View
    { field.view() } -> std::same_as<Kokkos::View<typename T::value_type**, typename T::device_type>>;

    
    { field(lid, component) } -> std::convertible_to<typename T::value_type>;

    { field.get_index_space() } -> LocalIndexSpace;
};

}
