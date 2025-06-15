#pragma once
#include "common.hpp"

namespace fem {

// forward declaration that implementers specialise for their map type
template<class Map> struct map_traits;

// ----------------------- IndexMap concept (host & device) ---------------------------------
template<class Map>
concept IndexMap = requires(const Map& m,
                            typename Map::global_index g,
                            typename Map::local_index  l)
{
    typename Map::global_index;
    typename Map::local_index;

    // host functions (-1 sentinel for invalid)
    { m.global_to_local(g) } -> std::same_as<typename Map::local_index>;
    { m.local_to_global(l) } -> std::same_as<typename Map::global_index>;
    { m.ownership(g) }      -> std::same_as<Ownership>;

    // device side (via traits to avoid duplicate member functions)
    { map_traits<Map>::device_global_to_local(m, g) } -> std::same_as<int>;
    { map_traits<Map>::device_local_to_global(m, l) } -> std::same_as<typename Map::global_index>;
};

} // namespace fem 

