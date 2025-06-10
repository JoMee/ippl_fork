#pragma once
#include "index_map.hpp"
#include "comm_pattern.hpp"

namespace fem {

// abstract DOF layout concept (built on top of IndexMap)
template<class L>
concept FieldLayout = IndexMap<L> && requires(const L& layout) {
    { layout.entity_dimension() } -> std::same_as<int>;
    { layout.pattern() }          -> std::same_as<const CommPattern&>;
};

} // namespace fem 

