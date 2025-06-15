#pragma once
#include "index_map.hpp"
#include "comm_pattern.hpp"

namespace fem {

// abstract DOF layout concept (built on top of IndexMap)
template<class Layout>
concept FieldLayout = requires(const Layout& layout)
{
    typename Layout::global_index;
    typename Layout::local_index;
    { layout.entity_dimension() } -> std::same_as<int>;
    { layout.pattern() }          -> std::same_as<const CommPattern&>;
    { layout.global_extent() }    -> requires std::ranges::range;
};
} // namespace fem 

