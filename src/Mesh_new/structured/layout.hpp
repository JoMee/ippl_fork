#pragma once
#include "Mesh_new/field_layout.hpp"

namespace fem {

// ---- IndexSpace concept for structured grids ---------------------------------------------
template<class S>
concept IndexSpace = requires(const S& s, int ax) {
    static_cast<int>(S::dimension);
    { s.lower(ax) }       -> std::same_as<int>;
    { s.upper(ax) }       -> std::same_as<int>;
    { s.is_periodic(ax) } -> std::same_as<bool>;
    { s.global_offset() } -> std::same_as<std::array<int, S::dimension>>;
};

// ---- Tensor-product DOF layout ------------------------------------------------------------
template<IndexSpace Space>
class StructuredLayout {
public:
    using global_index = default_gid_t;
    using local_index  = default_lid_t;

    StructuredLayout(int entity_dim, const CommPattern& p)
    : dim_{entity_dim}, pattern_{p} {}

    int entity_dimension() const            { return dim_; }
    const CommPattern& pattern() const      { return pattern_; }

private:
    int          dim_;
    CommPattern  pattern_;
};


} // namespace fem

