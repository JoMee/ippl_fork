#pragma once
#include "Mesh_new/common.hpp"
#include "Mesh_new/index_map.hpp"
#include "Mesh_new/structured/index_map_impl.hpp" // <- flatten/unflatten utils

namespace fem {

template<int Dim>
class StructuredIndexMap {
public:
    using global_index = std::uint64_t;
    using local_index  = int;

    StructuredIndexMap(std::array<int, Dim> global_shape,
                       std::array<int, Dim> global_offset,
                       std::array<int, Dim> local_shape,
                       int ghost_width)
    : global_shape_{global_shape},
      global_offset_{global_offset},
      local_shape_{local_shape},
      ghost_width_{ghost_width}
    {
        // precompute bounds
        for (int d = 0; d < Dim; ++d) {
            owned_lower_[d] = global_offset_[d];
            owned_upper_[d] = global_offset_[d] + local_shape_[d];

            ghost_lower_[d] = owned_lower_[d] - ghost_width_;
            ghost_upper_[d] = owned_upper_[d] + ghost_width_;
        }
    }

    local_index global_to_local(global_index g) const {
        auto idx = unflatten<Dim>(g, global_shape_);
        if (!within_bounds<Dim>(idx, ghost_lower_, ghost_upper_))
            return -1;

        for (int d = 0; d < Dim; ++d)
            idx[d] -= ghost_lower_[d];

        return flatten<Dim>(idx, extended_shape());
    }

    global_index local_to_global(local_index l) const {
        auto idx = unflatten<Dim>(l, extended_shape());
        for (int d = 0; d < Dim; ++d)
            idx[d] += ghost_lower_[d];
        return flatten<Dim>(idx, global_shape_);
    }

    Ownership ownership(global_index g) const {
        auto idx = unflatten<Dim>(g, global_shape_);
        if (within_bounds<Dim>(idx, owned_lower_, owned_upper_))
            return Ownership::Owned;
        else if (within_bounds<Dim>(idx, ghost_lower_, ghost_upper_))
            return Ownership::Ghost;
        else
            return Ownership::External;
    }

    const std::array<int, Dim>& global_shape() const { return global_shape_; }
    const std::array<int, Dim>& global_offset() const { return global_offset_; }
    const std::array<int, Dim>& local_shape() const { return local_shape_; }

private:
    std::array<int, Dim> global_shape_;
    std::array<int, Dim> global_offset_;
    std::array<int, Dim> local_shape_;
    int ghost_width_;

    std::array<int, Dim> owned_lower_, owned_upper_;
    std::array<int, Dim> ghost_lower_, ghost_upper_;

    std::array<int, Dim> extended_shape() const {
        std::array<int, Dim> ext;
        for (int d = 0; d < Dim; ++d)
            ext[d] = (ghost_upper_[d] - ghost_lower_[d]);
        return ext;
    }
};

// Device traits to satisfy concept
template<int Dim>
struct map_traits<StructuredIndexMap<Dim>> {
    KOKKOS_INLINE_FUNCTION
    static int device_global_to_local(const StructuredIndexMap<Dim>& map,
                                      std::uint64_t gid)
    {
        return map.global_to_local(gid);
    }

    KOKKOS_INLINE_FUNCTION
    static std::uint64_t device_local_to_global(const StructuredIndexMap<Dim>& map,
                                                int lid)
    {
        return map.local_to_global(lid);
    }
};

} // namespace fem

