#pragma once
#include <Kokkos_Core.hpp>
#include "common.hpp"

namespace fem {

// flattened CSR communication map (host & device readable)
struct CommPattern {
    Kokkos::View<const int*>           neighbor_ranks;          // size N
    Kokkos::View<const int*>           send_offset, recv_offset; // size N+1
    Kokkos::View<const default_gid_t*> send_ids,   recv_ids;    // flattened GIDs
};

} // namespace fem 

