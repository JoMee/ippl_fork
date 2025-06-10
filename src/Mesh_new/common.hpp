#pragma once
#include <cstdint>

namespace fem {

using default_gid_t = std::uint64_t;   // global entity id
using default_lid_t = std::int32_t;    // local  entity id

// ---- ownership flags (bit-maskable) -------------------------------------------------------
enum class Ownership : unsigned char { Owned = 1, Ghost = 2, External = 4 };

constexpr Ownership operator|(Ownership a, Ownership b)
{
    return static_cast<Ownership>(static_cast<unsigned>(a) |
                                  static_cast<unsigned>(b));
}
constexpr bool any(Ownership o, Ownership mask)
{
    return static_cast<unsigned>(o) & static_cast<unsigned>(mask);
}

// orientation sign (+1/-1) for edges / faces (FEEC-ready even if unused now)
using Orientation = signed char;

} // namespace fem 

