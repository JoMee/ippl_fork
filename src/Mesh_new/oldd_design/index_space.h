#include <concepts>
#include <Kokkos_Core.hpp>

namespace fem {

template<typename T>
concept LocalIndexSpace = requires(const T& space,
                                    typename T::local_id_type lid,
                                    typename T::global_id_type gid)
{
    typename T::local_id_type;
    typename T::global_id_type;

    { space.local_to_global(lid) } noexcept
        -> std::same_as<typename T::global_id_type>;

    /**
     * @brief Maps a global index to a local index, if it exists on this rank.
     * This is the key change from std::optional.
     * @return A Kokkos::pair where .first is true if lid corresponding to gid was found, 
     * and .second is the corresponding local_id. If .first is false, .second is undefined.
     */
    { space.global_to_local(gid) } noexcept
        -> std::same_as<Kokkos::pair<bool, typename T::local_id_type>>;

};

}
