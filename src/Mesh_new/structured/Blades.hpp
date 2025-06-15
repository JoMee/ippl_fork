#pragma once

#include <tuple>
#include <cstdint>

template<int... Dims>
struct Blade {
    static constexpr int Dim = sizeof...(Dims);
    static constexpr uint32_t Mask = ((1u << Dims) | ... | 0u);
};

// =================== Blade TMP Utilities ====================
namespace Detail {

    template<typename A, typename B>
    struct Concat;

    template<typename... T1, typename... T2>
    struct Concat<std::tuple<T1...>, std::tuple<T2...>> {
        using type = std::tuple<T1..., T2...>;
    };

    template<int D, typename Tuple>
    struct AddDimToAll;

    template<int D, typename... Blades>
    struct AddDimToAll<D, std::tuple<Blades...>> {
        template<typename B> struct Add;
        template<int... Ds> struct Add<Blade<Ds...>> {
            using type = Blade<Ds..., D>;
        };
        using type = std::tuple<typename Add<Blades>::type...>;
    };

    template<int N>
    struct GenerateBlades {
        using Prev = typename GenerateBlades<N - 1>::type;
        using Added = typename AddDimToAll<N - 1, Prev>::type;
        using type = typename Concat<Prev, Added>::type;
    };

    template<>
    struct GenerateBlades<0> {
        using type = std::tuple<Blade<>>;
    };

} // namespace Detail

