#pragma once

#include <tuple>
#include <cstdint>

namespace fem {

// =================== Blade TMP Utilities ====================
namespace Detail {

    template<int... Dims>
    struct Blade {
        static constexpr int Dim = sizeof...(Dims);
        static constexpr uint32_t Mask = ((1u << Dims) | ... | 0u);
        using DimsSequence = std::integer_sequence<int, Dims...>;
    };


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
    struct GenerateAllBlades {
        using Prev = typename GenerateAllBlades<N - 1>::type;
        using Added = typename AddDimToAll<N - 1, Prev>::type;
        using type = typename Concat<Prev, Added>::type;
    };

    template<>
    struct GenerateAllBlades<0> {
        using type = std::tuple<Blade<>>;
    };

    template<int k, typename TupleOfBlades, typename ResultTuple = std::tuple<>>
    struct FilterBladesByGrade;

    template<int k, typename... Result>
    struct FilterBladesByGrade<k, std::tuple<>, std::tuple<Result...>> {
        using type = std::tuple<Result...>;
    };

    template<int k, typename Head, typename... Tail, typename... Result>
    struct FilterBladesByGrade<k, std::tuple<Head, Tail...>, std::tuple<Result...>> {
        using next_tuple = std::conditional_t<
            (Head::Dim == k),
            std::tuple<Result..., Head>,
            std::tuple<Result...>
        >;
        using type = typename FilterBladesByGrade<k, std::tuple<Tail...>, next_tuple>::type;
    };


} // namespace Detail

template<int N, int k>
struct BladesForGrade {
    using AllBlades = typename Detail::GenerateAllBlades<N>::type;
    using type = typename Detail::FilterBladesByGrade<k, AllBlades>::type;
};

template<int... Dims>
using Blade = Detail::Blade<Dims...>;

} // namespace fem 
