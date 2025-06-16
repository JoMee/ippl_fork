#pragma once

#include <tuple>
#include <array>
#include <cstdint>
#include <stdexcept>
#include "Mesh_new/structured/Blades.hpp"

template<int N, typename Ret, typename... Args>
struct BladeDispatchTable {
    using FuncPtr = Ret(*)(Args...);
    std::array<FuncPtr, 1 << N> table{};

    template<template<typename> typename Impl>
    static constexpr BladeDispatchTable create() {
        BladeDispatchTable tbl;
        using AllBlades = typename Detail::GenerateBlades<N>::type;

        [&]<typename... Blades>(std::tuple<Blades...>) {
            ((tbl.table[Blades::Mask] = &Impl<Blades>::call), ...);
        }(AllBlades{});

        return tbl;
    }

    Ret operator()(uint32_t mask, Args... args) const {
        FuncPtr f = table[mask];
        if (f) return f(std::forward<Args>(args)...);
        throw std::runtime_error("Invalid or uninitialized blade mask");
    }
};

