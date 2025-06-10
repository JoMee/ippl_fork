#pragma once
#include "field.hpp"
#include "comm_pattern.hpp"

namespace fem {

struct CommsHandle { void* impl = nullptr; };

// user-facing concept
template<class E, class F>
concept HaloExchange =
    Field<F> &&
    requires(E ex, F& f, CommsHandle h)
{
    { ex.begin(f) } -> std::same_as<CommsHandle>;
    { ex.end(h) };
    { ex.poll(h) } -> std::same_as<bool>;
};

enum class CommPolicy { Persistent, OneShot };

// declaration of the default exchanger; heavy impl lives in halo_exchange.cpp
template<CommPolicy P = CommPolicy::Persistent>
class HaloExchanger {
public:
    explicit HaloExchanger(const CommPattern&);           // ctor in cpp

    CommsHandle begin(Field auto& f);                     // begin exchange
    void        end(CommsHandle);                         // wait
    bool        poll(CommsHandle);                        // progress
};

} // namespace fem 

