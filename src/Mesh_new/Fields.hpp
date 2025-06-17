#pragma once

#include "Mesh_new/SerialLayout.hpp"
#include <iostream>
#include <string>
#include <sstream>

namespace fem {

// --- Private Implementation Detail: ComponentField ---
namespace Detail {

// TMP helper to recursively generate a pointer type with N levels of indirection.
// e.g., NPtr<T, 3>::type is T***
template <typename T, int N>
struct NPtr {
    using type = typename NPtr<T, N - 1>::type*;
};
template <typename T>
struct NPtr<T, 1> {
    using type = T*;
};

template <typename T, typename BladeType, typename LayoutType>
class ComponentField {
public:
    static constexpr int Dim = LayoutType::DIM;
    // Use the general NPtr to create the data type for the Kokkos::View (e.g., double***)
    using view_data_type = typename NPtr<T, Dim>::type;
    // The final, N-dimensional Kokkos::View type
    using view_type = Kokkos::View<view_data_type>;

    explicit ComponentField(std::shared_ptr<const LayoutType> layout) : layout_(layout) {
        auto extent = layout->template get_alloc_extent<BladeType>();
        
        auto construct_view = [&]<size_t... Is>(std::index_sequence<Is...>) {
            return view_type("component", extent[Is]...);
        };
        view_ = construct_view(std::make_index_sequence<Dim>{});
        
        std::stringstream ss;
        for(int i=0; i < Dim; ++i) {
            ss << extent[i] << (i < Dim - 1 ? "x" : "");
        }
        std::cout << "    - ComponentField for Blade<" << typeid(BladeType).name() << "> allocated with size "
                  << ss.str() << std::endl;
    }

    ComponentField& operator=(T scalar) {
        Kokkos::deep_copy(view_, scalar);
        return *this;
    }

    void fillHalo() { layout_->fill_halo(); }
    view_type& view() { return view_; }
    const view_type& view() const { return view_; }
private:
    std::shared_ptr<const LayoutType> layout_;
    view_type view_;
};
} // namespace Detail


// --- Public-Facing Abstraction: Form ---
template <int k, typename T, typename LayoutType>
class Form {
private:
    // Helper to create a tuple of ComponentFields from a tuple of Blades
    template <typename> struct ComponentTupleFromBlades;
    template <typename... BladeTypes>
    struct ComponentTupleFromBlades<std::tuple<BladeTypes...>> {
        using type = std::tuple<Detail::ComponentField<T, BladeTypes, LayoutType>...>;
    };

    using BladeTuple = typename BladesForGrade<LayoutType::DIM, k>::type;
    using ComponentTuple = typename ComponentTupleFromBlades<BladeTuple>::type;

    std::shared_ptr<const LayoutType> layout_;
    ComponentTuple components_;

    static ComponentTuple create_components(std::shared_ptr<const LayoutType> layout) {
        return std::apply(
            [&](auto... blade_args) {
                return std::make_tuple(Detail::ComponentField<T, decltype(blade_args), LayoutType>(layout)...);
            }, BladeTuple{});
    }

public:
    explicit Form(std::shared_ptr<const LayoutType> layout)
        : layout_(layout), components_(create_components(layout)) {}

    void fillHalo() {
        std::cout << "Form<" << k << ">: Orchestrating halo exchange..." << std::endl;
        std::apply([](auto&... component) { (component.fillHalo(), ...); }, components_);
    }

    template <typename BladeType>
    const auto& get_component() const {
        return std::get<Detail::ComponentField<T, BladeType, LayoutType>>(components_);
    }
};

} // namespace fem
