#pragma once

#include "Mesh_new/Layout.hpp"
#include <iostream>
#include <string>
#include <sstream>

// Forward-declare the FunctionSpace, those are factories (via friending) 
// for Forms.
template <typename Family, int k, int r, typename LayoutType>
class FunctionSpace;

namespace fem {

// --- Private Implementation Detail: ComponentField ---
namespace Detail {

// The ComponentInitData struct is a public "Data Transfer Object".
// It is defined here so both FunctionSpace (the factory) and
// Form (the product) can use it.
template<typename GroupTag>
struct ComponentInitData {
    size_t num_entities;
    int num_coeffs;
};

template <typename T, typename GroupTag, typename LayoutType>
class ComponentField {
public:

  using ScalarType = T;
  using GroupTagType = GroupTag;

  using view_type = Kokkos::View<T**>;

  explicit ComponentField(
    std::shared_ptr<const LayoutType> layout,
    size_t num_entities,
    size_t num_coeffs)
    : layout_(layout)
  {
      view_ = view_type("component_field", num_entities, num_coeffs);

      std::cout << "    - ComponentField for GroupTag<...>"
                << " allocated with flat size " << num_entities << " x " << num_coeffs << " coeffs"
                << std::endl;
  }

  ComponentField& operator=(T scalar) {
      Kokkos::deep_copy(view_, scalar);
      return *this;
  }

  void fillHalo() {
      layout_->fillHalo();
  }

  view_type& view() { return view_; }
  const view_type& view() const { return view_; }

private:
    std::shared_ptr<const LayoutType> layout_;
    view_type view_;
};

} // namespace Detail


template <typename... ComponentFieldTypes>
class Form {
private:
    using ComponentTuple = std::tuple<ComponentFieldTypes...>;
    ComponentTuple components_;

    // The constructor is PRIVATE. Only a FunctionSpace can create a Form.
    // It takes a pre-constructed tuple of components.
    explicit Form(ComponentTuple&& components) : components_(std::move(components)) {}

    template <typename Family, int r, typename T, typename LayoutType>
    friend class FunctionSpace;

public:

    void fillHalo() {
        std::cout << "Form: Orchestrating halo exchange for all components..." << std::endl;
        std::apply(
            [](auto&... component) {
                (component.fillHalo(), ...);
            },
            components_
        );
    }

    template <typename ComponentFieldType>
    const auto& get_component() const {
        return std::get<ComponentFieldType>(components_);
    }

    template <typename ComponentFieldType>
    auto& get_component() {
        return std::get<ComponentFieldType>(components_);
    }
};

} // namespace fem
