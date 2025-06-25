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
  using view_type = Kokkos::View<T**>;

  explicit ComponentField(
    std::shared_ptr<const LayoutType> layout,
    size_t num_entities,
    int num_coeffs)
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

  // This method can now be specialized for different GroupTags in the future.
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
    template <typename> struct ComponentTupleFromGroups;
    template <typename... GroupTags>
    struct ComponentTupleFromGroups<std::tuple<GroupTags...>> {
        using type = std::tuple<Detail::ComponentField<T, GroupTags, LayoutType>...>;
    };

    // The policy defines what the GroupTags are for a given mesh type.
    using GroupTagTuple = typename LayoutType::Policy::template StorageModel<LayoutType::DIM, k>::GroupTagTuple;
    using ComponentTuple = typename ComponentTupleFromGroups<GroupTagTuple>::type;

    std::shared_ptr<const LayoutType> layout_;
    ComponentTuple components_;

    // The constructor is PRIVATE and takes the public DTO.
    template <typename... GroupTags>
    Form(std::shared_ptr<const LayoutType> layout,
         const std::tuple<Detail::ComponentInitData<GroupTags>...>& all_init_data)
        : layout_(layout),
          components_(std::make_tuple(
              Detail::ComponentField<T, GroupTags, LayoutType>(
                  layout,
                  std::get<Detail::ComponentInitData<GroupTags>>(all_init_data).num_entities,
                  std::get<Detail::ComponentInitData<GroupTags>>(all_init_data).num_coeffs
              )...
          ))
    {}

    // Declare all FunctionSpace instantiations as friends so they can call the private constructor.
    template <typename Family, int k_friend, int r, typename T_friend, typename LayoutType_friend>
    friend class FunctionSpace;

public:
    // The public interface for Form remains simple.
    void fillHalo() {
        std::cout << "Form<" << k << ">: Orchestrating halo exchange..." << std::endl;
        std::apply([](auto&... component) { (component.fillHalo(), ...); }, components_);
    }

    template <typename GroupTag>
    const auto& get_component() const {
        return std::get<Detail::ComponentField<T, GroupTag, LayoutType>>(components_);
    }
};

} // namespace fem
