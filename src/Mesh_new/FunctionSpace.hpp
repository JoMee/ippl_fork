#pragma once

#include "Mesh_new/SpaceTraits.hpp"
#include "Mesh_new/Layout.hpp"

namespace fem {

template <typename Family, int r, typename T, typename LayoutType>
class FunctionSpace {
private:
    using Traits = SpaceTraits<Family, r, T, LayoutType>;
    using ComponentTuple = typename Traits::ComponentTuple;

    template<typename... ComponentTypes>
    typename Traits::FormType create_form_impl(std::tuple<ComponentTypes...>* /*dummy_for_deduction*/) const {

        auto initialized_components = std::make_tuple(
            ComponentTypes(
                layout_,
                layout_->get_indexer().template get_num_allocated_entities<typename ComponentTypes::GroupTagType>(),
                Traits::template get_num_dofs_per_entity<typename ComponentTypes::GroupTagType>()
            )...
        );

        // Call the private Form constructor with the initialized tuple.
        return typename Traits::FormType(std::move(initialized_components));
    }
public:
    using FormType = typename Traits::FormType;


    explicit FunctionSpace(std::shared_ptr<const LayoutType> layout)
        : layout_(layout)
    {}

    FormType create_form() const {
      return create_form_impl(static_cast<ComponentTuple*>(nullptr));
    }

private:
    std::shared_ptr<const LayoutType> layout_;
};

} // namespace fem

