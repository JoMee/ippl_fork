#pragma once

#include "Mesh_new/Fields.hpp"
#include "Mesh_new/Layout.hpp"
#include "Mesh_new/PolynomialSpace.hpp"

namespace fem {

/**
 * @brief Represents a discrete function space over a mesh.
 *
 * This is the primary user-facing class for creating and interacting with
 * discrete fields (Forms). It bundles the mesh/layout information with the
 * mathematical properties of the finite element (family, degree, form grade).
 *
 * @tparam Family The element family (e.g., P_r_Lambda, Q_r_Family).
 * @tparam k The form grade (0-form, 1-form, etc.).
 * @tparam r The polynomial degree.
 * @tparam T The scalar type of the field (e.g., double).
 * @tparam LayoutType The type of the underlying mesh layout.
 */
template <typename Family, int k, int r, typename T, typename LayoutType>
class FunctionSpace {
public:
    using FormType = Form<k, T, LayoutType>;
    using PolySpaceType = PolynomialSpace<Family, r, LayoutType::DIM>;

    // The constructor takes the Layout, which holds all the mesh/indexer info.
    explicit FunctionSpace(std::shared_ptr<const LayoutType> layout)
        : layout_(layout)
    {}

    /**
     * @brief Factory method to create a new Form (a discrete field)
     * belonging to this function space.
     *
     * @return A new Form object, correctly initialized.
     */
    FormType create_form() const {
        using GroupTagTuple = typename LayoutType::Policy::template StorageModel<LayoutType::DIM, k>::GroupTagTuple;

        PolySpaceType poly_space;

        auto init_data_tuple = std::apply(
            [&](auto... tags) {
                return std::make_tuple(
                    Detail::ComponentInitData<decltype(tags)>{
                        .num_entities = layout_->get_indexer().template get_num_allocated_entities<decltype(tags)>(),
                        .num_coeffs = poly_space.get_num_coefficients(ElementShape::Line) // Placeholder shape
                    }...
                );
            }, GroupTagTuple{});

        return FormType(layout_, init_data_tuple);
    }

private:
    std::shared_ptr<const LayoutType> layout_;
};

} // namespace fem

