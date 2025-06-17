#include <Kokkos_Core.hpp>

namespace fem {

template<typename Derived, typename T, typename Blade>
struct FieldBase {
  using value_type = T;
  using blade_type = Blade;
  using index_type = Derived::index_type;
  
  static constexpr int dim = Blade::Dim;

  KOKKOS_INLINE_FUNCTION
  auto operator()(const index_type& id) const {
    return static_cast<const Derived&>(*this)(id);
  }
};

} // namespace fem
