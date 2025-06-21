
namespace fem {
// A strategy for a single-rank, serial layout.
struct SerialLayoutStrategy {
    template <typename LayoutType>
    static void fill_halo(const LayoutType& /*layout*/) {
        std::cout << "  -> (SerialLayoutStrategy) Performing trivial halo exchange." << std::endl;
    }
};
}
