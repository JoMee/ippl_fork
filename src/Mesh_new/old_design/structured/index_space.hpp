template<int Dim>
class StructuredIndexSpace {
public:
    constexpr static int dimension = Dim;

    StructuredIndexSpace(std::array<int, Dim> lower,
                         std::array<int, Dim> upper,
                         std::array<bool, Dim> periodic = {});

    int lower(int d) const;
    int upper(int d) const;
    bool is_periodic(int d) const;
    std::array<int, Dim> global_offset() const;
    int extent(int d) const;

private:
    std::array<int, Dim> lower_, upper_;
    std::array<bool, Dim> periodic_;
};

