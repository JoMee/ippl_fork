#pragma once
#include "mesh/concepts.hpp"
#include "mesh/field.hpp"

namespace mesh {

template<MeshModel M, typename T>
class HaloExchangePlan {
public:
    explicit HaloExchangePlan(const M& mesh) : mesh_(mesh) {}

    template<class PartitionT>
    void build(const PartitionT& part) {
        // build send/recv lists using part.mesh()
    }

    void exchange(Field<M,T>& fld) {
        // pack ‑> MPI ‑> unpack (CUDA‑aware or host bounce)
    }

private:
    const M& mesh_;
};

} // namespace mesh
