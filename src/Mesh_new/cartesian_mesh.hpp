#pragma once
#include "mesh/api.hpp"
#include <array>

namespace models {

struct CartesianMesh3D {
    // Required type aliases
    struct GlobalIndex { int i,j,k; };
    using LocalIndex = std::size_t;
    using Vec        = std::array<double,3>;
    using Scalar     = double;
    static constexpr int Dim = 3;

    // ctor
    CartesianMesh3D(Vec origin, Vec h, std::array<int,3> dims)
      : origin_(origin), h_(h), dims_(dims) {}

    // data
    Vec origin_, h_; std::array<int,3> dims_;
};

// -------- IndexMap specialisation -------------------------------------------

template<>
struct mesh::IndexMap<CartesianMesh3D> {
    using B = CartesianMesh3D; using G = B::GlobalIndex; using L = B::LocalIndex;
    const B* m;

    KOKKOS_FUNCTION std::optional<L> global_to_local(G g) const {
        if(g.i<0||g.i>=m->dims_[0]||g.j<0||g.j>=m->dims_[1]||g.k<0||g.k>=m->dims_[2])
            return std::nullopt;
        return g.k*m->dims_[0]*m->dims_[1] + g.j*m->dims_[0] + g.i;
    }
    KOKKOS_FUNCTION G local_to_global(L l) const {
        int nx=m->dims_[0], ny=m->dims_[1];
        int i =  l % nx;
        int j = (l/nx) % ny;
        int k =  l/(nx*ny);
        return {i,j,k};
    }
    KOKKOS_FUNCTION bool owns(G)   const { return true;  }
    KOKKOS_FUNCTION bool is_ghost(G)const{ return false; }
};

// -------- Geometry specialisation -------------------------------------------

template<>
struct mesh::Geometry<CartesianMesh3D> {
    using B = CartesianMesh3D; using G = B::GlobalIndex; using Vec = B::Vec;
    const B* m;
    KOKKOS_FUNCTION Vec    centroid(G g) const {
        return { m->origin_[0] + (g.i+0.5)*m->h_[0],
                 m->origin_[1] + (g.j+0.5)*m->h_[1],
                 m->origin_[2] + (g.k+0.5)*m->h_[2] };
    }
    KOKKOS_FUNCTION double measure(G) const {
        return m->h_[0]*m->h_[1]*m->h_[2];
    }
};

// -------- Topology specialisation (6‑pt neighbour) ---------------------------

template<>
struct mesh::Topology<CartesianMesh3D> {
    using B = CartesianMesh3D; using G = B::GlobalIndex;
    const B* m;
    KOKKOS_FUNCTION void incident(G c,int dim,G* out,std::size_t max,std::size_t& n) const {
        n = 0;
        if(dim!=3) return; // only cell neighbours here
        const int di[6]={1,-1,0,0,0,0};
        const int dj[6]={0,0,1,-1,0,0};
        const int dk[6]={0,0,0,0,1,-1};
        mesh::IndexMap<B> idx{m};
        for(int s=0;s<6 && n<max;++s){
            G nb{c.i+di[s],c.j+dj[s],c.k+dk[s]};
            if(idx.global_to_local(nb)) out[n++]=nb;
        }
    }
};

} // namespace models
