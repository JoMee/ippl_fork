#pragma once
#include "mesh/concepts.hpp"
#include <Kokkos_Core.hpp>

namespace mesh {

//---------- Index façade ------------------------------------------------------

template<MeshModel M>
struct IndexMap { static_assert(sizeof(M)==0, "IndexMap not specialised"); };

//---------- Geometry façade ---------------------------------------------------

template<MeshModel M>
struct Geometry  { static_assert(sizeof(M)==0, "Geometry not specialised"); };

//---------- Topology façade ---------------------------------------------------

template<MeshModel M>
struct Topology  { static_assert(sizeof(M)==0, "Topology not specialised"); };

// Helper aliases --------------------------------------------------------------

template<MeshModel M> using Global = typename M::GlobalIndex;
template<MeshModel M> using Local  = typename M::LocalIndex;
template<MeshModel M> using Vec    = typename M::Vec;

} // namespace mesh
