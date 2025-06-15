# -----------------------------------------------------------------------------
# src/Meshes/CMakeLists.txt
#
# Adds structured grid mesh headers to the IPPL target.
# -----------------------------------------------------------------------------

target_include_directories(ippl
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Install Meshes headers
install(FILES
    comm_pattern.hpp
    comm.hpp
    field.hpp
    field_layout.hpp
    halo_exchange.hpp
    incidence.hpp
    index_map.hpp
    mesh.hpp
    mesh_concepts.hpp
    structured/layout.hpp
    structured/index_map.hpp
    partition/rcb_partitioner.hpp
    DESTINATION include/Meshes_new
)

