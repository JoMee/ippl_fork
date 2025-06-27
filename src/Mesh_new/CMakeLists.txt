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
    api.hpp
    cartesian_mesh.hpp
    concepts.hpp
    field.hpp
    halo_exchange.hpp
    DESTINATION include/Meshes_new
)

