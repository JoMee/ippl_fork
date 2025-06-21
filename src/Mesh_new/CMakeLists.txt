
target_include_directories(ippl
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

install(FILES
    Blades.hpp
    Fields.hpp
    MeshPolicies.hpp
    Layout.hpp
    LayoutStrategy.hpp
    structured/IndexView.hpp
    structured/Connectivity.hpp
    structured/GrassmanIndex.hpp
    structured/StructuredMesh.hpp
    structured/Geometry.hpp
    DESTINATION include/Meshes_new
)
