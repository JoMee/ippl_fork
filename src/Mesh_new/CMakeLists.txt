
target_include_directories(ippl
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

install(FILES
    BladeDispatch.hpp
    Blades.hpp
    GrassmanIndex.hpp
    IndexView.hpp
    StructuredMesh.hpp
    Fields.hpp
    Mesh.hpp
    SerialLayout.hpp
    Connectivity.hpp
    DESTINATION include/Meshes_new
)
