
target_include_directories(ippl
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

install(FILES
    index_space.h
    structured/Blades.hpp
    structured/GrassmanIndex.hpp
    DESTINATION include/Meshes_new
)
