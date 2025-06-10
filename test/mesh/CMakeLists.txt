# -----------------------------------------------------------------------------
# test/mesh_new/CMakeLists.txt
#
# Integration tests for kokkos-related functionality.
# -----------------------------------------------------------------------------

file(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "🔧 Adding field integration tests from: ${_relPath}")

add_ippl_integration_test(TestMesh LABELS integration mesh)
