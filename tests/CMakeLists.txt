################################################################################
# Tests
################################################################################

add_executable(far_ccd_tests
  test_edge_edge_ccd.cpp
  test_point_triangle_ccd.cpp
  test_dataset.cpp
)

target_include_directories(far_ccd_tests PUBLIC ".")

################################################################################
# Required Libraries
################################################################################

target_link_libraries(far_ccd_tests PUBLIC far_ccd::far_ccd)

include(catch2)
target_link_libraries(far_ccd_tests PUBLIC Catch2::Catch2WithMain)

include(ccd_query_io)
target_link_libraries(far_ccd_tests PUBLIC ccd_io::ccd_io)

# Extra warnings (link last for highest priority)
include(far_ccd_warnings)
target_link_libraries(far_ccd_tests PRIVATE far_ccd::warnings)

################################################################################
# Compiler options
################################################################################

target_compile_definitions(far_ccd_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(far_ccd_tests)
