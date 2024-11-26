cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME global_optimization_methods)
project(${PROJECT_NAME})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin/release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin/debug)

set(PROJECT_LIBRARY "${PROJECT_NAME}")
set(PROJECT_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(PROJECT_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src")

include_directories("${PROJECT_INCLUDE}")

# OUTPUT DIR
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set(OUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})
else()
  set(OUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
endif()

# COPY SCRIPT
set(COPY_RESOURCE copy_resource)
add_custom_target(${COPY_RESOURCE} ALL)
add_custom_command(
  TARGET ${COPY_RESOURCE} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_SOURCE_DIR}/scripts/
  ${OUT_DIR}/scripts
)

# MAKE DIR
set(MAKE_DIR make_directory)
add_custom_target(${MAKE_DIR} ALL)
add_custom_command(
  TARGET ${MAKE_DIR} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory
  ${OUT_DIR}/output_graph

  ${OUT_DIR}/output_data

  ${OUT_DIR}/output_data/scanning_test
  ${OUT_DIR}/output_data/piyavsky_test
  ${OUT_DIR}/output_data/gsa_test
  ${OUT_DIR}/output_data/imgo_test
  ${OUT_DIR}/output_data/mggsa_test

  ${OUT_DIR}/output_data/onedimensional_operational_characteristics
  ${OUT_DIR}/output_data/onedimensional_operational_characteristics/scanning
  ${OUT_DIR}/output_data/onedimensional_operational_characteristics/piyavsky
  ${OUT_DIR}/output_data/onedimensional_operational_characteristics/gsa

  ${OUT_DIR}/output_data/multidimensional_constrained_operational_characteristics
  ${OUT_DIR}/output_data/multidimensional_constrained_operational_characteristics/mggsa

  ${OUT_DIR}/output_data/new_model_class/sample_fitting_problem
  ${OUT_DIR}/output_data/new_model_class/fitting_family_problems_operational_characteristics/mggsa
  ${OUT_DIR}/output_data/new_model_class/fitting_family_problems_plot
  ${OUT_DIR}/output_data/new_model_class/fitting_family_problems_search_solution
  ${OUT_DIR}/output_data/new_model_class/fitting_family_problems/slices
  ${OUT_DIR}/output_data/new_model_class/fitting_family_problems/experiment

  ${OUT_DIR}/output_data/mggsa_operational_characteristics
)

# GPROF
SET(GPROF FALSE)
if (GPROF)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

# GCGEN
add_subdirectory(GCGen)
include_directories(${GCGEN_INC_NAME})

# DIRECT
add_subdirectory(DiRect)
include_directories(${DIRECT_INC_NAME})

# BUILD
add_subdirectory(src)
add_subdirectory(samples)

# REPORT
message(STATUS "")
message(STATUS "General configuration for ${PROJECT_NAME}")
message(STATUS "======================================")
message(STATUS "")
message(STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
