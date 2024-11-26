file(GLOB_RECURSE samples_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

find_package(OpenMP)
if (NOT OpenMP_FOUND OR NOT OpenMP_CXX_FOUND)
  message(FATAL_ERROR "Package OpenMP required, but not found!")
endif()

find_package(NLopt)
if (NOT NLopt_FOUND)
  message(FATAL_ERROR "Package NLopt required, but not found!")
endif()

foreach(sample_filename ${samples_list})
  get_filename_component(sample ${sample_filename} NAME_WE)

  add_executable(${sample} ${sample_filename})
  target_link_libraries(${sample} ${PROJECT_LIBRARY} ${GCGEN_LIB_NAME} ${GCGEN_LIB_FAM_NAME} ${DIRECT_LIB_NAME})
  target_link_libraries(${sample} OpenMP::OpenMP_CXX nlopt)
  set_target_properties(${sample} PROPERTIES
    OUTPUT_NAME "${sample}"
    PROJECT_LABEL "${sample}"
    RUNTIME_OUTPUT_DIRECTORY "../"
  )

  # REPORT
  message( STATUS "")
  message( STATUS "General configuration for ${sample}")
  message( STATUS "OK")
  message( STATUS "")

endforeach()
