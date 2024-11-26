set(target "${PROJECT_LIBRARY}")

file(GLOB hdrs "${PROJECT_INCLUDE}/*.h*")
file(GLOB srcs "${PROJECT_SOURCE}/*.cpp")

find_package(Eigen3)
if (NOT Eigen3_FOUND)
  message(FATAL_ERROR "Package Eigen3 required, but not found!")
endif()

add_library(${target} STATIC ${srcs} ${hdrs})
target_link_libraries(${target} ${LIBRARY_DEPS} Eigen3::Eigen)

# REPORT
message( STATUS "")
message( STATUS "General configuration for lib_${target}")
message( STATUS "OK")
message( STATUS "")
