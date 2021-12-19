cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project (${PROJECT_LIBRARY_PROBLEMS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)


file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

add_library(${PROJECT_LIBRARY_PROBLEMS} STATIC ${SRC_FILES})
