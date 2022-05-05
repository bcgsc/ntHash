cmake_minimum_required(VERSION 3.22)
project(ntHash)

set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(include)
include_directories(vendor/btllib/include)
include_directories(vendor/argparse/include)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/nthash/*.cpp)

add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a
        COMMAND ${PROJECT_SOURCE_DIR}/vendor/btllib/compile --prefix ${PROJECT_BINARY_DIR}/btllib
)

add_library(nthash_lib STATIC ${SRC_FILES})
set_target_properties(nthash_lib PROPERTIES OUTPUT_NAME "nthash")

add_executable(nthash main.cpp ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a)
target_link_libraries(nthash ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a)
target_link_libraries(nthash ${PROJECT_BINARY_DIR}/libnthash.a)

add_executable(nthash_tests ${PROJECT_SOURCE_DIR}/tests/tests.cpp ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a)
target_link_libraries(nthash_tests ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a)
target_link_libraries(nthash_tests ${PROJECT_BINARY_DIR}/libnthash.a)
enable_testing()
add_test(NAME nthash_tests COMMAND $<TARGET_FILE:nthash_tests>)
