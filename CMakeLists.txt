cmake_minimum_required(VERSION 3.12)
project(alex)

set(CMAKE_CXX_STANDARD 14)

# Define the macro ‘DEBUG' in the debug mode
if(CMAKE_BUILD_TYPE STREQUAL Debug)        
    ADD_DEFINITIONS(-DDEBUG)               
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "/O2 /arch:AVX2 /W1 /EHsc")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-O3 -xHost")
else()
    # clang and gcc
    set(CMAKE_CXX_FLAGS "-march=native -Wall -Wextra")
endif()


include_directories(src/core)
include_directories(attack)

# directries where you place your Google OR-tool include and bianry
include_directories(/home/ubuntu/data/or-tools/include)
link_directories(/home/ubuntu/data/or-tools/lib)

add_executable(time_aca_blackbox src/benchmark/time_aca_blackbox.cpp)
add_executable(time_aca_graybox src/benchmark/time_aca_graybox.cpp)
add_executable(time_aca_whitebox src/benchmark/time_aca_whitebox.cpp)
add_executable(space_aca_in src/benchmark/space_aca_in.cpp)
add_executable(space_aca_dn src/benchmark/space_aca_dn.cpp)
add_executable(lis src/benchmark/lis.cpp)

target_link_libraries(time_aca_blackbox -lortools)
target_link_libraries(time_aca_graybox -lortools)
target_link_libraries(time_aca_whitebox -lortools)
target_link_libraries(space_aca_in -lortools)
target_link_libraries(space_aca_dn -lortools)
target_link_libraries(lis -lortools)


set(DOCTEST_DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/doctest)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
file(DOWNLOAD
    https://raw.githubusercontent.com/onqtam/doctest/2.4.6/doctest/doctest.h
    ${DOCTEST_DOWNLOAD_DIR}/doctest.h
    EXPECTED_HASH SHA1=40728d2bed7f074e80cb095844820114e7d16ce0
)

add_executable(test_alex test/unittest_main.cpp)
target_link_libraries(test_alex -lortools)
target_include_directories(test_alex PRIVATE ${DOCTEST_DOWNLOAD_DIR})

enable_testing()
add_test(test_alex test_alex)
