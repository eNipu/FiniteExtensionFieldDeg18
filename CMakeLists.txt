cmake_minimum_required(VERSION 3.10)
project(FiniteExtensionFieldDeg18 
        VERSION 1.0.0
        DESCRIPTION "Finite Extension Field Degree 18 Library" 
        LANGUAGES C)

# Set C standard to C11
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Add compiler flags
add_compile_options(-Wall -Wextra -O2)

# Find GMP library
find_library(GMP_LIBRARY gmp REQUIRED)
if(NOT GMP_LIBRARY)
    message(FATAL_ERROR "GMP library not found")
endif()

# Source files - exclude main.c for library target
file(GLOB_RECURSE SOURCES 
    "*.c"
)
list(FILTER SOURCES EXCLUDE REGEX ".*main.c$")

# Main executable
add_executable(fp18_arith main.c ${SOURCES})
target_link_libraries(fp18_arith ${GMP_LIBRARY} m)

# Library target
add_library(fp18 STATIC ${SOURCES})
target_link_libraries(fp18 ${GMP_LIBRARY} m)

# Test directory and executable
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
file(GLOB_RECURSE TEST_SOURCES "tests/*.c")
if(TEST_SOURCES)
    add_executable(run_tests ${TEST_SOURCES} ${SOURCES})
    target_link_libraries(run_tests ${GMP_LIBRARY} m)
    add_custom_target(test COMMAND run_tests)
endif()

# Benchmark directory and executable
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bench)
file(GLOB_RECURSE BENCH_SOURCES "bench/*.c")
if(BENCH_SOURCES)
    add_executable(run_benchmarks ${BENCH_SOURCES} ${SOURCES})
    target_link_libraries(run_benchmarks ${GMP_LIBRARY} m)
    add_custom_target(bench COMMAND run_benchmarks)
endif()

# Documentation via Doxygen (optional)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
    )
endif()

# Additional library and executables
add_library(finitefield STATIC
    fp.c fp3.c fp6.c fp18.c ec.c pairing.c parameters.c
)
target_include_directories(finitefield PUBLIC .)
target_link_libraries(finitefield PUBLIC gmp)

add_executable(main main.c)
target_link_libraries(main PRIVATE finitefield)

add_executable(test_fp tests/test_fp.c)
target_link_libraries(test_fp PRIVATE finitefield)

add_executable(test_main tests/test_main.c)
target_link_libraries(test_main PRIVATE finitefield)

add_executable(benchmark bench/benchmark.c)
target_link_libraries(benchmark PRIVATE finitefield)
