cmake_minimum_required(VERSION 3.10)

project(FiniteExtensionFieldDeg18 C)

# Find GMP library
find_package(GMP REQUIRED)

# Include directories
include_directories(${GMP_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

# Add source files
# Adjust this list based on actual dependencies and how files include each other.
# If Finite_Field.c and Elliptic_Curve.c are included directly by embedding_degree18.c,
# they might not need to be listed here.
add_executable(main_executable
    main.c
    Finite_Field.c
    Elliptic_Curve.c
    embedding_degree18.c
    BN.c
    degree18.c
)

# Link GMP library
target_link_libraries(main_executable ${GMP_LIBRARIES})

# Optional: Set C standard if needed
# set(CMAKE_C_STANDARD 99)
# set(CMAKE_C_STANDARD_REQUIRED ON)

# Optional: Enable tests using CTest (requires test code)
# enable_testing()
# add_test(NAME arithmetic_tests COMMAND test_executable) # Assuming a test_executable target exists
