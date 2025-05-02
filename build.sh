#!/bin/bash

# Function to display script usage
function show_usage {
    echo "Usage: $0 [OPTION]"
    echo "Build options:"
    echo "  --test        Build and run tests"
    echo "  --bench       Build and run benchmarks"
    echo "  --clean       Clean build files"
    echo "  --docs        Generate documentation"
    echo "  --cmake       Build using CMake"
    echo "  No argument   Build the main program"
}

# Clean build directory
function clean_build {
    echo "Cleaning build directory..."
    rm -rf build
    mkdir -p build
}

# Define colors for output
YELLOW='\033[0;33m'
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check for pkg-config
if ! command -v pkg-config &> /dev/null; then
    echo -e "${YELLOW}Warning: pkg-config is not installed. Assuming GMP is available.${NC}"
    CFLAGS="-lgmp -lm"
else
    CFLAGS="$(pkg-config --cflags --libs gmp) -lm"
fi

# Common compiler flags
CFLAGS="-Wall -Wextra -g -O2 $CFLAGS -DDEFINE_GLOBAL_VARIABLES"

# Process arguments
if [ "$1" == "--test" ]; then
    echo "Building and running tests..."
    mkdir -p build
    gcc $CFLAGS -o build/test_fp tests/test_fp.c parameters.c -I. && \
    gcc $CFLAGS -o build/test_main tests/test_main.c parameters.c -I. && \
    echo -e "${GREEN}Running tests...${NC}" && \
    ./build/test_fp && \
    ./build/test_main
    exit $?
elif [ "$1" == "--bench" ]; then
    echo "Building and running benchmarks..."
    mkdir -p build
    gcc $CFLAGS -o build/benchmark bench/benchmark.c parameters.c -I. && \
    echo -e "${GREEN}Running benchmarks...${NC}" && \
    ./build/benchmark
    exit $?
elif [ "$1" == "--clean" ]; then
    clean_build
    echo -e "${GREEN}Build directory cleaned.${NC}"
    exit 0
elif [ "$1" == "--docs" ]; then
    echo "Generating documentation..."
    doxygen Doxyfile
    echo -e "${GREEN}Documentation generated in docs/ directory.${NC}"
    exit $?
elif [ "$1" == "--cmake" ]; then
    echo "Building with CMake..."
    mkdir -p build
    cd build
    cmake .. && make
    exit $?
elif [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    show_usage
    exit 0
fi

# Main build process
echo "Building Finite Extension Field Degree 18 Library..."
clean_build

# Compile everything in a single command to prevent multiple definition errors
echo "Compiling main program..."
gcc $CFLAGS -o build/fp18_arith \
    parameters.c \
    main.c \
    BN.c \
    degree18.c \
    embedding_degree18.c \
    -I. 

status=$?
if [ $status -eq 0 ]; then
    echo -e "${GREEN}Build successful. Run ./build/fp18_arith to execute the program.${NC}"
else
    echo -e "${RED}Build failed with status $status.${NC}"
fi

exit $status