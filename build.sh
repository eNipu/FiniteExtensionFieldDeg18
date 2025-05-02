#!/bin/bash

# Build script for KSS degree 18 pairing implementation

set -e

# Default options
USE_CMAKE=false
CLEAN=false
BUILD_TESTS=false
BUILD_BENCH=false
BUILD_DOCS=false

# Parse command line options
for arg in "$@"
do
    case $arg in
        --cmake)
        USE_CMAKE=true
        ;;
        --clean)
        CLEAN=true
        ;;
        --test)
        BUILD_TESTS=true
        ;;
        --bench)
        BUILD_BENCH=true
        ;;
        --docs)
        BUILD_DOCS=true
        ;;
    esac
done

# Clean if requested
if [ "$CLEAN" = true ]; then
    echo "Cleaning build artifacts..."
    rm -f *.o pairing_deg18 test_runner benchmark
    rm -rf build
    exit 0
fi

# Build documentation if requested
if [ "$BUILD_DOCS" = true ]; then
    echo "Generating documentation..."
    doxygen Doxyfile
    exit 0
fi

# Build using CMake if requested
if [ "$USE_CMAKE" = true ]; then
    echo "Building with CMake..."
    mkdir -p build
    cd build
    cmake ..
    make
    cd ..
    exit 0
fi

# Standard build process
echo "Building KSS degree 18 pairing implementation..."

# Compile parameters module
echo "- Compiling parameters..."
gcc -c -g -Wall -Wextra parameters.c -lgmp

# Compile modular arithmetic modules
echo "- Compiling modular arithmetic modules..."
gcc -c -g -Wall -Wextra fp.c -lgmp
gcc -c -g -Wall -Wextra fp3.c -lgmp
gcc -c -g -Wall -Wextra fp6.c -lgmp
gcc -c -g -Wall -Wextra fp18.c -lgmp

# Compile elliptic curve module
echo "- Compiling elliptic curve module..."
gcc -c -g -Wall -Wextra ec.c -lgmp

# Compile pairing module
echo "- Compiling pairing module..."
gcc -c -g -Wall -Wextra pairing.c -lgmp

# Compile main program
echo "- Compiling main program..."
gcc -c -g -Wall -Wextra main.c -lgmp

# Link everything together
echo "- Linking..."
gcc -o pairing_deg18 main.o parameters.o fp.o fp3.o fp6.o fp18.o ec.o pairing.o -lgmp

# Build tests if requested
if [ "$BUILD_TESTS" = true ]; then
    echo "Building tests..."
    gcc -c -g -Wall -Wextra -I. tests/test_main.c -lgmp
    gcc -c -g -Wall -Wextra -I. tests/test_fp.c -lgmp
    gcc -o test_runner test_main.o test_fp.o parameters.o fp.o fp3.o fp6.o fp18.o ec.o pairing.o -lgmp
    echo "Tests built successfully. Run with: ./test_runner"
fi

# Build benchmark if requested
if [ "$BUILD_BENCH" = true ]; then
    echo "Building benchmarks..."
    gcc -c -g -Wall -Wextra bench/benchmark.c -lgmp
    gcc -o benchmark benchmark.o parameters.o fp.o fp3.o fp6.o fp18.o ec.o pairing.o -lgmp
    echo "Benchmark built successfully. Run with: ./benchmark"
fi

echo "Build completed successfully."