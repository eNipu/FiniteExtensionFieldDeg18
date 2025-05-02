# Finite Extension Field Degree 18 Library

A high-performance C library implementing arithmetic for finite extension fields of degree 18 (`𝔽p18`), elliptic curve operations, and pairing computations. This library is designed for cryptographic applications, particularly for pairing-based cryptography.

## Features

- Prime field (`𝔽p`) arithmetic
- Extension field arithmetic for `𝔽p3`, `𝔽p6`, and `𝔽p18`
- Elliptic curve operations over all field levels
- Optimal Ate pairing implementation
- Benchmarking tools for performance analysis
- Comprehensive test suite
- Efficient implementation using GMP for arbitrary-precision arithmetic

## Tower of Extension Fields

This library implements a tower of extension fields with the following structure:

1. Base prime field `𝔽p`
2. Cubic extension `𝔽p3 = 𝔽p[ω]/(ω³ - C₁)` where C₁ = 2
3. Quadratic extension `𝔽p6 = 𝔽p3[τ]/(τ² - ξ)` where ξ = ω
4. Cubic extension `𝔽p18 = 𝔽p6[v]/(v³ - τ)`

## Requirements

- GMP (GNU Multiple Precision Arithmetic Library)
- C compiler with C11 support (GCC or Clang recommended)
- CMake 3.10+ or Make for building

## Building

### Using Make

```bash
# Build the library and executable
make

# Run tests
make test

# Run benchmarks
make bench

# Generate documentation
make docs

# Clean build artifacts
make clean
```

### Using CMake

```bash
# Create build directory
mkdir -p build && cd build

# Configure and build
cmake ..
cmake --build .

# Run tests
ctest

# Run the executable
./fp18_arith
```

## Development Environment

This project includes a DevContainer configuration for Visual Studio Code, providing a consistent development environment with all necessary tools pre-installed.

To use the DevContainer:

1. Install [Visual Studio Code](https://code.visualstudio.com/)
2. Install the [Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension
3. Clone this repository
4. Open the repository in VS Code
5. When prompted, click "Reopen in Container"

## Project Structure

- `fp18_arith.h` - Main header file with all function declarations
- `Finite_Field.c` - Implementation of finite field arithmetic
- `Elliptic_Curve.c` - Implementation of elliptic curve operations
- `embedding_degree18.c` - Implementation of pairing operations
- `parameters.c` - Parameter management
- `main.c` - Example usage and demos
- `tests/` - Test suite
- `bench/` - Benchmarking tools

## Usage Example

```c
#include "fp18_arith.h"
#include <stdio.h>

int main() {
    // Initialize parameters
    init_parameters();
    
    // Set the generator value X
    mpz_set_str(X, "18446893747415302274", 10);
    
    // Generate curve parameters based on X
    generate_parameters();
    
    // Create two Fp elements
    struct Fp a, b, c;
    Fp_init(&a);
    Fp_init(&b);
    Fp_init(&c);
    
    // Set values
    Fp_set_ui(&a, 123);
    Fp_set_ui(&b, 456);
    
    // Perform addition
    Fp_add(&c, &a, &b);
    
    // Print result
    printf("Result of addition: ");
    Fp_printf(&c);
    printf("\n");
    
    // Clean up
    Fp_clear(&a);
    Fp_clear(&b);
    Fp_clear(&c);
    clear_parameters();
    
    return 0;
}
```

## Running Tests

The test suite ensures the correctness of the library's operations:

```bash
make test
```

## Performance Benchmarks

The library includes benchmarking tools to measure the performance of various operations:

```bash
make bench
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch: `git checkout -b feature/amazing-feature`
3. Commit your changes: `git commit -am 'Add some amazing feature'`
4. Push to the branch: `git push origin feature/amazing-feature`
5. Submit a pull request

## License

This project is licensed under the terms specified in the `LICENSE` file.

## References

- [The GMP Library](https://gmplib.org/)
- [Pairings for Beginners](https://www.craigcostello.com.au/s/PairingsForBeginners.pdf) by Craig Costello
- [BN Curves](https://tools.ietf.org/id/draft-kasamatsu-bncurves-01.html)
