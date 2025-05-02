# Finite Extension Field of Degree 18 (Fp^18)

This repository provides a C implementation of arithmetic operations within the finite extension field Fp^18, built upon the GMP (GNU Multiple Precision Arithmetic) library. It includes functionalities for base field (Fp), cubic extension (Fp^3), sextic extension (Fp^6), and the full Fp^18 extension. Elliptic curve operations and pairing functionalities might also be included or under development.

## Features

*   Arithmetic operations (add, sub, mul, inv, pow) for Fp, Fp^3, Fp^6, Fp^18.
*   Elliptic Curve Arithmetic (ECA, ECD, SCM) over relevant fields.
*   Pairing computations (specific algorithms may vary).
*   Utilizes GMP for all large integer calculations.
*   Build support via Makefile and CMake.
*   Basic test structure (under development).

## Prerequisites

*   **GCC** (or compatible C compiler)
*   **Make**
*   **CMake** (optional, alternative build system)
*   **GMP Library** (including development headers)
    *   On Debian/Ubuntu: `sudo apt-get update && sudo apt-get install libgmp-dev`
    *   On Fedora: `sudo dnf install gmp-devel`
    *   On macOS (using Homebrew): `brew install gmp`

## Building

### Using Makefile

1.  **Compile:**
    ```bash
    make
    ```
    This will create an executable named `main.out` (or similar, depending on the Makefile).

2.  **Clean:**
    ```bash
    make clean
    ```

### Using CMake

1.  **Configure:**
    ```bash
    cmake -S . -B build
    ```

2.  **Build:**
    ```bash
    cmake --build build
    ```
    This will create an executable in the `build/` directory.

3.  **Clean:**
    ```bash
    rm -rf build
    ```

## Running

After building, execute the main program:

```bash
./main.out
# or if using CMake
./build/main_executable # Adjust executable name based on CMakeLists.txt
```

## Testing

(Instructions for running tests will be added here once the test suite is implemented.)

```bash
# Example placeholder
./run_tests.sh
```

## Project Structure

```
.
├── BN.c                  # Big Number related functions (likely using GMP)
├── BN.h
├── degree18.c            # Potentially parameter generation or specific degree 18 functions
├── Elliptic_Curve.c      # Elliptic curve operations
├── embedding_degree18.c  # Main implementation file or pairing specific code
├── embedding_degree18.h  # Header for embedding_degree18.c
├── f18.h                 # Header potentially defining Fp18 structures/functions
├── Finite_Field.c        # Finite field arithmetic implementations (Fp, Fp3, Fp6, Fp18)
├── main.c                # Main executable source, example usage
├── README.md             # This file
├── LICENSE               # Project license (MIT)
├── Makefile              # Makefile build script
├── CMakeLists.txt        # CMake build script
├── .gitignore            # Git ignore rules
├── Project Requirement Document.md # Project requirements
├── .github/              # GitHub specific files
│   ├── copilot.md        # Instructions for GitHub Copilot
│   ├── ISSUE_TEMPLATE/
│   │   └── bug_report.md # Bug report template
│   └── workflows/
│       └── ci.yml        # GitHub Actions CI workflow
└── tests/                # Test suite directory (placeholder)
    └── test_main.c       # Placeholder test file
```

## Contributing

Please refer to the issue tracker and consider creating a bug report or feature request. Pull requests are welcome.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
