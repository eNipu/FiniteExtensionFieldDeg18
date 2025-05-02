# GitHub Copilot Instructions for FiniteExtensionFieldDeg18

## Project Overview

This library implements robust and high-performance arithmetic for Fp^18 using C and GMP. The codebase must be modular, well-documented, and optimized for production use.

## Key Technologies and Standards

* **Language:** C (C11)
* **Big-Integer Backend:** GMP (`mpz_` API)
* **Build Tools:** GNU Make and CMake
* **Testing:** Unit tests with a framework (e.g., Check) or custom test harness; benchmarks with Google Benchmark or custom timing.
* **CI/CD:** GitHub Actions
* **Containerization:** Dev container for development environment

## Coding Conventions

* **Naming:** snake\_case for functions/variables, UPPER\_SNAKE\_CASE for macros, StructTypes in PascalCase.
* **Headers:** Use include guards or `#pragma once`; group related declarations.
* **Documentation:** Write Doxygen comments for all public functions/types.
* **Error Handling:** Return error codes; assert on unrecoverable errors.
* **Memory:** `mpz_init`/`mpz_clear`; use RAII-like patterns where possible.
* **Performance:** Use Karatsuba, Montgomery reduction; inline critical functions with `static inline`.

## File Structure and Modules

* `src/`

  * `fp.c/h`, `fp3.c/h`, `fp6.c/h`, `fp18.c/h`
  * `ec.c/h` (elliptic curve)
  * `pairing.c/h` (optional)
* `include/` for public headers
* `tests/` for test cases
* `bench/` for benchmarks
* `examples/` for usage demos
* `tools/` for utility scripts
* `.devcontainer/` for VSCode dev environment

## Goals for Copilot Assistance

1. **Refactor and Rename**

   * Propose clear module and file names.
   * Standardize function and variable names per conventions.
2. **Optimize Implementation**

   * Suggest algorithmic improvements (e.g., Karatsuba splits, Montgomery arithmetic).
   * Inline and unroll critical loops.
3. **Build Configuration**

   * Generate Makefile and CMakeLists.txt with standard targets.
4. **Testing Setup**

   * Scaffold unit tests and benchmarks.
5. **Documentation**

   * Generate Doxygen comments and README sections.
