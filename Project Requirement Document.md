# Project Requirements: Finite Extension Field Degree 18 Library

## 1. Overview

This project implements efficient and maintainable arithmetic for the finite extension field Fp^18, leveraging the GMP (GNU Multiple Precision Arithmetic) library for large-integer computations. The library should follow professional C coding standards, with clear abstractions, modular design, and thorough testing.

## 2. Functional Requirements

1. **Field Arithmetic Modules**

   * **Fp (base field):** Addition, subtraction, multiplication, inversion, exponentiation modulo p.
   * **Fp3 (cubic extension):** Represent elements as triplets; support addition, subtraction, multiplication (Karatsuba), inversion, exponentiation.
   * **Fp6 (sextic extension):** Built over Fp3; implement optimal multiplication and inversion.
   * **Fp18 (degree-18 extension):** Built over Fp6; implement optimized multiplication, inversion, and powering.

2. **Elliptic Curve Operations**

   * Define a generic elliptic curve structure over Fp and Fp\[n] extensions.
   * Implement point addition, doubling, and scalar multiplication (using windowed methods) with clear interfaces.

3. **Pairing Computations** (Optional)

   * Support at least one pairing (e.g., Optimal Ate) on pairing-friendly curves defined over Fp18.

4. **Testing and Verification**

   * Comprehensive unit tests for each field and curve operation.
   * Property-based tests (e.g., associativity, distributivity, inversion correctness).
   * Provide example programs demonstrating usage.

## 3. Non-Functional Requirements

* **Performance:** Use algorithmic optimizations (e.g., Karatsuba, Montgomery representation) where appropriate. Benchmark critical routines.
* **Code Quality:**

  * Follow C11 standard; use `-Wall -Wextra -Werror` in build.
  * Use consistent naming conventions (snake\_case for functions and variables, uppercase for macros).
  * Proper memory management: initialize and clear GMP variables; avoid leaks.
  * Write clear comments and documentation blocks (Doxygen-compatible).
* **Modularity:** Organize code into separate modules and headers per field level (fp.c/h, fp3.c/h, fp6.c/h, fp18.c/h, ec.c/h, pairing.c/h).
* **Dependencies:** Only GMP; provide clear abstraction to swap out big-int backends if needed.
* **Portability:** Ensure code compiles on Linux and BSD systems; provide a dev-container setup.

## 4. Build and CI

* **Build Systems:**

  * Provide both a Makefile and CMakeLists.txt.
  * Targets: `all`, `clean`, `test`, `bench`.
* **Continuous Integration:** GitHub Actions workflow to build, run tests, and benchmarks on push and PR.
* **Dev Container:** Provide `.devcontainer/` configuration for VSCode.

## 5. Documentation and Maintainer Experience

* **README.md:** Project overview, installation, build instructions, examples.
* **CONTRIBUTING.md:** Guidelines for code style, testing, and pull requests.
* **LICENSE:** MIT.
* **Issue Templates:** Bug report and feature request.
* **GitHub Labels:** Use consistent labels for issues and PRs.
