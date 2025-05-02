# GitHub Copilot Instructions for FiniteExtensionFieldDeg18

## Project Overview

This project implements arithmetic for the finite extension field Fp¹⁸ using C and GMP. Key components include operations in Fp, Fp³, Fp⁶, Fp¹⁸, elliptic curve ops, and pairings.

## Key Technologies

- **Language:** C11  
- **Library:** GMP (`mpz_*` for big integers)  
- **Build:** Makefile, CMake  
- **Testing:** Custom C test suite + benchmarks  
- **CI/CD:** GitHub Actions  
- **Env:** Linux (dev container)

## Coding Conventions

- C11 with `-Wall -Wextra -Werror`.  
- `snake_case` for functions/variables, `UPPER_SNAKE` for macros.  
- `mpz_init`/`mpz_clear`; no leaks.  
- Header guards.  
- Doxygen comments for public APIs.

## GMP Usage

- Initialize `mpz_t` before use; clear when done.  
- Use `mpz_add`, `mpz_sub`, `mpz_mul`, `mpz_invert`, `mpz_powm`, `mpz_mod`, etc.

## Build and Testing

- Make targets: `make all`, `make test`, `make bench`, `make clean`.  
- CMake alias: `cmake --build . --target <name>`.  
- Link with `-lgmp`.

## Goals for Copilot Assistance

1. **Refactoring & Modularization**  
   - Help refactor existing C code into separate modules (fp, fp3, fp6, fp18, ec, pairing).  
   - Enforce naming and memory-management best practices.

2. **Parameter Script**  
   - Assist in writing a Python tool to compute and verify 128-bit KSS18 parameters:  
     - Compute *p* and *r* for u = 2^44 + 2^22 − 2^9 + 2.  
     - Primality tests.  
     - Emit C header (`kss18_params.h`) with `#define P …`, `#define R …`, `#define B 3`, `#define DELTA 3`.

3. **C Implementation & Verification**  
   - Generate pure-C code to:  
     1. Initialize GF(p) constants and do Fp, Fp³, Fp⁶, Fp¹⁸ arithmetic.  
     2. Define and operate on E: y² = x³ + 3 and its sextic twist.  
     3. Build GF(p¹⁸) via polynomial X¹⁸ − 3.  
     4. Implement Optimal Ate pairing.

4. **Test & Benchmark Generation**  
   - Draft C test vectors using outputs from the Python script.  
   - Scaffold `test/` files to compare C results against reference data.  
   - Create `bench/` harnesses for field, EC, and pairing routines.

5. **Documentation & Comments**  
   - Generate Doxygen skeletons for new modules.  
   - Write usage examples in README.

6. **CI Workflow**  
   - Define GitHub Actions steps to run `make test` and `make bench`.  
   - Fail on performance regressions (optional).

---

**Note:** Keep the existing tech stack (C, GMP, Make/CMake, GitHub Actions, dev container) unchanged.
