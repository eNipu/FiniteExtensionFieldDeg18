# Project Requirements: Finite Extension Field Degree 18 Library

## 1. Overview

This project implements efficient and maintainable arithmetic for the finite extension field Fp^18, leveraging the GMP (GNU Multiple Precision Arithmetic) library for large-integer computations. The library should follow professional C coding standards, with clear abstractions, modular design, and thorough testing.

## 2. Functional Requirements

1. **Field Arithmetic Modules**
   - **Fp (base field):** Addition, subtraction, multiplication, inversion, exponentiation modulo p.
   - **Fp³ (cubic extension):** Represent elements as triplets; support addition, subtraction, multiplication (Karatsuba), inversion, exponentiation.
   - **Fp⁶ (sextic extension):** Built over Fp³; implement optimal multiplication and inversion.
   - **Fp¹⁸ (degree-18 extension):** Built over Fp⁶; implement optimized multiplication, inversion, and powering.

2. **Elliptic Curve Operations**
   - Define a generic elliptic curve structure over Fp and Fpⁿ extensions.
   - Implement point addition, doubling, and scalar multiplication (using windowed methods) with clear interfaces.

3. **Pairing Computations**  
   - Support at least one pairing (e.g., Optimal Ate) on pairing-friendly curves defined over Fp¹⁸.

4. **Testing and Verification**
   - Comprehensive unit tests for each field and curve operation.
   - Property-based tests (e.g., associativity, distributivity, inversion correctness).
   - Provide example programs demonstrating usage.
   - **KSS18 Parameter Generation & C Implementation (128-bit):**
     - **Python parameter script:**  
       - In `/tools/`, write a small Python script to compute the 128-bit KSS18 parameters using  
         \[
           u = 2^{44} + 2^{22} - 2^9 + 2,\quad
           p = \frac{u^8 + 5u^7 + 7u^6 + 37u^5 + 188u^4 + 259u^3 + 343u^2 + 1763u + 2401}{21},\quad
           r = \frac{u^6 + 37u^3 + 343}{343},
         \]  
       - Verify that *p* and *r* are prime.
       - Output constants `p`, `r`, `B=3`, `δ=3` and sample field/curve values (e.g., a few random Fp elements, x³+3 points, twist–points, pairing outputs).
     - **C integration:**  
       - Incorporate the generated `p`, `r`, `B=3`, `δ=3` into the C codebase.
       - Implement pure-C modules to instantiate and test:
         1. **Field arithmetic** in Fp, Fp³, Fp⁶, Fp¹⁸ using the constant *p*.
         2. **Curve** E: y² = x³ + 3 over GF(p), and its sextic twist E′ (δ¹ᐟ³ in GF(p³)).
         3. **Extension field** GF(p¹⁸) defined by the irreducible polynomial X¹⁸ − 3.
         4. **Pairing** (e.g., Optimal Ate) on E/E′ producing elements in GF(p¹⁸).
     - **Tests & Benchmarks:**  
       - Under `make test`, write C-based tests that:
         - Compare each field operation against the Python reference outputs.
         - Check that #E(GF(p)) = r and that points on E′(GF(p³)) have order divisible by r.
         - Validate bilinearity and non-degeneracy of the pairing.
       - Under `make bench`, benchmark:
         - Fp, Fp³, Fp⁶, Fp¹⁸ operations.
         - Curve point operations (add/double/SCM).
         - Pairing runtime.
   - Provide benchmarks for the parameter generation, field, curve, and pairing routines; include these in the `bench` target.

## 3. Non-Functional Requirements

- **Performance:** Use algorithmic optimizations (e.g., Karatsuba, Montgomery representation) where appropriate. Benchmark critical routines.
- **Code Quality:**  
  - C11 standard; `-Wall -Wextra -Werror`.  
  - Snake_case for functions/vars, UPPER_SNAKE for macros.  
  - Initialize/clear GMP variables; no leaks.  
  - Doxygen-compatible comments.
- **Modularity:** Separate modules and headers per layer (fp.c/h, fp3.c/h, fp6.c/h, fp18.c/h, ec.c/h, pairing.c/h).
- **Dependencies:** Only GMP; abstract big-int backend.
- **Portability:** Linux & BSD; provide dev-container.

## 4. Build and CI

- **Build Systems:** Makefile & CMakeLists.txt with targets: `all`, `clean`, `test`, `bench`.
- **CI:** GitHub Actions to build, run tests, and benchmarks on push/PR.
- **Dev Container:** `.devcontainer/` for VSCode.

## 5. Documentation and Maintainer Experience

- **README.md:** Overview, install, build, examples.
- **CONTRIBUTING.md:** Style, testing, PRs.
- **LICENSE:** MIT.
- **Issue Templates:** Bug report & feature request.
- **Labels:** Consistent GitHub labels.
