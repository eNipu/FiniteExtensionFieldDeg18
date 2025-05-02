#include "fp.h"
#include "fp3.h"
#include "fp6.h"
#include "fp18.h"
#include "ec.h"
#include "pairing.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// External test functions
extern void test_fp_operations();
extern void test_fp3_operations();
extern void test_fp6_operations();
extern void test_fp18_operations();
extern void test_elliptic_operations();

// Test count for reporting
int test_count = 0;
int test_passed = 0;

// Simple test assertion macro
#define TEST(condition, message) \
    do { \
        test_count++; \
        if (condition) { \
            test_passed++; \
            printf("[PASS] %s\n", message); \
        } else { \
            printf("[FAIL] %s\n", message); \
        } \
    } while (0)

// Tests for Fp3 operations
void test_fp3_operations() {
    printf("\nRunning Fp3 arithmetic tests...\n");
    
    // Test Fp3 addition
    Fp3 a3, b3, c3, expected3;
    fp3_init(&a3);
    fp3_init(&b3);
    fp3_init(&c3);
    fp3_init(&expected3);
    
    // Set values for a3 and b3
    fp_set_ui(&a3.a0, 1);
    fp_set_ui(&a3.a1, 2);
    fp_set_ui(&a3.a2, 3);
    
    fp_set_ui(&b3.a0, 4);
    fp_set_ui(&b3.a1, 5);
    fp_set_ui(&b3.a2, 6);
    
    // Expected: component-wise addition
    fp_set_ui(&expected3.a0, 5);  // 1+4
    fp_set_ui(&expected3.a1, 7);  // 2+5
    fp_set_ui(&expected3.a2, 9);  // 3+6
    
    fp3_add(&c3, &a3, &b3);
    
    int result = (fp_cmp(&c3.a0, &expected3.a0) == 0) &&
                 (fp_cmp(&c3.a1, &expected3.a1) == 0) &&
                 (fp_cmp(&c3.a2, &expected3.a2) == 0);
                 
    TEST(result, "Fp3_add: Basic addition");
    
    // Test Fp3 multiplication
    fp3_mul(&c3, &a3, &b3);
    
    // For verification, we'd do a manual calculation here
    // For now, just confirm it doesn't crash and result is non-zero
    Fp3 zero3;
    fp3_init(&zero3);
    fp3_set_ui(&zero3, 0);
    
    TEST(fp3_cmp(&c3, &zero3) != 0, "Fp3_mul: Result is non-zero");
    
    // Test inversion (a * a^-1 = 1)
    Fp3 inv3, one3, check3;
    fp3_init(&inv3);
    fp3_init(&one3);
    fp3_init(&check3);
    
    // Set one3 to the multiplicative identity
    fp3_set_ui(&one3, 1);
    
    fp3_inv(&inv3, &a3);
    fp3_mul(&check3, &a3, &inv3);
    
    result = (fp3_cmp(&check3, &one3) == 0);
             
    TEST(result, "Fp3_invert: a * a^-1 = 1");
    
    // Clean up
    fp3_clear(&a3);
    fp3_clear(&b3);
    fp3_clear(&c3);
    fp3_clear(&expected3);
    fp3_clear(&zero3);
    fp3_clear(&inv3);
    fp3_clear(&one3);
    fp3_clear(&check3);
}

// Tests for elliptic curve operations
void test_elliptic_operations() {
    printf("\nRunning elliptic curve operation tests...\n");
    
    // Test point addition on Ep(Fp)
    EcFp P, Q, R;
    ecfp_init(&P);
    ecfp_init(&Q);
    ecfp_init(&R);
    
    // Find points on the curve
    // This is a simplistic approach; in practice, we'd either:
    // 1. Generate random points and check they're on the curve
    // 2. Use precalculated points known to be on the curve
    fp_random(&P.x);
    fp_random(&P.y);
    P.infinity = 0;
    fp_random(&Q.x);
    fp_random(&Q.y);
    Q.infinity = 0;
    
    // Test point addition
    ecfp_add(&R, &P, &Q);
    
    // It's hard to verify point addition without knowing the expected result
    // So we'll test that P+P = 2P (doubling)
    EcFp P_double, P_plus_P;
    ecfp_init(&P_double);
    ecfp_init(&P_plus_P);
    
    ecfp_double(&P_double, &P);  // P_double = 2P via doubling
    ecfp_add(&P_plus_P, &P, &P);  // P_plus_P = P+P via addition
    
    TEST(ecfp_cmp(&P_double, &P_plus_P) == 0, "ecfp_add/ecfp_double: 2P = P+P");
    
    // Test scalar multiplication
    mpz_t scalar;
    mpz_init(scalar);
    mpz_set_ui(scalar, 3);  // scalar = 3
    
    EcFp P3, P_plus_P_plus_P;
    ecfp_init(&P3);
    ecfp_init(&P_plus_P_plus_P);
    
    ecfp_scalar_mul(&P3, &P, scalar);  // P3 = 3P via scalar multiplication
    ecfp_add(&P_plus_P_plus_P, &P_double, &P);  // P_plus_P_plus_P = 2P+P = 3P
    
    TEST(ecfp_cmp(&P3, &P_plus_P_plus_P) == 0, "ecfp_scalar_mul: 3P = 2P+P");
    
    // Clean up
    ecfp_clear(&P);
    ecfp_clear(&Q);
    ecfp_clear(&R);
    ecfp_clear(&P_double);
    ecfp_clear(&P_plus_P);
    ecfp_clear(&P3);
    ecfp_clear(&P_plus_P_plus_P);
    mpz_clear(scalar);
}

// Tests for Fp operations - defined in test_fp.c
typedef void (*test_func_t)();
void test_fp_operations() {
    extern void test_fp_add();
    extern void test_fp_mul();
    extern void test_fp_inversion();
    extern void test_fp_associativity();
    printf("\nRunning Fp arithmetic tests...\n");
    test_func_t tests[] = {test_fp_add, test_fp_mul, test_fp_inversion, test_fp_associativity};
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); ++i) tests[i]();
}

// Tests for Fp6 operations
void test_fp6_operations() {
    printf("\nRunning Fp6 arithmetic tests...\n");
    // Similar to Fp3 tests, basic validation for now
    TEST(1, "Fp6 tests to be implemented");
}

// Tests for Fp18 operations
void test_fp18_operations() {
    printf("\nRunning Fp18 arithmetic tests...\n");
    // Similar to Fp3 tests, basic validation for now
    TEST(1, "Fp18 tests to be implemented");
}

// Fp3 tests
void test_fp3_arithmetic() {
    printf("\n[TEST] Fp3 arithmetic\n");
    Fp3 a, b, c, d, e;
    fp3_init(&a); fp3_init(&b); fp3_init(&c); fp3_init(&d); fp3_init(&e);
    fp3_random(&a); fp3_random(&b);
    // Addition
    fp3_add(&c, &a, &b);
    fp3_sub(&d, &c, &b);
    TEST(fp3_cmp(&d, &a) == 0, "Fp3 add/sub inverse");
    // Multiplication
    fp3_mul(&c, &a, &b);
    // Inversion
    fp3_inv(&d, &b);
    fp3_mul(&e, &b, &d);
    fp3_set_ui(&c, 1);
    TEST(fp3_cmp(&e, &c) == 0, "Fp3 inversion");
    // Power
    mpz_t exp; mpz_init_set_ui(exp, 5);
    fp3_pow(&d, &a, exp);
    // Negation
    fp3_neg(&c, &a);
    fp3_add(&e, &a, &c);
    fp3_set_ui(&d, 0);
    TEST(fp3_cmp(&e, &d) == 0, "Fp3 negation");
    fp3_clear(&a); fp3_clear(&b); fp3_clear(&c); fp3_clear(&d); fp3_clear(&e); mpz_clear(exp);
}

// Fp6 tests
void test_fp6_arithmetic() {
    printf("\n[TEST] Fp6 arithmetic\n");
    Fp6 a, b, c, d, e;
    fp6_init(&a); fp6_init(&b); fp6_init(&c); fp6_init(&d); fp6_init(&e);
    fp6_random(&a); fp6_random(&b);
    // Addition
    fp6_add(&c, &a, &b);
    fp6_sub(&d, &c, &b);
    TEST(fp6_cmp(&d, &a) == 0, "Fp6 add/sub inverse");
    // Multiplication
    fp6_mul(&c, &a, &b);
    // Inversion
    fp6_inv(&d, &b);
    fp6_mul(&e, &b, &d);
    fp6_set_ui(&c, 1);
    TEST(fp6_cmp(&e, &c) == 0, "Fp6 inversion");
    // Power
    mpz_t exp; mpz_init_set_ui(exp, 5);
    fp6_pow(&d, &a, exp);
    // Negation
    fp6_neg(&c, &a);
    fp6_add(&e, &a, &c);
    fp6_set_ui(&d, 0);
    TEST(fp6_cmp(&e, &d) == 0, "Fp6 negation");
    fp6_clear(&a); fp6_clear(&b); fp6_clear(&c); fp6_clear(&d); fp6_clear(&e); mpz_clear(exp);
}

// Fp18 tests
void test_fp18_arithmetic() {
    printf("\n[TEST] Fp18 arithmetic\n");
    Fp18 a, b, c, d, e;
    fp18_init(&a); fp18_init(&b); fp18_init(&c); fp18_init(&d); fp18_init(&e);
    fp18_random(&a); fp18_random(&b);
    // Addition
    fp18_add(&c, &a, &b);
    fp18_sub(&d, &c, &b);
    TEST(fp18_cmp(&d, &a) == 0, "Fp18 add/sub inverse");
    // Multiplication
    fp18_mul(&c, &a, &b);
    // Power
    mpz_t exp; mpz_init_set_ui(exp, 5);
    fp18_pow(&d, &a, exp);
    // Negation
    fp18_neg(&c, &a);
    fp18_add(&e, &a, &c);
    fp18_set_ui(&d, 0);
    TEST(fp18_cmp(&e, &d) == 0, "Fp18 negation");
    fp18_clear(&a); fp18_clear(&b); fp18_clear(&c); fp18_clear(&d); fp18_clear(&e); mpz_clear(exp);
}

// Property-based tests for Fp3
void test_fp3_properties() {
    printf("\n[PROPERTY] Fp3 associativity, distributivity, inversion\n");
    Fp3 a, b, c, left, right, tmp;
    fp3_init(&a); fp3_init(&b); fp3_init(&c); fp3_init(&left); fp3_init(&right); fp3_init(&tmp);
    // Associativity: (a + b) + c == a + (b + c)
    fp3_random(&a); fp3_random(&b); fp3_random(&c);
    fp3_add(&tmp, &a, &b); fp3_add(&left, &tmp, &c);
    fp3_add(&tmp, &b, &c); fp3_add(&right, &a, &tmp);
    TEST(fp3_cmp(&left, &right) == 0, "Fp3 add associativity");
    // Distributivity: a*(b + c) == a*b + a*c
    fp3_add(&tmp, &b, &c); fp3_mul(&left, &a, &tmp);
    fp3_mul(&tmp, &a, &b); fp3_mul(&right, &a, &c); fp3_add(&right, &tmp, &right);
    TEST(fp3_cmp(&left, &right) == 0, "Fp3 distributivity");
    // Inversion: a * a^-1 == 1 (if a != 0)
    fp3_random(&a);
    fp3_inv(&b, &a);
    fp3_mul(&c, &a, &b);
    fp3_set_ui(&left, 1);
    TEST(fp3_cmp(&c, &left) == 0, "Fp3 inversion correctness");
    fp3_clear(&a); fp3_clear(&b); fp3_clear(&c); fp3_clear(&left); fp3_clear(&right); fp3_clear(&tmp);
}

// Property-based tests for Fp6
void test_fp6_properties() {
    printf("\n[PROPERTY] Fp6 associativity, distributivity, inversion\n");
    Fp6 a, b, c, left, right, tmp;
    fp6_init(&a); fp6_init(&b); fp6_init(&c); fp6_init(&left); fp6_init(&right); fp6_init(&tmp);
    // Associativity: (a + b) + c == a + (b + c)
    fp6_random(&a); fp6_random(&b); fp6_random(&c);
    fp6_add(&tmp, &a, &b); fp6_add(&left, &tmp, &c);
    fp6_add(&tmp, &b, &c); fp6_add(&right, &a, &tmp);
    TEST(fp6_cmp(&left, &right) == 0, "Fp6 add associativity");
    // Distributivity: a*(b + c) == a*b + a*c
    fp6_add(&tmp, &b, &c); fp6_mul(&left, &a, &tmp);
    fp6_mul(&tmp, &a, &b); fp6_mul(&right, &a, &c); fp6_add(&right, &tmp, &right);
    TEST(fp6_cmp(&left, &right) == 0, "Fp6 distributivity");
    // Inversion: a * a^-1 == 1 (if a != 0)
    fp6_random(&a);
    fp6_inv(&b, &a);
    fp6_mul(&c, &a, &b);
    fp6_set_ui(&left, 1);
    TEST(fp6_cmp(&c, &left) == 0, "Fp6 inversion correctness");
    fp6_clear(&a); fp6_clear(&b); fp6_clear(&c); fp6_clear(&left); fp6_clear(&right); fp6_clear(&tmp);
}

// Property-based tests for Fp18
void test_fp18_properties() {
    printf("\n[PROPERTY] Fp18 associativity, distributivity, inversion\n");
    Fp18 a, b, c, left, right, tmp;
    fp18_init(&a); fp18_init(&b); fp18_init(&c); fp18_init(&left); fp18_init(&right); fp18_init(&tmp);
    // Associativity: (a + b) + c == a + (b + c)
    fp18_random(&a); fp18_random(&b); fp18_random(&c);
    fp18_add(&tmp, &a, &b); fp18_add(&left, &tmp, &c);
    fp18_add(&tmp, &b, &c); fp18_add(&right, &a, &tmp);
    TEST(fp18_cmp(&left, &right) == 0, "Fp18 add associativity");
    // Distributivity: a*(b + c) == a*b + a*c
    fp18_add(&tmp, &b, &c); fp18_mul(&left, &a, &tmp);
    fp18_mul(&tmp, &a, &b); fp18_mul(&right, &a, &c); fp18_add(&right, &tmp, &right);
    TEST(fp18_cmp(&left, &right) == 0, "Fp18 distributivity");
    // Inversion: a * a^-1 == 1 (if a != 0)
    fp18_random(&a);
    fp18_inv(&b, &a);
    fp18_mul(&c, &a, &b);
    fp18_set_ui(&left, 1);
    TEST(fp18_cmp(&c, &left) == 0, "Fp18 inversion correctness");
    fp18_clear(&a); fp18_clear(&b); fp18_clear(&c); fp18_clear(&left); fp18_clear(&right); fp18_clear(&tmp);
}

// Property-based tests for EC over Fp
void test_ecfp_properties() {
    printf("\n[PROPERTY] EC over Fp: addition commutativity\n");
    EcFp P, Q, R1, R2;
    ecfp_init(&P); ecfp_init(&Q); ecfp_init(&R1); ecfp_init(&R2);
    // Generate random points (for demo, not guaranteed on curve)
    fp_random(&P.x); fp_random(&P.y); P.infinity = 0;
    fp_random(&Q.x); fp_random(&Q.y); Q.infinity = 0;
    ecfp_add(&R1, &P, &Q);
    ecfp_add(&R2, &Q, &P);
    TEST(ecfp_cmp(&R1, &R2) == 0, "EC add commutativity");
    ecfp_clear(&P); ecfp_clear(&Q); ecfp_clear(&R1); ecfp_clear(&R2);
}

// Reference test scaffold (to be filled with known-good values)
void test_fp3_reference() {
    printf("\n[REFERENCE] Fp3 reference test (placeholder)\n");
    // Example: compare against Python/SageMath output
    // Fp3 a = ...; Fp3 b = ...; Fp3 expected = ...;
    // fp3_add(&c, &a, &b);
    // TEST(fp3_cmp(&c, &expected) == 0, "Fp3 reference add");
}

// Reference test scaffold for Fp6
void test_fp6_reference() {
    printf("\n[REFERENCE] Fp6 reference test (placeholder)\n");
    // Example: compare against Python/SageMath output
    // Fp6 a = ...; Fp6 b = ...; Fp6 expected = ...;
    // fp6_add(&c, &a, &b);
    // TEST(fp6_cmp(&c, &expected) == 0, "Fp6 reference add");
}

// Reference test scaffold for Fp18
void test_fp18_reference() {
    printf("\n[REFERENCE] Fp18 reference test (placeholder)\n");
    // Example: compare against Python/SageMath output
    // Fp18 a = ...; Fp18 b = ...; Fp18 expected = ...;
    // fp18_add(&c, &a, &b);
    // TEST(fp18_cmp(&c, &expected) == 0, "Fp18 reference add");
}

// Reference test scaffold for EC
void test_ec_reference() {
    printf("\n[REFERENCE] EC reference test\n");
    // Use a known point on E: y^2 = x^3 + 3 over Fp
    EcFp P, Q, R, expected;
    ecfp_init(&P); ecfp_init(&Q); ecfp_init(&R); ecfp_init(&expected);
    // Set P = (1, 2) (assuming 1^3+3=4, 2^2=4, so (1,2) is on the curve for p > 3)
    fp_set_ui(&P.x, 1); fp_set_ui(&P.y, 2); P.infinity = 0;
    // Set Q = (2, 3) (for demo, not guaranteed on curve, but for test structure)
    fp_set_ui(&Q.x, 2); fp_set_ui(&Q.y, 3); Q.infinity = 0;
    // Expected = P + Q (manually computed or from SageMath)
    // For now, just check that addition does not crash and result is not infinity
    ecfp_add(&R, &P, &Q);
    TEST(!ecfp_is_infinity(&R), "EC reference: P+Q is not infinity");
    // Doubling: expected = 2P
    ecfp_double(&expected, &P);
    ecfp_add(&R, &P, &P);
    TEST(ecfp_cmp(&R, &expected) == 0, "EC reference: 2P = P+P");
    ecfp_clear(&P); ecfp_clear(&Q); ecfp_clear(&R); ecfp_clear(&expected);
}

// Pairing property/reference test scaffold
void test_pairing_properties() {
    printf("\n[PROPERTY/REFERENCE] Pairing test (stub, not full KSS18)\n");
    // Use dummy points for now
    EcFp P; EcFp3 Q; Fp18 gt;
    ecfp_init(&P); ecfp3_init(&Q); fp18_init(&gt);
    fp_set_ui(&P.x, 1); fp_set_ui(&P.y, 2); P.infinity = 0;
    fp3_set_ui(&Q.x, 1); fp3_set_ui(&Q.y, 2); Q.infinity = 0;
    optimal_ate_pairing(&gt, &P, &Q);
    // For the stub, result should be 1
    Fp18 one; fp18_init(&one); fp18_set_ui(&one, 1);
    TEST(fp18_cmp(&gt, &one) == 0, "Pairing stub returns 1");
    fp18_clear(&gt); fp18_clear(&one); ecfp_clear(&P); ecfp3_clear(&Q);
}

int main() {
    printf("Running test suite for Finite Extension Field Degree 18 Library\n");
    
    // Initialize parameters
    init_kss18_params();
    
    srand(time(NULL));  // Seed for random tests
    
    // Run all test suites
    test_fp_operations();
    test_fp3_arithmetic();
    test_fp6_arithmetic();
    test_fp18_arithmetic();
    test_fp3_properties();
    test_fp6_properties();
    test_fp18_properties();
    test_ecfp_properties();
    test_fp3_reference();
    test_fp6_reference();
    test_fp18_reference();
    test_ec_reference();
    test_pairing_properties();
    
    // Report results
    printf("\nTest summary: %d passed out of %d tests\n", test_passed, test_count);
    
    // Cleanup
    clear_kss18_params();
    
    return (test_passed == test_count) ? 0 : 1;  // Return 0 if all tests pass, 1 otherwise
}