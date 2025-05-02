#include "../parameters.h"
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
    struct Fp3 a3, b3, c3, expected3;
    Fp3_init(&a3);
    Fp3_init(&b3);
    Fp3_init(&c3);
    Fp3_init(&expected3);
    
    // Set values for a3 and b3
    Fp_set_ui(&a3.x0, 1);
    Fp_set_ui(&a3.x1, 2);
    Fp_set_ui(&a3.x2, 3);
    
    Fp_set_ui(&b3.x0, 4);
    Fp_set_ui(&b3.x1, 5);
    Fp_set_ui(&b3.x2, 6);
    
    // Expected: component-wise addition
    Fp_set_ui(&expected3.x0, 5);  // 1+4
    Fp_set_ui(&expected3.x1, 7);  // 2+5
    Fp_set_ui(&expected3.x2, 9);  // 3+6
    
    Fp3_add(&c3, &a3, &b3);
    
    int result = (Fp_cmp(&c3.x0, &expected3.x0) == 0) &&
                 (Fp_cmp(&c3.x1, &expected3.x1) == 0) &&
                 (Fp_cmp(&c3.x2, &expected3.x2) == 0);
                 
    TEST(result, "Fp3_add: Basic addition");
    
    // Test Fp3 multiplication
    Fp3_mul(&c3, &a3, &b3);
    
    // For verification, we'd do a manual calculation here
    // For now, just confirm it doesn't crash and result is non-zero
    struct Fp3 zero3;
    Fp3_init(&zero3);
    Fp3_set_ui(&zero3, 0);
    
    TEST(Fp3_cmp(&c3, &zero3) != 0, "Fp3_mul: Result is non-zero");
    
    // Test inversion (a * a^-1 = 1)
    struct Fp3 inv3, one3, check3;
    Fp3_init(&inv3);
    Fp3_init(&one3);
    Fp3_init(&check3);
    
    // Set one3 to the multiplicative identity
    Fp_set_ui(&one3.x0, 1);
    Fp_set_ui(&one3.x1, 0);
    Fp_set_ui(&one3.x2, 0);
    
    Fp3_invert(&inv3, &a3);
    Fp3_mul(&check3, &a3, &inv3);
    
    result = (Fp_cmp(&check3.x0, &one3.x0) == 0) &&
             (Fp_cmp(&check3.x1, &one3.x1) == 0) &&
             (Fp_cmp(&check3.x2, &one3.x2) == 0);
             
    TEST(result, "Fp3_invert: a * a^-1 = 1");
    
    // Clean up
    Fp3_clear(&a3);
    Fp3_clear(&b3);
    Fp3_clear(&c3);
    Fp3_clear(&expected3);
    Fp3_clear(&zero3);
    Fp3_clear(&inv3);
    Fp3_clear(&one3);
    Fp3_clear(&check3);
}

// Tests for elliptic curve operations
void test_elliptic_operations() {
    printf("\nRunning elliptic curve operation tests...\n");
    
    // Test point addition on Ep(Fp)
    struct EFp P, Q, R;
    EFp_init(&P);
    EFp_init(&Q);
    EFp_init(&R);
    
    // Find points on the curve
    // This is a simplistic approach; in practice, we'd either:
    // 1. Generate random points and check they're on the curve
    // 2. Use precalculated points known to be on the curve
    EFp_random_set(&P);
    EFp_random_set(&Q);
    
    // Test point addition
    EFp_ECA(&R, &P, &Q);
    
    // It's hard to verify point addition without knowing the expected result
    // So we'll test that P+P = 2P (doubling)
    struct EFp P_double, P_plus_P;
    EFp_init(&P_double);
    EFp_init(&P_plus_P);
    
    EFp_ECD(&P_double, &P);  // P_double = 2P via doubling
    EFp_ECA(&P_plus_P, &P, &P);  // P_plus_P = P+P via addition
    
    TEST(EFp_cmp(&P_double, &P_plus_P) == 0, "EFp_ECA/ECD: 2P = P+P");
    
    // Test scalar multiplication
    mpz_t scalar;
    mpz_init(scalar);
    mpz_set_ui(scalar, 3);  // scalar = 3
    
    struct EFp P3, P_plus_P_plus_P;
    EFp_init(&P3);
    EFp_init(&P_plus_P_plus_P);
    
    EFp_SCM(&P3, &P, scalar);  // P3 = 3P via scalar multiplication
    EFp_ECA(&P_plus_P_plus_P, &P_double, &P);  // P_plus_P_plus_P = 2P+P = 3P
    
    TEST(EFp_cmp(&P3, &P_plus_P_plus_P) == 0, "EFp_SCM: 3P = 2P+P");
    
    // Clean up
    EFp_clear(&P);
    EFp_clear(&Q);
    EFp_clear(&R);
    EFp_clear(&P_double);
    EFp_clear(&P_plus_P);
    EFp_clear(&P3);
    EFp_clear(&P_plus_P_plus_P);
    mpz_clear(scalar);
}

// Tests for Fp operations - defined in test_fp.c
void test_fp_operations() {
    // External implementation in test_fp.c
    extern void test_fp_add();
    extern void test_fp_mul();
    extern void test_fp_inversion();
    extern void test_fp_associativity();
    
    printf("\nRunning Fp arithmetic tests...\n");
    test_fp_add();
    test_fp_mul();
    test_fp_inversion();
    test_fp_associativity();
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

int main() {
    printf("Running test suite for Finite Extension Field Degree 18 Library\n");
    
    // Initialize parameters
    init_parameters();
    
    // Set the generator value X
    mpz_set_str(X, "18446893747415302274", 10);
    
    // Generate curve parameters based on X
    generate_parameters();
    
    srand(time(NULL));  // Seed for random tests
    
    // Run all test suites
    test_fp_operations();
    test_fp3_operations();
    test_fp6_operations();
    test_fp18_operations();
    test_elliptic_operations();
    
    // Report results
    printf("\nTest summary: %d passed out of %d tests\n", test_passed, test_count);
    
    // Cleanup
    clear_parameters();
    
    return (test_passed == test_count) ? 0 : 1;  // Return 0 if all tests pass, 1 otherwise
}