#include "../fp18_arith.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

// Test addition in Fp
void test_fp_add() {
    struct Fp a, b, result, expected;
    Fp_init(&a);
    Fp_init(&b);
    Fp_init(&result);
    Fp_init(&expected);

    // Test case: a = 123, b = 456, expected = (123 + 456) mod p
    Fp_set_ui(&a, 123);
    Fp_set_ui(&b, 456);
    
    // Calculate expected result
    mpz_t temp;
    mpz_init(temp);
    mpz_add_ui(temp, a.x0, 456);
    mpz_mod(temp, temp, prime);
    
    Fp_set_mpz(&expected, temp);
    Fp_add(&result, &a, &b);
    
    TEST(Fp_cmp(&result, &expected) == 0, "Fp_add: 123 + 456");
    
    mpz_clear(temp);
    Fp_clear(&a);
    Fp_clear(&b);
    Fp_clear(&result);
    Fp_clear(&expected);
}

// Test multiplication in Fp
void test_fp_mul() {
    struct Fp a, b, result, expected;
    Fp_init(&a);
    Fp_init(&b);
    Fp_init(&result);
    Fp_init(&expected);

    // Test case: a = 123, b = 456, expected = (123 * 456) mod p
    Fp_set_ui(&a, 123);
    Fp_set_ui(&b, 456);
    
    // Calculate expected result
    mpz_t temp;
    mpz_init(temp);
    mpz_mul_ui(temp, a.x0, 456);
    mpz_mod(temp, temp, prime);
    
    Fp_set_mpz(&expected, temp);
    Fp_mul(&result, &a, &b);
    
    TEST(Fp_cmp(&result, &expected) == 0, "Fp_mul: 123 * 456");
    
    mpz_clear(temp);
    Fp_clear(&a);
    Fp_clear(&b);
    Fp_clear(&result);
    Fp_clear(&expected);
}

// Test inversion in Fp
void test_fp_inversion() {
    struct Fp a, inv_a, result;
    Fp_init(&a);
    Fp_init(&inv_a);
    Fp_init(&result);

    // Test a non-zero value
    Fp_set_ui(&a, 123);
    
    // Compute inverse of a
    Fp_invert(&inv_a, &a);
    
    // Multiply a * inv_a, should be 1
    Fp_mul(&result, &a, &inv_a);
    
    // Check if result is 1
    struct Fp one;
    Fp_init(&one);
    Fp_set_ui(&one, 1);
    
    TEST(Fp_cmp(&result, &one) == 0, "Fp_invert: a * a^(-1) = 1");
    
    Fp_clear(&a);
    Fp_clear(&inv_a);
    Fp_clear(&result);
    Fp_clear(&one);
}

// Test associativity of addition: (a + b) + c = a + (b + c)
void test_fp_associativity() {
    struct Fp a, b, c, left, right, temp1, temp2;
    Fp_init(&a);
    Fp_init(&b);
    Fp_init(&c);
    Fp_init(&left);
    Fp_init(&right);
    Fp_init(&temp1);
    Fp_init(&temp2);

    // Set random values
    Fp_random(&a);
    Fp_random(&b);
    Fp_random(&c);
    
    // Calculate (a + b) + c
    Fp_add(&temp1, &a, &b);
    Fp_add(&left, &temp1, &c);
    
    // Calculate a + (b + c)
    Fp_add(&temp2, &b, &c);
    Fp_add(&right, &a, &temp2);
    
    TEST(Fp_cmp(&left, &right) == 0, "Fp_add associativity: (a + b) + c = a + (b + c)");
    
    Fp_clear(&a);
    Fp_clear(&b);
    Fp_clear(&c);
    Fp_clear(&left);
    Fp_clear(&right);
    Fp_clear(&temp1);
    Fp_clear(&temp2);
}

int main() {
    printf("Running Fp arithmetic tests...\n");
    
    // Initialize parameters
    init_parameters();
    
    // Set the generator value X
    mpz_set_str(X, "18446893747415302274", 10);
    
    // Generate curve parameters based on X
    generate_parameters();
    
    srand(time(NULL));  // Seed for random tests
    
    // Run tests
    test_fp_add();
    test_fp_mul();
    test_fp_inversion();
    test_fp_associativity();
    
    // Report results
    printf("\nTest summary: %d passed out of %d tests\n", test_passed, test_count);
    
    // Cleanup
    clear_parameters();
    
    return (test_passed == test_count) ? 0 : 1;  // Return 0 if all tests pass, 1 otherwise
}