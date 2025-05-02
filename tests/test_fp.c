#include "fp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Test addition in Fp
void test_fp_add() {
    Fp a, b, result, expected;
    fp_init(&a); 
    fp_init(&b); 
    fp_init(&result); 
    fp_init(&expected);

    // Test case: a = 123, b = 456, expected = (123 + 456) mod p
    fp_set_ui(&a, 123);
    fp_set_ui(&b, 456);
    
    // Calculate expected result
    mpz_t temp;
    mpz_init(temp);
    mpz_add_ui(temp, a.x, 456);
    mpz_mod(temp, temp, kss18_p);
    
    fp_set_mpz(&expected, temp);
    fp_add(&result, &a, &b);
    
    mpz_clear(temp); 
    fp_clear(&a); 
    fp_clear(&b); 
    fp_clear(&result); 
    fp_clear(&expected);
}

// Test multiplication in Fp
void test_fp_mul() {
    Fp a, b, result, expected;
    fp_init(&a); 
    fp_init(&b); 
    fp_init(&result); 
    fp_init(&expected);

    // Test case: a = 123, b = 456, expected = (123 * 456) mod p
    fp_set_ui(&a, 123);
    fp_set_ui(&b, 456);
    
    // Calculate expected result
    mpz_t temp;
    mpz_init(temp);
    mpz_mul_ui(temp, a.x, 456);
    mpz_mod(temp, temp, kss18_p);
    
    fp_set_mpz(&expected, temp);
    fp_mul(&result, &a, &b);
    
    mpz_clear(temp); 
    fp_clear(&a); 
    fp_clear(&b); 
    fp_clear(&result); 
    fp_clear(&expected);
}

// Test inversion in Fp
void test_fp_inversion() {
    Fp a, inv_a, result, one;
    fp_init(&a); 
    fp_init(&inv_a); 
    fp_init(&result); 
    fp_init(&one);

    // Test a non-zero value
    fp_set_ui(&a, 123);
    
    // Compute inverse of a
    fp_inv(&inv_a, &a);
    
    // Multiply a * inv_a, should be 1
    fp_mul(&result, &a, &inv_a);
    
    // Check if result is 1
    fp_set_ui(&one, 1);
    
    fp_clear(&a); 
    fp_clear(&inv_a); 
    fp_clear(&result); 
    fp_clear(&one);
}

// Test associativity of addition: (a + b) + c = a + (b + c)
void test_fp_associativity() {
    Fp a, b, c, left, right, temp1, temp2;
    fp_init(&a); 
    fp_init(&b); 
    fp_init(&c); 
    fp_init(&left); 
    fp_init(&right); 
    fp_init(&temp1); 
    fp_init(&temp2);

    // Set random values
    fp_random(&a);
    fp_random(&b);
    fp_random(&c);
    
    // Calculate (a + b) + c
    fp_add(&temp1, &a, &b); 
    fp_add(&left, &temp1, &c);
    
    // Calculate a + (b + c)
    fp_add(&temp2, &b, &c); 
    fp_add(&right, &a, &temp2);
    
    fp_clear(&a); 
    fp_clear(&b); 
    fp_clear(&c); 
    fp_clear(&left); 
    fp_clear(&right); 
    fp_clear(&temp1); 
    fp_clear(&temp2);
}