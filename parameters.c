#define DEFINE_GLOBAL_VARIABLES
#include "parameters.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

// Global parameters for KSS degree 18 curve
mpz_t X;       // Parameter for the KSS curve
mpz_t prime;   // Field characteristic
mpz_t r_order; // Order of the subgroup
mpz_t t_trace; // Trace of Frobenius
mpz_t r_order_EFp; // Order of EFp
mpz_t b;       // Curve constant: y^2 = x^3 + b

// Constants used in the implementation
mpz_t c1_leg;
mpz_t c1_leg_bar;
mpz_t c1_omega;
mpz_t c1_omega_bar;

// Binary representation of X for efficient scalar multiplication
int *X_bit_binary = NULL;
int X_bit;

// Optimization counters
int add_count_miller = 0;
int add_count_finalexp = 0;
int sqr_count_miller = 0;
int sqr_count_finalexp = 0;
int inv_count = 0;

/**
 * Initialize parameters for KSS degree 18 curve
 */
void init_parameters(void) {
    // Initialize GMP variables
    mpz_init(X);
    mpz_init(prime);
    mpz_init(r_order);
    mpz_init(t_trace);
    mpz_init(r_order_EFp);
    mpz_init(b);
    
    // Initialize constants
    mpz_init(c1_leg);
    mpz_init(c1_leg_bar);
    mpz_init(c1_omega);
    mpz_init(c1_omega_bar);
    
    // Initialize counters
    add_count_miller = 0;
    add_count_finalexp = 0;
    sqr_count_miller = 0;
    sqr_count_finalexp = 0;
    inv_count = 0;
}

/**
 * Generate parameters for KSS degree 18 curve
 * 
 * This function sets the specific values for X and other parameters
 * that define the KSS degree 18 curve.
 */
void generate_parameters(void) {
    // Set X = 0x10000010100 (example - replace with actual value for KSS curve)
    mpz_set_str(X, "10000010100", 16);
    
    // Calculate prime p
    mpz_pow_ui(prime, X, 6);
    mpz_mul_ui(prime, prime, 81);
    mpz_add_ui(prime, prime, 3);
    mpz_mul_ui(prime, prime, 7);
    mpz_sub_ui(prime, prime, 7);
    
    // Calculate order r
    mpz_pow_ui(r_order, X, 4);
    
    // Fix the incorrect mpz_sub_ui call
    mpz_t temp;
    mpz_init(temp);
    mpz_pow_ui(temp, X, 2);
    mpz_sub(r_order, r_order, temp);
    mpz_clear(temp);
    
    mpz_add_ui(r_order, r_order, 1);
    
    // Calculate trace t
    mpz_pow_ui(t_trace, X, 3);
    mpz_add_ui(t_trace, t_trace, 1);
    
    // Set curve constant b
    mpz_set_ui(b, 3);
    
    // Calculate constants for arithmetic
    mpz_set_ui(c1_leg, 1);
    mpz_set_ui(c1_leg_bar, 1);
    mpz_set_ui(c1_omega, 1);
    mpz_set_ui(c1_omega_bar, 1);
    
    // Generate binary representation of X for scalar multiplication
    X_bit = mpz_sizeinbase(X, 2);
    X_bit_binary = (int*)malloc(X_bit * sizeof(int));
    
    if (X_bit_binary != NULL) {
        for (int i = 0; i < X_bit; i++) {
            X_bit_binary[i] = mpz_tstbit(X, i);
        }
    }
}
