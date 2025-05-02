#include "embedding_degree18.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef _Finite_Field_C_
#define _Finite_Field_C_

// Counter variables for performance evaluation
static int fp_add = 0;
static int fp_mul = 0;

// RNG state for random number generation
static gmp_randstate_t gmp_state;

// Initialize RNG if needed
static void init_rng_if_needed() {
    static int initialized = 0;
    if (!initialized) {
        gmp_randinit_default(gmp_state);
        gmp_randseed_ui(gmp_state, time(NULL));
        initialized = 1;
    }
}

// Global variables declared in embedding_degree18.h, implemented in parameters.c
// Do NOT redefine them here
unsigned long int c1 = 1;

// Fp arithmetic implementations
void Fp_init(struct Fp *A) {
    mpz_init(A->x_0);
}

void Fp_set(struct Fp *A, struct Fp *B) {
    mpz_set(A->x_0, B->x_0);
}

void Fp_set_ui(struct Fp *A, unsigned long int B) {
    mpz_set_ui(A->x_0, B);
    mpz_mod(A->x_0, A->x_0, prime);
}

void Fp_set_mpz(struct Fp *A, mpz_t B) {
    mpz_init_set(A->x_0, B);
    mpz_mod(A->x_0, A->x_0, prime);
}

void Fp_random(struct Fp *A) {
    mpz_t r;
    mpz_init(r);
    init_rng_if_needed();
    mpz_urandomm(r, gmp_state, prime);
    mpz_set(A->x_0, r);
    mpz_clear(r);
}

void Fp_clear(struct Fp *A) {
    mpz_clear(A->x_0);
}

void Fp_printf(struct Fp *A) {
    gmp_printf("%Zd", A->x_0);
}

void Fp_add(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_add(ANS->x_0, A->x_0, B->x_0);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    fp_add++;
}

void Fp_add_ui(struct Fp *ANS, struct Fp *A, unsigned long int B) {
    mpz_add_ui(ANS->x_0, A->x_0, B);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    fp_add++;
}

void Fp_sub(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_sub(ANS->x_0, A->x_0, B->x_0);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    fp_add++;
}

void Fp_sub_ui(struct Fp *ANS, struct Fp *A, unsigned long int B) {
    mpz_sub_ui(ANS->x_0, A->x_0, B);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    fp_add++;
}

void Fp_mul(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_mul(ANS->x_0, A->x_0, B->x_0);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    fp_mul++;
}

void Fp_mul_ui(struct Fp *ANS, struct Fp *A, unsigned long int B) {
    mpz_mul_ui(ANS->x_0, A->x_0, B);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    fp_mul++;
}

void Fp_invert(struct Fp *ANS, struct Fp *A) {
    mpz_t temp;
    mpz_init(temp);
    mpz_invert(temp, A->x_0, prime);
    mpz_set(ANS->x_0, temp);
    mpz_clear(temp);
    inv_count++;
}

void Fp_div(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_t temp;
    mpz_init(temp);
    mpz_invert(temp, B->x_0, prime);
    mpz_mul(ANS->x_0, A->x_0, temp);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
    mpz_clear(temp);
    fp_mul++;
}

void Fp_pow(struct Fp *ANS, struct Fp *A, mpz_t B) {
    mpz_powm(ANS->x_0, A->x_0, B, prime);
}

void Fp_sqrt(struct Fp *ANS, struct Fp *A) {
    mpz_t q, t, temp;
    mpz_init(q);
    mpz_init(t);
    mpz_init(temp);
    
    mpz_sub_ui(q, prime, 1);
    mpz_divexact_ui(q, q, 2);
    mpz_powm(temp, A->x_0, q, prime);
    
    if (mpz_cmp_ui(temp, 1) == 0) {
        mpz_add_ui(q, prime, 1);
        mpz_divexact_ui(q, q, 4);
        mpz_powm(ANS->x_0, A->x_0, q, prime);
    } else {
        // Handle non-quadratic residue
        gmp_printf("ERROR: sqrt is impossible");
        mpz_set_ui(ANS->x_0, 0);
    }
    
    mpz_clear(q);
    mpz_clear(t);
    mpz_clear(temp);
}

void Fp_neg(struct Fp *ANS, struct Fp *A) {
    mpz_sub(ANS->x_0, prime, A->x_0);
    mpz_mod(ANS->x_0, ANS->x_0, prime);
}

int Fp_cmp(struct Fp *A, struct Fp *B) {
    return mpz_cmp(A->x_0, B->x_0);
}

int Fp_cmp_mpz(struct Fp *A, mpz_t B) {
    return mpz_cmp(A->x_0, B);
}

// Implement the remaining Fp3, Fp6, and Fp18 functions as needed
// ...

// Implementation of Fp18_frobenius_map with power parameter
void Fp18_frobenius_map(struct Fp18 *ANS, struct Fp18 *A, int i) {
    // Apply the Frobenius map to each component with power i
    // For a field extension of degree 18, the Frobenius map raises elements to the power of p^i
    
    // Apply Frobenius map to each component based on the power i
    // This is a simplified implementation and would need to be replaced
    // with the actual KSS curve-specific implementation
    
    // Initialize temporary variables
    struct Fp6 t0, t1, t2;
    Fp6_init(&t0);
    Fp6_init(&t1);
    Fp6_init(&t2);
    
    // Apply Frobenius to each component
    // The actual implementation would use precomputed constants for efficiency
    
    // Copy A to ANS first
    Fp18_set(ANS, A);
    
    // Apply Frobenius map based on the power i
    // This is just a placeholder - actual implementation would be curve-specific
    for (int j = 0; j < i; j++) {
        // Apply one Frobenius map (i.e., raise to power p)
        Fp6_frobenius_map(&(ANS->m0), &(ANS->m0));
        Fp6_frobenius_map(&(ANS->m1), &(ANS->m1));
        Fp6_frobenius_map(&(ANS->m2), &(ANS->m2));
        
        // Apply additional twists as needed for the specific curve
        // This would involve multiplication by precomputed constants
    }
    
    // Clean up
    Fp6_clear(&t0);
    Fp6_clear(&t1);
    Fp6_clear(&t2);
}

#endif //Finite_Field_C_
