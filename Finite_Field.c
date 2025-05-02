#include "embedding_degree18.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef _Finite_Field_C_
#define _Finite_Field_C_

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
    mpz_mod(A->x_0, A->x_0, kss18_p);
}

void Fp_set_mpz(struct Fp *A, mpz_t B) {
    mpz_init_set(A->x_0, B);
    mpz_mod(A->x_0, A->x_0, kss18_p);
}

void Fp_random(struct Fp *A) {
    mpz_t r;
    mpz_init(r);
    init_rng_if_needed();
    mpz_urandomm(r, gmp_state, kss18_p);
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
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
}

void Fp_add_ui(struct Fp *ANS, struct Fp *A, unsigned long int B) {
    mpz_add_ui(ANS->x_0, A->x_0, B);
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
}

void Fp_sub(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_sub(ANS->x_0, A->x_0, B->x_0);
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
}

void Fp_sub_ui(struct Fp *ANS, struct Fp *A, unsigned long int B) {
    mpz_sub_ui(ANS->x_0, A->x_0, B);
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
}

void Fp_mul(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_mul(ANS->x_0, A->x_0, B->x_0);
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
}

void Fp_mul_ui(struct Fp *ANS, struct Fp *A, unsigned long int B) {
    mpz_mul_ui(ANS->x_0, A->x_0, B);
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
}

void Fp_invert(struct Fp *ANS, struct Fp *A) {
    mpz_t temp;
    mpz_init(temp);
    mpz_invert(temp, A->x_0, kss18_p);
    mpz_set(ANS->x_0, temp);
    mpz_clear(temp);
}

void Fp_div(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_t temp;
    mpz_init(temp);
    mpz_invert(temp, B->x_0, kss18_p);
    mpz_mul(ANS->x_0, A->x_0, temp);
    mpz_mod(ANS->x_0, ANS->x_0, kss18_p);
    mpz_clear(temp);
}

void Fp_pow(struct Fp *ANS, struct Fp *A, mpz_t B) {
    mpz_powm(ANS->x_0, A->x_0, B, kss18_p);
}

void Fp_sqrt(struct Fp *ANS, struct Fp *A) {
    mpz_t q, t, temp;
    mpz_init(q);
    mpz_init(t);
    mpz_init(temp);

    mpz_sub_ui(q, kss18_p, 1);
    mpz_divexact_ui(q, q, 2);
    mpz_powm(temp, A->x_0, q, kss18_p);

    if (mpz_cmp_ui(temp, 1) == 0) {
        mpz_t p_mod_4;
        mpz_init(p_mod_4);
        mpz_mod_ui(p_mod_4, kss18_p, 4);
        if (mpz_cmp_ui(p_mod_4, 3) == 0) {
            mpz_add_ui(q, kss18_p, 1);
            mpz_divexact_ui(q, q, 4);
            mpz_powm(ANS->x_0, A->x_0, q, kss18_p);
        } else {
            gmp_printf("ERROR: Fp_sqrt currently only supports p = 3 mod 4\n");
            mpz_set_ui(ANS->x_0, 0);
        }
        mpz_clear(p_mod_4);
    } else if (mpz_cmp_ui(temp, 0) == 0) {
        mpz_set_ui(ANS->x_0, 0);
    } else {
        gmp_printf("ERROR: sqrt is impossible (non-quadratic residue)\n");
        mpz_set_ui(ANS->x_0, 0);
    }

    mpz_clear(q);
    mpz_clear(t);
    mpz_clear(temp);
}

void Fp_neg(struct Fp *ANS, struct Fp *A) {
    if (mpz_sgn(A->x_0) == 0) {
        mpz_set_ui(ANS->x_0, 0);
    } else {
        mpz_sub(ANS->x_0, kss18_p, A->x_0);
    }
}

int Fp_cmp(struct Fp *A, struct Fp *B) {
    return mpz_cmp(A->x_0, B->x_0);
}

int Fp_cmp_mpz(struct Fp *A, mpz_t B) {
    return mpz_cmp(A->x_0, B);
}

// Update Fp3 functions (example for Fp3_add)
void Fp3_add(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B) {
    Fp_add(&ANS->a0, &A->a0, &B->a0);
    Fp_add(&ANS->a1, &A->a1, &B->a1);
    Fp_add(&ANS->a2, &A->a2, &B->a2);
}

// Update Fp6 functions (example for Fp6_add)
void Fp6_add(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B) {
    Fp3_add(&ANS->a0, &A->a0, &B->a0);
    Fp3_add(&ANS->a1, &A->a1, &B->a1);
}

// Update Fp18 functions (example for Fp18_add)
void Fp18_add(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B) {
    Fp6_add(&ANS->m0, &A->m0, &B->m0);
    Fp6_add(&ANS->m1, &A->m1, &B->m1);
    Fp6_add(&ANS->m2, &A->m2, &B->m2);
}

// Update Fp18_frobenius_map - This needs careful review based on KSS18 specifics
void Fp18_frobenius_map(struct Fp18 *ANS, struct Fp18 *A, int i) {
    fprintf(stderr, "WARNING: Fp18_frobenius_map is a placeholder and likely incorrect for KSS18.\n");

    Fp18_set(ANS, A);
    for (int j = 0; j < i; ++j) {
        Fp6_frobenius_map(&ANS->m0, &ANS->m0);
        Fp6_frobenius_map(&ANS->m1, &ANS->m1);
        Fp6_frobenius_map(&ANS->m2, &ANS->m2);
    }
}

#endif //Finite_Field_C_
