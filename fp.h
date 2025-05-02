#ifndef FP_H
#define FP_H

#include <gmp.h>
#include "parameters.h"

/**
 * @file fp.h
 * @brief Base field Fp arithmetic for KSS18 (C11, GMP)
 */

/**
 * @struct Fp
 * @brief Element of the base field Fp
 */
typedef struct {
    mpz_t x;
} Fp;

/** Initialize Fp element */
void fp_init(Fp *a);
/** Clear Fp element */
void fp_clear(Fp *a);
/** Set Fp element (copy) */
void fp_set(Fp *rop, const Fp *op);
/** Set Fp element from unsigned long */
void fp_set_ui(Fp *rop, unsigned long v);
/** Set Fp element from mpz_t */
void fp_set_mpz(Fp *rop, const mpz_t v);
/** Random Fp element */
void fp_random(Fp *rop);
/** Print Fp element */
void fp_print(const Fp *a);
/** Addition: rop = a + b mod p */
void fp_add(Fp *rop, const Fp *a, const Fp *b);
/** Subtraction: rop = a - b mod p */
void fp_sub(Fp *rop, const Fp *a, const Fp *b);
/** Multiplication: rop = a * b mod p */
void fp_mul(Fp *rop, const Fp *a, const Fp *b);
/** Inversion: rop = a^{-1} mod p */
void fp_inv(Fp *rop, const Fp *a);
/** Exponentiation: rop = a^e mod p */
void fp_pow(Fp *rop, const Fp *a, const mpz_t e);
/** Negation: rop = -a mod p */
void fp_neg(Fp *rop, const Fp *a);
/** Compare: returns 0 if equal */
int fp_cmp(const Fp *a, const Fp *b);

#endif // FP_H
