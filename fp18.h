#ifndef FP18_H
#define FP18_H

#include "fp6.h"

/**
 * @file fp18.h
 * @brief Degree-18 extension field Fp18 arithmetic for KSS18
 */

typedef struct {
    Fp6 a0, a1, a2;
} Fp18;

void fp18_init(Fp18 *a);
void fp18_clear(Fp18 *a);
void fp18_set(Fp18 *rop, const Fp18 *op);
void fp18_set_ui(Fp18 *rop, unsigned long v);
void fp18_random(Fp18 *rop);
void fp18_print(const Fp18 *a);
void fp18_add(Fp18 *rop, const Fp18 *a, const Fp18 *b);
void fp18_sub(Fp18 *rop, const Fp18 *a, const Fp18 *b);
void fp18_mul(Fp18 *rop, const Fp18 *a, const Fp18 *b);
void fp18_inv(Fp18 *rop, const Fp18 *a);
void fp18_pow(Fp18 *rop, const Fp18 *a, const mpz_t e);
void fp18_neg(Fp18 *rop, const Fp18 *a);
int fp18_cmp(const Fp18 *a, const Fp18 *b);

#endif // FP18_H
