#ifndef FP6_H
#define FP6_H

#include "fp3.h"

/**
 * @file fp6.h
 * @brief Sextic extension field Fp6 arithmetic for KSS18
 */

typedef struct {
    Fp3 a0, a1;
} Fp6;

void fp6_init(Fp6 *a);
void fp6_clear(Fp6 *a);
void fp6_set(Fp6 *rop, const Fp6 *op);
void fp6_set_ui(Fp6 *rop, unsigned long v);
void fp6_random(Fp6 *rop);
void fp6_print(const Fp6 *a);
void fp6_add(Fp6 *rop, const Fp6 *a, const Fp6 *b);
void fp6_sub(Fp6 *rop, const Fp6 *a, const Fp6 *b);
void fp6_mul(Fp6 *rop, const Fp6 *a, const Fp6 *b);
void fp6_inv(Fp6 *rop, const Fp6 *a);
void fp6_pow(Fp6 *rop, const Fp6 *a, const mpz_t e);
void fp6_neg(Fp6 *rop, const Fp6 *a);
int fp6_cmp(const Fp6 *a, const Fp6 *b);

#endif // FP6_H
