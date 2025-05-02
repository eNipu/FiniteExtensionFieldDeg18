#ifndef FP3_H
#define FP3_H

#include "fp.h"

/**
 * @file fp3.h
 * @brief Cubic extension field Fp3 arithmetic for KSS18
 */

typedef struct {
    Fp a0, a1, a2;
} Fp3;

void fp3_init(Fp3 *a);
void fp3_clear(Fp3 *a);
void fp3_set(Fp3 *rop, const Fp3 *op);
void fp3_set_ui(Fp3 *rop, unsigned long v);
void fp3_random(Fp3 *rop);
void fp3_print(const Fp3 *a);
void fp3_add(Fp3 *rop, const Fp3 *a, const Fp3 *b);
void fp3_sub(Fp3 *rop, const Fp3 *a, const Fp3 *b);
void fp3_mul(Fp3 *rop, const Fp3 *a, const Fp3 *b);
void fp3_inv(Fp3 *rop, const Fp3 *a);
void fp3_pow(Fp3 *rop, const Fp3 *a, const mpz_t e);
void fp3_neg(Fp3 *rop, const Fp3 *a);
int fp3_cmp(const Fp3 *a, const Fp3 *b);

#endif // FP3_H
