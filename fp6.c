#include "fp6.h"
#include <stdio.h>

void fp6_init(Fp6 *a) {
    fp3_init(&a->a0);
    fp3_init(&a->a1);
}

void fp6_clear(Fp6 *a) {
    fp3_clear(&a->a0);
    fp3_clear(&a->a1);
}

void fp6_set(Fp6 *rop, const Fp6 *op) {
    fp3_set(&rop->a0, &op->a0);
    fp3_set(&rop->a1, &op->a1);
}

void fp6_set_ui(Fp6 *rop, unsigned long v) {
    fp3_set_ui(&rop->a0, v);
    fp3_set_ui(&rop->a1, 0);
}

void fp6_random(Fp6 *rop) {
    fp3_random(&rop->a0);
    fp3_random(&rop->a1);
}

void fp6_print(const Fp6 *a) {
    printf("(");
    fp3_print(&a->a0);
    printf(", ");
    fp3_print(&a->a1);
    printf(")");
}

void fp6_add(Fp6 *rop, const Fp6 *a, const Fp6 *b) {
    fp3_add(&rop->a0, &a->a0, &b->a0);
    fp3_add(&rop->a1, &a->a1, &b->a1);
}

void fp6_sub(Fp6 *rop, const Fp6 *a, const Fp6 *b) {
    fp3_sub(&rop->a0, &a->a0, &b->a0);
    fp3_sub(&rop->a1, &a->a1, &b->a1);
}

// Multiplication for Fp6: (a0 + a1*v)*(b0 + b1*v), v^2 = xi, xi = 3 for KSS18
void fp6_mul(Fp6 *rop, const Fp6 *a, const Fp6 *b) {
    Fp3 t0, t1, t2, t3;
    fp3_init(&t0); fp3_init(&t1); fp3_init(&t2); fp3_init(&t3);
    // t0 = a0 * b0
    fp3_mul(&t0, &a->a0, &b->a0);
    // t1 = a1 * b1
    fp3_mul(&t1, &a->a1, &b->a1);
    // rop->a0 = t0 + xi * t1
    fp3_set_ui(&t2, 3); // xi = 3
    fp3_mul(&t3, &t1, &t2);
    fp3_add(&rop->a0, &t0, &t3);
    // rop->a1 = (a0 + a1)*(b0 + b1) - t0 - t1
    fp3_add(&t2, &a->a0, &a->a1);
    fp3_add(&t3, &b->a0, &b->a1);
    fp3_mul(&t2, &t2, &t3);
    fp3_sub(&t2, &t2, &t0);
    fp3_sub(&rop->a1, &t2, &t1);
    fp3_clear(&t0); fp3_clear(&t1); fp3_clear(&t2); fp3_clear(&t3);
}

void fp6_inv(Fp6 *rop, const Fp6 *a) {
    // Inversion in Fp6: (a0 + a1*v)^{-1} = (a0 - a1*v) / (a0^2 - xi*a1^2)
    Fp3 t0, t1, t2, t3;
    fp3_init(&t0); fp3_init(&t1); fp3_init(&t2); fp3_init(&t3);
    // t0 = a0^2
    fp3_mul(&t0, &a->a0, &a->a0);
    // t1 = a1^2
    fp3_mul(&t1, &a->a1, &a->a1);
    // t2 = xi * t1
    fp3_set_ui(&t2, 3);
    fp3_mul(&t2, &t2, &t1);
    // denom = t0 - t2
    fp3_sub(&t3, &t0, &t2);
    fp3_inv(&t3, &t3);
    // rop->a0 = a0 * denom
    fp3_mul(&rop->a0, &a->a0, &t3);
    // rop->a1 = -a1 * denom
    fp3_neg(&t0, &a->a1);
    fp3_mul(&rop->a1, &t0, &t3);
    fp3_clear(&t0); fp3_clear(&t1); fp3_clear(&t2); fp3_clear(&t3);
}

void fp6_pow(Fp6 *rop, const Fp6 *a, const mpz_t e) {
    Fp6 res, base;
    fp6_init(&res); fp6_init(&base);
    fp6_set(&base, a);
    fp6_set_ui(&res, 1);
    mpz_t exp; mpz_init_set(exp, e);
    while (mpz_sgn(exp) > 0) {
        if (mpz_odd_p(exp)) {
            fp6_mul(&res, &res, &base);
        }
        fp6_mul(&base, &base, &base);
        mpz_fdiv_q_2exp(exp, exp, 1);
    }
    fp6_set(rop, &res);
    fp6_clear(&res); fp6_clear(&base); mpz_clear(exp);
}

void fp6_neg(Fp6 *rop, const Fp6 *a) {
    fp3_neg(&rop->a0, &a->a0);
    fp3_neg(&rop->a1, &a->a1);
}

int fp6_cmp(const Fp6 *a, const Fp6 *b) {
    return fp3_cmp(&a->a0, &b->a0) || fp3_cmp(&a->a1, &b->a1);
}
