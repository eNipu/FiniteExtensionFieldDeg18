#include "fp18.h"
#include <stdio.h>

void fp18_init(Fp18 *a) {
    fp6_init(&a->a0);
    fp6_init(&a->a1);
    fp6_init(&a->a2);
}

void fp18_clear(Fp18 *a) {
    fp6_clear(&a->a0);
    fp6_clear(&a->a1);
    fp6_clear(&a->a2);
}

void fp18_set(Fp18 *rop, const Fp18 *op) {
    fp6_set(&rop->a0, &op->a0);
    fp6_set(&rop->a1, &op->a1);
    fp6_set(&rop->a2, &op->a2);
}

void fp18_set_ui(Fp18 *rop, unsigned long v) {
    fp6_set_ui(&rop->a0, v);
    fp6_set_ui(&rop->a1, 0);
    fp6_set_ui(&rop->a2, 0);
}

void fp18_random(Fp18 *rop) {
    fp6_random(&rop->a0);
    fp6_random(&rop->a1);
    fp6_random(&rop->a2);
}

void fp18_print(const Fp18 *a) {
    printf("(");
    fp6_print(&a->a0);
    printf(", ");
    fp6_print(&a->a1);
    printf(", ");
    fp6_print(&a->a2);
    printf(")");
}

void fp18_add(Fp18 *rop, const Fp18 *a, const Fp18 *b) {
    fp6_add(&rop->a0, &a->a0, &b->a0);
    fp6_add(&rop->a1, &a->a1, &b->a1);
    fp6_add(&rop->a2, &a->a2, &b->a2);
}

void fp18_sub(Fp18 *rop, const Fp18 *a, const Fp18 *b) {
    fp6_sub(&rop->a0, &a->a0, &b->a0);
    fp6_sub(&rop->a1, &a->a1, &b->a1);
    fp6_sub(&rop->a2, &a->a2, &b->a2);
}

// Multiplication for Fp18: (a0 + a1*w + a2*w^2)*(b0 + b1*w + b2*w^2), w^3 = delta, delta = 3 for KSS18
void fp18_mul(Fp18 *rop, const Fp18 *a, const Fp18 *b) {
    Fp6 t0, t1, t2, t3, t4, t5, t6;
    fp6_init(&t0); fp6_init(&t1); fp6_init(&t2); fp6_init(&t3); fp6_init(&t4); fp6_init(&t5); fp6_init(&t6);
    // v0 = a0*b0
    fp6_mul(&t0, &a->a0, &b->a0);
    // v1 = a1*b1
    fp6_mul(&t1, &a->a1, &b->a1);
    // v2 = a2*b2
    fp6_mul(&t2, &a->a2, &b->a2);
    // s1 = (a1 + a2)*(b1 + b2) - v1 - v2
    fp6_add(&t3, &a->a1, &a->a2);
    fp6_add(&t4, &b->a1, &b->a2);
    fp6_mul(&t5, &t3, &t4);
    fp6_sub(&t5, &t5, &t1);
    fp6_sub(&t5, &t5, &t2);
    // s2 = (a0 + a1)*(b0 + b1) - v0 - v1
    fp6_add(&t3, &a->a0, &a->a1);
    fp6_add(&t4, &b->a0, &b->a1);
    fp6_mul(&t6, &t3, &t4);
    fp6_sub(&t6, &t6, &t0);
    fp6_sub(&t6, &t6, &t1);
    // s3 = (a0 + a2)*(b0 + b2) - v0 - v2
    fp6_add(&t3, &a->a0, &a->a2);
    fp6_add(&t4, &b->a0, &b->a2);
    fp6_mul(&t3, &t3, &t4);
    fp6_sub(&t3, &t3, &t0);
    fp6_sub(&t3, &t3, &t2);
    // rop->a0 = v0 + 3*s1
    fp6_set_ui(&t4, 3);
    fp6_mul(&t5, &t5, &t4);
    fp6_add(&rop->a0, &t0, &t5);
    // rop->a1 = s2 + 3*v2
    fp6_mul(&t4, &t2, &t4);
    fp6_add(&rop->a1, &t6, &t4);
    // rop->a2 = s3 + t1
    fp6_add(&rop->a2, &t3, &t1);
    fp6_clear(&t0); fp6_clear(&t1); fp6_clear(&t2); fp6_clear(&t3); fp6_clear(&t4); fp6_clear(&t5); fp6_clear(&t6);
}

void fp18_inv(Fp18 *rop, const Fp18 *a) {
    // Inversion in Fp18 is non-trivial; use standard formula for cubic extensions
    // TODO: Optimize for KSS18 if needed
    // Placeholder: set to zero
    fp18_set_ui(rop, 0);
}

void fp18_pow(Fp18 *rop, const Fp18 *a, const mpz_t e) {
    Fp18 res, base;
    fp18_init(&res); fp18_init(&base);
    fp18_set(&base, a);
    fp18_set_ui(&res, 1);
    mpz_t exp; mpz_init_set(exp, e);
    while (mpz_sgn(exp) > 0) {
        if (mpz_odd_p(exp)) {
            fp18_mul(&res, &res, &base);
        }
        fp18_mul(&base, &base, &base);
        mpz_fdiv_q_2exp(exp, exp, 1);
    }
    fp18_set(rop, &res);
    fp18_clear(&res); fp18_clear(&base); mpz_clear(exp);
}

void fp18_neg(Fp18 *rop, const Fp18 *a) {
    fp6_neg(&rop->a0, &a->a0);
    fp6_neg(&rop->a1, &a->a1);
    fp6_neg(&rop->a2, &a->a2);
}

int fp18_cmp(const Fp18 *a, const Fp18 *b) {
    return fp6_cmp(&a->a0, &b->a0) || fp6_cmp(&a->a1, &b->a1) || fp6_cmp(&a->a2, &b->a2);
}
