#include "fp3.h"
#include <stdio.h>

void fp3_init(Fp3 *a) {
    fp_init(&a->a0);
    fp_init(&a->a1);
    fp_init(&a->a2);
}

void fp3_clear(Fp3 *a) {
    fp_clear(&a->a0);
    fp_clear(&a->a1);
    fp_clear(&a->a2);
}

void fp3_set(Fp3 *rop, const Fp3 *op) {
    fp_set(&rop->a0, &op->a0);
    fp_set(&rop->a1, &op->a1);
    fp_set(&rop->a2, &op->a2);
}

void fp3_set_ui(Fp3 *rop, unsigned long v) {
    fp_set_ui(&rop->a0, v);
    fp_set_ui(&rop->a1, 0);
    fp_set_ui(&rop->a2, 0);
}

void fp3_random(Fp3 *rop) {
    fp_random(&rop->a0);
    fp_random(&rop->a1);
    fp_random(&rop->a2);
}

void fp3_print(const Fp3 *a) {
    printf("(");
    fp_print(&a->a0);
    printf(", ");
    fp_print(&a->a1);
    printf(", ");
    fp_print(&a->a2);
    printf(")");
}

void fp3_add(Fp3 *rop, const Fp3 *a, const Fp3 *b) {
    fp_add(&rop->a0, &a->a0, &b->a0);
    fp_add(&rop->a1, &a->a1, &b->a1);
    fp_add(&rop->a2, &a->a2, &b->a2);
}

void fp3_sub(Fp3 *rop, const Fp3 *a, const Fp3 *b) {
    fp_sub(&rop->a0, &a->a0, &b->a0);
    fp_sub(&rop->a1, &a->a1, &b->a1);
    fp_sub(&rop->a2, &a->a2, &b->a2);
}

// Karatsuba multiplication for Fp3: (a0 + a1*x + a2*x^2)*(b0 + b1*x + b2*x^2)
// with x^3 = xi, where xi is a non-residue in Fp (for KSS18, xi = 3)
void fp3_mul(Fp3 *rop, const Fp3 *a, const Fp3 *b) {
    Fp t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    fp_init(&t0); fp_init(&t1); fp_init(&t2); fp_init(&t3); fp_init(&t4);
    fp_init(&t5); fp_init(&t6); fp_init(&t7); fp_init(&t8); fp_init(&t9);
    // Karatsuba-like method
    // v0 = a0*b0
    fp_mul(&t0, &a->a0, &b->a0);
    // v1 = a1*b1
    fp_mul(&t1, &a->a1, &b->a1);
    // v2 = a2*b2
    fp_mul(&t2, &a->a2, &b->a2);
    // s1 = (a1 + a2)*(b1 + b2) - v1 - v2
    fp_add(&t3, &a->a1, &a->a2);
    fp_add(&t4, &b->a1, &b->a2);
    fp_mul(&t5, &t3, &t4);
    fp_sub(&t5, &t5, &t1);
    fp_sub(&t5, &t5, &t2);
    // s2 = (a0 + a1)*(b0 + b1) - v0 - v1
    fp_add(&t6, &a->a0, &a->a1);
    fp_add(&t7, &b->a0, &b->a1);
    fp_mul(&t8, &t6, &t7);
    fp_sub(&t8, &t8, &t0);
    fp_sub(&t8, &t8, &t1);
    // s3 = (a0 + a2)*(b0 + b2) - v0 - v2
    fp_add(&t6, &a->a0, &a->a2);
    fp_add(&t7, &b->a0, &b->a2);
    fp_mul(&t9, &t6, &t7);
    fp_sub(&t9, &t9, &t0);
    fp_sub(&t9, &t9, &t2);
    // rop->a0 = v0 + 3*s1
    fp_set_ui(&t6, 3);
    fp_mul(&t5, &t5, &t6);
    fp_add(&rop->a0, &t0, &t5);
    // rop->a1 = s2 + 3*v2
    fp_mul(&t6, &t2, &t6);
    fp_add(&rop->a1, &t8, &t6);
    // rop->a2 = s3 + t1
    fp_add(&rop->a2, &t9, &t1);
    fp_clear(&t0); fp_clear(&t1); fp_clear(&t2); fp_clear(&t3); fp_clear(&t4);
    fp_clear(&t5); fp_clear(&t6); fp_clear(&t7); fp_clear(&t8); fp_clear(&t9);
}

// Inversion in Fp3 using the method for cubic extensions
void fp3_inv(Fp3 *rop, const Fp3 *a) {
    // Let x = a0 + a1*x + a2*x^2, x^3 = xi
    // Compute the norm N(x) = a0^3 + xi*a1^3 + xi^2*a2^3 - 3*xi*a0*a1*a2
    Fp t0, t1, t2, t3, t4, norm, norm_inv;
    fp_init(&t0); fp_init(&t1); fp_init(&t2); fp_init(&t3); fp_init(&t4);
    fp_init(&norm); fp_init(&norm_inv);
    // t0 = a0^3
    fp_mul(&t0, &a->a0, &a->a0);
    fp_mul(&t0, &t0, &a->a0);
    // t1 = a1^3
    fp_mul(&t1, &a->a1, &a->a1);
    fp_mul(&t1, &t1, &a->a1);
    // t2 = a2^3
    fp_mul(&t2, &a->a2, &a->a2);
    fp_mul(&t2, &t2, &a->a2);
    // t3 = a0*a1*a2
    fp_mul(&t3, &a->a0, &a->a1);
    fp_mul(&t3, &t3, &a->a2);
    // norm = t0 + 3*t1 + 9*t2 - 9*t3
    fp_set_ui(&t4, 3);
    fp_mul(&t1, &t1, &t4); // 3*a1^3
    fp_set_ui(&t4, 9);
    fp_mul(&t2, &t2, &t4); // 9*a2^3
    fp_mul(&t3, &t3, &t4); // 9*a0*a1*a2
    fp_add(&norm, &t0, &t1);
    fp_add(&norm, &norm, &t2);
    fp_sub(&norm, &norm, &t3);
    // norm_inv = norm^{-1}
    fp_inv(&norm_inv, &norm);
    // Compute adjugate
    // rop->a0 = a0^2 - 3*a1*a2
    fp_mul(&t0, &a->a0, &a->a0);
    fp_mul(&t1, &a->a1, &a->a2);
    fp_set_ui(&t2, 3);
    fp_mul(&t1, &t1, &t2);
    fp_sub(&rop->a0, &t0, &t1);
    fp_mul(&rop->a0, &rop->a0, &norm_inv);
    // rop->a1 = a2^2 - 3*a0*a1
    fp_mul(&t0, &a->a2, &a->a2);
    fp_mul(&t1, &a->a0, &a->a1);
    fp_mul(&t1, &t1, &t2);
    fp_sub(&rop->a1, &t0, &t1);
    fp_mul(&rop->a1, &rop->a1, &norm_inv);
    // rop->a2 = a1^2 - 3*a0*a2
    fp_mul(&t0, &a->a1, &a->a1);
    fp_mul(&t1, &a->a0, &a->a2);
    fp_mul(&t1, &t1, &t2);
    fp_sub(&rop->a2, &t0, &t1);
    fp_mul(&rop->a2, &rop->a2, &norm_inv);
    fp_clear(&t0); fp_clear(&t1); fp_clear(&t2); fp_clear(&t3); fp_clear(&t4); fp_clear(&norm); fp_clear(&norm_inv);
}

void fp3_pow(Fp3 *rop, const Fp3 *a, const mpz_t e) {
    Fp3 res, base;
    fp3_init(&res); fp3_init(&base);
    fp3_set(&base, a);
    fp3_set_ui(&res, 1);
    mpz_t exp; mpz_init_set(exp, e);
    while (mpz_sgn(exp) > 0) {
        if (mpz_odd_p(exp)) {
            fp3_mul(&res, &res, &base);
        }
        fp3_mul(&base, &base, &base);
        mpz_fdiv_q_2exp(exp, exp, 1);
    }
    fp3_set(rop, &res);
    fp3_clear(&res); fp3_clear(&base); mpz_clear(exp);
}

void fp3_neg(Fp3 *rop, const Fp3 *a) {
    fp_neg(&rop->a0, &a->a0);
    fp_neg(&rop->a1, &a->a1);
    fp_neg(&rop->a2, &a->a2);
}

int fp3_cmp(const Fp3 *a, const Fp3 *b) {
    return fp_cmp(&a->a0, &b->a0) || fp_cmp(&a->a1, &b->a1) || fp_cmp(&a->a2, &b->a2);
}
