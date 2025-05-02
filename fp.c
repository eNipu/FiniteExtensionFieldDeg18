#include "fp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static gmp_randstate_t fp_rand_state;
static int fp_rand_initialized = 0;

static void fp_init_rand() {
    if (!fp_rand_initialized) {
        gmp_randinit_default(fp_rand_state);
        gmp_randseed_ui(fp_rand_state, (unsigned long)time(NULL));
        fp_rand_initialized = 1;
    }
}

void fp_init(Fp *a) {
    mpz_init(a->x);
}

void fp_clear(Fp *a) {
    mpz_clear(a->x);
}

void fp_set(Fp *rop, const Fp *op) {
    mpz_set(rop->x, op->x);
}

void fp_set_ui(Fp *rop, unsigned long v) {
    mpz_set_ui(rop->x, v);
    mpz_mod(rop->x, rop->x, kss18_p);
}

void fp_set_mpz(Fp *rop, const mpz_t v) {
    mpz_set(rop->x, v);
    mpz_mod(rop->x, rop->x, kss18_p);
}

void fp_random(Fp *rop) {
    fp_init_rand();
    mpz_urandomm(rop->x, fp_rand_state, kss18_p);
}

void fp_print(const Fp *a) {
    gmp_printf("%Zd", a->x);
}

void fp_add(Fp *rop, const Fp *a, const Fp *b) {
    mpz_add(rop->x, a->x, b->x);
    mpz_mod(rop->x, rop->x, kss18_p);
}

void fp_sub(Fp *rop, const Fp *a, const Fp *b) {
    mpz_sub(rop->x, a->x, b->x);
    mpz_mod(rop->x, rop->x, kss18_p);
}

void fp_mul(Fp *rop, const Fp *a, const Fp *b) {
    mpz_mul(rop->x, a->x, b->x);
    mpz_mod(rop->x, rop->x, kss18_p);
}

void fp_inv(Fp *rop, const Fp *a) {
    mpz_invert(rop->x, a->x, kss18_p);
}

void fp_pow(Fp *rop, const Fp *a, const mpz_t e) {
    mpz_powm(rop->x, a->x, e, kss18_p);
}

void fp_neg(Fp *rop, const Fp *a) {
    if (mpz_sgn(a->x) == 0) {
        mpz_set_ui(rop->x, 0);
    } else {
        mpz_sub(rop->x, kss18_p, a->x);
    }
}

int fp_cmp(const Fp *a, const Fp *b) {
    return mpz_cmp(a->x, b->x);
}
