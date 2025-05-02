#include "ec.h"
#include <string.h>

// --- Fp curve ops ---
void ecfp_init(EcFp *p) {
    fp_init(&p->x);
    fp_init(&p->y);
    p->infinity = 1;
}
void ecfp_clear(EcFp *p) {
    fp_clear(&p->x);
    fp_clear(&p->y);
}
void ecfp_set(EcFp *rop, const EcFp *op) {
    fp_set(&rop->x, &op->x);
    fp_set(&rop->y, &op->y);
    rop->infinity = op->infinity;
}
void ecfp_set_infinity(EcFp *p) {
    fp_set_ui(&p->x, 0);
    fp_set_ui(&p->y, 0);
    p->infinity = 1;
}
int ecfp_is_infinity(const EcFp *p) {
    return p->infinity;
}
int ecfp_cmp(const EcFp *p, const EcFp *q) {
    return fp_cmp(&p->x, &q->x) || fp_cmp(&p->y, &q->y) || (p->infinity != q->infinity);
}
void ecfp_add(EcFp *rop, const EcFp *p, const EcFp *q) {
    if (ecfp_is_infinity(p)) {
        ecfp_set(rop, q);
        return;
    }
    if (ecfp_is_infinity(q)) {
        ecfp_set(rop, p);
        return;
    }
    if (fp_cmp(&p->x, &q->x) == 0) {
        if (fp_cmp(&p->y, &q->y) != 0 || fp_cmp(&p->y, &(Fp){0}) == 0) {
            ecfp_set_infinity(rop);
            return;
        } else {
            ecfp_double(rop, p);
            return;
        }
    }
    Fp lambda, t;
    fp_init(&lambda); fp_init(&t);
    fp_sub(&t, &q->y, &p->y);
    Fp denom;
    fp_init(&denom);
    fp_sub(&denom, &q->x, &p->x);
    fp_inv(&denom, &denom);
    fp_mul(&lambda, &t, &denom);
    fp_mul(&t, &lambda, &lambda);
    fp_sub(&t, &t, &p->x);
    fp_sub(&t, &t, &q->x);
    fp_set(&rop->x, &t);
    fp_sub(&t, &p->x, &rop->x);
    fp_mul(&t, &t, &lambda);
    fp_sub(&t, &t, &p->y);
    fp_set(&rop->y, &t);
    rop->infinity = 0;
    fp_clear(&lambda); fp_clear(&t); fp_clear(&denom);
}

void ecfp_double(EcFp *rop, const EcFp *p) {
    if (ecfp_is_infinity(p) || fp_cmp(&p->y, &(Fp){0}) == 0) {
        ecfp_set_infinity(rop);
        return;
    }
    Fp lambda, t;
    fp_init(&lambda); fp_init(&t);
    Fp denom;
    fp_init(&denom);
    fp_set_ui(&t, 3);
    fp_mul(&t, &t, &p->x);
    fp_mul(&t, &t, &p->x);
    fp_set_ui(&denom, 2);
    fp_mul(&denom, &denom, &p->y);
    fp_inv(&denom, &denom);
    fp_mul(&lambda, &t, &denom);
    fp_mul(&t, &lambda, &lambda);
    fp_sub(&t, &t, &p->x);
    fp_sub(&t, &t, &p->x);
    fp_set(&rop->x, &t);
    fp_sub(&t, &p->x, &rop->x);
    fp_mul(&t, &t, &lambda);
    fp_sub(&t, &t, &p->y);
    fp_set(&rop->y, &t);
    rop->infinity = 0;
    fp_clear(&lambda); fp_clear(&t); fp_clear(&denom);
}

void ecfp_scalar_mul(EcFp *rop, const EcFp *p, const mpz_t scalar) {
    EcFp res, tmp;
    ecfp_init(&res); ecfp_init(&tmp);
    ecfp_set_infinity(&res);
    size_t n = mpz_sizeinbase(scalar, 2);
    for (ssize_t i = n - 1; i >= 0; --i) {
        ecfp_double(&res, &res);
        if (mpz_tstbit(scalar, i)) {
            ecfp_add(&res, &res, p);
        }
    }
    ecfp_set(rop, &res);
    ecfp_clear(&res); ecfp_clear(&tmp);
}

// --- Fp3 curve ops ---
void ecfp3_init(EcFp3 *p) {
    fp3_init(&p->x);
    fp3_init(&p->y);
    p->infinity = 1;
}
void ecfp3_clear(EcFp3 *p) {
    fp3_clear(&p->x);
    fp3_clear(&p->y);
}
void ecfp3_set(EcFp3 *rop, const EcFp3 *op) {
    fp3_set(&rop->x, &op->x);
    fp3_set(&rop->y, &op->y);
    rop->infinity = op->infinity;
}
void ecfp3_set_infinity(EcFp3 *p) {
    fp3_set_ui(&p->x, 0);
    fp3_set_ui(&p->y, 0);
    p->infinity = 1;
}
int ecfp3_is_infinity(const EcFp3 *p) {
    return p->infinity;
}
int ecfp3_cmp(const EcFp3 *p, const EcFp3 *q) {
    return fp3_cmp(&p->x, &q->x) || fp3_cmp(&p->y, &q->y) || (p->infinity != q->infinity);
}
void ecfp3_add(EcFp3 *rop, const EcFp3 *p, const EcFp3 *q) {
    if (ecfp3_is_infinity(p)) {
        ecfp3_set(rop, q);
        return;
    }
    if (ecfp3_is_infinity(q)) {
        ecfp3_set(rop, p);
        return;
    }
    if (fp3_cmp(&p->x, &q->x) == 0) {
        if (fp3_cmp(&p->y, &q->y) != 0) {
            ecfp3_set_infinity(rop);
            return;
        } else {
            ecfp3_double(rop, p);
            return;
        }
    }
    Fp3 lambda, t;
    fp3_init(&lambda); fp3_init(&t);
    fp3_sub(&t, &q->y, &p->y);
    Fp3 denom;
    fp3_init(&denom);
    fp3_sub(&denom, &q->x, &p->x);
    fp3_inv(&denom, &denom);
    fp3_mul(&lambda, &t, &denom);
    fp3_mul(&t, &lambda, &lambda);
    fp3_sub(&t, &t, &p->x);
    fp3_sub(&t, &t, &q->x);
    fp3_set(&rop->x, &t);
    fp3_sub(&t, &p->x, &rop->x);
    fp3_mul(&t, &t, &lambda);
    fp3_sub(&t, &t, &p->y);
    fp3_set(&rop->y, &t);
    rop->infinity = 0;
    fp3_clear(&lambda); fp3_clear(&t); fp3_clear(&denom);
}
void ecfp3_double(EcFp3 *rop, const EcFp3 *p) {
    if (ecfp3_is_infinity(p)) {
        ecfp3_set_infinity(rop);
        return;
    }
    Fp3 lambda, t;
    fp3_init(&lambda); fp3_init(&t);
    Fp3 denom;
    fp3_init(&denom);
    fp3_set_ui(&t, 3);
    fp3_mul(&t, &t, &p->x);
    fp3_mul(&t, &t, &p->x);
    fp3_set_ui(&denom, 2);
    fp3_mul(&denom, &denom, &p->y);
    fp3_inv(&denom, &denom);
    fp3_mul(&lambda, &t, &denom);
    fp3_mul(&t, &lambda, &lambda);
    fp3_sub(&t, &t, &p->x);
    fp3_sub(&t, &t, &p->x);
    fp3_set(&rop->x, &t);
    fp3_sub(&t, &p->x, &rop->x);
    fp3_mul(&t, &t, &lambda);
    fp3_sub(&t, &t, &p->y);
    fp3_set(&rop->y, &t);
    rop->infinity = 0;
    fp3_clear(&lambda); fp3_clear(&t); fp3_clear(&denom);
}
void ecfp3_scalar_mul(EcFp3 *rop, const EcFp3 *p, const mpz_t scalar) {
    EcFp3 res;
    ecfp3_init(&res);
    ecfp3_set_infinity(&res);
    size_t n = mpz_sizeinbase(scalar, 2);
    for (ssize_t i = n - 1; i >= 0; --i) {
        ecfp3_double(&res, &res);
        if (mpz_tstbit(scalar, i)) {
            ecfp3_add(&res, &res, p);
        }
    }
    ecfp3_set(rop, &res);
    ecfp3_clear(&res);
}

// --- Fp18 curve ops ---
void ecfp18_init(EcFp18 *p) {
    fp18_init(&p->x);
    fp18_init(&p->y);
    p->infinity = 1;
}
void ecfp18_clear(EcFp18 *p) {
    fp18_clear(&p->x);
    fp18_clear(&p->y);
}
void ecfp18_set(EcFp18 *rop, const EcFp18 *op) {
    fp18_set(&rop->x, &op->x);
    fp18_set(&rop->y, &op->y);
    rop->infinity = op->infinity;
}
void ecfp18_set_infinity(EcFp18 *p) {
    fp18_set_ui(&p->x, 0);
    fp18_set_ui(&p->y, 0);
    p->infinity = 1;
}
int ecfp18_is_infinity(const EcFp18 *p) {
    return p->infinity;
}
int ecfp18_cmp(const EcFp18 *p, const EcFp18 *q) {
    return fp18_cmp(&p->x, &q->x) || fp18_cmp(&p->y, &q->y) || (p->infinity != q->infinity);
}
void ecfp18_add(EcFp18 *rop, const EcFp18 *p, const EcFp18 *q) {
    if (ecfp18_is_infinity(p)) {
        ecfp18_set(rop, q);
        return;
    }
    if (ecfp18_is_infinity(q)) {
        ecfp18_set(rop, p);
        return;
    }
    if (fp18_cmp(&p->x, &q->x) == 0) {
        if (fp18_cmp(&p->y, &q->y) != 0) {
            ecfp18_set_infinity(rop);
            return;
        } else {
            ecfp18_double(rop, p);
            return;
        }
    }
    Fp18 lambda, t;
    fp18_init(&lambda); fp18_init(&t);
    fp18_sub(&t, &q->y, &p->y);
    Fp18 denom;
    fp18_init(&denom);
    fp18_sub(&denom, &q->x, &p->x);
    fp18_inv(&denom, &denom);
    fp18_mul(&lambda, &t, &denom);
    fp18_mul(&t, &lambda, &lambda);
    fp18_sub(&t, &t, &p->x);
    fp18_sub(&t, &t, &q->x);
    fp18_set(&rop->x, &t);
    fp18_sub(&t, &p->x, &rop->x);
    fp18_mul(&t, &t, &lambda);
    fp18_sub(&t, &t, &p->y);
    fp18_set(&rop->y, &t);
    rop->infinity = 0;
    fp18_clear(&lambda); fp18_clear(&t); fp18_clear(&denom);
}
void ecfp18_double(EcFp18 *rop, const EcFp18 *p) {
    if (ecfp18_is_infinity(p)) {
        ecfp18_set_infinity(rop);
        return;
    }
    Fp18 lambda, t;
    fp18_init(&lambda); fp18_init(&t);
    Fp18 denom;
    fp18_init(&denom);
    fp18_set_ui(&t, 3);
    fp18_mul(&t, &t, &p->x);
    fp18_mul(&t, &t, &p->x);
    fp18_set_ui(&denom, 2);
    fp18_mul(&denom, &denom, &p->y);
    fp18_inv(&denom, &denom);
    fp18_mul(&lambda, &t, &denom);
    fp18_mul(&t, &lambda, &lambda);
    fp18_sub(&t, &t, &p->x);
    fp18_sub(&t, &t, &p->x);
    fp18_set(&rop->x, &t);
    fp18_sub(&t, &p->x, &rop->x);
    fp18_mul(&t, &t, &lambda);
    fp18_sub(&t, &t, &p->y);
    fp18_set(&rop->y, &t);
    rop->infinity = 0;
    fp18_clear(&lambda); fp18_clear(&t); fp18_clear(&denom);
}
void ecfp18_scalar_mul(EcFp18 *rop, const EcFp18 *p, const mpz_t scalar) {
    EcFp18 res;
    ecfp18_init(&res);
    ecfp18_set_infinity(&res);
    size_t n = mpz_sizeinbase(scalar, 2);
    for (ssize_t i = n - 1; i >= 0; --i) {
        ecfp18_double(&res, &res);
        if (mpz_tstbit(scalar, i)) {
            ecfp18_add(&res, &res, p);
        }
    }
    ecfp18_set(rop, &res);
    ecfp18_clear(&res);
}
