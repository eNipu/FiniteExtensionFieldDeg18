#ifndef EC_H
#define EC_H

#include "fp.h"
#include "fp3.h"
#include "fp6.h"
#include "fp18.h"
#include <gmp.h>

/**
 * @file ec.h
 * @brief Elliptic curve operations over Fp and extensions for KSS18
 */

/**
 * @struct EcFp
 * @brief Point on E: y^2 = x^3 + 3 over Fp
 */
typedef struct {
    Fp x, y;
    int infinity;
} EcFp;

/**
 * @struct EcFp3
 * @brief Point on E' (twist) over Fp3
 */
typedef struct {
    Fp3 x, y;
    int infinity;
} EcFp3;

/**
 * @struct EcFp6
 * @brief Point on E over Fp6 (optional, for completeness)
 */
typedef struct {
    Fp6 x, y;
    int infinity;
} EcFp6;

/**
 * @struct EcFp18
 * @brief Point on E over Fp18 (for pairing target group)
 */
typedef struct {
    Fp18 x, y;
    int infinity;
} EcFp18;

// --- Fp curve ops ---
void ecfp_init(EcFp *p);
void ecfp_clear(EcFp *p);
void ecfp_set(EcFp *rop, const EcFp *op);
void ecfp_set_infinity(EcFp *p);
int ecfp_is_infinity(const EcFp *p);
void ecfp_add(EcFp *rop, const EcFp *p, const EcFp *q);
void ecfp_double(EcFp *rop, const EcFp *p);
void ecfp_scalar_mul(EcFp *rop, const EcFp *p, const mpz_t scalar);
int ecfp_cmp(const EcFp *p, const EcFp *q);

// --- Fp3 curve ops ---
void ecfp3_init(EcFp3 *p);
void ecfp3_clear(EcFp3 *p);
void ecfp3_set(EcFp3 *rop, const EcFp3 *op);
void ecfp3_set_infinity(EcFp3 *p);
int ecfp3_is_infinity(const EcFp3 *p);
void ecfp3_add(EcFp3 *rop, const EcFp3 *p, const EcFp3 *q);
void ecfp3_double(EcFp3 *rop, const EcFp3 *p);
void ecfp3_scalar_mul(EcFp3 *rop, const EcFp3 *p, const mpz_t scalar);
int ecfp3_cmp(const EcFp3 *p, const EcFp3 *q);

// --- Fp18 curve ops (for completeness) ---
void ecfp18_init(EcFp18 *p);
void ecfp18_clear(EcFp18 *p);
void ecfp18_set(EcFp18 *rop, const EcFp18 *op);
void ecfp18_set_infinity(EcFp18 *p);
int ecfp18_is_infinity(const EcFp18 *p);
void ecfp18_add(EcFp18 *rop, const EcFp18 *p, const EcFp18 *q);
void ecfp18_double(EcFp18 *rop, const EcFp18 *p);
void ecfp18_scalar_mul(EcFp18 *rop, const EcFp18 *p, const mpz_t scalar);
int ecfp18_cmp(const EcFp18 *p, const EcFp18 *q);

#endif // EC_H
