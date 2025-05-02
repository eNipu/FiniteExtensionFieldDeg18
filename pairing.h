#ifndef PAIRING_H
#define PAIRING_H

#include "ec.h"
#include "fp18.h"
#include <gmp.h>

/**
 * @file pairing.h
 * @brief Pairing computations (e.g., Optimal Ate) for KSS18
 */

/**
 * Compute the Optimal Ate pairing e(P, Q) for P in E(Fp)[r], Q in E'(Fp^3)[r].
 * @param result Output: Fp18 element (GT)
 * @param p Point in E(Fp)[r]
 * @param q Point in E'(Fp3)[r]
 */
void optimal_ate_pairing(Fp18 *result, const EcFp *p, const EcFp3 *q);

#endif // PAIRING_H
