#define DEFINE_GLOBAL_VARIABLES
#include "parameters.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

// Define global GMP variables declared in parameters.h
mpz_t kss18_p;
mpz_t kss18_r;
mpz_t kss18_b;
mpz_t kss18_t; // Trace of Frobenius (calculated)
mpz_t kss18_order_efp; // Order of E(Fp) = p + 1 - t (calculated)

// Optimization counters (keep for now, might be refactored)
int fp_add_count = 0;
int fp_mul_count = 0;
int inv_count = 0;

/**
 * Initialize global parameters for KSS18 curve.
 * Reads p, r, B from defines in parameters.h and calculates t.
 */
void init_kss18_params(void) {
    // Initialize p, r, B using the inline function from the header
    kss18_init_globals();

    // Calculate u from r (or p) to find t
    // r = (u^6 + 37u^3 + 343) / 343 => 343r = u^6 + 37u^3 + 343
    // p = (u^8 + 5u^7 + ... + 2401) / 21
    // t = (u^4 + 16u + 7) / 7
    // For simplicity, let's hardcode u for now, as deriving it is complex.
    // u = 2^44 + 2^22 - 2^9 + 2 = 17592190238210
    mpz_t u, u_pow_4, tmp1, tmp2;
    mpz_init_set_str(u, "17592190238210", 10);
    mpz_init(u_pow_4);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(kss18_t);
    mpz_init(kss18_order_efp);

    // Calculate t = (u^4 + 16u + 7) / 7
    mpz_pow_ui(u_pow_4, u, 4);
    mpz_mul_ui(tmp1, u, 16);
    mpz_add(tmp2, u_pow_4, tmp1);
    mpz_add_ui(tmp2, tmp2, 7);
    mpz_fdiv_q_ui(kss18_t, tmp2, 7); // Use fdiv for integer division

    // Calculate #E(Fp) = p + 1 - t
    mpz_add_ui(kss18_order_efp, kss18_p, 1);
    mpz_sub(kss18_order_efp, kss18_order_efp, kss18_t);

    // Clear temporary variables
    mpz_clear(u);
    mpz_clear(u_pow_4);
    mpz_clear(tmp1);
    mpz_clear(tmp2);

    // Reset counters
    fp_add_count = 0;
    fp_mul_count = 0;
    inv_count = 0;
}

/**
 * Clear global parameters.
 */
void clear_kss18_params(void) {
    kss18_clear_globals(); // Clears p, r, b using inline function
    mpz_clear(kss18_t);
    mpz_clear(kss18_order_efp);
}
