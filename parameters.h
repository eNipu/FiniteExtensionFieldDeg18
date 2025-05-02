#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <gmp.h>

// Declare global parameters as extern
// These will be defined in parameters.c
extern mpz_t X; // Parameter for generating p and r
extern mpz_t prime; // Finite field prime p
extern mpz_t r_order; // Order r of the subgroup G1/G2
extern mpz_t t_trace; // Trace of Frobenius
extern mpz_t r_order_EFp; // Order of the elliptic curve group E(Fp) = p + 1 - t
extern mpz_t b; // Elliptic curve equation y^2 = x^3 + b parameter

// Frobenius calculation related parameters (if needed, otherwise remove)
// extern mpz_t c1_leg, c1_leg_bar, c1_omega, c1_omega_bar;

// Function prototypes for parameter management
void init_parameters(); // Initialize GMP variables
void generate_parameters(); // Generate curve parameters based on X
void clear_parameters(); // Clear GMP variables

#endif // PARAMETERS_H
