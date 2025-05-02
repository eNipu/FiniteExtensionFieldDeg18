#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

// Include the main header file for structure definitions
#include "embedding_degree18.h"

// Declaration of global parameters for KSS degree 18 curve
extern mpz_t X;       // Parameter for the KSS curve
extern mpz_t prime;   // Field characteristic
extern mpz_t r_order; // Order of the subgroup
extern mpz_t t_trace; // Trace of Frobenius
extern mpz_t r_order_EFp; // Order of EFp
extern mpz_t b;       // Curve constant: y^2 = x^3 + b

// Constants used in the implementation
extern mpz_t c1_leg;
extern mpz_t c1_leg_bar;
extern mpz_t c1_omega;
extern mpz_t c1_omega_bar;

// Binary representation of X for efficient scalar multiplication
extern int *X_bit_binary;
extern int X_bit;

// Optimization counters
extern int add_count_miller;
extern int add_count_finalexp;
extern int sqr_count_miller;
extern int sqr_count_finalexp;
extern int inv_count;

// Function to initialize the parameters
void init_parameters(void);

// Function to generate curve parameters
void generate_parameters(void);

#endif // PARAMETERS_H
