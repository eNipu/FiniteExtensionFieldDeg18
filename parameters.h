#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <gmp.h>

// Define the EXTERN macro to handle global variable declarations
#ifdef DEFINE_GLOBAL_VARIABLES
  #define EXTERN
#else
  #define EXTERN extern
#endif

// Global variables for curve parameters
EXTERN mpz_t X;
EXTERN mpz_t prime;
EXTERN mpz_t order;
EXTERN mpz_t trace;
EXTERN mpz_t order_EFp;
EXTERN mpz_t b;

// Common constants
EXTERN unsigned long int c1;

// Function declarations
void EFp_set_EC_parameter(void);
void generate_X(void);

#endif // PARAMETERS_H
