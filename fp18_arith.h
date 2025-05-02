#ifndef FP18_ARITH_H
#define FP18_ARITH_H

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// Define constants (Consider moving these to parameters.h or defining appropriately)
#define C1 2
#define C2 -1
#define C1_sq 4

// Define boolean type if not available
#ifndef __cplusplus
typedef enum { FALSE, TRUE } bool;
#endif

// Structure Definitions

// Base Field Fp
struct Fp{
    mpz_t x0; // Using x0 consistent with Finite_Field.c
};

// Extension Field Fp3 (Fp[omega]/(omega^3 - c1))
struct Fp3{
    struct Fp x0, x1, x2; // Using x0, x1, x2 consistent with Finite_Field.c
};

// Extension Field Fp6 (Fp3[tau]/(tau^2 - xi)) where xi is omega
struct Fp6{
    struct Fp3 x0, x1; // Using x0, x1 consistent with Finite_Field.c
};

// Extension Field Fp18 (Fp6[w]/(w^3 - v)) where v is tau
struct Fp18{
    struct Fp6 x0, x1, x2; // Using x0, x1, x2 consistent with Finite_Field.c
};

// Elliptic Curve Point over Fp: y^2 = x^3 + b
struct EFp{
    struct Fp x, y; // Using x, y consistent with Elliptic_Curve.c & BN.c
    bool infity; // Using infity consistent with Elliptic_Curve.c
};

// Elliptic Curve Point over Fp3
struct EFp3{
    struct Fp3 x, y; // Using x, y consistent with Elliptic_Curve.c
    bool infity; // Using infity consistent with Elliptic_Curve.c
};

// Elliptic Curve Point over Fp6
struct EFp6{
    struct Fp6 x, y; // Using x, y consistent with Elliptic_Curve.c & BN.c
    bool infity; // Using infity consistent with Elliptic_Curve.c (PaI in BN.c) - Standardizing to infity
};

// Elliptic Curve Point over Fp18
struct EFp18{
    struct Fp18 x, y; // Using x, y consistent with embedding_degree18.c (p18x/y in main.c) - Standardizing to x, y
    bool infity; // Using infity consistent with embedding_degree18.c (isInfinity in main.c) - Standardizing to infity
};


// Global Parameters (Declare as extern, define in parameters.c)
extern mpz_t X;
extern mpz_t prime, r_order, t_trace, r_order_EFp, b;
extern mpz_t c1_leg, c1_leg_bar, c1_omega, c1_omega_bar; // Keep if used, otherwise remove

// Function Prototypes

// Parameter Generation
void generate_parameters();
void init_parameters();
void clear_parameters();

// Fp Arithmetic (from Finite_Field.c / main.c)
void Fp_init(struct Fp *A);
void Fp_set(struct Fp *ANS, struct Fp *A);
void Fp_set_ui(struct Fp *A, signed long int B);
void Fp_set_mpz(struct Fp *RES, mpz_t a); // Added from main.c
void Fp_random(struct Fp *A);
void Fp_clear(struct Fp *A);
void Fp_printf(struct Fp *A);
void Fp_add(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_add_ui(struct Fp *ANS, struct Fp *A, unsigned long int B);
void Fp_sub(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_sub_ui(struct Fp *ANS, struct Fp *A, unsigned long int B);
void Fp_mul(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_mul_ui(struct Fp *ANS, struct Fp *A, unsigned long int B);
void Fp_mul_c1(struct Fp *RES, struct Fp *a); // Added from main.c (specific to c1=2)
void Fp_inv(struct Fp *RES, struct Fp *a); // Added from main.c (uses mpz_invert)
void Fp_div(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_pow(struct Fp *ANS, struct Fp *A, mpz_t B);
void Fp_sqrt(struct Fp *ANS, struct Fp *A);
void Fp_neg(struct Fp *ANS, struct Fp *A);
int Fp_cmp(struct Fp *A, struct Fp *B);
int Fp_cmp_mpz(struct Fp *A, mpz_t B);

// Fp3 Arithmetic (from Finite_Field.c / main.c)
void Fp3_init(struct Fp3 *A);
void Fp3_set(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_set_ui(struct Fp3 *A, unsigned long int B); // Changed param type to match Finite_Field.c
void Fp3_random(struct Fp3 *A);
void Fp3_clear(struct Fp3 *A);
void Fp3_printf(struct Fp3 *A);
void Fp3_add(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_add_ui(struct Fp3 *ANS, struct Fp3 *A, unsigned long int B);
void Fp3_sub(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_mul(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_mul_Fp(struct Fp3 *ANS, struct Fp3 *A, struct Fp *B);
void Fp3_mul_ui(struct Fp3 *ANS, struct Fp3 *A, unsigned long int B);
void Fp3_mul_omega(struct Fp3 *RES, struct Fp3 *A); // Added from main.c (renamed from Fp3_mul_xi) - specific to Fp3 construction
void Fp3_invert(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_div(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_pow(struct Fp3 *ANS, struct Fp3 *A, mpz_t B);
void Fp3_sqrt(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_neg(struct Fp3 *ANS, struct Fp3 *A);
int Fp3_cmp(struct Fp3 *A, struct Fp3 *B);
int Fp3_cmp_mpz(struct Fp3 *A, mpz_t B);
int Fp3_legendre(struct Fp3 *a);
void Fp3_frobenius_map(struct Fp3 *ANS, struct Fp3 *A); // Needs review/consolidation (different implementations in main.c and Finite_Field.c)

// Fp6 Arithmetic (from Finite_Field.c / main.c)
void Fp6_init(struct Fp6 *A);
void Fp6_set(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_set_ui(struct Fp6 *A, signed long int B);
void Fp6_random(struct Fp6 *A);
void Fp6_clear(struct Fp6 *A);
void Fp6_printf(struct Fp6 *A);
void Fp6_add(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int B);
void Fp6_sub(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_mul(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_mul_Fp(struct Fp6 *ANS, struct Fp6 *A, struct Fp *B);
void Fp6_mul_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int B);
void Fp6_mul_tau(struct Fp6 *RES, struct Fp6 *A); // Added from main.c - specific to Fp6 construction (tau^2 = -1)
void Fp6_mul_omega(struct Fp6 *RES, struct Fp6 *A); // Added from main.c - multiplication by Fp3 generator
void Fp6_mul_v(struct Fp6 *ANS, struct Fp6 *A); // Added from Finite_Field.c (renamed from Fp6_mul_xi) - specific to Fp6 construction (tau^2 = xi)
void Fp6_invert(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_div(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_pow(struct Fp6 *ANS, struct Fp6 *A, mpz_t B);
void Fp6_sqrt(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_neg(struct Fp6 *ANS, struct Fp6 *A);
int Fp6_cmp(struct Fp6 *A, struct Fp6 *B);
int Fp6_cmp_mpz(struct Fp6 *A, mpz_t B);
int Fp6_legendre(struct Fp6 *A);
void Fp6_frobenius_map(struct Fp6 *ANS, struct Fp6 *A); // Added from Finite_Field.c

// Fp18 Arithmetic (from Finite_Field.c / main.c)
void Fp18_init(struct Fp18 *A);
void Fp18_set(struct Fp18 *ANS, struct Fp18 *A);
void Fp18_set_ui(struct Fp18 *A, signed long int B);
void Fp18_random(struct Fp18 *A);
void Fp18_clear(struct Fp18 *A);
void Fp18_printf(struct Fp18 *A);
void Fp18_add(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);
void Fp18_add_ui(struct Fp18 *ANS, struct Fp18 *A, unsigned long int B);
void Fp18_sub(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);
void Fp18_mul(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);
void Fp18_mul_Fp(struct Fp18 *ANS, struct Fp18 *A, struct Fp *B);
void Fp18_mul_ui(struct Fp18 *ANS, struct Fp18 *A, unsigned long int B);
void Fp18_invert(struct Fp18 *ANS, struct Fp18 *A);
void Fp18_div(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);
void Fp18_pow(struct Fp18 *ANS, struct Fp18 *A, mpz_t B);
void Fp18_sqrt(struct Fp18 *ANS, struct Fp18 *A);
void Fp18_neg(struct Fp18 *ANS, struct Fp18 *A);
int Fp18_cmp(struct Fp18 *A, struct Fp18 *B);
int Fp18_cmp_mpz(struct Fp18 *A, mpz_t B);
int Fp18_legendre(struct Fp18 *A);
void Fp18_frobenius_map(struct Fp18 *ANS, struct Fp18 *A); // Needs review/consolidation (different implementations in main.c and Finite_Field.c)

// Elliptic Curve Arithmetic (Prototypes from f18.h, Elliptic_Curve.c, BN.h, embedding_degree18.h, main.c)
// EFp
void EFp_init(struct EFp *A);
void EFp_set(struct EFp *A, struct EFp *B);
void EFp_set_infity(struct EFp *A); // Standardized name
void EFp_clear(struct EFp *A);
void EFp_printf(struct EFp *A);
void EFp_SCM(struct EFp *ANS, struct EFp *P, mpz_t scalar);
void EFp_ECD(struct EFp *ANS, struct EFp *P); // Point Doubling
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2); // Point Addition
int EFp_cmp(struct EFp *A, struct EFp *B);
void EFp_random_set(struct EFp *ANS); // Find random point on curve
void EFp_neg(struct EFp *RES, struct EFp *A); // Added, usually needed

// EFp3
void EFp3_init(struct EFp3 *A);
void EFp3_set(struct EFp3 *A, struct EFp3 *B);
void EFp3_set_infity(struct EFp3 *A);
void EFp3_clear(struct EFp3 *A);
void EFp3_printf(struct EFp3 *A);
void EFp3_SCM(struct EFp3 *ANS, struct EFp3 *P, mpz_t scalar);
void EFp3_ECD(struct EFp3 *ANS, struct EFp3 *P);
void EFp3_ECA(struct EFp3 *ANS, struct EFp3 *P1, struct EFp3 *P2);
int EFp3_cmp(struct EFp3 *A, struct EFp3 *B);
void EFp3_random_set(struct EFp3 *ANS); // Find random point on curve
void EFp3_neg(struct EFp3 *RES, struct EFp3 *A);
void EFp3_set_EFp(struct EFp3 *A, struct EFp *B); // Map point from base field

// EFp6
void EFp6_init(struct EFp6 *A);
void EFp6_set(struct EFp6 *A, struct EFp6 *B);
void EFp6_set_infity(struct EFp6 *A); // Standardized name
void EFp6_clear(struct EFp6 *A);
void EFp6_printf(struct EFp6 *A);
void EFp6_SCM(struct EFp6 *ANS, struct EFp6 *P, mpz_t scalar);
void EFp6_ECD(struct EFp6 *ANS, struct EFp6 *P);
void EFp6_ECA(struct EFp6 *ANS, struct EFp6 *P1, struct EFp6 *P2);
int EFp6_cmp(struct EFp6 *A, struct EFp6 *B);
void EFp6_random_set(struct EFp6 *ANS); // Find random point on curve
void EFp6_neg(struct EFp6 *RES, struct EFp6 *A); // Added, usually needed

// EFp18
void EFp18_init(struct EFp18 *A);
void EFp18_set(struct EFp18 *A, struct EFp18 *B);
void EFp18_set_infity(struct EFp18 *A); // Standardized name
void EFp18_clear(struct EFp18 *A);
void EFp18_printf(struct EFp18 *A);
void EFp18_SCM(struct EFp18 *ANS, struct EFp18 *P, mpz_t scalar);
void EFp18_ECD(struct EFp18 *ANS, struct EFp18 *P);
void EFp18_ECA(struct EFp18 *ANS, struct EFp18 *P1, struct EFp18 *P2);
int EFp18_cmp(struct EFp18 *A, struct EFp18 *B);
void EFp18_random_set(struct EFp18 *ANS); // Find random point on curve
void EFp18_random_set_G2(struct EFp18 *ANS); // Specific G2 point generation (from embedding_degree18.c)
void EFp18_neg(struct EFp18 *RES, struct EFp18 *A); // Added, usually needed

// Pairing Functions (Prototypes from embedding_degree18.c, BN.h)
// Need to decide which pairing (Tate, Ate) and which implementation to keep/refactor
void Miller_algo(struct Fp18 *ANS, struct EFp18 *P, struct EFp18 *Q, mpz_t loop); // From embedding_degree18.c (Fp18 result)
// void Miller_algo(struct Fp12 *ANS,struct EFp12 *P,struct EFp12 *Q,mpz_t roop); // From BN.c (Fp12 result) - Decide which degree/version
void Final_Exp(struct Fp18 *ANS, struct Fp18 *A); // From embedding_degree18.c (Fp18 result)
// void Final_Exp(struct Fp12 *ANS,struct Fp12 *A); // From BN.c (Fp12 result) - Decide which degree/version
void Tate_Pairing(struct Fp18 *ANS, struct EFp18 *P, struct EFp18 *Q); // From embedding_degree18.c (Fp18 result)
// void Tate_Pairing(struct Fp12 *ANS,struct EFp12 *A,struct EFp12 *B); // From BN.c (Fp12 result) - Decide which degree/version
void Optimal_Ate_Pairing(struct Fp18 *ANS, struct EFp *P, struct EFp18 *Q); // From embedding_degree18.c

// Input Functions (from main.c) - Consider moving to a test/example utility file
void Fp_take_input(struct Fp *RES);
void Fp3_take_input(struct Fp3 *RES);
void Fp6_take_input(struct Fp6 *RES);
void Fp18_take_input(struct Fp18 *RES);


#endif // FP18_ARITH_H
