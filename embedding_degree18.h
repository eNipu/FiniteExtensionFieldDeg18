#ifndef EMBEDDING_DEGREE18_H
#define EMBEDDING_DEGREE18_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

// Define TRUE and FALSE constants if not already defined
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

// Structure definitions for finite field extensions
struct Fp {
    mpz_t x_0;  // Using the naming convention from f18.h
};

struct Fp3 {
    struct Fp a0, a1, a2;  // Using the naming convention from f18.h
};

struct Fp6 {
    struct Fp3 a0, a1;  // Using the naming convention from f18.h
};

struct Fp18 {
    struct Fp6 m0, m1, m2;  // Using the naming convention from f18.h
};

// Elliptic curve structures
struct EFp {
    struct Fp px, py;
    int isInfinity;  // Using the naming convention from f18.h
};

struct EFp3 {
    struct Fp3 p3x, p3y;  // Using the naming convention from f18.h
    int isInfinity;
};

struct EFp6 {
    struct Fp6 p6x, p6y;  // Using the naming convention from f18.h
    int isInfinity;
};

struct EFp18 {
    struct Fp18 p18x, p18y;  // Using the naming convention from f18.h
    int isInfinity;
};

// Global variables for parameters
extern mpz_t X;  // variable to find p,r
extern mpz_t prime, r_order, t_trace, r_order_EFp, b;
extern mpz_t c1_leg, c1_leg_bar, c1_omega, c1_omega_bar;
extern int *X_bit_binary;
extern int X_bit;

// Function declarations for Fp arithmetic
void Fp_init(struct Fp *A);
void Fp_set(struct Fp *A, struct Fp *B);
void Fp_set_ui(struct Fp *A, unsigned long int B);
void Fp_set_mpz(struct Fp *A, mpz_t B);
void Fp_random(struct Fp *A);
void Fp_clear(struct Fp *A);
void Fp_printf(struct Fp *A);
void Fp_add(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_add_ui(struct Fp *ANS, struct Fp *A, unsigned long int B);
void Fp_sub(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_sub_ui(struct Fp *ANS, struct Fp *A, unsigned long int B);
void Fp_mul(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_mul_ui(struct Fp *ANS, struct Fp *A, unsigned long int B);
void Fp_div(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_pow(struct Fp *ANS, struct Fp *A, mpz_t B);
void Fp_sqrt(struct Fp *ANS, struct Fp *A);
void Fp_neg(struct Fp *ANS, struct Fp *A);
void Fp_invert(struct Fp *ANS, struct Fp *A);
int Fp_cmp(struct Fp *A, struct Fp *B);
int Fp_cmp_mpz(struct Fp *A, mpz_t B);

// Function declarations for Fp3 arithmetic
void Fp3_init(struct Fp3 *A);
void Fp3_set(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_set_ui(struct Fp3 *A, unsigned long int B);
void Fp3_random(struct Fp3 *A);
void Fp3_clear(struct Fp3 *A);
void Fp3_printf(struct Fp3 *A);
void Fp3_add(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_add_ui(struct Fp3 *ANS, struct Fp3 *A, unsigned long int B);
void Fp3_sub(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_mul(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_mul_omega(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_mul_Fp(struct Fp3 *ANS, struct Fp3 *A, struct Fp *B);
void Fp3_mul_ui(struct Fp3 *ANS, struct Fp3 *A, unsigned long int B);
void Fp3_invert(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_div(struct Fp3 *ANS, struct Fp3 *A, struct Fp3 *B);
void Fp3_pow(struct Fp3 *ANS, struct Fp3 *A, mpz_t B);
int Fp3_cmp(struct Fp3 *A, struct Fp3 *B);
int Fp3_cmp_mpz(struct Fp3 *A, mpz_t B);
int Fp3_legendre(struct Fp3 *A);
void Fp3_sqrt(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_neg(struct Fp3 *ANS, struct Fp3 *A);
void Fp3_frobenius_map(struct Fp3 *ANS, struct Fp3 *A);

// Function declarations for Fp6 arithmetic
void Fp6_init(struct Fp6 *A);
void Fp6_set(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_set_ui(struct Fp6 *A, unsigned long int B);
void Fp6_random(struct Fp6 *A);
void Fp6_clear(struct Fp6 *A);
void Fp6_printf(struct Fp6 *A);
void Fp6_add(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int B);
void Fp6_sub(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_mul(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_mul_tau(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_mul_Fp(struct Fp6 *ANS, struct Fp6 *A, struct Fp *B);
void Fp6_mul_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int B);
void Fp6_invert(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_div(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_pow(struct Fp6 *ANS, struct Fp6 *A, mpz_t B);
int Fp6_cmp(struct Fp6 *A, struct Fp6 *B);
int Fp6_cmp_mpz(struct Fp6 *A, mpz_t B);
int Fp6_legendre(struct Fp6 *A);
void Fp6_sqrt(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_neg(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_map(struct Fp6 *ANS, struct Fp6 *A);

// Function declarations for Fp18 arithmetic
void Fp18_init(struct Fp18 *A);
void Fp18_set(struct Fp18 *ANS, struct Fp18 *A);
void Fp18_set_ui(struct Fp18 *A, unsigned long int B);
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
int Fp18_cmp(struct Fp18 *A, struct Fp18 *B);
int Fp18_cmp_mpz(struct Fp18 *A, mpz_t B);
void Fp18_neg(struct Fp18 *ANS, struct Fp18 *A);
void Fp18_frobenius_map(struct Fp18 *ANS, struct Fp18 *A, int i);

// Function declarations for elliptic curve operations
void EFp_init(struct EFp *P);
void EFp_set(struct EFp *P, struct EFp *Q);
void EFp_set_infity(struct EFp *P);
void EFp_clear(struct EFp *P);
void EFp_printf(struct EFp *P);
void EFp_SCM(struct EFp *P, struct EFp *Q, mpz_t scalar);
void EFp_ECD(struct EFp *P, struct EFp *Q);
void EFp_ECA(struct EFp *ANS, struct EFp *P, struct EFp *Q);

// Elliptic curve over Fp3
void EFp3_init(struct EFp3 *P);
void EFp3_set(struct EFp3 *P, struct EFp3 *Q);
void EFp3_set_infity(struct EFp3 *P);
void EFp3_clear(struct EFp3 *P);
void EFp3_printf(struct EFp3 *P);
void EFp3_ECD(struct EFp3 *ANS, struct EFp3 *P);
void EFp3_ECA(struct EFp3 *ANS, struct EFp3 *P, struct EFp3 *Q);
void EFp3_SCM(struct EFp3 *ANS, struct EFp3 *P, mpz_t scalar);
void EFp3_type2_SCM(struct EFp3 *ANS, struct EFp3 *P, mpz_t scalar);
void EFp3_set_EFp(struct EFp3 *ANS, struct EFp *A);
int EFp3_set_EFp18_Sparse(struct EFp3 *ANS, struct EFp18 *A);
void EFp3_neg(struct EFp3 *ANS, struct EFp3 *P);

// Elliptic curve over Fp18
void EFp18_init(struct EFp18 *P);
void EFp18_set(struct EFp18 *P, struct EFp18 *Q);
void EFp18_set_infity(struct EFp18 *P);
void EFp18_clear(struct EFp18 *P);
void EFp18_printf(struct EFp18 *P);
void EFp18_ECD(struct EFp18 *ANS, struct EFp18 *P);
void EFp18_ECA(struct EFp18 *ANS, struct EFp18 *P, struct EFp18 *Q);

// KSS Pairing functions
void Sparse_type1_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Sparse_type1_Optimal_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Sparse_type1_ADD_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *P, struct EFp3 *Q, struct Fp3 *Qx_neg);
void Sparse_type1_DBL_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *Q, struct Fp3 *Qx_neg);
void Sparse_type1_mul(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);

void Sparse_type2_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Sparse_type2_Optimal_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Sparse_type2_ADD_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *P, struct EFp3 *Q, struct Fp3 *Qx_neg);
void Sparse_type2_DBL_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *Q, struct Fp3 *Qx_neg);
void Sparse_type2_mul(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);

void Pseudo_type1_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Pseudo_type1_Optimal_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Pseudo_type1_ADD_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *P, struct EFp3 *Q, struct Fp3 *L);
void Pseudo_type1_DBL_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *Q, struct Fp3 *L);
void Pseudo_type1_mul(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);

void Pseudo_type2_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Pseudo_type2_Optimal_Miller(struct Fp18 *ANS, struct EFp3 *P, struct EFp3 *Q, mpz_t loop);
void Pseudo_type2_ADD_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *P, struct EFp3 *Q, struct Fp3 *L);
void Pseudo_type2_DBL_LINE(struct Fp18 *l_ANS, struct EFp3 *T_ANS, struct EFp3 *T, struct EFp3 *Q, struct Fp3 *L);
void Pseudo_type2_mul(struct Fp18 *ANS, struct Fp18 *A, struct Fp18 *B);

void Pseudo_Sparse_Ate_Pairing(struct Fp18 *ANS, struct EFp *G1, struct EFp18 *G2);
void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp18 *ANS, struct EFp *G1, struct EFp18 *G2);

void Final_Exp(struct Fp18 *ANS, struct Fp18 *A);

// Parameter generation and utility functions
void generate_parameters(void);
void check_Pairing(void);
void Masure_pairing_time(void);

#endif // EMBEDDING_DEGREE18_H
