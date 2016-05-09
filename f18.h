//
//  f18.h
//  Fp18_Arith
//
//  Created by Khandaker Md. Al-Amin on 4/29/16.
//  Copyright © 2016 Khandaker Md. Al-Amin. All rights reserved.
//

#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define TRUE 1
#define FALSE 0

//======Fp18 Arith Variables========
struct Fp{
    mpz_t x_0;
};
struct Fp3{
    struct Fp a0,a1,a2;
};
struct Fp6{
    struct Fp3 a0,a1;
};
struct Fp18{
    struct Fp6 m0,m1,m2;
};

//======ECC Variables========
struct EFp{
    struct Fp px,py;
    int isInfinity;
};
struct EFp3{
    struct Fp3 p3x,p3y;
    int isInfinity;
};
struct EFp6{
    struct Fp6 p6x,p6y;
    int isInfinity;
};
struct EFp18{
    struct Fp18 p18x,p18y;
    int isInfinity;
};

#pragma mark  Parameters Variables
mpz_t X; // variable to find p,r
mpz_t prime,r_order,t_trace,r_order_EFp,b;
mpz_t c1_leg,c1_leg_bar,c1_omega,c1_omega_bar;

#pragma mark Fp methods
void Fp_take_input(struct Fp *RES);
void Fp_init(struct Fp *A);
void Fp_set(struct Fp *RES,struct Fp *a);
void Fp_set_mpz(struct Fp *RES,mpz_t a);
void Fp_random(struct Fp *a);
void Fp_clear(struct Fp *A);
void Fp_printf(struct Fp *A);

void Fp_add(struct Fp *RES,struct Fp *a,struct Fp *b);
void Fp_add_ui(struct Fp *RES,struct Fp *a,unsigned long int b);
void Fp_sub(struct Fp *RES,struct Fp *a,struct Fp *b);
void Fp_sub_ui(struct Fp *RES,struct Fp *a,unsigned long int b);
void Fp_mul(struct Fp *RES,struct Fp *a,struct Fp *b);
void Fp_mul_c1(struct Fp *RES,struct Fp *a);
void Fp_inv(struct Fp *RES,struct Fp *a);
void Fp_div(struct Fp *RES,struct Fp *a,struct Fp *b);
void Fp_neg(struct Fp *RES,struct Fp *a);

void Fp_pow(struct Fp *RES,struct Fp *a,mpz_t b);
void Fp_sqrt(struct Fp *RES,struct Fp *a);
int  Fp_cmp_mpz(struct Fp *a,mpz_t b);


#pragma mark Fp3 methods
void Fp3_take_input(struct Fp3 *RES);
void Fp3_init(struct Fp3 *A);
void Fp3_set(struct Fp3 *RES,struct Fp3 *A);
void Fp3_set_ui(struct Fp3 *A,signed long int B);
void Fp3_random(struct Fp3 *A);
void Fp3_clear(struct Fp3 *A);
void Fp3_printf(struct Fp3 *A);

void Fp3_add(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B);
void Fp3_add_ui(struct Fp3 *RES,struct Fp3 *A,unsigned long int B);
void Fp3_sub(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B);
void Fp3_mul(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B);
void Fp3_mul_omega(struct Fp3 *RES,struct Fp3 *A);
void Fp3_mul_ui(struct Fp3 *RES,struct Fp3 *A,unsigned long int B);
void Fp3_mul_Fp(struct Fp3 *RES,struct Fp3 *A,struct Fp *B);
void Fp3_neg(struct Fp3 *RES,struct Fp3 *A);
void Fp3_frobenius_map(struct Fp3 *RES,struct Fp3 *A);
void Fp3_invert(struct Fp3 *RES,struct Fp3 *A);
void Fp3_div(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B);

void Fp3_pow(struct Fp3 *RES,struct Fp3 *A,mpz_t B);
void Fp3_sqrt(struct Fp3 *RES,struct Fp3 *A);//x^2=a mod p
int  Fp3_legendre(struct Fp3 *A);
int  Fp3_cmp(struct Fp3 *A,struct Fp3 *B);
int  Fp3_cmp_mpz(struct Fp3 *A,mpz_t B);

#pragma mark Fp6 methods
void Fp6_take_input(struct Fp6 *RES);
void Fp6_init(struct Fp6 *A);
void Fp6_set(struct Fp6 *RES,struct Fp6 *A);
void Fp6_set_ui(struct Fp6 *A,signed long int B);
void Fp6_random(struct Fp6 *A);
void Fp6_clear(struct Fp6 *A);
void Fp6_printf(struct Fp6 *A);

void Fp6_add(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *RES,struct Fp6 *A,unsigned long int B);
void Fp6_sub(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B);
void Fp6_neg(struct Fp6 *RES,struct Fp6 *A);
void Fp6_mul(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B);
void Fp6_mul_tau(struct Fp6 *RES,struct Fp6 *A);
void Fp6_mul_omega(struct Fp6 *RES,struct Fp6 *A);
void Fp6_mul_ui(struct Fp6 *RES,struct Fp6 *A,unsigned long int B);
void Fp6_mul_Fp(struct Fp6 *RES,struct Fp6 *A,struct Fp *B);
void Fp6_invert(struct Fp6 *RES,struct Fp6 *A);
void Fp6_div(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B);

void Fp6_pow(struct Fp6 *RES,struct Fp6 *A,mpz_t B);
void Fp6_sqrt(struct Fp6 *RES,struct Fp6 *A);//x^2=a mod p
int  Fp6_legendre(struct Fp6 *A);
int  Fp6_cmp(struct Fp6 *A,struct Fp6 *B);
int  Fp6_cmp_mpz(struct Fp6 *A,mpz_t B);


#pragma mark Fp18 methods
void Fp18_take_input(struct Fp18 *RES);
void Fp18_init(struct Fp18 *A);
void Fp18_set(struct Fp18 *RES,struct Fp18 *A);
void Fp18_set_ui(struct Fp18 *A,signed long int B);
void Fp18_random(struct Fp18 *A);
void Fp18_clear(struct Fp18 *A);
void Fp18_printf(struct Fp18 *A);

void Fp18_add(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B);
void Fp18_add_ui(struct Fp18 *RES,struct Fp18 *A,unsigned long int B);
void Fp18_sub(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B);
void Fp18_mul(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B);
void Fp18_mul_ui(struct Fp18 *RES,struct Fp18 *A,unsigned long int B);
void Fp18_mul_Fp(struct Fp18 *RES,struct Fp18 *A,struct Fp *B);//TODO: not required
void Fp18_neg(struct Fp18 *RES,struct Fp18 *A);
void Fp18_frobenius_map(struct Fp18 *RES,struct Fp18 *A);
void Fp18_invert(struct Fp18 *RES,struct Fp18 *A);
void Fp18_div(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B);

void Fp18_pow(struct Fp18 *RES,struct Fp18 *A,mpz_t B);
void Fp18_sqrt(struct Fp18 *RES,struct Fp18 *A);//x^2=a mod p
int  Fp18_legendre(struct Fp18 *A);
int  Fp18_cmp(struct Fp18 *A,struct Fp18 *B);
int  Fp18_cmp_mpz(struct Fp18 *A,mpz_t B);


#pragma mark EFp methods
void EFp_init(struct EFp *A);
void EFp_set(struct EFp *A,struct EFp *B);
void EFp_set_infity(struct EFp *A);
void EFp_clear(struct EFp *A);
void EFp_printf(struct EFp *A);
void EFp_SCM(struct EFp *RES, struct EFp *P,mpz_t scalar);
void EFp_ECD(struct EFp *RES, struct EFp *P);//RES=2*P
void EFp_ECA(struct EFp *RES, struct EFp *P, struct EFp *Q);//RES=P1+P2
int  EFp_cmp(struct EFp *A,struct EFp *B);


#pragma mark EFp3 methods
void EFp3_init(struct EFp3 *A);
void EFp3_set(struct EFp3 *A,struct EFp3 *B);
void EFp3_set_infity(struct EFp3 *A);
void EFp3_clear(struct EFp3 *A);
void EFp3_printf(struct EFp3 *A);
void EFp3_ECD(struct EFp3 *RES, struct EFp3 *P);
void EFp3_ECA(struct EFp3 *RES, struct EFp3 *P, struct EFp3 *Q);
int  EFp3_cmp(struct EFp3 *A,struct EFp3 *B);
void EFp3_set_EFp(struct EFp3 *A,struct EFp *B);
void EFp3_SCM(struct EFp3 *RES, struct EFp3 *P, mpz_t scalar);
void EFp3_neg(struct EFp3 *RES, struct EFp3 *A);


#pragma mark EFp6 methods
void EFp6_init(struct EFp6 *A);
void EFp6_set(struct EFp6 *A,struct EFp6 *B);
void EFp6_set_infity(struct EFp6 *A);
void EFp6_clear(struct EFp6 *A);
void EFp6_printf(struct EFp6 *A);
void EFp6_ECD(struct EFp6 *RES, struct EFp6 *P);
void EFp6_ECA(struct EFp6 *RES, struct EFp6 *P1, struct EFp6 *P2);
void EFp6_SCM(struct EFp6 *RES, struct EFp6 *P, mpz_t S);
int  EFp6_cmp(struct EFp6 *A,struct EFp6 *B);

#pragma mark EFp18 methods
void EFp18_init(struct EFp18 *A);
void EFp18_set(struct EFp18 *A,struct EFp18 *B);
void EFp18_set_infity(struct EFp18 *A);
void EFp18_clear(struct EFp18 *A);
void EFp18_printf(struct EFp18 *A);
void EFp18_ECD(struct EFp18 *RES, struct EFp18 *P);
void EFp18_ECA(struct EFp18 *RES, struct EFp18 *P1, struct EFp18 *P2);
void EFp18_SCM(struct EFp18 *RES,struct EFp18 *P,mpz_t j);
int  EFp18_cmp(struct EFp18 *A,struct EFp18 *B);

#pragma mark Parameters methods
void generate_parameters(void);
void get_C1_C1bar();
void get_C1omega_C1omegabar();

#pragma mark Previous methods
void Fp_pow_bin(struct Fp *RES,struct Fp *A,mpz_t B);
void Fp_mul_ui(struct Fp *RES,struct Fp *a,unsigned long int b);
void Fp3_mul_prev(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B);
void Fp3_invert_prev(struct Fp3 *ANS, struct Fp3 *A);
void Fp6_pow_prev(struct Fp6 *RES,struct Fp6 *A,mpz_t B);
void Fp6_mul_prev(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B);
void Fp6_invert_prev(struct Fp6 *ANS, struct Fp6 *A);
void Fp18_invert_prev(struct Fp18 *RES, struct Fp18 *A);
