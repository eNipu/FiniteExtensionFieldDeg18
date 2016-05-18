//gcc -I/usr/local/include -L/usr/local/lib -Wall -O3 -o BN2.out BN2.c -lgmp


#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define TRUE 1
#define FALSE 0
struct Fp{
	mpz_t x0;
};
struct Fp2{
	struct Fp x0,x1;
};
struct Fp6{
	struct Fp2 x0,x1,x2;
};
struct Fp12{
	struct Fp6 x0,x1;
};
struct EFp{
	struct Fp x,y;
	int PaI;//Point at Infinity
};
struct EFp2{
	struct Fp2 x,y;
	int PaI;
};
struct EFp6{
	struct Fp6 x,y;
	int PaI;
};
struct EFp12{
	struct Fp12 x,y;
	int PaI;
};
mpz_t X;

mpz_t p,t,r;//r=p-t+1
mpz_t b;
void Fp_init(struct Fp *A);
void Fp_set(struct Fp *ANS,struct Fp *A);
void Fp_set_si(struct Fp *A,signed long int B);
void Fp_random(struct Fp *A);
void Fp_clear(struct Fp *A);
void Fp_printf(struct Fp *A);
void Fp_add(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a+b mod p
void Fp_add_ui(struct Fp *ans,struct Fp *a,unsigned long int b);//ans=a+b mod p
void Fp_sub(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a-b mod p
void Fp_sub_ui(struct Fp *ans,struct Fp *a,unsigned long int b);//ans=a+b mod p
void Fp_mul(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a*b mod p
void Fp_mul_ui(struct Fp *ans,struct Fp *a,unsigned long int b);//ans=a*b mod p
void Fp_div(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a/b mod p
void Fp_pow(struct Fp *ans,struct Fp *a,mpz_t b);
void Fp_sqrt(struct Fp *ans,struct Fp *a);//x^2=a mod p
int Fp_cmp_mpz(struct Fp *A,mpz_t B);
//-----------------------------------------------------------------------------------------

void Fp2_init(struct Fp2 *A);
void Fp2_set(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_set_si(struct Fp2 *A,signed long int B);
void Fp2_random(struct Fp2 *A);
void Fp2_clear(struct Fp2 *A);
void Fp2_printf(struct Fp2 *A);
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B);
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_mul_i(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B);
void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B);
void Fp2_invert(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A);//x^2=a mod p
int Fp2_legendre(struct Fp2 *A);
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B);
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t B);
void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);
//-----------------------------------------------------------------------------------------

void Fp6_init(struct Fp6 *A);
void Fp6_set(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_set_si(struct Fp6 *A,signed long int B);
void Fp6_random(struct Fp6 *A);
void Fp6_clear(struct Fp6 *A);
void Fp6_printf(struct Fp6 *A);
void Fp6_add(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B);
void Fp6_sub(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_mul(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_mul_v(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_mul_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B);
void Fp6_mul_Fp(struct Fp6 *ANS,struct Fp6 *A,struct Fp *B);
void Fp6_neg(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_invert(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_div(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_pow(struct Fp6 *ANS,struct Fp6 *A,mpz_t B);
void Fp6_sqrt(struct Fp6 *ANS,struct Fp6 *A);//x^2=a mod p
int Fp6_legendre(struct Fp6 *A);
int Fp6_cmp(struct Fp6 *A,struct Fp6 *B);
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t B);
void Fp6_neg(struct Fp6 *ANS,struct Fp6 *A);
//-----------------------------------------------------------------------------------------

void Fp12_init(struct Fp12 *A);
void Fp12_set(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_set_si(struct Fp12 *A,signed long int B);
void Fp12_random(struct Fp12 *A);
void Fp12_clear(struct Fp12 *A);
void Fp12_printf(struct Fp12 *A);
void Fp12_add(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_add_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int B);
void Fp12_sub(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_mul_Fp(struct Fp12 *ANS,struct Fp12 *A,struct Fp *B);
void Fp12_invert(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_div(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_pow(struct Fp12 *ANS,struct Fp12 *A,mpz_t B);
void Fp12_sqrt(struct Fp12 *ANS,struct Fp12 *A);//x^2=a mod p
int Fp12_legendre(struct Fp12 *A);
int Fp12_cmp(struct Fp12 *A,struct Fp12 *B);
int Fp12_cmp_mpz(struct Fp12 *A,mpz_t B);
void Fp12_neg(struct Fp12 *ANS,struct Fp12 *A);
//-----------------------------------------------------------------------------------------

void EFp_init(struct EFp *A);
void EFp_set(struct EFp *A,struct EFp *B);
void EFp_set_PaI(struct EFp *A);
void EFp_clear(struct EFp *A);
void EFp_printf(struct EFp *A);
void EFp_SCM(struct EFp *ANS, struct EFp *P,mpz_t j);
void EFp_ECD(struct EFp *ANS, struct EFp *P);//ANS=2*P
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2);//ANS=P1+P2
int EFp_cmp(struct EFp *A,struct EFp *B);
void EFp_random_set(struct EFp *ANS);//random set EFp on curve
void EFp_to_EFp12(struct EFp12 *A,struct EFp *B);
//-----------------------------------------------------------------------------------------

void EFp2_init(struct EFp2 *A);
void EFp2_set(struct EFp2 *A,struct EFp2 *B);
void EFp2_set_PaI(struct EFp2 *A);
void EFp2_clear(struct EFp2 *A);
void EFp2_printf(struct EFp2 *A);
void EFp2_SCM(struct EFp2 *ANS,struct EFp2 *P,mpz_t j);
void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P);//ANS=2*P
void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2);//ANS=P1+P2
int EFp2_cmp(struct EFp2 *A,struct EFp2 *B);
void EFp2_random_set(struct EFp2 *ANS);

//-----------------------------------------------------------------------------------------

void EFp6_init(struct EFp6 *A);
void EFp6_set(struct EFp6 *A,struct EFp6 *B);
void EFp6_set_PaI(struct EFp6 *A);
void EFp6_clear(struct EFp6 *A);
void EFp6_printf(struct EFp6 *A);
void EFp6_ECD(struct EFp6 *ANS, struct EFp6 *P);//ANS=2*P
void EFp6_ECA(struct EFp6 *ANS, struct EFp6 *P1, struct EFp6 *P2);//ANS=P1+P2
void EFp6_SCM(struct EFp6 *ANS,struct EFp6 *P,mpz_t j);
int EFp6_cmp(struct EFp6 *A,struct EFp6 *B);
void EFp6_random_set(struct EFp6 *ANS);
//-----------------------------------------------------------------------------------------

void EFp12_init(struct EFp12 *A);
void EFp12_set(struct EFp12 *A,struct EFp12 *B);
void EFp12_set_PaI(struct EFp12 *A);
void EFp12_clear(struct EFp12 *A);
void EFp12_printf(struct EFp12 *A);
void EFp12_ECD(struct EFp12 *ANS, struct EFp12 *P);//ANS=2*P
void EFp12_ECA(struct EFp12 *ANS, struct EFp12 *P1, struct EFp12 *P2);//ANS=P1+P2
int EFp12_cmp(struct EFp12 *A,struct EFp12 *B);
void EFp12_SCM(struct EFp12 *ANS, struct EFp12 *P, mpz_t j);
void EFp12_random_set(struct EFp12 *ANS);
void EFp12_random_set_for_Ate(struct EFp12 *ANS);
void EFp12_frobenius_map(struct EFp12 *ANS,struct EFp12 *A);
//-----------------------------------------------------------------------------------------

void set_BN_parameter(struct EFp *ANS);

//-----------------------------------------------------------------------------------------
void Miller_algo(struct Fp12 *ANS,struct EFp12 *P,struct EFp12 *Q,mpz_t roop);
void Final_Exp(struct Fp12 *ANS,struct Fp12 *A);
void Tate_Pairing(struct Fp12 *ANS,struct EFp12 *A,struct EFp12 *B);
void Ate_Pairing(struct Fp12 *ANS,struct EFp12 *A,struct EFp12 *B);
void Optimal_Ate_Pairing(struct Fp12 *ANS,struct EFp12 *A,struct EFp12 *B);
void check_Pairing(void);
void ltt_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *Q);
void v2t_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *Q);
void ltp_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *P,struct EFp12 *Q);
void vtp_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *Q);
