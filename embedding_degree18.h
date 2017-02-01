#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#define TRUE 1
#define FALSE 0


struct Fp{
	mpz_t x0;
};
struct Fp3{
	struct Fp x0,x1,x2;
};
struct Fp6{
	struct Fp3 x0,x1;
};
struct Fp18{
	struct Fp6 x0,x1,x2;
};
struct EFp{
	struct Fp x,y;
	int infity;
};
struct EFp3{
	struct Fp3 x,y;
	int infity;
};
struct EFp6{
	struct Fp6 x,y;
	int infity;
};
struct EFp18{
	struct Fp18 x,y;
	int infity;
};
mpz_t X;

mpz_t prime,order,trace,order_EFp,b;


void test();


void Fp_init(struct Fp *A);
void Fp_set(struct Fp *ANS,struct Fp *A);
void Fp_set_ui(struct Fp *A,signed long int B);
void Fp_random(struct Fp *A);
void Fp_clear(struct Fp *A);
void Fp_printf(struct Fp *A);
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B);//ANS=A+B mod p
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int B);//ANS=A+B mod p
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B);//ANS=A-B mod p
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int B);//ANS=A+B mod p
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B);//ANS=A*B mod p
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int B);//ANS=A*B mod p
void Fp_div(struct Fp *ANS,struct Fp *A,struct Fp *B);//ANS=A/B mod p
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t B);//ANS=A^B mod p
void Fp_sqrt(struct Fp *ANS,struct Fp *A);//ANS^2=A
int Fp_cmp_mpz(struct Fp *A,mpz_t B);//A==B->0, A!=B->1
void Fp_neg(struct Fp *ANS,struct Fp *A);

//-----------------------------------------------------------------------------------------
void Fp3_init(struct Fp3 *A);
void Fp3_set(struct Fp3 *ANS,struct Fp3 *A);
void Fp3_set_ui(struct Fp3 *A,unsigned long int B);
void Fp3_random(struct Fp3 *A);
void Fp3_clear(struct Fp3 *A);
void Fp3_printf(struct Fp3 *A);
void Fp3_add(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B);
void Fp3_add_ui(struct Fp3 *ANS,struct Fp3 *A,unsigned long int B);
void Fp3_sub(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B);
void Fp3_mul_old(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B);
void Fp3_mul(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B);
void Fp3_square(struct Fp3 *ANS,struct Fp3 *A);
void Fp3_mul_xi(struct Fp3 *ANS,struct Fp3 *A);//ANS=A*(basis i)
void Fp3_mul_xi_inv(struct Fp3 *ANS,struct Fp3 *A);//ANS=A*(basis i)
void Fp3_mul_ui(struct Fp3 *ANS,struct Fp3 *A,unsigned long int B);
void Fp3_neg(struct Fp3 *ANS,struct Fp3 *A);//ANS=-A
void Fp3_invert(struct Fp3 *ANS,struct Fp3 *A);//ANS=A^-1
void Fp3_div(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B);
void Fp3_pow(struct Fp3 *ANS,struct Fp3 *A,mpz_t B);
void Fp3_sqrt(struct Fp3 *ANS,struct Fp3 *A);//ANS^2=A mod p
int Fp3_legendre(struct Fp3 *A);
int Fp3_cmp(struct Fp3 *A,struct Fp3 *B);
int Fp3_cmp_mpz(struct Fp3 *A,mpz_t B);
void Fp3_neg(struct Fp3 *ANS,struct Fp3 *A);
void Fp3_frobenius_map(struct Fp3 *ANS, struct Fp3 *A);

//-----------------------------------------------------------------------------------------
void Fp6_init(struct Fp6 *A);
void Fp6_set(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_set_ui(struct Fp6 *A,signed long int B);
void Fp6_random(struct Fp6 *A);
void Fp6_clear(struct Fp6 *A);
void Fp6_printf(struct Fp6 *A);
void Fp6_add(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B);
void Fp6_sub(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_mul(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_mul_v(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_mul_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B);
void Fp6_neg(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_invert(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_div(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_pow(struct Fp6 *ANS,struct Fp6 *A,mpz_t B);
void Fp6_sqrt(struct Fp6 *ANS,struct Fp6 *A);
int Fp6_legendre(struct Fp6 *A);
int Fp6_cmp(struct Fp6 *A,struct Fp6 *B);
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t B);
void Fp6_neg(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_map(struct Fp6 *ANS, struct Fp6 *A);

//-----------------------------------------------------------------------------------------
void Fp18_init(struct Fp18 *A);
void Fp18_set(struct Fp18 *ANS,struct Fp18 *A);
void Fp18_set_ui(struct Fp18 *A,signed long int B);
void Fp18_random(struct Fp18 *A);
void Fp18_clear(struct Fp18 *A);
void Fp18_printf(struct Fp18 *A);
void Fp18_add(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);
void Fp18_add_ui(struct Fp18 *ANS,struct Fp18 *A,unsigned long int B);
void Fp18_sub(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);
void Fp18_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);
void Fp18_mul_ui(struct Fp18 *ANS,struct Fp18 *A,unsigned long int B);
void Fp18_neg(struct Fp18 *ANS,struct Fp18 *A);
void Fp18_invert(struct Fp18 *ANS,struct Fp18 *A);
void Fp18_div(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);
void Fp18_pow(struct Fp18 *ANS,struct Fp18 *A,mpz_t B);
void Fp18_sqrt(struct Fp18 *ANS,struct Fp18 *A);//x^2=a mod p
int Fp18_legendre(struct Fp18 *A);
int Fp18_cmp(struct Fp18 *A,struct Fp18 *B);
int Fp18_cmp_mpz(struct Fp18 *A,mpz_t B);
void Fp18_neg(struct Fp18 *ANS,struct Fp18 *A);
void Fp18_frobenius_map(struct Fp18 *ANS, struct Fp18 *A,int times);

//-----------------------------------------------------------------------------------------

void EFp_init(struct EFp *A);
void EFp_set(struct EFp *A,struct EFp *B);
void EFp_set_infity(struct EFp *A);
void EFp_clear(struct EFp *A);
void EFp_printf(struct EFp *A);
void EFp_SCM(struct EFp *ANS, struct EFp *P,mpz_t j);//ANS=[j]P
void EFp_ECD(struct EFp *ANS, struct EFp *P);//ANS=2*P
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2);//ANS=P1+P2
int EFp_cmp(struct EFp *A,struct EFp *B);//A==B->0
void EFp_random_set(struct EFp *ANS);//generate random rational point

//-----------------------------------------------------------------------------------------
void EFp3_init(struct EFp3 *A);
void EFp3_set(struct EFp3 *A,struct EFp3 *B);
void EFp3_set_infity(struct EFp3 *A);
void EFp3_clear(struct EFp3 *A);
void EFp3_printf(struct EFp3 *A);
void EFp3_ECD(struct EFp3 *ANS, struct EFp3 *P);//ANS=2*P
void EFp3_ECA(struct EFp3 *ANS, struct EFp3 *P1, struct EFp3 *P2);//ANS=P1+P2
int EFp3_cmp(struct EFp3 *A,struct EFp3 *B);
void EFp3_set_EFp(struct EFp3 *A,struct EFp *B);//set EFp rational point for EFp3
int EFp3_set_EFp18_Sparse(struct EFp3 *ANS,struct EFp18 *A);//set EFp18 sparse point for EFp3
void EFp3_random_set(struct EFp3 *ANS);
void EFp3_SCM(struct EFp3 *ANS, struct EFp3 *P, mpz_t j);
void EFp3_neg(struct EFp3 *ANS, struct EFp3 *A);//ANS=-A

//-----------------------------------------------------------------------------------------
void EFp6_init(struct EFp6 *A);
void EFp6_set(struct EFp6 *A,struct EFp6 *B);
void EFp6_set_infity(struct EFp6 *A);
void EFp6_clear(struct EFp6 *A);
void EFp6_printf(struct EFp6 *A);
void EFp6_ECD(struct EFp6 *ANS, struct EFp6 *P);//ANS=2*P
void EFp6_ECA(struct EFp6 *ANS, struct EFp6 *P1, struct EFp6 *P2);//ANS=P1+P2
int EFp6_cmp(struct EFp6 *A,struct EFp6 *B);
void EFp6_random_set(struct EFp6 *ANS);
void EFp6_SCM(struct EFp6 *ANS, struct EFp6 *P, mpz_t j);

//-----------------------------------------------------------------------------------------
void EFp18_init(struct EFp18 *A);
void EFp18_set(struct EFp18 *A,struct EFp18 *B);
void EFp18_set_infity(struct EFp18 *A);
void EFp18_clear(struct EFp18 *A);
void EFp18_printf(struct EFp18 *A);
void EFp18_ECD(struct EFp18 *ANS, struct EFp18 *P);//ANS=2*P
void EFp18_ECA(struct EFp18 *ANS, struct EFp18 *P1, struct EFp18 *P2);//ANS=P1+P2
int EFp18_cmp(struct EFp18 *A,struct EFp18 *B);
void EFp18_set_EFp(struct EFp18 *A,struct EFp *B);
void EFp18_random_set(struct EFp18 *ANS);
void EFp18_random_set_G2(struct EFp18 *ANS);
void EFp18_SCM(struct EFp18 *ANS, struct EFp18 *P, mpz_t j);
void EFp18_frobenius_map(struct EFp18 *ANS,struct EFp18 *A,int i);

//-----------------------------------------------------------------------------------------

void EFp_set_EC_parameter(void);
//-----------------------------------------------------------------------------------------

void Miller_algo(struct Fp18 *ANS,struct EFp18 *P,struct EFp18 *Q,mpz_t loop);
void Optimal_Miller(struct Fp18 *ANS,struct EFp18 *P,struct EFp18 *Q,mpz_t loop);
void Tate_Pairing(struct Fp18 *ANS,struct EFp18 *P,struct EFp18 *Q);
void Ate_Pairing(struct Fp18 *ANS,struct EFp18 *P,struct EFp18 *Q);
void Optimal_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2);
void ADD_LINE(struct Fp18 *l_ANS,struct EFp18 *T_ANS,struct EFp18 *T,struct EFp18 *P,struct EFp18 *Q,struct Fp18 *Px_neg);
void DBL_LINE(struct Fp18 *l_ANS,struct EFp18 *T_ANS,struct EFp18 *T,struct EFp18 *Q,struct Fp18 *Px_neg);

void Sparse_Ate_Pairing(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q);
void Sparse_Optimal_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2);

void Sparse_type1_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Sparse_type1_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Sparse_type1_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *Px_neg);
void Sparse_type1_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *Px_neg);
void Sparse_type1_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);

void Sparse_type2_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Sparse_type2_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Sparse_type2_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *Px_neg);
void Sparse_type2_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *Px_neg);
void Sparse_type2_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);

void Pseudo_Sparse_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2);
void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2);
void Pseudo_type1_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Pseudo_type1_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Pseudo_type1_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *L);
void Pseudo_type1_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *L);
void Pseudo_type1_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);

void Pseudo_type2_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Pseudo_type2_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop);
void Pseudo_type2_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *L);
void Pseudo_type2_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *L);
void Pseudo_type2_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B);

void ltt_q(struct Fp18 *ANS,struct EFp18 *Q,struct EFp18 *R);
void v2t_q(struct Fp18 *ANS,struct EFp18 *Q,struct EFp18 *R);
void ltp_q(struct Fp18 *ANS,struct EFp18 *Q,struct EFp18 *R,struct EFp18 *P);
void vtp_q(struct Fp18 *ANS,struct EFp18 *Q,struct EFp18 *R);
void Final_Exp(struct Fp18 *ANS,struct Fp18 *A);


void check_Pairing(void);
void Masure_pairing_time(void);
