#include"embedding_degree18.h"
#include"Finite_Field.c"

#ifndef _Elliptic_Curve_C_
#define _Elliptic_Curve_C_

void EFp_init(struct EFp *A){
	Fp_init(&A->x);
	Fp_init(&A->y);
	A->infity=FALSE;
}
void EFp_set(struct EFp *A,struct EFp *B){
	Fp_set(&A->x,&B->x);
	Fp_set(&A->y,&B->y);
	A->infity=B->infity;
}
void EFp_set_infity(struct EFp *A){
	Fp_set_ui(&A->x,0);
	Fp_set_ui(&A->y,0);
	A->infity=TRUE;
}
void EFp_clear(struct EFp *A){
	Fp_clear(&A->x);
	Fp_clear(&A->y);
}
void EFp_printf(struct EFp *A){
	gmp_printf("(%Zd,%Zd)\n",A->x.x0,A->y.x0);
}
void EFp_SCM(struct EFp *ANS, struct EFp *P,mpz_t j){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(j,2);

	struct EFp Q;
	EFp_init(&Q);
	EFp_set(&Q,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp_ECD(&Q,&Q);
			EFp_ECA(&Q,&Q,P);
		}else{
			EFp_ECD(&Q,&Q);
		}
	}

	EFp_set(ANS,&Q);
	EFp_clear(&Q);
	return;
}
void EFp_ECD(struct EFp *ANS, struct EFp *P){
	if(P->infity==TRUE){
		EFp_set(ANS,P);
		return;
	}
	if(mpz_sgn(P->y.x0)==0){//P.y==0
		EFp_set_infity(ANS);
		return;
	}

	struct Fp x,y,lambda,tmp;
	struct EFp t_ans;
	Fp_init(&x);
	Fp_init(&lambda);
	Fp_init(&tmp);
	Fp_init(&y);
	EFp_init(&t_ans);

	Fp_mul(&x,&P->x,&P->x);
	Fp_add(&tmp,&x,&x);
	Fp_add(&x,&tmp,&x);//3x^2
	// gmp_printf("%Zd\n",tmp_a);
	Fp_add(&y,&P->y,&P->y);//2y

	Fp_div(&lambda,&x,&y);
	Fp_mul(&tmp,&lambda,&lambda);
	Fp_add(&x,&P->x,&P->x);
	Fp_sub(&x,&tmp,&x);
	Fp_sub(&tmp,&P->x,&x);

	Fp_set(&t_ans.x,&x);

	Fp_mul(&tmp,&tmp,&lambda);
	Fp_sub(&t_ans.y,&tmp,&P->y);

	EFp_set(ANS,&t_ans);

	Fp_clear(&x);
	Fp_clear(&lambda);
	Fp_clear(&y);
	Fp_clear(&tmp);
	EFp_clear(&t_ans);
}
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2){
	if(P2->infity==TRUE){//if P2==inf
		EFp_set(ANS,P1);
		return;
	}
	else if(P1->infity==TRUE){//if P1==inf
		EFp_set(ANS,P2);
		return;
	}
	else if(Fp_cmp(&P1->x,&P2->x)==0&&Fp_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp_set_infity(ANS);
		return;
	}
	else if(EFp_cmp(P1,P2)==0){ // P=Q
		EFp_ECD(ANS,P1);
		return;
	}

	struct Fp x,y,lambda,tmp;
	struct EFp t_ans;

	Fp_init(&x);
	Fp_init(&y);
	Fp_init(&lambda);
	Fp_init(&tmp);
	EFp_init(&t_ans);

	Fp_sub(&x,&P2->x,&P1->x);
	Fp_sub(&y,&P2->y,&P1->y);
	Fp_div(&lambda,&y,&x);
	Fp_mul(&tmp,&lambda,&lambda);
	Fp_add(&x,&P1->x,&P2->x);
	Fp_sub(&x,&tmp,&x);
	Fp_sub(&tmp,&P1->x,&x);
	Fp_set(&t_ans.x,&x);
	Fp_mul(&tmp,&tmp,&lambda);
	Fp_sub(&t_ans.y,&tmp,&P1->y);

	EFp_set(ANS,&t_ans);

	Fp_clear(&x);
	Fp_clear(&y);
	Fp_clear(&lambda);
	Fp_clear(&tmp);
	EFp_clear(&t_ans);
}
int EFp_cmp(struct EFp *A,struct EFp *B){
	if(Fp_cmp(&A->x,&B->x)==0 && Fp_cmp(&A->y,&B->y)==0){
		return 0;
	}
	return 1;
}
void EFp_random_set(struct EFp *ANS){
	struct Fp a,x;
	Fp_init(&a);
	Fp_init(&x);

	struct EFp P,Q;
	EFp_init(&P);
	EFp_init(&Q);

	mpz_t r_div_2;
	mpz_init(r_div_2);
	mpz_div_ui(r_div_2,order,2);
	// gmp_printf("%Zd\n",r_div_2);
	do{
		do{
			Fp_random(&x);
			Fp_mul(&a,&x,&x);
			Fp_mul(&a,&a,&x);
			mpz_add(a.x0,a.x0,b);
		}while(mpz_legendre(a.x0,prime)!=1);
		Q.infity=0;
		Fp_sqrt(&P.y,&a);
		Fp_set(&P.x,&x);

		EFp_SCM(&Q,&P,r_div_2);
	}while(Q.infity==TRUE);
	EFp_set(ANS,&P);


	Fp_clear(&a);
	Fp_clear(&x);
	EFp_clear(&P);
	EFp_clear(&Q);
}

//-----------------------------------------------------------------------------------------

void EFp3_init(struct EFp3 *A){
	Fp3_init(&A->x);
	Fp3_init(&A->y);
	A->infity=FALSE;
}
void EFp3_set(struct EFp3 *A,struct EFp3 *B){
	Fp3_set(&A->x,&B->x);
	Fp3_set(&A->y,&B->y);
	A->infity=B->infity;
}
void EFp3_set_infity(struct EFp3 *A){
	Fp3_set_ui(&A->x,0);
	Fp3_set_ui(&A->y,0);
	A->infity=TRUE;
}
void EFp3_clear(struct EFp3 *A){
	Fp3_clear(&A->x);
	Fp3_clear(&A->y);
}
void EFp3_printf(struct EFp3 *A){
	gmp_printf("((%Zd,%Zd,%Zd),(%Zd,%Zd,%Zd))\n\n",A->x.x0.x0,A->x.x1.x0,A->x.x2.x0,A->y.x0.x0,A->y.x1.x0,A->y.x2.x0);
}
void EFp3_ECD(struct EFp3 *ANS, struct EFp3 *P){
	if(P->infity==TRUE){
		EFp3_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp3_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp3_set_infity(ANS);
		return;
	}

	struct Fp3 x,y,lambda,tmp;
	struct EFp3 t_ans;
	Fp3_init(&x);
	Fp3_init(&lambda);
	Fp3_init(&tmp);
	Fp3_init(&y);
	EFp3_init(&t_ans);

	Fp3_mul(&x,&P->x,&P->x);
	Fp3_add(&tmp,&x,&x);
	Fp3_add(&x,&tmp,&x);
	Fp3_add(&y,&P->y,&P->y);
	Fp3_div(&lambda,&x,&y);
	Fp3_mul(&tmp,&lambda,&lambda);
	Fp3_add(&x,&P->x,&P->x);
	Fp3_sub(&x,&tmp,&x);
	Fp3_sub(&tmp,&P->x,&x);
	Fp3_set(&t_ans.x,&x);
	Fp3_mul(&tmp,&tmp,&lambda);
	Fp3_sub(&t_ans.y,&tmp,&P->y);

	EFp3_set(ANS,&t_ans);

	Fp3_clear(&x);
	Fp3_clear(&lambda);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&t_ans);
}
void EFp3_ECA(struct EFp3 *ANS, struct EFp3 *P1, struct EFp3 *P2){
	if(P2->infity==TRUE){//if P2==inf
		EFp3_set(ANS,P1);
		return;
	}
	else if(P1->infity==TRUE){//if P1==inf
		EFp3_set(ANS,P2);
		return;
	}
	else if(Fp3_cmp(&P1->x,&P2->x)==0&&Fp3_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp3_set_infity(ANS);
		return;
	}
	else if(EFp3_cmp(P1,P2)==0){ // P=Q
		EFp3_ECD(ANS,P1);
		return;
	}

	struct Fp3 x,y,lambda,tmp;
	struct EFp3 t_ans;

	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&lambda);
	Fp3_init(&tmp);
	EFp3_init(&t_ans);

	Fp3_sub(&x,&P2->x,&P1->x);
	Fp3_sub(&y,&P2->y,&P1->y);
	Fp3_div(&lambda,&y,&x);
	Fp3_mul(&tmp,&lambda,&lambda);
	Fp3_add(&x,&P1->x,&P2->x);
	Fp3_sub(&x,&tmp,&x);
	Fp3_sub(&tmp,&P1->x,&x);
	Fp3_set(&t_ans.x,&x);
	Fp3_mul(&tmp,&tmp,&lambda);
	Fp3_sub(&t_ans.y,&tmp,&P1->y);

	EFp3_set(ANS,&t_ans);

	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&lambda);
	Fp3_clear(&tmp);
	EFp3_clear(&t_ans);
}
void EFp3_SCM(struct EFp3 *ANS, struct EFp3 *P,mpz_t j){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(j,2);

	struct EFp3 Q;
	EFp3_init(&Q);
	EFp3_set(&Q,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp3_ECD(&Q,&Q);
			// EFp3_printf(&Q);
			EFp3_ECA(&Q,&Q,P);
		}else{
			EFp3_ECD(&Q,&Q);
		}
	}

	EFp3_set(ANS,&Q);
	EFp3_clear(&Q);
	return;
}
int EFp3_cmp(struct EFp3 *A,struct EFp3 *B){
	if(Fp3_cmp(&A->x,&B->x)==0 && Fp3_cmp(&A->y,&B->y)==0){
		return 0;
	}
	return 1;
}
void EFp3_set_EFp(struct EFp3 *ANS,struct EFp *A){
	Fp3_set_ui(&ANS->x,0);
	Fp3_set_ui(&ANS->y,0);

	Fp_set(&ANS->x.x0,&A->x);
	Fp_set(&ANS->y.x0,&A->y);
	ANS->infity=A->infity;
}
int EFp3_set_EFp18_Sparse(struct EFp3 *ANS,struct EFp18 *A){
	Fp3_set_ui(&ANS->x,0);
	Fp3_set_ui(&ANS->y,0);
	struct Fp3 cmp;
	Fp3_init(&cmp);
	if(Fp3_cmp(&A->x.x2.x0,&cmp) && Fp3_cmp(&A->y.x0.x1,&cmp)){
		Fp3_set(&ANS->x,&A->x.x2.x0);
		Fp3_set(&ANS->y,&A->y.x0.x1);
		ANS->infity=A->infity;
		return 1;
	}else if(Fp3_cmp(&A->x.x1.x1,&cmp) && Fp3_cmp(&A->y.x0.x1,&cmp)){
		Fp3_set(&ANS->x,&A->x.x1.x1);
		Fp3_set(&ANS->y,&A->y.x0.x1);
		ANS->infity=A->infity;
		return 2;
	}else{
		return 0;
	}
}
void EFp3_random_set(struct EFp3 *ANS){
	struct Fp3 a,x;
	Fp3_init(&a);
	Fp3_init(&x);

	struct EFp3 P,Q;
	EFp3_init(&P);
	EFp3_init(&Q);

	mpz_t p3,t3,r3,tmp;
	mpz_init(p3);
	mpz_init(t3);
	mpz_init(r3);
	mpz_init(tmp);

	mpz_mul(p3,prime,prime);
	mpz_mul(p3,p3,prime);//p3=p^3

	mpz_mul(t3,trace,trace);
	mpz_mul(t3,t3,trace);//t3=t^3

	mpz_mul(tmp,trace,prime);
	mpz_mul_ui(tmp,tmp,3);//tmp=3tp

	mpz_sub(r3,p3,t3);//r3=p3-t3
	mpz_add(r3,r3,tmp);//r3=p3-t3+3tp
	mpz_add_ui(r3,r3,1);//r3=p3-t3+3tp+1

	mpz_tdiv_q(r3,r3,order);
	mpz_tdiv_q(r3,r3,order);
	do{
		Fp3_random(&x);
		Fp3_mul(&a,&x,&x);
		Fp3_mul(&a,&a,&x);
		mpz_add(a.x0.x0,a.x0.x0,b);
	}while(Fp3_legendre(&a)!=1);
	P.infity=FALSE;
	Fp3_sqrt(&P.y,&a);
	Fp3_set(&P.x,&x);

	EFp3_SCM(&P,&P,r3);
	EFp3_set(ANS,&P);


	Fp3_clear(&a);
	Fp3_clear(&x);
	EFp3_clear(&P);
	mpz_clear(p3);
}
void EFp3_neg(struct EFp3 *ANS, struct EFp3 *A){
	struct EFp3 tmp;
	EFp3_init(&tmp);
	Fp3_neg(&tmp.y,&A->y);
	Fp3_set(&tmp.x,&A->x);

	EFp3_set(ANS,&tmp);
	EFp3_clear(&tmp);
}

//-----------------------------------------------------------------------------------------

void EFp6_init(struct EFp6 *A){
	Fp6_init(&A->x);
	Fp6_init(&A->y);
	A->infity=FALSE;
}
void EFp6_set(struct EFp6 *A,struct EFp6 *B){
	Fp6_set(&A->x,&B->x);
	Fp6_set(&A->y,&B->y);
	A->infity=B->infity;
}
void EFp6_set_infity(struct EFp6 *A){
	Fp6_set_ui(&A->x,0);
	Fp6_set_ui(&A->y,0);
	A->infity=TRUE;
}
void EFp6_clear(struct EFp6 *A){
	Fp6_clear(&A->x);
	Fp6_clear(&A->y);
}
void EFp6_printf(struct EFp6 *A){
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x.x0.x0.x0,A->x.x0.x1.x0,A->x.x0.x2.x0,A->x.x1.x0.x0,A->x.x1.x1.x0,A->x.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->y.x0.x0.x0,A->y.x0.x1.x0,A->y.x0.x2.x0,A->y.x1.x0.x0,A->y.x1.x1.x0,A->y.x1.x2.x0);
}
void EFp6_ECD(struct EFp6 *ANS, struct EFp6 *P){
	if(P->infity==TRUE){
		EFp6_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp6_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp6_set_infity(ANS);
		return;
	}

	struct Fp6 x,y,lambda,tmp;
	struct EFp6 t_ans;
	Fp6_init(&x);
	Fp6_init(&lambda);
	Fp6_init(&tmp);
	Fp6_init(&y);
	EFp6_init(&t_ans);

	Fp6_mul(&x,&P->x,&P->x);
	Fp6_add(&tmp,&x,&x);
	Fp6_add(&x,&tmp,&x);
	Fp6_add(&y,&P->y,&P->y);
	Fp6_div(&lambda,&x,&y);
	Fp6_mul(&tmp,&lambda,&lambda);
	Fp6_add(&x,&P->x,&P->x);
	Fp6_sub(&x,&tmp,&x);
	Fp6_sub(&tmp,&P->x,&x);
	Fp6_set(&t_ans.x,&x);
	Fp6_mul(&tmp,&tmp,&lambda);
	Fp6_sub(&t_ans.y,&tmp,&P->y);

	EFp6_set(ANS,&t_ans);

	Fp6_clear(&x);
	Fp6_clear(&lambda);
	Fp6_clear(&y);
	Fp6_clear(&tmp);
	EFp6_clear(&t_ans);
}
void EFp6_ECA(struct EFp6 *ANS, struct EFp6 *P1, struct EFp6 *P2){
	if(P2->infity==TRUE){//if P2==inf
		EFp6_set(ANS,P1);
		return;
	}
	else if(P1->infity==TRUE){//if P1==inf
		EFp6_set(ANS,P2);
		return;
	}
	else if(Fp6_cmp(&P1->x,&P2->x)==0&&Fp6_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp6_set_infity(ANS);
		return;
	}
	else if(EFp6_cmp(P1,P2)==0){ // P=Q
		EFp6_ECD(ANS,P1);
		return;
	}

	struct Fp6 x,y,lambda,tmp;
	struct EFp6 t_ans;

	Fp6_init(&x);
	Fp6_init(&y);
	Fp6_init(&lambda);
	Fp6_init(&tmp);
	EFp6_init(&t_ans);

	Fp6_sub(&x,&P2->x,&P1->x);
	Fp6_sub(&y,&P2->y,&P1->y);
	Fp6_div(&lambda,&y,&x);
	Fp6_mul(&tmp,&lambda,&lambda);
	Fp6_add(&x,&P1->x,&P2->x);
	Fp6_sub(&x,&tmp,&x);
	Fp6_sub(&tmp,&P1->x,&x);
	Fp6_set(&t_ans.x,&x);
	Fp6_mul(&tmp,&tmp,&lambda);
	Fp6_sub(&t_ans.y,&tmp,&P1->y);

	EFp6_set(ANS,&t_ans);

	Fp6_clear(&x);
	Fp6_clear(&y);
	Fp6_clear(&lambda);
	Fp6_clear(&tmp);
	EFp6_clear(&t_ans);
}
void EFp6_SCM(struct EFp6 *ANS, struct EFp6 *P,mpz_t j){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(j,2);

	struct EFp6 Q;
	EFp6_init(&Q);
	EFp6_set(&Q,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp6_ECD(&Q,&Q);
			EFp6_ECA(&Q,&Q,P);
		}else{
			EFp6_ECD(&Q,&Q);
		}
	}

	EFp6_set(ANS,&Q);
	EFp6_clear(&Q);
	return;
}
int EFp6_cmp(struct EFp6 *A,struct EFp6 *B){
	if(Fp6_cmp(&A->x,&B->x)==0 && Fp6_cmp(&A->y,&B->y)==0){
		return 0;
	}
	return 1;
}
void EFp6_random_set(struct EFp6 *ANS){
	struct EFp6 P;
	EFp6_init(&P);

	struct Fp6 x,a;
	Fp6_init(&a);
	Fp6_init(&x);

	//t12=a^12+b^12={(t^3-3pt)^2-2p^3}^2-2p^6
	mpz_t tmp1,tmp2,p_pow;
	mpz_t r6,t6;
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(p_pow);
	mpz_init(r6);
	mpz_init(t6);

	mpz_pow_ui(tmp1,trace,3);
	mpz_mul(p_pow,prime,trace);
	mpz_mul_ui(p_pow,p_pow,3);
	mpz_sub(tmp1,tmp1,p_pow);
	mpz_mul(tmp1,tmp1,tmp1);

	mpz_pow_ui(p_pow,prime,3);
	mpz_mul_ui(tmp2,p_pow,2);
	mpz_sub(t6,tmp1,tmp2);
	// mpz_pow_ui(tmp1,tmp1,2);

	mpz_pow_ui(p_pow,p_pow,2);
	// mpz_mul_ui(tmp2,p_pow,2);
	// mpz_sub(t6,tmp1,tmp2);

	// mpz_pow_ui(p_pow,p_pow,2);
	mpz_add_ui(tmp1,p_pow,1);
	mpz_sub(r6,tmp1,t6);

	do{
		Fp6_random(&x);
		Fp6_mul(&a,&x,&x);
		Fp6_mul(&a,&a,&x);
		mpz_add(a.x0.x0.x0,a.x0.x0.x0,b);
	}while(Fp6_legendre(&a)!=1);
	Fp6_sqrt(&P.y,&a);
	Fp6_set(&P.x,&x);

	mpz_t r6_div_r2;
	mpz_init(r6_div_r2);
	mpz_div(r6_div_r2,r6,order);
	mpz_div(r6_div_r2,r6_div_r2,order);

	EFp6_SCM(ANS,&P,r6_div_r2);

	EFp6_clear(&P);
	Fp6_clear(&a);
	Fp6_clear(&x);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(p_pow);
	mpz_clear(r6);
	mpz_clear(t6);
	mpz_clear(r6_div_r2);
}

void EFp3_type2_ECD(struct EFp3 *ANS, struct EFp3 *P){
	if(P->infity==TRUE){
		EFp3_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp3_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp3_set_infity(ANS);
		return;
	}

	struct Fp3 x,y,lambda,tmp;
	struct EFp3 t_ans;
	Fp3_init(&x);
	Fp3_init(&lambda);
	Fp3_init(&tmp);
	Fp3_init(&y);
	EFp3_init(&t_ans);

	Fp3_mul(&x,&P->x,&P->x);
	Fp3_add(&tmp,&x,&x);
	Fp3_add(&x,&tmp,&x);
	Fp3_add(&y,&P->y,&P->y);
	Fp3_div(&lambda,&x,&y);
	Fp3_mul(&tmp,&lambda,&lambda);
	Fp3_mul_xi(&tmp,&tmp);
	Fp3_add(&x,&P->x,&P->x);
	Fp3_sub(&x,&tmp,&x);
	Fp3_sub(&tmp,&P->x,&x);
	Fp3_set(&t_ans.x,&x);
	Fp3_mul(&tmp,&tmp,&lambda);
	Fp3_mul_xi(&tmp,&tmp);
	Fp3_sub(&t_ans.y,&tmp,&P->y);

	EFp3_set(ANS,&t_ans);

	Fp3_clear(&x);
	Fp3_clear(&lambda);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&t_ans);
}
void EFp3_type2_ECA(struct EFp3 *ANS, struct EFp3 *P1, struct EFp3 *P2){
	if(P2->infity==TRUE){//if P2==inf
		EFp3_set(ANS,P1);
		return;
	}
	else if(P1->infity==TRUE){//if P1==inf
		EFp3_set(ANS,P2);
		return;
	}
	else if(Fp3_cmp(&P1->x,&P2->x)==0&&Fp3_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp3_set_infity(ANS);
		return;
	}
	else if(EFp3_cmp(P1,P2)==0){ // P=Q
		EFp3_type2_ECD(ANS,P1);
		return;
	}

	struct Fp3 x,y,lambda,tmp;
	struct EFp3 t_ans;

	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&lambda);
	Fp3_init(&tmp);
	EFp3_init(&t_ans);

	Fp3_sub(&x,&P2->x,&P1->x);
	Fp3_sub(&y,&P2->y,&P1->y);
	Fp3_div(&lambda,&y,&x);
	Fp3_mul(&tmp,&lambda,&lambda);
	Fp3_mul_xi_inv(&tmp,&tmp);
	Fp3_add(&x,&P1->x,&P2->x);
	Fp3_sub(&x,&tmp,&x);
	Fp3_sub(&tmp,&P1->x,&x);
	Fp3_set(&t_ans.x,&x);
	Fp3_mul(&tmp,&tmp,&lambda);
	Fp3_sub(&t_ans.y,&tmp,&P1->y);

	EFp3_set(ANS,&t_ans);

	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&lambda);
	Fp3_clear(&tmp);
	EFp3_clear(&t_ans);
}
void EFp3_type2_SCM(struct EFp3 *ANS, struct EFp3 *P,mpz_t j){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(j,2);

	struct EFp3 Q;
	EFp3_init(&Q);
	EFp3_set(&Q,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp3_type2_ECD(&Q,&Q);
			// EFp3_printf(&Q);
			EFp3_type2_ECA(&Q,&Q,P);
		}else{
			EFp3_type2_ECD(&Q,&Q);
		}
	}

	EFp3_set(ANS,&Q);
	EFp3_clear(&Q);
	return;
}

//-----------------------------------------------------------------------------------------

void EFp18_init(struct EFp18 *A){
	Fp18_init(&A->x);
	Fp18_init(&A->y);
	A->infity=FALSE;
}
void EFp18_set(struct EFp18 *A,struct EFp18 *B){
	Fp18_set(&A->x,&B->x);
	Fp18_set(&A->y,&B->y);
	A->infity=B->infity;
}
void EFp18_set_infity(struct EFp18 *A){
	Fp18_set_ui(&A->x,0);
	Fp18_set_ui(&A->y,0);
	A->infity=TRUE;
}
void EFp18_clear(struct EFp18 *A){
	Fp18_clear(&A->x);
	Fp18_clear(&A->y);
}
void EFp18_printf(struct EFp18 *A){
	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x.x0.x0.x0.x0,A->x.x0.x0.x1.x0,A->x.x0.x0.x2.x0,A->x.x0.x1.x0.x0,A->x.x0.x1.x1.x0,A->x.x0.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x.x1.x0.x0.x0,A->x.x1.x0.x1.x0,A->x.x1.x0.x2.x0,A->x.x1.x1.x0.x0,A->x.x1.x1.x1.x0,A->x.x1.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",A->x.x2.x0.x0.x0,A->x.x2.x0.x1.x0,A->x.x2.x0.x2.x0,A->x.x2.x1.x0.x0,A->x.x2.x1.x1.x0,A->x.x2.x1.x2.x0);

	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->y.x0.x0.x0.x0,A->y.x0.x0.x1.x0,A->y.x0.x0.x2.x0,A->y.x0.x1.x0.x0,A->y.x0.x1.x1.x0,A->y.x0.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->y.x1.x0.x0.x0,A->y.x1.x0.x1.x0,A->y.x1.x0.x2.x0,A->y.x1.x1.x0.x0,A->y.x1.x1.x1.x0,A->y.x1.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n\n",A->y.x2.x0.x0.x0,A->y.x2.x0.x1.x0,A->y.x2.x0.x2.x0,A->y.x2.x1.x0.x0,A->y.x2.x1.x1.x0,A->y.x2.x1.x2.x0);
}
void EFp18_ECD(struct EFp18 *ANS, struct EFp18 *P){
	if(P->infity==TRUE){
		EFp18_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp18_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp18_set_infity(ANS);
		return;
	}

	struct Fp18 x,y,lambda,tmp;
	struct EFp18 t_ans;
	Fp18_init(&x);
	Fp18_init(&lambda);
	Fp18_init(&tmp);
	Fp18_init(&y);
	EFp18_init(&t_ans);

	Fp18_mul(&x,&P->x,&P->x);
	Fp18_add(&tmp,&x,&x);
	Fp18_add(&x,&tmp,&x);
	Fp18_add(&y,&P->y,&P->y);
	Fp18_div(&lambda,&x,&y);
	Fp18_mul(&tmp,&lambda,&lambda);
	Fp18_add(&x,&P->x,&P->x);
	Fp18_sub(&x,&tmp,&x);
	Fp18_sub(&tmp,&P->x,&x);
	Fp18_set(&t_ans.x,&x);
	Fp18_mul(&tmp,&tmp,&lambda);
	Fp18_sub(&t_ans.y,&tmp,&P->y);

	EFp18_set(ANS,&t_ans);

	Fp18_clear(&x);
	Fp18_clear(&lambda);
	Fp18_clear(&y);
	Fp18_clear(&tmp);
	EFp18_clear(&t_ans);
}
void EFp18_ECA(struct EFp18 *ANS, struct EFp18 *P1, struct EFp18 *P2){
	if(P2->infity==TRUE){//if P2==inf
		EFp18_set(ANS,P1);
		return;
	}
	else if(P1->infity==TRUE){//if P1==inf
		EFp18_set(ANS,P2);
		return;
	}
	else if(Fp18_cmp(&P1->x,&P2->x)==0&&Fp18_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp18_set_infity(ANS);
		return;
	}
	else if(EFp18_cmp(P1,P2)==0){ // P=Q
		EFp18_ECD(ANS,P1);
		return;
	}

	struct Fp18 x,y,lambda,tmp;
	struct EFp18 t_ans;

	Fp18_init(&x);
	Fp18_init(&y);
	Fp18_init(&lambda);
	Fp18_init(&tmp);
	EFp18_init(&t_ans);

	Fp18_sub(&x,&P2->x,&P1->x);
	Fp18_sub(&y,&P2->y,&P1->y);
	Fp18_div(&lambda,&y,&x);
	Fp18_mul(&tmp,&lambda,&lambda);
	Fp18_add(&x,&P1->x,&P2->x);
	Fp18_sub(&x,&tmp,&x);
	Fp18_sub(&tmp,&P1->x,&x);
	Fp18_set(&t_ans.x,&x);
	Fp18_mul(&tmp,&tmp,&lambda);
	Fp18_sub(&t_ans.y,&tmp,&P1->y);

	EFp18_set(ANS,&t_ans);

	Fp18_clear(&x);
	Fp18_clear(&y);
	Fp18_clear(&lambda);
	Fp18_clear(&tmp);
	EFp18_clear(&t_ans);
}
void EFp18_SCM(struct EFp18 *ANS, struct EFp18 *P,mpz_t j){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(j,2);

	struct EFp18 Q;
	EFp18_init(&Q);
	EFp18_set(&Q,P);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			EFp18_ECD(&Q,&Q);
			// EFp18_printf(&Q);
			EFp18_ECA(&Q,&Q,P);
		}else{
			EFp18_ECD(&Q,&Q);
		}
	}

	EFp18_set(ANS,&Q);
	EFp18_clear(&Q);
	return;
}
int EFp18_cmp(struct EFp18 *A,struct EFp18 *B){
	if(Fp18_cmp(&A->x,&B->x)==0 && Fp18_cmp(&A->y,&B->y)==0){
		return 0;
	}
	return 1;
}
void EFp18_random_set(struct EFp18 *ANS){
	int i=1,j;
	struct EFp18 P;
	EFp18_init(&P);

	struct Fp18 x,a;
	mpz_t r18_q,r18_r;
	mpz_t tmp1,tmp2,tmp3,tmp,t18,r18;
	struct EFp18 t[256];
	struct EFp18 t_ans;

	Fp18_init(&a);
	Fp18_init(&x);

	mpz_init(r18_q);
	mpz_init(r18_r);

	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(tmp3);
	mpz_init(tmp);
	mpz_init(t18);
	mpz_init(r18);
	for(i=0;i<256;i++){
		EFp18_init(&t[i]);
	}
	EFp18_init(&t_ans);
	//t18=a^18+b^18={(t^3-3pt)^3-3p^3(t^3-3pt)}^2-2p^9

	mpz_mul(tmp1,prime,trace);
	mpz_mul_ui(tmp1,tmp1,3);

	mpz_mul(tmp2,trace,trace);
	mpz_mul(tmp2,tmp2,trace);

	mpz_sub(tmp1,tmp2,tmp1);//tmp1=t^3-3pt

	mpz_mul(tmp,tmp1,tmp1);
	mpz_mul(tmp,tmp,tmp1);//tmp=(t^3-3pt)^3

	mpz_mul(tmp2,prime,prime);
	mpz_mul(tmp2,tmp2,prime);//tmp2=p^3

	mpz_mul(tmp3,tmp1,tmp2);
	mpz_mul_ui(tmp3,tmp3,3);//tmp3=3p^3(t^3-3pt)

	mpz_sub(tmp,tmp,tmp3);//tmp=(t^3-3pt)^3-3p^3(t^3-3pt)
	mpz_mul(tmp,tmp,tmp);//tmp=tmp^2

	mpz_pow_ui(tmp2,tmp2,3);//tmp2=p^9
	mpz_add(tmp3,tmp2,tmp2);//tmp3=2*tmp2

	mpz_sub(t18,tmp,tmp3);//t18=tmp-tmp3

	mpz_mul(tmp3,tmp2,tmp2);//tmp3=p^18

	mpz_add_ui(tmp3,tmp3,1);
	mpz_sub(r18,tmp3,t18);
	mpz_div(r18,r18,order);
	mpz_div(r18,r18,order);

	do{
		Fp18_random(&x);
		Fp18_mul(&a,&x,&x);
		Fp18_mul(&a,&a,&x);
		mpz_add(a.x0.x0.x0.x0,a.x0.x0.x0.x0,b);
	}while(Fp18_legendre(&a)!=1);
	Fp18_sqrt(&P.y,&a);
	Fp18_set(&P.x,&x);

	mpz_set(r18_q,r18);
	int r18_length=(int)(mpz_sizeinbase(r18_q,2)-1)/8+1;
	int r18_bit_separate[r18_length];

	i=0;
	//256-adic representation
	while(mpz_cmp_ui(r18_q,0)!=0){
		mpz_tdiv_qr_ui(r18_q,r18_r,r18_q,256);
		r18_bit_separate[i]=(unsigned int)mpz_get_ui(r18_r);
		// printf("%d\n",r18_bit_separate[i]);
		i++;
	}
	EFp18_set(&t[1],&P);
	EFp18_ECD(&t[2],&t[1]);
	EFp18_ECD(&t[4],&t[2]);
	EFp18_ECD(&t[8],&t[4]);
	EFp18_ECD(&t[16],&t[8]);
	EFp18_ECD(&t[32],&t[16]);
	EFp18_ECD(&t[64],&t[32]);
	EFp18_ECD(&t[128],&t[64]);

	for(i=1;i<256;i=i*2){
		for(j=1;j<i;j++){
			EFp18_ECA(&t[i+j],&t[i],&t[j]);
		}
	}
	EFp18_set(&t_ans,&t[r18_bit_separate[r18_length-1]]);
	for(i=r18_length-2;i>=0;i--){
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		EFp18_ECD(&t_ans,&t_ans);
		if(r18_bit_separate[i]!=0){
			EFp18_ECA(&t_ans,&t_ans,&t[r18_bit_separate[i]]);
		}
	}

	EFp18_set(ANS,&t_ans);

	EFp18_clear(&P);
	Fp18_clear(&a);
	Fp18_clear(&x);

	mpz_clear(r18_q);
	mpz_clear(r18_r);

	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	mpz_clear(tmp);
	mpz_clear(t18);
	mpz_clear(r18);
	for(i=0;i<256;i++){
		EFp18_clear(&t[i]);
	}
	EFp18_clear(&t_ans);

}
void EFp18_set_EFp(struct EFp18 *A,struct EFp *B){
	Fp18_set_ui(&A->x,0);
	Fp18_set_ui(&A->y,0);

	Fp_set(&A->x.x0.x0.x0,&B->x);
	Fp_set(&A->y.x0.x0.x0,&B->y);
	A->infity=B->infity;
}
void EFp18_frobenius_map(struct EFp18 *ANS,struct EFp18 *A,int i){
	struct EFp18 tmp;
	EFp18_init(&tmp);

	Fp18_frobenius_map(&tmp.x,&A->x,i);
	Fp18_frobenius_map(&tmp.y,&A->y,i);

	EFp18_set(ANS,&tmp);

	EFp18_clear(&tmp);
}
void EFp18_random_set_G2(struct EFp18 *ANS){
	struct EFp18 P,P_frobenius,tmp_EFp18;
	EFp18_init(&P);
	EFp18_init(&P_frobenius);
	EFp18_init(&tmp_EFp18);

	EFp18_random_set(&P);

	EFp18_frobenius_map(&P_frobenius,&P,1);
	Fp18_neg(&tmp_EFp18.y,&P.y);
	Fp18_set(&tmp_EFp18.x,&P.x);

	EFp18_ECA(&tmp_EFp18,&tmp_EFp18,&P_frobenius);

	EFp18_set(ANS,&tmp_EFp18);

	EFp18_clear(&P);
	EFp18_clear(&P_frobenius);
	EFp18_clear(&tmp_EFp18);
}

#endif //_Elliptic_Curve_C_
