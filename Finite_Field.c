#include"embedding_degree18.h"
#ifndef _Finite_Field_C_
#define _Finite_Field_C_

void Fp_init(struct Fp *A){
	mpz_init(A->x0);
}
void Fp_set(struct Fp *ANS,struct Fp *A){
	mpz_set(ANS->x0,A->x0);
}
void Fp_set_ui(struct Fp *A,signed long int B){
	mpz_set_ui(A->x0,B);
}
void Fp_random(struct Fp *A){
	mpz_random(A->x0,10);
	mpz_mod(A->x0,A->x0,prime);
}
void Fp_clear(struct Fp *A){
	mpz_clear(A->x0);
}
void Fp_printf(struct Fp *A){
	gmp_printf("%Zd\n",A->x0);
}
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_add(tmp.x0,A->x0,B->x0);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_add_ui(tmp.x0,A->x0,B);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_sub(tmp.x0,A->x0,B->x0);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_sub_ui(tmp.x0,A->x0,B);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_mul(tmp.x0,A->x0,B->x0);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_mul_ui(tmp.x0,A->x0,B);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_div(struct Fp *ANS,struct Fp *A,struct Fp *B){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_invert(tmp.x0,B->x0,prime);
	mpz_mul(tmp.x0,A->x0,tmp.x0);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t j){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(j,2);
	//printf("r= %d\n",r);

	struct Fp answer_tmp;
	Fp_init(&answer_tmp);
	Fp_set(&answer_tmp,A);

	struct Fp in_tmp;
	Fp_init(&in_tmp);
	Fp_set(&in_tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(j,i)==1){
			Fp_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
			Fp_mul(&answer_tmp,&answer_tmp,&in_tmp);//*a
		}else{
			Fp_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
		}
	}

	Fp_set(ANS,&answer_tmp);

	Fp_clear(&answer_tmp);
	Fp_clear(&in_tmp);
}
void Fp_sqrt(struct Fp *ANS,struct Fp *A){
	struct Fp n_tmp,y_tmp,x_tmp,b_tmp,t_tmp,tmp_Fp;
	Fp_init(&n_tmp);
	Fp_init(&y_tmp);
	Fp_init(&x_tmp);
	Fp_init(&b_tmp);
	Fp_init(&t_tmp);
	Fp_init(&tmp_Fp);

	Fp_set(&n_tmp,A);

	mpz_t tmp_mpz,q_tmp,e_tmp,r_tmp,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q_tmp);
	mpz_init(e_tmp);
	mpz_init(r_tmp);
	mpz_init(set_1);
	mpz_init(set_2);

	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(mpz_legendre(n_tmp.x0,prime)!=-1){
		Fp_add_ui(&n_tmp,&n_tmp,1);
	}

	mpz_set(q_tmp,prime);
	mpz_sub_ui(q_tmp,q_tmp,1);
	mpz_set_ui(e_tmp,0);

	while(mpz_odd_p(q_tmp)==0){
		mpz_add_ui(e_tmp,e_tmp,1);
		mpz_div_ui(q_tmp,q_tmp,2);
	}

	Fp_pow(&y_tmp,&n_tmp,q_tmp);

	mpz_set(r_tmp,e_tmp);

	mpz_sub_ui(tmp_mpz,q_tmp,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);
	Fp_pow(&x_tmp,A,tmp_mpz);
	Fp_pow(&tmp_Fp,&x_tmp,set_2);
	Fp_mul(&b_tmp,&tmp_Fp,A);
	Fp_mul(&x_tmp,&x_tmp,A);

	int m;

	while(Fp_cmp_mpz(&b_tmp,set_1)==1){
		m=-1;
		while(Fp_cmp_mpz(&tmp_Fp,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp_pow(&tmp_Fp,&b_tmp,tmp_mpz);
		}
		// gmp_printf("%Zd,%Zd\n",x.x0.x0,x.x1.x0);
		mpz_sub_ui(tmp_mpz,r_tmp,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,prime);
		Fp_pow(&t_tmp,&y_tmp,tmp_mpz);
		Fp_pow(&y_tmp,&t_tmp,set_2);
		mpz_set_ui(r_tmp,m);
		Fp_mul(&x_tmp,&x_tmp,&t_tmp);
		Fp_mul(&b_tmp,&b_tmp,&y_tmp);
		Fp_set(&tmp_Fp,&b_tmp);
	}

	Fp_set(ANS,&x_tmp);

	Fp_clear(&n_tmp);
	Fp_clear(&y_tmp);
	Fp_clear(&x_tmp);
	Fp_clear(&b_tmp);
	Fp_clear(&t_tmp);
	Fp_clear(&tmp_Fp);
	mpz_clear(tmp_mpz);
	mpz_clear(q_tmp);
	mpz_clear(e_tmp);
	mpz_clear(r_tmp);
	mpz_clear(set_1);
	mpz_clear(set_2);
}
void Fp_neg(struct Fp *ANS,struct Fp *A){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_sub(tmp.x0,prime,A->x0);
	mpz_mod(tmp.x0,tmp.x0,prime);

	Fp_set(ANS,&tmp);

	Fp_clear(&tmp);
}
int Fp_cmp(struct Fp *A,struct Fp *B){
	if(mpz_cmp(A->x0,B->x0)==0){
		return 0;
	}
	return 1;
}
int Fp_cmp_mpz(struct Fp *A,mpz_t B){
	if(mpz_cmp(A->x0,B)==0){
		return 0;
	}
	return 1;
}

//-----------------------------------------------------------------------------------------

void Fp3_init(struct Fp3 *A){
	Fp_init(&A->x0);
	Fp_init(&A->x1);
	Fp_init(&A->x2);
}
void Fp3_set(struct Fp3 *ANS,struct Fp3 *A){
	Fp_set(&ANS->x0,&A->x0);
	Fp_set(&ANS->x1,&A->x1);
	Fp_set(&ANS->x2,&A->x2);
}
void Fp3_set_ui(struct Fp3 *A,unsigned long int B){
	Fp_set_ui(&A->x0,B);
	Fp_set_ui(&A->x1,B);
	Fp_set_ui(&A->x2,B);
}
void Fp3_random(struct Fp3 *A){
	Fp_random(&A->x0);
	Fp_random(&A->x1);
	Fp_random(&A->x2);
}
void Fp3_clear(struct Fp3 *A){
	Fp_clear(&A->x0);
	Fp_clear(&A->x1);
	Fp_clear(&A->x2);
}
void Fp3_printf(struct Fp3 *A){
	gmp_printf("%Zd,%Zd,%Zd\n",A->x0.x0,A->x1.x0,A->x2.x0);
}
void Fp3_add(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B){
	struct Fp3 tmp;
	Fp3_init(&tmp);

	Fp_add(&tmp.x0,&A->x0,&B->x0);
	Fp_add(&tmp.x1,&A->x1,&B->x1);
	Fp_add(&tmp.x2,&A->x2,&B->x2);

	Fp3_set(ANS,&tmp);

	Fp3_clear(&tmp);
}
void Fp3_add_ui(struct Fp3 *ANS,struct Fp3 *A,unsigned long int B){
	struct Fp3 tmp;
	Fp3_init(&tmp);

	Fp_add_ui(&tmp.x0,&A->x0,B);
	Fp_add_ui(&tmp.x1,&A->x1,B);
	Fp_add_ui(&tmp.x2,&A->x2,B);

	Fp3_set(ANS,&tmp);

	Fp3_clear(&tmp);
}
void Fp3_sub(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B){
	struct Fp3 tmp;
	Fp3_init(&tmp);

	Fp_sub(&tmp.x0,&A->x0,&B->x0);
	Fp_sub(&tmp.x1,&A->x1,&B->x1);
	Fp_sub(&tmp.x2,&A->x2,&B->x2);

	Fp3_set(ANS,&tmp);

	Fp3_clear(&tmp);
}
void Fp3_mul(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B){
	//(x0,x1,x2)*(y0,y1,y2)=(x0y0+xi((x1+x2)(y1+y2)-x1y1-x2y2),xix2y2+(x0+x1)(y0+y1)-x0y0-x1y1,x1y1+(x0+x2)(y0+y2)-x0y0-x2y2)
	struct Fp tmp00,tmp11,tmp22,tmpx01,tmpx12,tmpx20,tmpy01,tmpy12,tmpy20,t0,t1,t2,tmp;
	struct Fp3 t_ans;
	Fp_init(&tmp00);
	Fp_init(&tmp11);
	Fp_init(&tmp22);
	Fp_init(&tmpx01);
	Fp_init(&tmpx12);
	Fp_init(&tmpx20);
	Fp_init(&tmpy01);
	Fp_init(&tmpy12);
	Fp_init(&tmpy20);
	Fp_init(&t0);
	Fp_init(&t1);
	Fp_init(&t2);
	Fp_init(&tmp);
	Fp3_init(&t_ans);

	Fp_mul(&tmp00,&A->x0,&B->x0);//x0*y0
	Fp_mul(&tmp11,&A->x1,&B->x1);//x1*y1
	Fp_mul(&tmp22,&A->x2,&B->x2);//x2*y2

	Fp_add(&tmpx01,&A->x0,&A->x1);//x0+x1
	Fp_add(&tmpx12,&A->x1,&A->x2);//x1+x2
	Fp_add(&tmpx20,&A->x0,&A->x2);//x2+x0
	Fp_add(&tmpy01,&B->x0,&B->x1);//y0+y1
	Fp_add(&tmpy12,&B->x1,&B->x2);//y1+y2
	Fp_add(&tmpy20,&B->x0,&B->x2);//y2+y0

	Fp_mul(&t0,&tmpx01,&tmpy01);//(x0+x1)(y0+y1)
	Fp_mul(&t1,&tmpx12,&tmpy12);//(x1+x2)(y1+y2)
	Fp_mul(&t2,&tmpx20,&tmpy20);//(x2+x0)(y2+y0)
	//(x0,x1,x2)*(y0,y1,y2)=(x0y0+xi((x1+x2)(y1+y2)-x1y1-x2y2),xix2y2+(x0+x1)(y0+y1)-x0y0-x1y1,x1y1+(x0+x2)(y0+y2)-x0y0-x2y2)
	Fp_sub(&t1,&t1,&tmp11);
	Fp_sub(&t1,&t1,&tmp22);//(x1+x2)(y1+y2)-x1y1-x2y2
	Fp_mul_ui(&tmp,&t1,c1);
	Fp_add(&t_ans.x0,&tmp00,&tmp);

	Fp_sub(&t0,&t0,&tmp00);
	Fp_sub(&t0,&t0,&tmp11);
	Fp_mul_ui(&tmp,&tmp22,c1);
	Fp_add(&t_ans.x1,&tmp,&t0);

	Fp_sub(&t2,&t2,&tmp00);
	Fp_sub(&t2,&t2,&tmp22);
	Fp_add(&t_ans.x2,&tmp11,&t2);

	Fp3_set(ANS,&t_ans);

	Fp_clear(&tmp00);
	Fp_clear(&tmp11);
	Fp_clear(&tmp22);
	Fp_clear(&tmpx01);
	Fp_clear(&tmpx12);
	Fp_clear(&tmpx20);
	Fp_clear(&tmpy01);
	Fp_clear(&tmpy12);
	Fp_clear(&tmpy20);
	Fp_clear(&t0);
	Fp_clear(&t1);
	Fp_clear(&t2);
	Fp_clear(&tmp);
	Fp3_clear(&t_ans);
}
void Fp3_mul_old(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B){
	struct Fp tmp1,tmp2,s1,s2,s3;
	struct Fp a1b1,a2b2,a3b3;
	struct Fp3 t_ans;
	Fp_init(&tmp1);
	Fp_init(&tmp2);
	Fp_init(&s1);
	Fp_init(&s2);
	Fp_init(&s3);
	Fp_init(&a1b1);
	Fp_init(&a2b2);
	Fp_init(&a3b3);
	Fp3_init(&t_ans);

	Fp_sub(&tmp1,&A->x0,&A->x1);
	Fp_sub(&tmp2,&B->x1,&B->x0);
	Fp_mul(&s1,&tmp1,&tmp2);

	Fp_sub(&tmp1,&A->x1,&A->x2);
	Fp_sub(&tmp2,&B->x2,&B->x1);
	Fp_mul(&s2,&tmp1,&tmp2);

	Fp_sub(&tmp1,&A->x0,&A->x2);
	Fp_sub(&tmp2,&B->x2,&B->x0);
	Fp_mul(&s3,&tmp1,&tmp2);

	Fp_mul(&a1b1,&A->x0,&B->x0);
	Fp_mul(&a2b2,&A->x1,&B->x1);
	Fp_mul(&a3b3,&A->x2,&B->x2);

	Fp_add(&tmp1,&s1,&s2);
	Fp_sub(&t_ans.x0,&tmp1,&a1b1);
	Fp_add(&tmp1,&s2,&s3);
	Fp_sub(&t_ans.x1,&tmp1,&a2b2);
	Fp_add(&tmp1,&s3,&s1);
	Fp_sub(&t_ans.x2,&tmp1,&a3b3);
	Fp3_set(ANS,&t_ans);
}
void Fp3_mul_Fp(struct Fp3 *ANS,struct Fp3 *A,struct Fp *B){
	struct Fp3 tmp;
	Fp3_init(&tmp);

	Fp_mul(&tmp.x0,&A->x0,B);
	Fp_mul(&tmp.x1,&A->x1,B);
	Fp_mul(&tmp.x2,&A->x2,B);

	Fp3_set(ANS,&tmp);

	Fp3_clear(&tmp);
}
void Fp3_mul_xi(struct Fp3 *ANS,struct Fp3 *A){
	struct Fp3 tmp;
	Fp3_init(&tmp);

	Fp_mul_ui(&tmp.x0,&A->x2,c1);
	Fp_set(&tmp.x1,&A->x0);
	Fp_set(&tmp.x2,&A->x1);

	Fp3_set(ANS,&tmp);

	Fp3_clear(&tmp);
}
void Fp3_mul_xi_inv(struct Fp3 *ANS,struct Fp3 *A){
	struct Fp3 tmp;
	Fp3_init(&tmp);
	struct Fp Fp_c1;
	Fp_init(&Fp_c1);
	Fp_set_ui(&Fp_c1,c1);

	Fp_div(&tmp.x2,&A->x0,&Fp_c1);
	Fp_set(&tmp.x0,&A->x1);
	Fp_set(&tmp.x1,&A->x2);

	Fp3_set(ANS,&tmp);

	Fp3_clear(&tmp);
	Fp_clear(&Fp_c1);
}
void Fp3_mul_ui(struct Fp3 *ANS,struct Fp3 *A,unsigned long int B){
	struct Fp3 tmp;
	Fp3_init(&tmp);

	Fp_mul_ui(&tmp.x0,&A->x0,B);
	Fp_mul_ui(&tmp.x1,&A->x1,B);
	Fp_mul_ui(&tmp.x2,&A->x2,B);

	Fp3_set(ANS,&tmp);
	Fp3_clear(&tmp);
}
void Fp3_invert(struct Fp3 *ANS, struct Fp3 *A){
	struct Fp3 t_ans;
	Fp3_init(&t_ans);

	struct Fp T0,T1,t0,t1,t2,t3;
	Fp_init(&T0);
	Fp_init(&T1);
	Fp_init(&t0);
	Fp_init(&t1);
	Fp_init(&t2);
	Fp_init(&t3);

	// An optimized version of Grewal's Algo. 3   (a,b,c)
	Fp_mul(&T0,&A->x0,&A->x0);
	Fp_mul_ui(&t0,&A->x1,c1);

	Fp_mul(&T1,&t0,&A->x2);
	Fp_sub(&t1,&T0,&T1); // t1=(a^2-bci) mod q

	Fp_mul(&T0,&A->x2,&A->x2);
	Fp_mul_ui(&T0,&T0,c1);
	Fp_mul(&T1,&A->x0,&A->x1);
	Fp_sub(&t2,&T0,&T1); // t2=(c^2i-ab) mod q

	Fp_mul(&T0,&A->x1,&A->x1);
	Fp_mul(&T1,&A->x0,&A->x2);
	Fp_sub(&t3,&T0,&T1); // t3=(b^2-ac) mod q

	Fp_mul(&T0,&t0,&t3);
	Fp_mul(&T1,&A->x0,&t1);
	Fp_add(&T0,&T0,&T1); // T0={bi(b^2-ac)+a(a^2-bci)} mod q

	Fp_mul_ui(&t0,&A->x2,c1);
	Fp_mul(&T1,&t0,&t2);
	Fp_add(&t0,&T0,&T1); // t0={ci(c^2i-ab)+{bi(b^2-ac)+a(a^2-bci)}} mod q .0

	mpz_invert(t0.x0,t0.x0,prime);

	Fp_mul(&t_ans.x0,&t1,&t0);
	Fp_mul(&t_ans.x1,&t2,&t0);
	Fp_mul(&t_ans.x2,&t3,&t0);

	Fp3_set(ANS,&t_ans);

	Fp3_clear(&t_ans);
	Fp_clear(&T0);
	Fp_clear(&T1);
	Fp_clear(&t0);
	Fp_clear(&t1);
	Fp_clear(&t2);
	Fp_clear(&t3);
}
void Fp3_div(struct Fp3 *ANS,struct Fp3 *A,struct Fp3 *B){
	struct Fp3 tmp,t_ans;
	Fp3_init(&tmp);
	Fp3_init(&t_ans);

	Fp3_invert(&tmp,B);
	Fp3_mul(&t_ans,A,&tmp);

	Fp3_set(ANS,&t_ans);

	Fp3_clear(&tmp);
	Fp3_clear(&t_ans);
}
void Fp3_pow(struct Fp3 *ANS,struct Fp3 *A,mpz_t B){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(B,2);
	//printf("r= %d\n",r);

	struct Fp3 answer_tmp;
	Fp3_init(&answer_tmp);
	Fp3_set(&answer_tmp,A);

	struct Fp3 in_tmp;
	Fp3_init(&in_tmp);
	Fp3_set(&in_tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			Fp3_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
			Fp3_mul(&answer_tmp,&answer_tmp,&in_tmp);//*a
		}else{
			Fp3_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
		}
	}

	Fp3_set(ANS,&answer_tmp);

	Fp3_clear(&answer_tmp);
	Fp3_clear(&in_tmp);
}
void Fp3_sqrt(struct Fp3 *ANS,struct Fp3 *A){
	struct Fp3 n,y,x,b,t,tmp_Fp3;
	Fp3_init(&n);
	Fp3_init(&y);
	Fp3_init(&x);
	Fp3_init(&b);
	Fp3_init(&t);
	Fp3_init(&tmp_Fp3);
	Fp3_set(&n,A);

	mpz_t tmp_mpz,q,e,r,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q);
	mpz_init(e);
	mpz_init(r);
	mpz_init(set_1);
	mpz_init(set_2);
	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(Fp3_legendre(&n)!=-1){
		Fp3_random(&n);
	}

	mpz_pow_ui(q,prime,3);
	mpz_sub_ui(q,q,1);
	mpz_set_ui(e,0);
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}
	Fp3_pow(&y,&n,q);
	mpz_set(r,e);
	mpz_sub_ui(tmp_mpz,q,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);
	Fp3_pow(&x,A,tmp_mpz);
	Fp3_pow(&tmp_Fp3,&x,set_2);
	Fp3_mul(&b,&tmp_Fp3,A);
	Fp3_mul(&x,&x,A);
	int m;
	while(Fp3_cmp_mpz(&b,set_1)==1){
		m=-1;
		Fp3_set(&tmp_Fp3,&b);
		while(Fp3_cmp_mpz(&tmp_Fp3,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp3_pow(&tmp_Fp3,&b,tmp_mpz);
		}
		mpz_sub_ui(tmp_mpz,r,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,prime);
		Fp3_pow(&t,&y,tmp_mpz);
		Fp3_pow(&y,&t,set_2);
		mpz_set_ui(r,m);
		Fp3_mul(&x,&x,&t);
		Fp3_mul(&b,&b,&y);
		Fp3_set(&tmp_Fp3,&b);
	}

	Fp3_set(ANS,&x);

	Fp3_clear(&n);
	Fp3_clear(&y);
	Fp3_clear(&x);
	Fp3_clear(&b);
	Fp3_clear(&t);
	Fp3_clear(&tmp_Fp3);
	mpz_clear(tmp_mpz);
	mpz_clear(q);
	mpz_clear(e);
	mpz_clear(r);
	mpz_clear(set_1);
	mpz_clear(set_2);
}
int Fp3_cmp(struct Fp3 *A,struct Fp3 *B){
	if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0 && Fp_cmp(&A->x2,&B->x2)==0){
		return 0;
	}
	return 1;
}
int Fp3_cmp_mpz(struct Fp3 *A,mpz_t B){
	struct Fp3 tmp;
	Fp3_init(&tmp);
	if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp(&A->x1,&tmp.x1)==0 && Fp_cmp(&A->x2,&tmp.x2)==0){
		Fp3_clear(&tmp);
		return 0;
	}
	Fp3_clear(&tmp);
	return 1;
}
int Fp3_legendre(struct Fp3 *a){
	mpz_t i,cmp;
	mpz_init(i);
	mpz_init(cmp);
	mpz_set_ui(cmp,1);
	struct Fp3 tmp;
	Fp3_init(&tmp);
	mpz_pow_ui(i,prime,3);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp3_pow(&tmp,a,i);

	if((Fp3_cmp_mpz(&tmp,cmp))==0){
		Fp3_clear(&tmp);
		mpz_clear(cmp);
		mpz_clear(i);
		return 1;
	}else{
		Fp3_clear(&tmp);
		mpz_clear(cmp);
		mpz_clear(i);
		return -1;
	}
}
void Fp3_neg(struct Fp3 *ans,struct Fp3 *a){
	struct Fp3 tmp;
	Fp3_init(&tmp);
	Fp_neg(&tmp.x0,&a->x0);
	Fp_neg(&tmp.x1,&a->x1);
	Fp_neg(&tmp.x2,&a->x2);
	Fp3_set(ans,&tmp);
	Fp3_clear(&tmp);
}
void Fp3_frobenius_map(struct Fp3 *ANS, struct Fp3 *A){
	struct Fp pm1d3,pm1d3p2,tmp,set_c1;
	struct Fp3 t_ans;
	Fp_init(&pm1d3);
	Fp_init(&pm1d3p2);
	Fp_init(&tmp);
	Fp_init(&set_c1);
	Fp3_init(&t_ans);

	mpz_tdiv_q_ui(pm1d3.x0,prime,3);

	Fp_set_ui(&set_c1,c1);
	Fp_pow(&pm1d3,&set_c1,pm1d3.x0);
	Fp_mul(&pm1d3p2,&pm1d3,&pm1d3);

	Fp_set(&t_ans.x0,&A->x0);
	Fp_mul(&t_ans.x1,&A->x1,&pm1d3);
	Fp_mul(&t_ans.x2,&A->x2,&pm1d3p2);

	Fp3_set(ANS,&t_ans);

	Fp_clear(&pm1d3);
	Fp_clear(&pm1d3p2);
	Fp_clear(&tmp);
	Fp_clear(&set_c1);
	Fp3_clear(&t_ans);

}

//-----------------------------------------------------------------------------------------

void Fp6_init(struct Fp6 *A){
	Fp3_init(&A->x0);
	Fp3_init(&A->x1);
}
void Fp6_set(struct Fp6 *ANS,struct Fp6 *A){
	Fp3_set(&ANS->x0,&A->x0);
	Fp3_set(&ANS->x1,&A->x1);
}
void Fp6_set_ui(struct Fp6 *A,signed long int B){
	Fp3_set_ui(&A->x0,B);
	Fp3_set_ui(&A->x1,B);
}
void Fp6_random(struct Fp6 *A){
	Fp3_random(&A->x0);
	Fp3_random(&A->x1);
}
void Fp6_clear(struct Fp6 *A){
	Fp3_clear(&A->x0);
	Fp3_clear(&A->x1);
}
void Fp6_printf(struct Fp6 *A){
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0,A->x0.x1.x0,A->x0.x2.x0,A->x1.x0.x0,A->x1.x1.x0,A->x1.x2.x0);
}
void Fp6_add(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp3_add(&tmp.x0,&A->x0,&B->x0);
	Fp3_add(&tmp.x1,&A->x1,&B->x1);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_add_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp3_add_ui(&tmp.x0,&A->x0,B);
	Fp3_add_ui(&tmp.x1,&A->x1,B);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_sub(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp3_sub(&tmp.x0,&A->x0,&B->x0);
	Fp3_sub(&tmp.x1,&A->x1,&B->x1);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_mul(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	//x^2-v=0
	struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&tmp5);
	Fp3_init(&tmp6);

	struct Fp6 t_ans;
	Fp6_init(&t_ans);

	Fp3_mul(&tmp1,&A->x0,&B->x0);//a*c
	Fp3_mul(&tmp2,&A->x1,&B->x1);//b*d
	Fp3_mul_xi(&tmp3,&tmp2);//b*d*v
	Fp3_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
	Fp3_add(&tmp4,&A->x0,&A->x1);//a+b
	Fp3_add(&tmp5,&B->x0,&B->x1);//c+d
	Fp3_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
	Fp3_sub(&t_ans.x1,&tmp6,&tmp1);
	Fp3_sub(&t_ans.x1,&t_ans.x1,&tmp2);

	Fp6_set(ANS,&t_ans);

	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&tmp5);
	Fp3_clear(&tmp6);
	Fp6_clear(&t_ans);
}
void Fp6_mul_Fp(struct Fp6 *ANS,struct Fp6 *A,struct Fp *B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp3_mul_Fp(&tmp.x0,&A->x0,B);
	Fp3_mul_Fp(&tmp.x1,&A->x1,B);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_mul_v(struct Fp6 *ANS,struct Fp6 *A){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp3_mul_xi(&tmp.x0,&A->x1);
	Fp3_set(&tmp.x1,&A->x0);

	Fp6_set(ANS,&tmp);
	Fp6_clear(&tmp);
}
void Fp6_mul_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp3_mul_ui(&tmp.x0,&A->x0,B);
	Fp3_mul_ui(&tmp.x1,&A->x1,B);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_invert(struct Fp6 *ANS, struct Fp6 *A){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	// tmp=A^(q^6)=(x0,-x1)
	Fp3_set(&tmp.x0,&A->x0);
	Fp3_neg(&tmp.x1,&A->x1);

	struct Fp3 c,a,b;
	Fp3_init(&c);
	Fp3_init(&a);
	Fp3_init(&b);

	Fp3_mul(&a,&A->x0,&A->x0); // a=x0^2
	Fp3_mul(&b,&A->x1,&A->x1); // b=x1^2
	Fp3_mul_xi(&b,&b); // b=x1^2*v
	Fp3_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q

	Fp3_invert(&c,&c);

	// ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp6-element and tmp is a vector A Fp6
	Fp3_mul(&tmp.x0,&tmp.x0,&c);
	Fp3_mul(&tmp.x1,&tmp.x1,&c);

	Fp6_set(ANS,&tmp);

	Fp3_clear(&c);
	Fp3_clear(&a);
	Fp3_clear(&b);
	Fp6_clear(&tmp);
}
void Fp6_div(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	struct Fp6 tmp,t_ans;
	Fp6_init(&tmp);
	Fp6_init(&t_ans);

	Fp6_invert(&tmp,B);
	Fp6_mul(&t_ans,A,&tmp);

	Fp6_set(ANS,&t_ans);

	Fp6_clear(&tmp);
	Fp6_clear(&t_ans);
}
void Fp6_pow(struct Fp6 *ANS,struct Fp6 *A,mpz_t B){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(B,2);
	//printf("r= %d\n",r);

	struct Fp6 answer_tmp;
	Fp6_init(&answer_tmp);
	Fp6_set(&answer_tmp,A);

	struct Fp6 in_tmp;
	Fp6_init(&in_tmp);
	Fp6_set(&in_tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			Fp6_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
			Fp6_mul(&answer_tmp,&answer_tmp,&in_tmp);//*a
		}else{
			Fp6_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
		}
	}

	Fp6_set(ANS,&answer_tmp);

	Fp6_clear(&answer_tmp);
	Fp6_clear(&in_tmp);
}
void Fp6_sqrt(struct Fp6 *ANS,struct Fp6 *A){
	struct Fp6 n,y,x,b,t,tmp_Fp6;
	Fp6_init(&n);
	Fp6_init(&y);
	Fp6_init(&x);
	Fp6_init(&b);
	Fp6_init(&t);
	Fp6_init(&tmp_Fp6);
	Fp6_set(&n,A);

	mpz_t tmp_mpz,q,e,r,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q);
	mpz_init(e);
	mpz_init(r);
	mpz_init(set_1);
	mpz_init(set_2);
	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(Fp6_legendre(&n)!=-1){
		Fp6_random(&n);
	}
	mpz_pow_ui(q,prime,6);
	mpz_sub_ui(q,q,1);
	mpz_set_ui(e,0);
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}
	Fp6_pow(&y,&n,q);

	mpz_set(r,e);

	mpz_sub_ui(tmp_mpz,q,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);

	Fp6_pow(&x,A,tmp_mpz);
	Fp6_pow(&tmp_Fp6,&x,set_2);
	Fp6_mul(&b,&tmp_Fp6,A);
	Fp6_mul(&x,&x,A);

	int m;

	while(Fp6_cmp_mpz(&b,set_1)==1){
		m=-1;
		Fp6_set(&tmp_Fp6,&b);
		while(Fp6_cmp_mpz(&tmp_Fp6,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp6_pow(&tmp_Fp6,&b,tmp_mpz);
		}
		mpz_sub_ui(tmp_mpz,r,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,prime);
		// gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
		Fp6_pow(&t,&y,tmp_mpz);
		Fp6_pow(&y,&t,set_2);
		// gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
		mpz_set_ui(r,m);
		Fp6_mul(&x,&x,&t);
		Fp6_mul(&b,&b,&y);
		Fp6_set(&tmp_Fp6,&b);
	}

	Fp6_set(ANS,&x);

	Fp6_clear(&n);
	Fp6_clear(&y);
	Fp6_clear(&x);
	Fp6_clear(&b);
	Fp6_clear(&t);
	Fp6_clear(&tmp_Fp6);
	mpz_clear(tmp_mpz);
	mpz_clear(q);
	mpz_clear(e);
	mpz_clear(r);
	mpz_clear(set_1);
	mpz_clear(set_2);
}
int Fp6_legendre(struct Fp6 *A){
	mpz_t i,cmp;
	struct Fp6 tmp;
	Fp6_init(&tmp);
	mpz_init(i);
	mpz_init(cmp);
	mpz_set_ui(cmp,1);
	mpz_pow_ui(i,prime,6);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp6_pow(&tmp,A,i);

	if((Fp6_cmp_mpz(&tmp,cmp))==0){
		Fp6_clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return 1;
	}else{
		Fp6_clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return -1;
	}
}
int Fp6_cmp(struct Fp6 *A,struct Fp6 *B){
	if(Fp3_cmp(&A->x0,&B->x0)==0 && Fp3_cmp(&A->x1,&B->x1)==0){
		return 0;
	}
	return 1;
}
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t B){
	struct Fp6 tmp;
	Fp6_init(&tmp);
	if(Fp3_cmp_mpz(&A->x0,B)==0 && Fp3_cmp(&A->x1,&tmp.x1)==0){
		Fp6_clear(&tmp);
		return 0;
	}
	Fp6_clear(&tmp);
	return 1;
}
void Fp6_neg(struct Fp6 *ANS,struct Fp6 *A){
	struct Fp6 tmp;
	Fp6_init(&tmp);
	Fp3_neg(&tmp.x0,&A->x0);
	Fp3_neg(&tmp.x1,&A->x1);
	Fp6_set(ANS,&tmp);
	Fp6_clear(&tmp);
}
void Fp6_frobenius_map(struct Fp6 *ANS, struct Fp6 *A){
	struct Fp pm1d3,pm1d3p2,pm1d6,tmp,set_c1;
	struct Fp6 t_ans;
	Fp_init(&pm1d3);
	Fp_init(&pm1d3p2);
	Fp_init(&pm1d6);
	Fp_init(&tmp);
	Fp_init(&set_c1);
	Fp6_init(&t_ans);

	mpz_tdiv_q_ui(pm1d6.x0,prime,6);

	Fp_set_ui(&set_c1,c1);
	Fp_pow(&pm1d6,&set_c1,pm1d6.x0);
	Fp_mul(&pm1d3,&pm1d6,&pm1d6);
	Fp_mul(&pm1d3p2,&pm1d3,&pm1d3);

	Fp_set(&t_ans.x0.x0,&A->x0.x0);
	Fp_mul(&t_ans.x0.x1,&A->x0.x1,&pm1d3);
	Fp_mul(&t_ans.x0.x2,&A->x0.x2,&pm1d3p2);

	Fp_mul(&t_ans.x1.x0,&A->x1.x0,&pm1d6);
	Fp_mul(&t_ans.x1.x1,&A->x1.x1,&pm1d6);
	Fp_mul(&t_ans.x1.x1,&t_ans.x1.x1,&pm1d3);
	Fp_mul(&t_ans.x1.x2,&A->x1.x2,&pm1d6);
	Fp_mul(&t_ans.x1.x2,&t_ans.x1.x2,&pm1d3p2);

	Fp6_set(ANS,&t_ans);

	Fp_clear(&pm1d3);
	Fp_clear(&pm1d3p2);
	Fp_clear(&pm1d6);
	Fp_clear(&tmp);
	Fp_clear(&set_c1);
	Fp6_clear(&t_ans);
}

//-----------------------------------------------------------------------------------------

void Fp18_init(struct Fp18 *A){
	Fp6_init(&A->x0);
	Fp6_init(&A->x1);
	Fp6_init(&A->x2);
}
void Fp18_set(struct Fp18 *ANS,struct Fp18 *A){
	Fp6_set(&ANS->x0,&A->x0);
	Fp6_set(&ANS->x1,&A->x1);
	Fp6_set(&ANS->x2,&A->x2);
}
void Fp18_set_ui(struct Fp18 *A,signed long int B){
	Fp6_set_ui(&A->x0,B);
	Fp6_set_ui(&A->x1,B);
	Fp6_set_ui(&A->x2,B);
}
void Fp18_random(struct Fp18 *A){
	Fp6_random(&A->x0);
	Fp6_random(&A->x1);
	Fp6_random(&A->x2);
}
void Fp18_clear(struct Fp18 *A){
	Fp6_clear(&A->x0);
	Fp6_clear(&A->x1);
	Fp6_clear(&A->x2);
}
void Fp18_printf(struct Fp18 *A){
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0.x0,A->x0.x0.x1.x0,A->x0.x0.x2.x0,A->x0.x1.x0.x0,A->x0.x1.x1.x0,A->x0.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x1.x0.x0.x0,A->x1.x0.x1.x0,A->x1.x0.x2.x0,A->x1.x1.x0.x0,A->x1.x1.x1.x0,A->x1.x1.x2.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n\n",A->x2.x0.x0.x0,A->x2.x0.x1.x0,A->x2.x0.x2.x0,A->x2.x1.x0.x0,A->x2.x1.x1.x0,A->x2.x1.x2.x0);
}
void Fp18_add(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp18 tmp;
	Fp18_init(&tmp);

	Fp6_add(&tmp.x0,&A->x0,&B->x0);
	Fp6_add(&tmp.x1,&A->x1,&B->x1);
	Fp6_add(&tmp.x2,&A->x2,&B->x2);

	Fp18_set(ANS,&tmp);

	Fp18_clear(&tmp);
}
void Fp18_add_ui(struct Fp18 *ANS,struct Fp18 *A,unsigned long int B){
	struct Fp18 tmp;
	Fp18_init(&tmp);

	Fp6_add_ui(&tmp.x0,&A->x0,B);
	Fp6_add_ui(&tmp.x1,&A->x1,B);
	Fp6_add_ui(&tmp.x2,&A->x2,B);

	Fp18_set(ANS,&tmp);

	Fp18_clear(&tmp);
}
void Fp18_sub(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp18 tmp;
	Fp18_init(&tmp);

	Fp6_sub(&tmp.x0,&A->x0,&B->x0);
	Fp6_sub(&tmp.x1,&A->x1,&B->x1);
	Fp6_sub(&tmp.x2,&A->x2,&B->x2);

	Fp18_set(ANS,&tmp);

	Fp18_clear(&tmp);
}
void Fp18_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	//(x0,x1,x2)*(y0,y1,y2)=(x0y0+xi((x1+x2)(y1+y2)-x1y1-x2y2),xix2y2+(x0+x1)(y0+y1)-x0y0-x1y1,x1y1+(x0+x2)(y0+y2)-x0y0-x2y2)
	struct Fp6 tmp00,tmp11,tmp22,tmpx01,tmpx12,tmpx20,tmpy01,tmpy12,tmpy20,t0,t1,t2,tmp;
	struct Fp18 t_ans;
	Fp6_init(&tmp00);
	Fp6_init(&tmp11);
	Fp6_init(&tmp22);
	Fp6_init(&tmpx01);
	Fp6_init(&tmpx12);
	Fp6_init(&tmpx20);
	Fp6_init(&tmpy01);
	Fp6_init(&tmpy12);
	Fp6_init(&tmpy20);
	Fp6_init(&t0);
	Fp6_init(&t1);
	Fp6_init(&t2);
	Fp6_init(&tmp);
	Fp18_init(&t_ans);

	Fp6_mul(&tmp00,&A->x0,&B->x0);//x0*y0
	Fp6_mul(&tmp11,&A->x1,&B->x1);//x1*y1
	Fp6_mul(&tmp22,&A->x2,&B->x2);//x2*y2

	Fp6_add(&tmpx01,&A->x0,&A->x1);//x0+x1
	Fp6_add(&tmpx12,&A->x1,&A->x2);//x1+x2
	Fp6_add(&tmpx20,&A->x0,&A->x2);//x2+x0
	Fp6_add(&tmpy01,&B->x0,&B->x1);//y0+y1
	Fp6_add(&tmpy12,&B->x1,&B->x2);//y1+y2
	Fp6_add(&tmpy20,&B->x0,&B->x2);//y2+y0

	Fp6_mul(&t0,&tmpx01,&tmpy01);//(x0+x1)(y0+y1)
	Fp6_mul(&t1,&tmpx12,&tmpy12);//(x1+x2)(y1+y2)
	Fp6_mul(&t2,&tmpx20,&tmpy20);//(x2+x0)(y2+y0)

	Fp6_sub(&t1,&t1,&tmp11);
	Fp6_sub(&t1,&t1,&tmp22);//(x1+x2)(y1+y2)-x1y1-x2y2
	Fp6_mul_v(&tmp,&t1);
	Fp6_add(&t_ans.x0,&tmp00,&tmp);

	Fp6_sub(&t0,&t0,&tmp00);
	Fp6_sub(&t0,&t0,&tmp11);
	Fp6_mul_v(&tmp,&tmp22);
	Fp6_add(&t_ans.x1,&tmp,&t0);

	Fp6_sub(&t2,&t2,&tmp00);
	Fp6_sub(&t2,&t2,&tmp22);
	Fp6_add(&t_ans.x2,&tmp11,&t2);

	Fp18_set(ANS,&t_ans);

	Fp6_clear(&tmp00);
	Fp6_clear(&tmp11);
	Fp6_clear(&tmp22);
	Fp6_clear(&tmpx01);
	Fp6_clear(&tmpx12);
	Fp6_clear(&tmpx20);
	Fp6_clear(&tmpy01);
	Fp6_clear(&tmpy12);
	Fp6_clear(&tmpy20);
	Fp6_clear(&t0);
	Fp6_clear(&t1);
	Fp6_clear(&t2);
	Fp6_clear(&tmp);
	Fp18_clear(&t_ans);
}
void Fp18_mul_Fp(struct Fp18 *ANS,struct Fp18 *A,struct Fp *B){
	struct Fp18 tmp;
	Fp18_init(&tmp);

	Fp6_mul_Fp(&tmp.x0,&A->x0,B);
	Fp6_mul_Fp(&tmp.x1,&A->x1,B);
	Fp6_mul_Fp(&tmp.x2,&A->x2,B);

	Fp18_set(ANS,&tmp);

	Fp18_clear(&tmp);
}
void Fp18_mul_ui(struct Fp18 *ANS,struct Fp18 *A,unsigned long int B){
	struct Fp18 tmp;
	Fp18_init(&tmp);

	Fp6_mul_ui(&tmp.x0,&A->x0,B);
	Fp6_mul_ui(&tmp.x1,&A->x1,B);
	Fp6_mul_ui(&tmp.x2,&A->x2,B);

	Fp18_set(ANS,&tmp);
	Fp18_clear(&tmp);
}
void Fp18_invert(struct Fp18 *ANS, struct Fp18 *A){
	struct Fp18 t_ans;
	Fp18_init(&t_ans);

	struct Fp6 T0,T1,t0,t1,t2,t3;
	Fp6_init(&T0);
	Fp6_init(&T1);
	Fp6_init(&t0);
	Fp6_init(&t1);
	Fp6_init(&t2);
	Fp6_init(&t3);

	// An optimized version of Grewal's Algo. 3   (a,b,c)
	Fp6_mul(&T0,&A->x0,&A->x0);
	Fp6_mul_v(&t0,&A->x1);

	Fp6_mul(&T1,&t0,&A->x2);
	Fp6_sub(&t1,&T0,&T1); // t1=(a^2-bci) mod q

	Fp6_mul(&T0,&A->x2,&A->x2);
	Fp6_mul_v(&T0,&T0);
	Fp6_mul(&T1,&A->x0,&A->x1);
	Fp6_sub(&t2,&T0,&T1); // t2=(c^2i-ab) mod q

	Fp6_mul(&T0,&A->x1,&A->x1);
	Fp6_mul(&T1,&A->x0,&A->x2);
	Fp6_sub(&t3,&T0,&T1); // t3=(b^2-ac) mod q

	Fp6_mul(&T0,&t0,&t3);
	Fp6_mul(&T1,&A->x0,&t1);
	Fp6_add(&T0,&T0,&T1); // T0={bi(b^2-ac)+a(a^2-bci)} mod q

	Fp6_mul_v(&t0,&A->x2);
	Fp6_mul(&T1,&t0,&t2);
	Fp6_add(&t0,&T0,&T1); // t0={ci(c^2i-ab)+{bi(b^2-ac)+a(a^2-bci)}} mod q .0

	Fp6_invert(&t0,&t0);

	Fp6_mul(&t_ans.x0,&t1,&t0);
	Fp6_mul(&t_ans.x1,&t2,&t0);
	Fp6_mul(&t_ans.x2,&t3,&t0);

	Fp18_set(ANS,&t_ans);

	Fp18_clear(&t_ans);
	Fp6_clear(&T0);
	Fp6_clear(&T1);
	Fp6_clear(&t0);
	Fp6_clear(&t1);
	Fp6_clear(&t2);
	Fp6_clear(&t3);
}
void Fp18_div(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp18 tmp,t_ans;
	Fp18_init(&tmp);
	Fp18_init(&t_ans);

	Fp18_invert(&tmp,B);
	Fp18_mul(&t_ans,A,&tmp);

	Fp18_set(ANS,&t_ans);

	Fp18_clear(&tmp);
	Fp18_clear(&t_ans);
}
void Fp18_pow(struct Fp18 *ANS,struct Fp18 *A,mpz_t B){
	int i;
	int r;//bit数
	r= (int)mpz_sizeinbase(B,2);
	//printf("r= %d\n",r);

	struct Fp18 answer_tmp;
	Fp18_init(&answer_tmp);
	Fp18_set(&answer_tmp,A);

	struct Fp18 in_tmp;
	Fp18_init(&in_tmp);
	Fp18_set(&in_tmp,A);

	for(i=r-2;i>=0;i--){
		if(mpz_tstbit(B,i)==1){
			Fp18_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
			Fp18_mul(&answer_tmp,&answer_tmp,&in_tmp);//*a
		}else{
			Fp18_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
		}
	}

	Fp18_set(ANS,&answer_tmp);

	Fp18_clear(&answer_tmp);
	Fp18_clear(&in_tmp);
}
void Fp18_sqrt(struct Fp18 *ANS,struct Fp18 *A){
	struct Fp18 n,y,x,b,t,tmp_Fp18;
	Fp18_init(&n);
	Fp18_init(&y);
	Fp18_init(&x);
	Fp18_init(&b);
	Fp18_init(&t);
	Fp18_init(&tmp_Fp18);
	Fp18_set(&n,A);

	mpz_t tmp_mpz,q,e,r,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q);
	mpz_init(e);
	mpz_init(r);
	mpz_init(set_1);
	mpz_init(set_2);
	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(Fp18_legendre(&n)!=-1){
		Fp18_random(&n);
	}

	mpz_pow_ui(q,prime,18);
	mpz_sub_ui(q,q,1);
	mpz_set_ui(e,0);

	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	Fp18_pow(&y,&n,q);
	mpz_set(r,e);
	mpz_sub_ui(tmp_mpz,q,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);
	Fp18_pow(&x,A,tmp_mpz);
	Fp18_pow(&tmp_Fp18,&x,set_2);
	Fp18_mul(&b,&tmp_Fp18,A);
	Fp18_mul(&x,&x,A);

	int m;

	while(Fp18_cmp_mpz(&b,set_1)==1){
		m=-1;
		Fp18_set(&tmp_Fp18,&b);
		while(Fp18_cmp_mpz(&tmp_Fp18,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp18_pow(&tmp_Fp18,&b,tmp_mpz);
		}
		mpz_sub_ui(tmp_mpz,r,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,prime);
		Fp18_pow(&t,&y,tmp_mpz);
		Fp18_pow(&y,&t,set_2);
		mpz_set_ui(r,m);
		Fp18_mul(&x,&x,&t);
		Fp18_mul(&b,&b,&y);
		Fp18_set(&tmp_Fp18,&b);
	}

	Fp18_set(ANS,&x);

	Fp18_clear(&n);
	Fp18_clear(&y);
	Fp18_clear(&x);
	Fp18_clear(&b);
	Fp18_clear(&t);
	Fp18_clear(&tmp_Fp18);
	mpz_clear(tmp_mpz);
	mpz_clear(q);
	mpz_clear(e);
	mpz_clear(r);
	mpz_clear(set_1);
	mpz_clear(set_2);
}
int Fp18_cmp(struct Fp18 *A,struct Fp18 *B){
	if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0 && Fp6_cmp(&A->x2,&B->x2)==0){
		return 0;
	}
	return 1;
}
int Fp18_cmp_mpz(struct Fp18 *A,mpz_t B){
	struct Fp18 tmp;
	Fp18_init(&tmp);
	if(Fp6_cmp_mpz(&A->x0,B)==0 && Fp6_cmp(&A->x1,&tmp.x1)==0 && Fp6_cmp(&A->x2,&tmp.x2)==0){
		Fp18_clear(&tmp);
		return 0;
	}
	Fp18_clear(&tmp);
	return 1;
}
int Fp18_legendre(struct Fp18 *A){
	mpz_t cmp,pow;
	struct Fp18 tmp;
	struct Fp18 frobenius;
	int i;

	mpz_init(cmp);
	mpz_init(pow);
	Fp18_init(&tmp);
	Fp18_init(&frobenius);

	mpz_set_ui(cmp,1);
	Fp18_frobenius_map(&tmp,A,9);
	Fp18_mul(&tmp,&tmp,A);

	Fp18_frobenius_map(&frobenius,&tmp,1);
	Fp18_mul(&tmp,&tmp,&frobenius);
	for(i=1;i<8;i++){
		Fp18_frobenius_map(&frobenius,&frobenius,1);
		Fp18_mul(&tmp,&tmp,&frobenius);
	}

	mpz_sub_ui(pow,prime,1);
	mpz_div_ui(pow,pow,2);
	Fp18_pow(&tmp,&tmp,pow);

	if((Fp18_cmp_mpz(&tmp,cmp))==0){
		Fp18_clear(&tmp);
		mpz_clear(cmp);
		return 1;
	}else{
		Fp18_clear(&tmp);
		mpz_clear(cmp);
		return -1;
	}
}
void Fp18_neg(struct Fp18 *ANS,struct Fp18 *A){
	struct Fp18 tmp;
	Fp18_init(&tmp);
	Fp6_neg(&tmp.x0,&A->x0);
	Fp6_neg(&tmp.x1,&A->x1);
	Fp6_neg(&tmp.x2,&A->x2);
	Fp18_set(ANS,&tmp);
	Fp18_clear(&tmp);
}
void Fp18_frobenius_map(struct Fp18 *ANS, struct Fp18 *A,int times){
	struct Fp i,ip2,v,w,wp2,set_c1;
	struct Fp iv,ip2v;
	struct Fp18 t_ans,t_Fp18;
	Fp_init(&i);
	Fp_init(&ip2);
	Fp_init(&v);
	Fp_init(&w);
	Fp_init(&wp2);
	Fp_init(&set_c1);
	Fp18_init(&t_ans);
	Fp18_init(&t_Fp18);

	Fp_init(&iv);
	Fp_init(&ip2v);

	struct Fp Switch,tmp,tmp2;
	Fp_init(&Switch);
	Fp_init(&tmp);
	Fp_init(&tmp2);

	mpz_tdiv_q_ui(v.x0,prime,6);

	mpz_tdiv_r_ui(Switch.x0,v.x0,3);

	Fp_set_ui(&set_c1,c1);
	Fp_pow(&v,&set_c1,v.x0);//c1^p-1/6
	Fp_mul(&i,&v,&v);//c1^p-1/3
	Fp_mul(&ip2,&i,&i);//(c1^p-1/3)^2

	Fp_mul(&iv,&i,&v);
	Fp_mul(&ip2v,&ip2,&v);
	int loop;
	Fp18_set(&t_Fp18,A);
	for(loop=0;loop<times;loop++){
		//1
		Fp_set(&t_ans.x0.x0.x0,&t_Fp18.x0.x0.x0);
		Fp_mul(&t_ans.x0.x0.x1,&t_Fp18.x0.x0.x1,&i);
		Fp_mul(&t_ans.x0.x0.x2,&t_Fp18.x0.x0.x2,&ip2);

		Fp_mul(&t_ans.x0.x1.x0,&t_Fp18.x0.x1.x0,&v);
		Fp_mul(&t_ans.x0.x1.x1,&t_Fp18.x0.x1.x1,&iv);
		Fp_mul(&t_ans.x0.x1.x2,&t_Fp18.x0.x1.x2,&ip2v);
		switch (mpz_get_ui(Switch.x0)){
			case 0:
			mpz_tdiv_q_ui(w.x0,prime,18);
			Fp_pow(&w,&set_c1,w.x0);//c1^p-1/18
			Fp_mul(&wp2,&w,&w);//(c1^p-1/3)^2

			//w
			Fp_mul(&t_ans.x1.x0.x0,&t_Fp18.x1.x0.x0,&w);
			Fp_mul(&t_ans.x1.x0.x1,&t_Fp18.x1.x0.x1,&w);
			Fp_mul(&t_ans.x1.x0.x1,&t_ans.x1.x0.x1,&i);
			Fp_mul(&t_ans.x1.x0.x2,&t_Fp18.x1.x0.x2,&w);
			Fp_mul(&t_ans.x1.x0.x2,&t_ans.x1.x0.x2,&ip2);

			Fp_mul(&t_ans.x1.x1.x0,&t_Fp18.x1.x1.x0,&w);
			Fp_mul(&t_ans.x1.x1.x0,&t_ans.x1.x1.x0,&v);
			Fp_mul(&t_ans.x1.x1.x1,&t_Fp18.x1.x1.x1,&w);
			Fp_mul(&t_ans.x1.x1.x1,&t_ans.x1.x1.x1,&iv);
			Fp_mul(&t_ans.x1.x1.x2,&t_Fp18.x1.x1.x2,&w);
			Fp_mul(&t_ans.x1.x1.x2,&t_ans.x1.x1.x2,&ip2v);

			//w^2
			Fp_mul(&t_ans.x2.x0.x0,&t_Fp18.x2.x0.x0,&wp2);
			Fp_mul(&t_ans.x2.x0.x1,&t_Fp18.x2.x0.x1,&wp2);
			Fp_mul(&t_ans.x2.x0.x1,&t_ans.x2.x0.x1,&i);
			Fp_mul(&t_ans.x2.x0.x2,&t_Fp18.x2.x0.x2,&wp2);
			Fp_mul(&t_ans.x2.x0.x2,&t_ans.x2.x0.x2,&ip2);

			Fp_mul(&t_ans.x2.x1.x0,&t_Fp18.x2.x1.x0,&wp2);
			Fp_mul(&t_ans.x2.x1.x0,&t_ans.x2.x1.x0,&v);
			Fp_mul(&t_ans.x2.x1.x1,&t_Fp18.x2.x1.x1,&wp2);
			Fp_mul(&t_ans.x2.x1.x1,&t_ans.x2.x1.x1,&iv);
			Fp_mul(&t_ans.x2.x1.x2,&t_Fp18.x2.x1.x2,&wp2);
			Fp_mul(&t_ans.x2.x1.x2,&t_ans.x2.x1.x2,&ip2v);
			printf("frobenius map error\n");
			break;

			case 1:
			mpz_sub_ui(w.x0,prime,7);
			mpz_tdiv_q_ui(w.x0,w.x0,18);
			Fp_pow(&w,&set_c1,w.x0);//c1^p-1/18
			Fp_mul(&wp2,&w,&w);//(c1^p-1/3)^2
			Fp_mul(&tmp,&w,&set_c1);

			//w
			Fp_mul(&t_ans.x1.x0.x1,&t_Fp18.x1.x0.x0,&w);
			Fp_mul(&t_ans.x1.x0.x2,&t_Fp18.x1.x0.x1,&w);
			Fp_mul(&t_ans.x1.x0.x2,&t_ans.x1.x0.x2,&i);
			Fp_mul(&t_ans.x1.x0.x0,&t_Fp18.x1.x0.x2,&tmp);
			Fp_mul(&t_ans.x1.x0.x0,&t_ans.x1.x0.x0,&ip2);

			Fp_mul(&t_ans.x1.x1.x1,&t_Fp18.x1.x1.x0,&w);
			Fp_mul(&t_ans.x1.x1.x1,&t_ans.x1.x1.x1,&v);
			Fp_mul(&t_ans.x1.x1.x2,&t_Fp18.x1.x1.x1,&w);
			Fp_mul(&t_ans.x1.x1.x2,&t_ans.x1.x1.x2,&iv);
			Fp_mul(&t_ans.x1.x1.x0,&t_Fp18.x1.x1.x2,&tmp);
			Fp_mul(&t_ans.x1.x1.x0,&t_ans.x1.x1.x0,&ip2v);

			Fp_mul(&tmp,&wp2,&set_c1);
			//w^2
			Fp_mul(&t_ans.x2.x0.x2,&t_Fp18.x2.x0.x0,&wp2);
			Fp_mul(&t_ans.x2.x0.x0,&t_Fp18.x2.x0.x1,&tmp);
			Fp_mul(&t_ans.x2.x0.x0,&t_ans.x2.x0.x0,&i);
			Fp_mul(&t_ans.x2.x0.x1,&t_Fp18.x2.x0.x2,&tmp);
			Fp_mul(&t_ans.x2.x0.x1,&t_ans.x2.x0.x1,&ip2);

			Fp_mul(&t_ans.x2.x1.x2,&t_Fp18.x2.x1.x0,&wp2);
			Fp_mul(&t_ans.x2.x1.x2,&t_ans.x2.x1.x2,&v);
			Fp_mul(&t_ans.x2.x1.x0,&t_Fp18.x2.x1.x1,&tmp);
			Fp_mul(&t_ans.x2.x1.x0,&t_ans.x2.x1.x0,&iv);
			Fp_mul(&t_ans.x2.x1.x1,&t_Fp18.x2.x1.x2,&tmp);
			Fp_mul(&t_ans.x2.x1.x1,&t_ans.x2.x1.x1,&ip2v);
			// printf("1\n");
			break;

			case 2:
			mpz_sub_ui(w.x0,prime,13);
			mpz_tdiv_q_ui(w.x0,w.x0,18);
			Fp_pow(&w,&set_c1,w.x0);//c1^p-1/18
			Fp_mul(&wp2,&w,&w);//(c1^p-1/3)^2
			Fp_mul(&tmp,&w,&set_c1);

			//w
			Fp_mul(&t_ans.x1.x0.x2,&t_Fp18.x1.x0.x0,&w);
			Fp_mul(&t_ans.x1.x0.x0,&t_Fp18.x1.x0.x1,&tmp);
			Fp_mul(&t_ans.x1.x0.x0,&t_ans.x1.x0.x0,&i);
			Fp_mul(&t_ans.x1.x0.x1,&t_Fp18.x1.x0.x2,&tmp);
			Fp_mul(&t_ans.x1.x0.x1,&t_ans.x1.x0.x1,&ip2);

			Fp_mul(&t_ans.x1.x1.x2,&t_Fp18.x1.x1.x0,&w);
			Fp_mul(&t_ans.x1.x1.x2,&t_ans.x1.x1.x2,&v);
			Fp_mul(&t_ans.x1.x1.x0,&t_Fp18.x1.x1.x1,&tmp);
			Fp_mul(&t_ans.x1.x1.x0,&t_ans.x1.x1.x0,&iv);
			Fp_mul(&t_ans.x1.x1.x1,&t_Fp18.x1.x1.x2,&tmp);
			Fp_mul(&t_ans.x1.x1.x1,&t_ans.x1.x1.x1,&ip2v);

			Fp_mul(&tmp,&wp2,&set_c1);
			Fp_mul(&tmp2,&tmp,&set_c1);
			//w^2
			Fp_mul(&t_ans.x2.x0.x1,&t_Fp18.x2.x0.x0,&tmp);
			Fp_mul(&t_ans.x2.x0.x2,&t_Fp18.x2.x0.x1,&tmp);
			Fp_mul(&t_ans.x2.x0.x2,&t_ans.x2.x0.x2,&i);
			Fp_mul(&t_ans.x2.x0.x0,&t_Fp18.x2.x0.x2,&tmp2);
			Fp_mul(&t_ans.x2.x0.x0,&t_ans.x2.x0.x0,&ip2);

			Fp_mul(&t_ans.x2.x1.x1,&t_Fp18.x2.x1.x0,&tmp);
			Fp_mul(&t_ans.x2.x1.x1,&t_ans.x2.x1.x1,&v);
			Fp_mul(&t_ans.x2.x1.x2,&t_Fp18.x2.x1.x1,&tmp);
			Fp_mul(&t_ans.x2.x1.x2,&t_ans.x2.x1.x2,&iv);
			Fp_mul(&t_ans.x2.x1.x0,&t_Fp18.x2.x1.x2,&tmp2);
			Fp_mul(&t_ans.x2.x1.x0,&t_ans.x2.x1.x0,&ip2v);
			// printf("2\n");
			break;
		}
		Fp18_set(&t_Fp18,&t_ans);
	}


	Fp18_set(ANS,&t_ans);

	Fp_clear(&i);
	Fp_clear(&ip2);
	Fp_clear(&v);
	Fp_clear(&w);
	Fp_clear(&wp2);
	Fp_clear(&set_c1);
	Fp18_clear(&t_ans);
	Fp18_clear(&t_Fp18);
	Fp_clear(&iv);
	Fp_clear(&ip2v);
	Fp_clear(&tmp);
	Fp_clear(&tmp2);
}

#endif //Finite_Field_C_
