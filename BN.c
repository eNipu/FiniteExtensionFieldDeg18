//gcc -I/usr/local/include -L/usr/local/lib -Wall -O3 -o BN.out BN.c -lgmp


#include"BN.h"

int main(void){
	mpz_init(X);
	// mpz_set_str(X,"45756572214084355210000000873",10);
	// mpz_set_str(X,"4575657221408435521",10);
	mpz_set_str(X,"100997",10);
	// mpz_set_str(X,"41",10);
	// mpz_set_str(X,"1",10);
	mpz_inits(p,t,r);
	mpz_init(b);
	struct EFp P;
	EFp_init(&P);

	set_BN_parameter(&P);
	gmp_printf("p=%Zd\n",p);
	gmp_printf("t=%Zd\n",t);
	gmp_printf("r=%Zd\n",r);

	printf("p = %dbit\n",(int)mpz_sizeinbase(r,2));

	gmp_printf("y^2=x^3+%Zd\n",b);

	// struct EFp12 AA,BB;
	// EFp12_init(&AA);
	// EFp12_init(&BB);
	//
	// EFp12_random_set(&AA);
	// EFp12_printf(&AA);
	// EFp12_SCM(&AA,&AA,r);
	// EFp12_printf(&AA);

	check_Pairing();


	EFp_clear(&P);
	mpz_clear(p);
	mpz_clear(t);
	mpz_clear(r);
	mpz_clear(b);
	return 0;
}
//-----------------------------------------------------------------------------------------

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
	mpz_mod(A->x0,A->x0,p);
}
void Fp_clear(struct Fp *A){
	mpz_clear(A->x0);
}
void Fp_printf(struct Fp *A){
	gmp_printf("%Zd\n",A->x0);
}
void Fp_add(struct Fp *ans,struct Fp *a,struct Fp *b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_add(tmp.x0,a->x0,b->x0);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_add_ui(struct Fp *ans,struct Fp *a,unsigned long int b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_add_ui(tmp.x0,a->x0,b);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_sub_ui(struct Fp *ans,struct Fp *a,unsigned long int b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_sub_ui(tmp.x0,a->x0,b);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_sub(struct Fp *ans,struct Fp *a,struct Fp *b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_sub(tmp.x0,a->x0,b->x0);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_mul(struct Fp *ans,struct Fp *a,struct Fp *b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_mul(tmp.x0,a->x0,b->x0);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_mul_ui(struct Fp *ans,struct Fp *a,unsigned long int b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_mul_ui(tmp.x0,a->x0,b);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_div(struct Fp *ans,struct Fp *a,struct Fp *b){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_invert(tmp.x0,b->x0,p);
	mpz_mul(tmp.x0,a->x0,tmp.x0);
	mpz_mod(tmp.x0,tmp.x0,p);

	Fp_set(ans,&tmp);

	Fp_clear(&tmp);
}
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t B){
	int i;
	char B_binary[512];
	mpz_get_str(B_binary,2,B);
	struct Fp tmp;
	Fp_init(&tmp);
	Fp_set(&tmp,A);
	for(i=1;B_binary[i]!='\0';i++){
		Fp_mul(&tmp,&tmp,&tmp);
		if(B_binary[i]=='1'){
			Fp_mul(&tmp,&tmp,A);
		}
	}
	Fp_set(ANS,&tmp);
	Fp_clear(&tmp);
}
void Fp_sqrt(struct Fp *ANS,struct Fp *A){
	struct Fp n,y,x,b,t,tmp_Fp;
	Fp_init(&n);
	Fp_init(&y);
	Fp_init(&x);
	Fp_init(&b);
	Fp_init(&t);
	Fp_init(&tmp_Fp);

	Fp_set(&n,A);

	mpz_t tmp_mpz,q,e,r,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q);
	mpz_init(e);
	mpz_init(r);
	mpz_init(set_1);
	mpz_init(set_2);

	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(mpz_legendre(n.x0,p)!=-1){
		Fp_add_ui(&n,&n,1);
	}

	mpz_set(q,p);
	mpz_sub_ui(q,q,1);
	mpz_set_ui(e,0);

	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	Fp_pow(&y,&n,q);

	mpz_set(r,e);

	mpz_sub_ui(tmp_mpz,q,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);
	Fp_pow(&x,A,tmp_mpz);
	Fp_pow(&tmp_Fp,&x,set_2);
	Fp_mul(&b,&tmp_Fp,A);
	Fp_mul(&x,&x,A);

	int m;

	while(Fp_cmp_mpz(&b,set_1)==1){
		m=-1;
		while(Fp_cmp_mpz(&tmp_Fp,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp_pow(&tmp_Fp,&b,tmp_mpz);
		}
			// gmp_printf("%Zd,%Zd\n",x.x0.x0,x.x1.x0);
		mpz_sub_ui(tmp_mpz,r,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,p);
		Fp_pow(&t,&y,tmp_mpz);
		Fp_pow(&y,&t,set_2);
		mpz_set_ui(r,m);
		Fp_mul(&x,&x,&t);
		Fp_mul(&b,&b,&y);
	}

	Fp_set(ANS,&x);

	Fp_clear(&n);
	Fp_clear(&y);
	Fp_clear(&x);
	Fp_clear(&b);
	Fp_clear(&t);
	Fp_clear(&tmp_Fp);
	mpz_clear(tmp_mpz);
	mpz_clear(q);
	mpz_clear(e);
	mpz_clear(r);
	mpz_clear(set_1);
}
void Fp_neg(struct Fp *ans,struct Fp *a){
	struct Fp tmp;
	Fp_init(&tmp);

	mpz_sub(tmp.x0,p,a->x0);

	Fp_set(ans,&tmp);

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

void Fp2_init(struct Fp2 *A){
	Fp_init(&A->x0);
	Fp_init(&A->x1);
}
void Fp2_set(struct Fp2 *ANS,struct Fp2 *A){
	Fp_set(&ANS->x0,&A->x0);
	Fp_set(&ANS->x1,&A->x1);
}
void Fp2_set_ui(struct Fp2 *A,signed long int B){
	Fp_set_ui(&A->x0,B);
	Fp_set_ui(&A->x1,B);
}
void Fp2_random(struct Fp2 *A){
	Fp_random(&A->x0);
	Fp_random(&A->x1);
}
void Fp2_clear(struct Fp2 *A){
	Fp_clear(&A->x0);
	Fp_clear(&A->x1);
}
void Fp2_printf(struct Fp2 *A){
	gmp_printf("%Zd,%Zd\n",A->x0.x0,A->x1.x0);
}
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
	struct Fp2 tmp;
	Fp2_init(&tmp);

	Fp_add(&tmp.x0,&A->x0,&B->x0);
	Fp_add(&tmp.x1,&A->x1,&B->x1);

	Fp2_set(ANS,&tmp);

	Fp2_clear(&tmp);
}
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
	struct Fp2 tmp;
	Fp2_init(&tmp);

	Fp_add_ui(&tmp.x0,&A->x0,B);
	Fp_add_ui(&tmp.x1,&A->x1,B);

	Fp2_set(ANS,&tmp);

	Fp2_clear(&tmp);
}
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
	struct Fp2 tmp;
	Fp2_init(&tmp);

	Fp_sub(&tmp.x0,&A->x0,&B->x0);
	Fp_sub(&tmp.x1,&A->x1,&B->x1);

	Fp2_set(ANS,&tmp);

	Fp2_clear(&tmp);
}
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
	// (A0,A1)(B0,B1)=(A0B0-A1B1,A0B1+A1B0)	x^2+1=0
	struct Fp tmp1,tmp2,tmp3,tmp4,tmp5;
	struct Fp2 t_ans;
	Fp_init(&tmp1);
	Fp_init(&tmp2);
	Fp_init(&tmp3);
	Fp_init(&tmp4);
	Fp_init(&tmp5);
	Fp2_init(&t_ans);

	Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
	Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
	Fp_sub(&t_ans.x0,&tmp1,&tmp2);//a*c-b*d
	Fp_add(&tmp3,&A->x0,&A->x1);//a+b
	Fp_add(&tmp4,&B->x0,&B->x1);//c+d
	Fp_mul(&tmp5,&tmp3,&tmp4);//(a+b)(c+d)
	Fp_sub(&t_ans.x1,&tmp5,&tmp1);
	Fp_sub(&t_ans.x1,&t_ans.x1,&tmp2);

	Fp2_set(ANS,&t_ans);

	Fp_clear(&tmp1);
	Fp_clear(&tmp2);
	Fp_clear(&tmp3);
	Fp_clear(&tmp4);
	Fp_clear(&tmp5);
	Fp2_clear(&t_ans);
}
void Fp2_mul_i(struct Fp2 *ANS,struct Fp2 *A){
	//(a,b)(1,1)=(a-b,a+b)
	struct Fp2 tmp;
	Fp2_init(&tmp);
	Fp_sub(&tmp.x0,&A->x0,&A->x1);
	Fp_add(&tmp.x1,&A->x0,&A->x1);
	Fp2_set(ANS,&tmp);
	Fp2_clear(&tmp);
}
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
	struct Fp2 tmp;
	Fp2_init(&tmp);
	Fp_mul_ui(&tmp.x0,&A->x0,B);
	Fp_mul_ui(&tmp.x1,&A->x1,B);
	Fp2_set(ANS,&tmp);
	Fp2_clear(&tmp);
}
void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B){
	struct Fp2 tmp;
	Fp2_init(&tmp);

	Fp_mul(&tmp.x0,&A->x0,B);
	Fp_mul(&tmp.x1,&A->x1,B);

	Fp2_set(ANS,&tmp);

	Fp2_clear(&tmp);
}
void Fp2_invert(struct Fp2 *ANS,struct Fp2 *A){
	struct Fp c,c_inv,tmp1,tmp2;
	struct Fp2 A_asterisk,t_ans;
	Fp_init(&c);
	Fp_init(&c_inv);
	Fp_init(&tmp1);
	Fp_init(&tmp2);
	Fp2_init(&A_asterisk);
	Fp2_init(&t_ans);

	Fp_set(&A_asterisk.x0,&A->x0);
	Fp_neg(&A_asterisk.x1,&A->x1);//(a0,-a1)
	//(a,b)(a,-b)=(a^2+b^2,0)

	Fp_mul(&tmp1,&A->x0,&A->x0);
	Fp_mul(&tmp2,&A->x1,&A->x1);
	Fp_add(&c,&tmp1,&tmp2);
	mpz_invert(c_inv.x0,c.x0,p);

	Fp_mul(&t_ans.x0,&A_asterisk.x0,&c_inv);
	Fp_mul(&t_ans.x1,&A_asterisk.x1,&c_inv);

	Fp2_set(ANS,&t_ans);

	Fp_clear(&c);
	Fp_clear(&c_inv);
	Fp_clear(&tmp1);
	Fp_clear(&tmp2);
	Fp2_clear(&A_asterisk);
	Fp2_clear(&t_ans);
}
void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
	struct Fp2 tmp,t_ans;
	Fp2_init(&tmp);
	Fp2_init(&t_ans);

	Fp2_invert(&tmp,B);
	Fp2_mul(&t_ans,A,&tmp);

	Fp2_set(ANS,&t_ans);

	Fp2_clear(&tmp);
	Fp2_clear(&t_ans);
}
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
	int i;
	char B_binary[512];
	mpz_get_str(B_binary,2,B);
	struct Fp2 tmp;
	Fp2_init(&tmp);
	Fp2_set(&tmp,A);

	for(i=1;B_binary[i]!='\0';i++){
		Fp2_mul(&tmp,&tmp,&tmp);
		if(B_binary[i]=='1'){
			Fp2_mul(&tmp,&tmp,A);
		}
	}

	Fp2_set(ANS,&tmp);

	Fp2_clear(&tmp);
}
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A){
	struct Fp2 n,y,x,b,t,tmp_Fp2;
	Fp2_init(&n);
	Fp2_init(&y);
	Fp2_init(&x);
	Fp2_init(&b);
	Fp2_init(&t);
	Fp2_init(&tmp_Fp2);
	Fp2_set(&n,A);

	mpz_t tmp_mpz,q,e,r,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q);
	mpz_init(e);
	mpz_init(r);
	mpz_init(set_1);
	mpz_init(set_2);
	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(Fp2_legendre(&n)!=-1){
		Fp2_random(&n);
	}

	mpz_pow_ui(q,p,2);
	mpz_sub_ui(q,q,1);
	mpz_set_ui(e,0);

	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}

	Fp2_pow(&y,&n,q);
	mpz_set(r,e);
	mpz_sub_ui(tmp_mpz,q,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);
	Fp2_pow(&x,A,tmp_mpz);
	Fp2_pow(&tmp_Fp2,&x,set_2);
	Fp2_mul(&b,&tmp_Fp2,A);
	Fp2_mul(&x,&x,A);

	int m;

	while(Fp2_cmp_mpz(&b,set_1)==1){
		m=-1;
		Fp2_set(&tmp_Fp2,&b);
		while(Fp2_cmp_mpz(&tmp_Fp2,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp2_pow(&tmp_Fp2,&b,tmp_mpz);
		}
		mpz_sub_ui(tmp_mpz,r,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,p);
		Fp2_pow(&t,&y,tmp_mpz);
		Fp2_pow(&y,&t,set_2);
		mpz_set_ui(r,m);
		Fp2_mul(&x,&x,&t);
		Fp2_mul(&b,&b,&y);
	}

	Fp2_set(ANS,&x);

	Fp2_clear(&n);
	Fp2_clear(&y);
	Fp2_clear(&x);
	Fp2_clear(&b);
	Fp2_clear(&t);
	Fp2_clear(&tmp_Fp2);
	mpz_clear(tmp_mpz);
	mpz_clear(q);
	mpz_clear(e);
	mpz_clear(r);
	mpz_clear(set_1);
}
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B){
	if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
		return 0;
	}
	return 1;
}
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t B){
	struct Fp2 tmp;
	Fp2_init(&tmp);
	if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp(&A->x1,&tmp.x1)==0){
		Fp2_clear(&tmp);
		return 0;
	}
	Fp2_clear(&tmp);
	return 1;
}
int Fp2_legendre(struct Fp2 *a){
	mpz_t i;
	struct Fp2 tmp;
	Fp2_init(&tmp);
	mpz_init(i);

	mpz_pow_ui(i,p,2);
	mpz_sub_ui(i,i,1);
	mpz_div_ui(i,i,2);
	Fp2_pow(&tmp,a,i);

	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,1);

	if((Fp2_cmp_mpz(&tmp,cmp))==0){
		Fp2_clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return 1;
	}else{
		Fp2_clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return -1;
	}
}
void Fp2_neg(struct Fp2 *ans,struct Fp2 *a){
	struct Fp2 tmp;
	Fp2_init(&tmp);

	Fp_neg(&tmp.x0,&a->x0);
	Fp_neg(&tmp.x1,&a->x1);

	Fp2_set(ans,&tmp);

	Fp2_clear(&tmp);
}
//-----------------------------------------------------------------------------------------

void Fp6_init(struct Fp6 *A){
	Fp2_init(&A->x0);
	Fp2_init(&A->x1);
	Fp2_init(&A->x2);
}
void Fp6_set(struct Fp6 *ANS,struct Fp6 *A){
	Fp2_set(&ANS->x0,&A->x0);
	Fp2_set(&ANS->x1,&A->x1);
	Fp2_set(&ANS->x2,&A->x2);
}
void Fp6_set_ui(struct Fp6 *A,signed long int B){
	Fp2_set_ui(&A->x0,B);
	Fp2_set_ui(&A->x1,B);
	Fp2_set_ui(&A->x2,B);
}
void Fp6_random(struct Fp6 *A){
	Fp2_random(&A->x0);
	Fp2_random(&A->x1);
	Fp2_random(&A->x2);
}
void Fp6_clear(struct Fp6 *A){
	Fp2_clear(&A->x0);
	Fp2_clear(&A->x1);
	Fp2_clear(&A->x2);
}
void Fp6_printf(struct Fp6 *A){
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0,A->x0.x1.x0,A->x1.x0.x0,A->x1.x1.x0,A->x2.x0.x0,A->x2.x1.x0);
}
void Fp6_add(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp2_add(&tmp.x0,&A->x0,&B->x0);
	Fp2_add(&tmp.x1,&A->x1,&B->x1);
	Fp2_add(&tmp.x2,&A->x2,&B->x2);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_add_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp2_add_ui(&tmp.x0,&A->x0,B);
	Fp2_add_ui(&tmp.x1,&A->x1,B);
	Fp2_add_ui(&tmp.x2,&A->x2,B);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_sub(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp2_sub(&tmp.x0,&A->x0,&B->x0);
	Fp2_sub(&tmp.x1,&A->x1,&B->x1);
	Fp2_sub(&tmp.x2,&A->x2,&B->x2);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_mul(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
	//(x0,x1,x2)*(y0,y1,y2)=(x0y0+xi((x1+x2)(y1+y2)-x1y1-x2y2),xix2y2+(x0+x1)(y0+y1)-x0y0-x1y1,x1y1+(x0+x2)(y0+y2)-x0y0-x2y2)
	struct Fp2 tmp00,tmp11,tmp22,tmpx01,tmpx12,tmpx20,tmpy01,tmpy12,tmpy20,t0,t1,t2,tmp;
	struct Fp6 t_ans;
	Fp2_init(&tmp00);
	Fp2_init(&tmp11);
	Fp2_init(&tmp22);
	Fp2_init(&tmpx01);
	Fp2_init(&tmpx12);
	Fp2_init(&tmpx20);
	Fp2_init(&tmpy01);
	Fp2_init(&tmpy12);
	Fp2_init(&tmpy20);
	Fp2_init(&t0);
	Fp2_init(&t1);
	Fp2_init(&t2);
	Fp2_init(&tmp);
	Fp6_init(&t_ans);

	Fp2_mul(&tmp00,&A->x0,&B->x0);//x0*y0
	Fp2_mul(&tmp11,&A->x1,&B->x1);//x1*y1
	Fp2_mul(&tmp22,&A->x2,&B->x2);//x2*y2

	Fp2_add(&tmpx01,&A->x0,&A->x1);//x0+x1
	Fp2_add(&tmpx12,&A->x1,&A->x2);//x1+x2
	Fp2_add(&tmpx20,&A->x0,&A->x2);//x2+x0
	Fp2_add(&tmpy01,&B->x0,&B->x1);//y0+y1
	Fp2_add(&tmpy12,&B->x1,&B->x2);//y1+y2
	Fp2_add(&tmpy20,&B->x0,&B->x2);//y2+y0

	Fp2_mul(&t0,&tmpx01,&tmpy01);//(x0+x1)(y0+y1)
	Fp2_mul(&t1,&tmpx12,&tmpy12);//(x1+x2)(y1+y2)
	Fp2_mul(&t2,&tmpx20,&tmpy20);//(x2+x0)(y2+y0)

	Fp2_sub(&t1,&t1,&tmp11);
	Fp2_sub(&t1,&t1,&tmp22);//(x1+x2)(y1+y2)-x1y1-x2y2
	Fp2_mul_i(&tmp,&t1);
	Fp2_add(&t_ans.x0,&tmp00,&tmp);

	Fp2_sub(&t0,&t0,&tmp00);
	Fp2_sub(&t0,&t0,&tmp11);
	Fp2_mul_i(&tmp,&tmp22);
	Fp2_add(&t_ans.x1,&tmp,&t0);

	Fp2_sub(&t2,&t2,&tmp00);
	Fp2_sub(&t2,&t2,&tmp22);
	Fp2_add(&t_ans.x2,&tmp11,&t2);

	Fp6_set(ANS,&t_ans);

	Fp2_clear(&tmp00);
	Fp2_clear(&tmp11);
	Fp2_clear(&tmp22);
	Fp2_clear(&tmpx01);
	Fp2_clear(&tmpx12);
	Fp2_clear(&tmpx20);
	Fp2_clear(&tmpy01);
	Fp2_clear(&tmpy12);
	Fp2_clear(&tmpy20);
	Fp2_clear(&t0);
	Fp2_clear(&t1);
	Fp2_clear(&t2);
	Fp2_clear(&tmp);
}
void Fp6_mul_v(struct Fp6 *ANS,struct Fp6 *A){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp2_mul_i(&tmp.x0,&A->x2);
	Fp2_set(&tmp.x1,&A->x0);
	Fp2_set(&tmp.x2,&A->x1);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_mul_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp2_mul_ui(&tmp.x0,&A->x0,B);
	Fp2_mul_ui(&tmp.x1,&A->x1,B);
	Fp2_mul_ui(&tmp.x2,&A->x2,B);

	Fp6_set(ANS,&tmp);
	Fp6_clear(&tmp);
}
void Fp6_mul_Fp(struct Fp6 *ANS,struct Fp6 *A,struct Fp *B){
	struct Fp6 tmp;
	Fp6_init(&tmp);

	Fp2_mul_Fp(&tmp.x0,&A->x0,B);
	Fp2_mul_Fp(&tmp.x1,&A->x1,B);
	Fp2_mul_Fp(&tmp.x2,&A->x2,B);

	Fp6_set(ANS,&tmp);

	Fp6_clear(&tmp);
}
void Fp6_invert(struct Fp6 *ANS, struct Fp6 *A){
	struct Fp6 t_ans;
	Fp6_init(&t_ans);

	struct Fp2 T0,T1,t0,t1,t2,t3;
	Fp2_init(&T0);
	Fp2_init(&T1);
	Fp2_init(&t0);
	Fp2_init(&t1);
	Fp2_init(&t2);
	Fp2_init(&t3);

	// An optimized version of Grewal's Algo. 3   (a,b,c)
	Fp2_mul(&T0,&A->x0,&A->x0);
	Fp2_mul_i(&t0,&A->x1);

	Fp2_mul(&T1,&t0,&A->x2);
	Fp2_sub(&t1,&T0,&T1); // t1=(a^2-bci) mod q

	Fp2_mul(&T0,&A->x2,&A->x2);
	Fp2_mul_i(&T0,&T0);
	Fp2_mul(&T1,&A->x0,&A->x1);
	Fp2_sub(&t2,&T0,&T1); // t2=(c^2i-ab) mod q

	Fp2_mul(&T0,&A->x1,&A->x1);
	Fp2_mul(&T1,&A->x0,&A->x2);
	Fp2_sub(&t3,&T0,&T1); // t3=(b^2-ac) mod q

	Fp2_mul(&T0,&t0,&t3);
	Fp2_mul(&T1,&A->x0,&t1);
	Fp2_add(&T0,&T0,&T1); // T0={bi(b^2-ac)+a(a^2-bci)} mod q

	Fp2_mul_i(&t0,&A->x2);
	Fp2_mul(&T1,&t0,&t2);
	Fp2_add(&t0,&T0,&T1); // t0={ci(c^2i-ab)+{bi(b^2-ac)+a(a^2-bci)}} mod q .0

	Fp2_invert(&t0,&t0);

	Fp2_mul(&t_ans.x0,&t1,&t0);
	Fp2_mul(&t_ans.x1,&t2,&t0);
	Fp2_mul(&t_ans.x2,&t3,&t0);

	Fp6_set(ANS,&t_ans);

	Fp6_clear(&t_ans);
	Fp2_clear(&T0);
	Fp2_clear(&T1);
	Fp2_clear(&t0);
	Fp2_clear(&t1);
	Fp2_clear(&t2);
	Fp2_clear(&t3);
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
	// int i;
	// char B_binary[512];
	// mpz_get_str(B_binary,2,B);
	// struct Fp6 tmp;
	// Fp6_init(&tmp);
	// Fp6_set(&tmp,A);
	// for(i=1;B_binary[i]!='\0';i++){
	// 	Fp6_mul(&tmp,&tmp,&tmp);
	// 	if(B_binary[i]=='1'){
	// 		Fp6_mul(&tmp,&tmp,A);
	// 	}
	// }
	// Fp6_set(ANS,&tmp);
	// Fp6_clear(&tmp);
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

	mpz_pow_ui(q,p,6);
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
		mpz_powm(tmp_mpz,set_2,tmp_mpz,p);
		Fp6_pow(&t,&y,tmp_mpz);
		Fp6_pow(&y,&t,set_2);
		mpz_set_ui(r,m);
		Fp6_mul(&x,&x,&t);
		Fp6_mul(&b,&b,&y);
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
}
int Fp6_cmp(struct Fp6 *A,struct Fp6 *B){
	if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0 && Fp2_cmp(&A->x2,&B->x2)==0){
		return 0;
	}
	return 1;
}
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t B){
	struct Fp6 tmp;
	Fp6_init(&tmp);
	if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp(&A->x1,&tmp.x1)==0 && Fp2_cmp(&A->x2,&tmp.x2)==0){
		Fp6_clear(&tmp);
		return 0;
	}
	Fp6_clear(&tmp);
	return 1;
}
int Fp6_legendre(struct Fp6 *a){
	mpz_t i,cmp;
	mpz_init(i);
	mpz_init(cmp);
	mpz_set_ui(cmp,1);
	struct Fp6 tmp;
	Fp6_init(&tmp);
	mpz_pow_ui(i,p,6);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp6_pow(&tmp,a,i);

	if((Fp6_cmp_mpz(&tmp,cmp))==0){
		Fp6_clear(&tmp);
		mpz_clear(i);
		return 1;
	}else{
		Fp6_clear(&tmp);
		mpz_clear(i);
		return -1;
	}
}
void Fp6_neg(struct Fp6 *ans,struct Fp6 *a){
	struct Fp6 tmp;
	Fp6_init(&tmp);
	Fp2_neg(&tmp.x0,&a->x0);
	Fp2_neg(&tmp.x1,&a->x1);
	Fp2_neg(&tmp.x2,&a->x2);
	Fp6_set(ans,&tmp);
	Fp6_clear(&tmp);
}
//-----------------------------------------------------------------------------------------

void Fp12_init(struct Fp12 *A){
	Fp6_init(&A->x0);
	Fp6_init(&A->x1);
}
void Fp12_set(struct Fp12 *ANS,struct Fp12 *A){
	Fp6_set(&ANS->x0,&A->x0);
	Fp6_set(&ANS->x1,&A->x1);
}
void Fp12_set_ui(struct Fp12 *A,signed long int B){
	Fp6_set_ui(&A->x0,B);
	Fp6_set_ui(&A->x1,B);
}
void Fp12_random(struct Fp12 *A){
	Fp6_random(&A->x0);
	Fp6_random(&A->x1);
}
void Fp12_clear(struct Fp12 *A){
	Fp6_clear(&A->x0);
	Fp6_clear(&A->x1);
}
void Fp12_printf(struct Fp12 *A){
	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,",A->x0.x0.x0.x0,A->x0.x0.x1.x0,A->x0.x1.x0.x0,A->x0.x1.x1.x0,A->x0.x2.x0.x0,A->x0.x2.x1.x0);
	gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",A->x1.x0.x0.x0,A->x1.x0.x1.x0,A->x1.x1.x0.x0,A->x1.x1.x1.x0,A->x1.x2.x0.x0,A->x1.x2.x1.x0);
}
void Fp12_add(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
	struct Fp12 tmp;
	Fp12_init(&tmp);

	Fp6_add(&tmp.x0,&A->x0,&B->x0);
	Fp6_add(&tmp.x1,&A->x1,&B->x1);

	Fp12_set(ANS,&tmp);

	Fp12_clear(&tmp);
}
void Fp12_add_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int B){
	struct Fp12 tmp;
	Fp12_init(&tmp);

	Fp6_add_ui(&tmp.x0,&A->x0,B);
	Fp6_add_ui(&tmp.x1,&A->x1,B);

	Fp12_set(ANS,&tmp);

	Fp12_clear(&tmp);
}
void Fp12_sub(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
	struct Fp12 tmp;
	Fp12_init(&tmp);

	Fp6_sub(&tmp.x0,&A->x0,&B->x0);
	Fp6_sub(&tmp.x1,&A->x1,&B->x1);

	Fp12_set(ANS,&tmp);

	Fp12_clear(&tmp);
}
void Fp12_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
	//x^2-v=0
	struct Fp6 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	Fp6_init(&tmp1);
	Fp6_init(&tmp2);
	Fp6_init(&tmp3);
	Fp6_init(&tmp4);
	Fp6_init(&tmp5);
	Fp6_init(&tmp6);

	struct Fp12 t_ans;
	Fp12_init(&t_ans);

	Fp6_mul(&tmp1,&A->x0,&B->x0);//a*c
	Fp6_mul(&tmp2,&A->x1,&B->x1);//b*d
	Fp6_mul_v(&tmp3,&tmp2);//b*d*v
	Fp6_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
	Fp6_add(&tmp4,&A->x0,&A->x1);//a+b
	Fp6_add(&tmp5,&B->x0,&B->x1);//c+d
	Fp6_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
	Fp6_sub(&t_ans.x1,&tmp6,&tmp1);
	Fp6_sub(&t_ans.x1,&t_ans.x1,&tmp2);

	Fp12_set(ANS,&t_ans);

	Fp6_clear(&tmp1);
	Fp6_clear(&tmp2);
	Fp6_clear(&tmp3);
	Fp6_clear(&tmp4);
	Fp6_clear(&tmp5);
	Fp6_clear(&tmp6);
	Fp12_clear(&t_ans);
}
void Fp12_mul_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int B){
	struct Fp12 tmp;
	Fp12_init(&tmp);

	Fp6_mul_ui(&tmp.x0,&A->x0,B);
	Fp6_mul_ui(&tmp.x1,&A->x1,B);

	Fp12_set(ANS,&tmp);

	Fp12_clear(&tmp);
}
void Fp12_mul_Fp(struct Fp12 *ANS,struct Fp12 *A,struct Fp *B){
	struct Fp12 tmp;
	Fp12_init(&tmp);

	Fp6_mul_Fp(&tmp.x0,&A->x0,B);
	Fp6_mul_Fp(&tmp.x1,&A->x1,B);

	Fp12_set(ANS,&tmp);

	Fp12_clear(&tmp);
}
void Fp12_invert(struct Fp12 *ANS, struct Fp12 *A){
    struct Fp12 tmp;
    Fp12_init(&tmp);

    // tmp=A^(q^6)=(x0,-x1)
    Fp6_set(&tmp.x0,&A->x0);
    Fp6_neg(&tmp.x1,&A->x1);

    struct Fp6 c,a,b;
    Fp6_init(&c);
    Fp6_init(&a);
    Fp6_init(&b);

    Fp6_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp6_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp6_mul_v(&b,&b); // b=x1^2*v
    Fp6_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q

    Fp6_invert(&c,&c);

    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp6-element and tmp is a vector A Fp12
    Fp6_mul(&tmp.x0,&tmp.x0,&c);
    Fp6_mul(&tmp.x1,&tmp.x1,&c);

    Fp12_set(ANS,&tmp);

    Fp12_clear(&tmp);
}
void Fp12_div(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
	struct Fp12 tmp,t_ans;
	Fp12_init(&tmp);
	Fp12_init(&t_ans);

	Fp12_invert(&tmp,B);
	Fp12_mul(&t_ans,A,&tmp);

	Fp12_set(ANS,&t_ans);

	Fp12_clear(&tmp);
	Fp12_clear(&t_ans);
}
void Fp12_pow(struct Fp12 *ANS,struct Fp12 *A,mpz_t B){
	int i,length;
	length= (int)mpz_sizeinbase(B,2);
	char B_binary[length];
	mpz_get_str(B_binary,2,B);
	struct Fp12 tmp;
	Fp12_init(&tmp);
	Fp12_set(&tmp,A);
	for(i=1;B_binary[i]!='\0';i++){
		Fp12_mul(&tmp,&tmp,&tmp);
		if(B_binary[i]=='1'){
			Fp12_mul(&tmp,&tmp,A);
		}
	}
	Fp12_set(ANS,&tmp);
	Fp12_clear(&tmp);
}
void Fp12_sqrt(struct Fp12 *ANS,struct Fp12 *A){
	struct Fp12 n,y,x,b,t,tmp_Fp12;
	Fp12_init(&n);
	Fp12_init(&y);
	Fp12_init(&x);
	Fp12_init(&b);
	Fp12_init(&t);
	Fp12_init(&tmp_Fp12);
	Fp12_set(&n,A);

	mpz_t tmp_mpz,q,e,r,set_1,set_2;
	mpz_init(tmp_mpz);
	mpz_init(q);
	mpz_init(e);
	mpz_init(r);
	mpz_init(set_1);
	mpz_init(set_2);
	mpz_set_ui(set_1,1);
	mpz_set_ui(set_2,2);

	while(Fp12_legendre(&n)!=-1){
		Fp12_random(&n);
	}
	mpz_pow_ui(q,p,12);
	mpz_sub_ui(q,q,1);
	mpz_set_ui(e,0);
	while(mpz_odd_p(q)==0){
		mpz_add_ui(e,e,1);
		mpz_div_ui(q,q,2);
	}
	Fp12_pow(&y,&n,q);

	mpz_set(r,e);

	mpz_sub_ui(tmp_mpz,q,1);
	mpz_div_ui(tmp_mpz,tmp_mpz,2);

	Fp12_pow(&x,A,tmp_mpz);
	Fp12_pow(&tmp_Fp12,&x,set_2);
	Fp12_mul(&b,&tmp_Fp12,A);
	Fp12_mul(&x,&x,A);

	int m;

	while(Fp12_cmp_mpz(&b,set_1)==1){
		m=-1;
		Fp12_set(&tmp_Fp12,&b);
		while(Fp12_cmp_mpz(&tmp_Fp12,set_1)==1){
			m++;
			mpz_pow_ui(tmp_mpz,set_2,m);
			Fp12_pow(&tmp_Fp12,&b,tmp_mpz);
		}
		mpz_sub_ui(tmp_mpz,r,m);
		mpz_sub_ui(tmp_mpz,tmp_mpz,1);
		mpz_powm(tmp_mpz,set_2,tmp_mpz,p);
		// gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
		Fp12_pow(&t,&y,tmp_mpz);
		Fp12_pow(&y,&t,set_2);
		// gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
		mpz_set_ui(r,m);
		Fp12_mul(&x,&x,&t);
		Fp12_mul(&b,&b,&y);
	}

	Fp12_set(ANS,&x);

	Fp12_clear(&n);
	Fp12_clear(&y);
	Fp12_clear(&x);
	Fp12_clear(&b);
	Fp12_clear(&t);
	Fp12_clear(&tmp_Fp12);
	mpz_clear(tmp_mpz);
	mpz_clear(q);
	mpz_clear(e);
	mpz_clear(r);
	mpz_clear(set_1);
}
int Fp12_legendre(struct Fp12 *a){
	mpz_t i,cmp;
	struct Fp12 tmp;
	Fp12_init(&tmp);
	mpz_init(i);
	mpz_init(cmp);
	mpz_set_ui(cmp,1);
	mpz_pow_ui(i,p,12);
	mpz_sub_ui(i,i,1);
	mpz_tdiv_q_ui(i,i,2);
	Fp12_pow(&tmp,a,i);

	if((Fp12_cmp_mpz(&tmp,cmp))==0){
		Fp12_clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return 1;
	}else{
		Fp12_clear(&tmp);
		mpz_clear(i);
		mpz_clear(cmp);
		return -1;
	}
}
int Fp12_cmp(struct Fp12 *A,struct Fp12 *B){
	if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
		return 0;
	}
	return 1;
}
int Fp12_cmp_mpz(struct Fp12 *A,mpz_t B){
	struct Fp12 tmp;
	Fp12_init(&tmp);
	if(Fp6_cmp_mpz(&A->x0,B)==0 && Fp6_cmp(&A->x1,&tmp.x1)==0){
		Fp12_clear(&tmp);
		return 0;
	}
	Fp12_clear(&tmp);
	return 1;
}
void Fp12_neg(struct Fp12 *ans,struct Fp12 *a){
	struct Fp12 tmp;
	Fp12_init(&tmp);
	Fp6_neg(&tmp.x0,&a->x0);
	Fp6_neg(&tmp.x1,&a->x1);
	Fp12_set(ans,&tmp);
	Fp12_clear(&tmp);
}
//-----------------------------------------------------------------------------------------

void EFp_init(struct EFp *A){
	Fp_init(&A->x);
	Fp_init(&A->y);
	A->PaI=FALSE;
}
void EFp_set(struct EFp *A,struct EFp *B){
	Fp_set(&A->x,&B->x);
	Fp_set(&A->y,&B->y);
	A->PaI=B->PaI;
}
void EFp_set_PaI(struct EFp *A){
	Fp_set_ui(&A->x,0);
	Fp_set_ui(&A->y,0);
	A->PaI=TRUE;
}
void EFp_clear(struct EFp *A){
	Fp_clear(&A->x);
	Fp_clear(&A->y);
}
void EFp_printf(struct EFp *A){
	gmp_printf("(%Zd,%Zd)\n",A->x.x0,A->y.x0);
}
void EFp_SCM(struct EFp *ANS, struct EFp *P,mpz_t j){
	int i,length;
	length= (int)mpz_sizeinbase(j,2);
	// int eca=0,ecd=0;
	char r_binary[length];
	mpz_get_str(r_binary,2,j);
	struct EFp Q,R;
	EFp_init(&Q);
	EFp_set(&Q,P);
	EFp_init(&R);
	for(i=1;r_binary[i]!='\0';i++){
		EFp_ECD(&Q,&Q);
		// ecd++;
		if(r_binary[i]=='1'){
			EFp_ECA(&Q,&Q,P);
			// eca++;
		}
	}
	EFp_set(ANS,&Q);
	// printf("%d,%d\n",eca,ecd);
	EFp_clear(&Q);
	EFp_clear(&R);
	return;
}
void EFp_ECD(struct EFp *ANS, struct EFp *P){
	if(P->PaI==TRUE){
		EFp_set(ANS,P);
		return;
	}
	if(mpz_sgn(P->y.x0)==0){//P.y==0
		EFp_set_PaI(ANS);
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
	if(P2->PaI==TRUE){//if P2==inf
		EFp_set(ANS,P1);
		return;
	}
	else if(P1->PaI==TRUE){//if P1==inf
		EFp_set(ANS,P2);
		return;
	}
	else if(Fp_cmp(&P1->x,&P2->x)==0&&Fp_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp_set_PaI(ANS);
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
	// int legendle;
	struct EFp P;
	EFp_init(&P);

	struct Fp x,a;
	Fp_init(&a);

	Fp_init(&x);
	mpz_t set_3;
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	do{
		Fp_random(&x);
		Fp_pow(&a,&x,set_3);
		mpz_add(a.x0,a.x0,b);
	}while(mpz_legendre(a.x0,p)!=1);

	Fp_sqrt(&P.y,&a);
	Fp_set(&P.x,&x);
	EFp_set(ANS,&P);

	EFp_clear(&P);
	Fp_clear(&a);
	Fp_clear(&x);
	mpz_clear(set_3);
}
void EFp12_set_EFp(struct EFp12 *A,struct EFp *B){
	Fp12_set_ui(&A->x,0);
	Fp12_set_ui(&A->y,0);

	Fp_set(&A->x.x0.x0.x0,&B->x);
	Fp_set(&A->y.x0.x0.x0,&B->y);
	A->PaI=B->PaI;
}
//-----------------------------------------------------------------------------------------

void EFp2_init(struct EFp2 *A){
	Fp2_init(&A->x);
	Fp2_init(&A->y);
	A->PaI=FALSE;
}
void EFp2_set(struct EFp2 *A,struct EFp2 *B){
	Fp2_set(&A->x,&B->x);
	Fp2_set(&A->y,&B->y);
	A->PaI=B->PaI;
}
void EFp2_set_PaI(struct EFp2 *A){
	Fp2_set_ui(&A->x,0);
	Fp2_set_ui(&A->y,0);
	A->PaI=TRUE;
}
void EFp2_clear(struct EFp2 *A){
	Fp2_clear(&A->x);
	Fp2_clear(&A->y);
}
void EFp2_printf(struct EFp2 *A){
	gmp_printf("(%Zd,%Zd)(%Zd,%Zd)\n",A->x.x0.x0,A->x.x1.x0,A->y.x0.x0,A->y.x1.x0);
}
void EFp2_SCM(struct EFp2 *ANS,struct EFp2 *P,mpz_t j){
	int i,length;
	length= (int)mpz_sizeinbase(j,2);
	char j_binary[length];
	mpz_get_str(j_binary,2,j);
	struct EFp2 Q,R;
	EFp2_init(&Q);
	EFp2_set(&Q,P);
	EFp2_init(&R);
	for(i=1;j_binary[i]!='\0';i++){
		EFp2_ECD(&Q,&Q);
		if(j_binary[i]=='1'){
			EFp2_ECA(&Q,&Q,P);
		}
	}
	EFp2_set(ANS,&Q);

	EFp2_clear(&Q);
	EFp2_clear(&R);
	return;
}
void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P){
	if(P->PaI==TRUE){
		EFp2_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp2_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp2_set_PaI(ANS);
		return;
	}

	struct Fp2 x,y,lambda,tmp;
	struct EFp2 t_ans;
	Fp2_init(&x);
	Fp2_init(&lambda);
	Fp2_init(&tmp);
	Fp2_init(&y);
	EFp2_init(&t_ans);

	Fp2_mul(&x,&P->x,&P->x);
	Fp2_add(&tmp,&x,&x);
	Fp2_add(&x,&tmp,&x);
	Fp2_add(&y,&P->y,&P->y);
	Fp2_div(&lambda,&x,&y);
	Fp2_mul(&tmp,&lambda,&lambda);
	Fp2_add(&x,&P->x,&P->x);
	Fp2_sub(&x,&tmp,&x);
	Fp2_sub(&tmp,&P->x,&x);
	Fp2_set(&t_ans.x,&x);
	Fp2_mul(&tmp,&tmp,&lambda);
	Fp2_sub(&t_ans.y,&tmp,&P->y);

	EFp2_set(ANS,&t_ans);

	Fp2_clear(&x);
	Fp2_clear(&lambda);
	Fp2_clear(&y);
	Fp2_clear(&tmp);
	EFp2_clear(&t_ans);
}
void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2){
	if(P2->PaI==TRUE){//if P2==inf
		EFp2_set(ANS,P1);
		return;
	}
	else if(P1->PaI==TRUE){//if P1==inf
		EFp2_set(ANS,P2);
		return;
	}
	else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp2_set_PaI(ANS);
		return;
	}
	else if(EFp2_cmp(P1,P2)==0){ // P=Q
		EFp2_ECD(ANS,P1);
		return;
	}

	struct Fp2 x,y,lambda,tmp;
	struct EFp2 t_ans;

	Fp2_init(&x);
	Fp2_init(&y);
	Fp2_init(&lambda);
	Fp2_init(&tmp);
	EFp2_init(&t_ans);

	Fp2_sub(&x,&P2->x,&P1->x);
	Fp2_sub(&y,&P2->y,&P1->y);
	Fp2_div(&lambda,&y,&x);
	Fp2_mul(&tmp,&lambda,&lambda);
	Fp2_add(&x,&P1->x,&P2->x);
	Fp2_sub(&x,&tmp,&x);
	Fp2_sub(&tmp,&P1->x,&x);
	Fp2_set(&t_ans.x,&x);
	Fp2_mul(&tmp,&tmp,&lambda);
	Fp2_sub(&t_ans.y,&tmp,&P1->y);

	EFp2_set(ANS,&t_ans);

	Fp2_clear(&x);
	Fp2_clear(&y);
	Fp2_clear(&lambda);
	Fp2_clear(&tmp);
	EFp2_clear(&t_ans);
}
int EFp2_cmp(struct EFp2 *A,struct EFp2 *B){
	if(Fp2_cmp(&A->x,&B->x)==0 && Fp2_cmp(&A->y,&B->y)==0){
		return 0;
	}
	return 1;
}
void EFp2_random_set(struct EFp2 *ANS){
	struct EFp2 P;
	EFp2_init(&P);

	struct Fp2 x,a;
	Fp2_init(&a);
	Fp2_init(&x);

	mpz_t t2,p2,p22,tmp,r2;
	mpz_t set_3;
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	mpz_init(t2);
	mpz_init(p2);
	mpz_init(p22);
	mpz_init(tmp);
	mpz_init(r2);

	mpz_add(p22,p,p);
	mpz_mul(p2,p,p);
	mpz_mul(t2,t,t);
	mpz_sub(tmp,t2,p22);
	mpz_add_ui(r2,p2,1);
	mpz_sub(r2,r2,tmp);

	do{
		Fp2_random(&x);
		Fp2_pow(&a,&x,set_3);
		mpz_add(a.x0.x0,a.x0.x0,b);
	}while(Fp2_legendre(&a)!=1);
	Fp2_sqrt(&P.y,&a);
	Fp2_set(&P.x,&x);

	mpz_t r12_div_r2;
	mpz_init(r12_div_r2);
	mpz_div(r12_div_r2,r2,r);
	mpz_div(r12_div_r2,r12_div_r2,r);

	EFp2_SCM(ANS,&P,r12_div_r2);

	EFp2_clear(&P);
	Fp2_clear(&a);
	Fp2_clear(&x);
}
//-----------------------------------------------------------------------------------------

void EFp6_init(struct EFp6 *A){
	Fp6_init(&A->x);
	Fp6_init(&A->y);
	A->PaI=FALSE;
}
void EFp6_set(struct EFp6 *A,struct EFp6 *B){
	Fp6_set(&A->x,&B->x);
	Fp6_set(&A->y,&B->y);
	A->PaI=B->PaI;
}
void EFp6_set_PaI(struct EFp6 *A){
	Fp6_set_ui(&A->x,0);
	Fp6_set_ui(&A->y,0);
	A->PaI=TRUE;
}
void EFp6_clear(struct EFp6 *A){
	Fp6_clear(&A->x);
	Fp6_clear(&A->y);
}
void EFp6_printf(struct EFp6 *A){
	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)",A->x.x0.x0.x0,A->x.x0.x1.x0,A->x.x1.x0.x0,A->x.x1.x1.x0,A->x.x2.x0.x0,A->x.x2.x1.x0);
	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",A->y.x0.x0.x0,A->y.x0.x1.x0,A->y.x1.x0.x0,A->y.x1.x1.x0,A->y.x2.x0.x0,A->y.x2.x1.x0);
}
void EFp6_ECD(struct EFp6 *ANS, struct EFp6 *P){
	if(P->PaI==TRUE){
		EFp6_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp6_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp6_set_PaI(ANS);
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
	if(P2->PaI==TRUE){//if P2==inf
		EFp6_set(ANS,P1);
		return;
	}
	else if(P1->PaI==TRUE){//if P1==inf
		EFp6_set(ANS,P2);
		return;
	}
	else if(Fp6_cmp(&P1->x,&P2->x)==0&&Fp6_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp6_set_PaI(ANS);
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
void EFp6_SCM(struct EFp6 *ANS,struct EFp6 *P,mpz_t j){
	int i,length;
	length= (int)mpz_sizeinbase(j,2);
	char j_binary[length];
	mpz_get_str(j_binary,2,j);
	struct EFp6 Q,R;
	EFp6_init(&Q);
	EFp6_set(&Q,P);
	EFp6_init(&R);
	for(i=1;j_binary[i]!='\0';i++){
		EFp6_ECD(&Q,&Q);
		if(j_binary[i]=='1'){
			EFp6_ECA(&Q,&Q,P);
		}
	}
	EFp6_set(ANS,&Q);

	EFp6_clear(&Q);
	EFp6_clear(&R);
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

	mpz_t set_3;
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	//t12=a^12+b^12={(t^3-3pt)^2-2p^3}^2-2p^6
	mpz_t tmp1,tmp2,p_pow;
	mpz_t r12,t12;
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(p_pow);
	mpz_init(r12);
	mpz_init(t12);

	mpz_pow_ui(tmp1,t,3);
	mpz_mul(p_pow,p,t);
	mpz_mul_ui(p_pow,p_pow,3);
	mpz_sub(tmp1,tmp1,p_pow);
	mpz_pow_ui(tmp1,tmp1,2);

	mpz_pow_ui(p_pow,p,3);
	mpz_mul_ui(tmp2,p_pow,2);
	mpz_sub(t12,tmp1,tmp2);
	// mpz_pow_ui(tmp1,tmp1,2);

	mpz_pow_ui(p_pow,p_pow,2);
	// mpz_mul_ui(tmp2,p_pow,2);
	// mpz_sub(t12,tmp1,tmp2);

	// mpz_pow_ui(p_pow,p_pow,2);
	mpz_add_ui(tmp1,p_pow,1);
	mpz_sub(r12,tmp1,t12);

	do{
		Fp6_random(&x);
		Fp6_pow(&a,&x,set_3);
		mpz_add(a.x0.x0.x0,a.x0.x0.x0,b);
	}while(Fp6_legendre(&a)!=1);
	Fp6_sqrt(&P.y,&a);
	Fp6_set(&P.x,&x);

	mpz_t r12_div_r2;
	mpz_init(r12_div_r2);
	mpz_div(r12_div_r2,r12,r);
	mpz_div(r12_div_r2,r12_div_r2,r);

	EFp6_SCM(ANS,&P,r12_div_r2);

	EFp6_clear(&P);
	Fp6_clear(&a);
	Fp6_clear(&x);
	mpz_clear(set_3);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(p_pow);
	mpz_clear(r12);
	mpz_clear(t12);
	mpz_clear(r12_div_r2);
}
//-----------------------------------------------------------------------------------------

void EFp12_init(struct EFp12 *A){
	Fp12_init(&A->x);
	Fp12_init(&A->y);
	A->PaI=FALSE;
}
void EFp12_set(struct EFp12 *A,struct EFp12 *B){
	Fp12_set(&A->x,&B->x);
	Fp12_set(&A->y,&B->y);
	A->PaI=B->PaI;
}
void EFp12_set_PaI(struct EFp12 *A){
	Fp12_set_ui(&A->x,0);
	Fp12_set_ui(&A->y,0);
	A->PaI=TRUE;
}
void EFp12_clear(struct EFp12 *A){
	Fp12_clear(&A->x);
	Fp12_clear(&A->y);
}
void EFp12_printf(struct EFp12 *A){
	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)",A->x.x0.x0.x0.x0,A->x.x0.x0.x1.x0,A->x.x0.x1.x0.x0,A->x.x0.x1.x1.x0,A->x.x0.x2.x0.x0,A->x.x0.x2.x1.x0,A->x.x1.x0.x0.x0,A->x.x1.x0.x1.x0,A->x.x1.x1.x0.x0,A->x.x1.x1.x1.x0,A->x.x1.x2.x0.x0,A->x.x1.x2.x1.x0);
	gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",A->y.x0.x0.x0.x0,A->y.x0.x0.x1.x0,A->y.x0.x1.x0.x0,A->y.x0.x1.x1.x0,A->y.x0.x2.x0.x0,A->y.x0.x2.x1.x0,A->y.x1.x0.x0.x0,A->y.x1.x0.x1.x0,A->y.x1.x1.x0.x0,A->y.x1.x1.x1.x0,A->y.x1.x2.x0.x0,A->y.x1.x2.x1.x0);
}
void EFp12_ECD(struct EFp12 *ANS, struct EFp12 *P){
	if(P->PaI==TRUE){
		EFp12_set(ANS,P);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp12_cmp_mpz(&P->y,cmp)==0){//P.y==0
		EFp12_set_PaI(ANS);
		return;
	}

	struct Fp12 x,y,lambda,tmp;
	struct EFp12 t_ans;
	Fp12_init(&x);
	Fp12_init(&lambda);
	Fp12_init(&tmp);
	Fp12_init(&y);
	EFp12_init(&t_ans);

	Fp12_mul(&x,&P->x,&P->x);
	Fp12_add(&tmp,&x,&x);
	Fp12_add(&x,&tmp,&x);
	Fp12_add(&y,&P->y,&P->y);

	Fp12_div(&lambda,&x,&y);
	Fp12_mul(&tmp,&lambda,&lambda);




	// Fp12_printf(&P->x);
	Fp12_add(&x,&P->x,&P->x);
	Fp12_sub(&x,&tmp,&x);
	Fp12_sub(&tmp,&P->x,&x);
	Fp12_set(&t_ans.x,&x);
	Fp12_mul(&tmp,&tmp,&lambda);
	Fp12_sub(&t_ans.y,&tmp,&P->y);

	EFp12_set(ANS,&t_ans);

	Fp12_clear(&x);
	Fp12_clear(&lambda);
	Fp12_clear(&y);
	Fp12_clear(&tmp);
	EFp12_clear(&t_ans);
}
void EFp12_ECA(struct EFp12 *ANS, struct EFp12 *P1, struct EFp12 *P2){
	if(P2->PaI==TRUE){//if P2==inf
		EFp12_set(ANS,P1);
		return;
	}
	else if(P1->PaI==TRUE){//if P1==inf
		EFp12_set(ANS,P2);
		return;
	}
	else if(Fp12_cmp(&P1->x,&P2->x)==0&&Fp12_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
		EFp12_set_PaI(ANS);
		return;
	}
	else if(EFp12_cmp(P1,P2)==0){ // P=Q
		EFp12_ECD(ANS,P1);
		return;
	}

	struct Fp12 x,y,lambda,tmp;
	struct EFp12 t_ans;

	Fp12_init(&x);
	Fp12_init(&y);
	Fp12_init(&lambda);
	Fp12_init(&tmp);
	EFp12_init(&t_ans);

	Fp12_sub(&x,&P2->x,&P1->x);
	Fp12_sub(&y,&P2->y,&P1->y);
	Fp12_div(&lambda,&y,&x);
	Fp12_mul(&tmp,&lambda,&lambda);
	Fp12_add(&x,&P1->x,&P2->x);
	Fp12_sub(&x,&tmp,&x);
	Fp12_sub(&tmp,&P1->x,&x);
	Fp12_set(&t_ans.x,&x);
	Fp12_mul(&tmp,&tmp,&lambda);
	Fp12_sub(&t_ans.y,&tmp,&P1->y);

	EFp12_set(ANS,&t_ans);

	Fp12_clear(&x);
	Fp12_clear(&y);
	Fp12_clear(&lambda);
	Fp12_clear(&tmp);
	EFp12_clear(&t_ans);
}
int EFp12_cmp(struct EFp12 *A,struct EFp12 *B){
	if(Fp12_cmp(&A->x,&B->x)==0 && Fp12_cmp(&A->y,&B->y)==0){
		return 0;
	}
	return 1;
}
void EFp12_SCM(struct EFp12 *ANS, struct EFp12 *P, mpz_t j){
	int i,length;
	length= (int)mpz_sizeinbase(j,2);
	char r_binary[length];
	// printf("%d\n",length);
	mpz_get_str(r_binary,2,j);
	struct EFp12 Q,R;
	EFp12_init(&Q);
	EFp12_set(&Q,P);
	EFp12_init(&R);
	for(i=1;r_binary[i]!='\0';i++){
		EFp12_ECD(&Q,&Q);
		if(r_binary[i]=='1'){
			EFp12_ECA(&Q,&Q,P);
		}
	}
	EFp12_set(ANS,&Q);

	EFp12_clear(&Q);
	EFp12_clear(&R);
	return;
}
void EFp12_frobenius_map(struct EFp12 *ANS,struct EFp12 *A){
	struct EFp12 tmp;
	EFp12_init(&tmp);

	Fp12_pow(&tmp.x,&A->x,p);
	Fp12_pow(&tmp.y,&A->y,p);

	EFp12_set(ANS,&tmp);

	EFp12_clear(&tmp);
}
void EFp12_random_set(struct EFp12 *ANS){
	struct EFp12 P;
	EFp12_init(&P);

	struct Fp12 x,a;
	Fp12_init(&a);
	Fp12_init(&x);

	mpz_t set_3;
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	//t12=a^12+b^12={(t^3-3pt)^2-2p^3}^2-2p^6
	mpz_t tmp1,tmp2,p_pow;
	mpz_t r12,t12;
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(p_pow);
	mpz_init(r12);
	mpz_init(t12);

	mpz_pow_ui(tmp1,t,3);
	mpz_mul(p_pow,p,t);
	mpz_mul_ui(p_pow,p_pow,3);
	mpz_sub(tmp1,tmp1,p_pow);
	mpz_pow_ui(tmp1,tmp1,2);

	mpz_pow_ui(p_pow,p,3);
	mpz_mul_ui(tmp2,p_pow,2);
	mpz_sub(tmp1,tmp1,tmp2);
	mpz_pow_ui(tmp1,tmp1,2);

	mpz_pow_ui(p_pow,p_pow,2);
	mpz_mul_ui(tmp2,p_pow,2);
	mpz_sub(t12,tmp1,tmp2);

	mpz_pow_ui(p_pow,p_pow,2);
	mpz_add_ui(tmp1,p_pow,1);
	mpz_sub(r12,tmp1,t12);

	do{
		Fp12_random(&x);
		Fp12_pow(&a,&x,set_3);
		mpz_add(a.x0.x0.x0.x0,a.x0.x0.x0.x0,b);
	}while(Fp12_legendre(&a)!=1);
	Fp12_sqrt(&P.y,&a);
	Fp12_set(&P.x,&x);

	mpz_t r12_div_r2;
	mpz_init(r12_div_r2);
	mpz_div(r12_div_r2,r12,r);
	mpz_div(r12_div_r2,r12_div_r2,r);
	EFp12_SCM(ANS,&P,r12_div_r2);

	EFp12_clear(&P);
	Fp12_clear(&a);
	Fp12_clear(&x);
	mpz_clear(set_3);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(p_pow);
	mpz_clear(r12);
	mpz_clear(t12);
	mpz_clear(r12_div_r2);
}
void EFp12_random_set_for_Ate(struct EFp12 *ANS){
	struct EFp12 P,P_frobenius,tmp_EFp12;
	EFp12_init(&P);
	EFp12_init(&P_frobenius);
	EFp12_init(&tmp_EFp12);

	EFp12_random_set(&P);

	EFp12_frobenius_map(&P_frobenius,&P);
	Fp12_neg(&tmp_EFp12.y,&P.y);
	Fp12_set(&tmp_EFp12.x,&P.x);

	EFp12_ECA(&tmp_EFp12,&tmp_EFp12,&P_frobenius);

	EFp12_set(ANS,&tmp_EFp12);

	EFp12_clear(&P);
	EFp12_clear(&P_frobenius);
	EFp12_clear(&tmp_EFp12);
}
//-----------------------------------------------------------------------------------------

void set_BN_parameter(struct EFp *ANS){
	//set p,r,t
	mpz_t x1,x2,x3,x4;
	int i,j,k;
	mpz_init(x1);
	mpz_init(x2);
	mpz_init(x3);
	mpz_init(x4);

	mpz_set(x1,X);

	mpz_mul(x2,x1,x1);
	mpz_mul(x3,x2,x1);
	mpz_mul(x4,x2,x2);
	mpz_mul_si(x4,x4,36);
	mpz_mul_si(x3,x3,-36);
	mpz_mul_si(x2,x2,24);
	mpz_mul_si(x1,x1,-6);

	mpz_add(p,x4,x3);
	mpz_add(p,p,x2);
	mpz_add(p,p,x1);
	mpz_add_ui(p,p,1);

	mpz_set(x1,X);
	mpz_mul(x2,x1,x1);
	mpz_mul_si(x2,x2,6);
	mpz_add_ui(t,x2,1);
	mpz_sub(r,p,t);
	mpz_add_ui(r,r,1);

	i=mpz_probab_prime_p(p,25);
	mpz_t mod,tmp;
	mpz_init(mod);
	mpz_init(tmp);
	mpz_set_ui(mod,4);
	mpz_mod(tmp,p,mod);
	j=mpz_get_ui(tmp);
	mpz_set_ui(mod,8);
	mpz_mod(tmp,p,mod);
	k=mpz_get_ui(tmp);
	if(i==0){
		gmp_printf("not prime number\n");
		mpz_clear(x4);
		mpz_clear(x3);
		mpz_clear(x2);
		mpz_clear(x1);
		mpz_clear(tmp);
		mpz_clear(mod);
		exit(0);
	}else if(j==1){
		gmp_printf("p-1 is multiples of 4\n");
		mpz_clear(x4);
		mpz_clear(x3);
		mpz_clear(x2);
		mpz_clear(x1);
		mpz_clear(tmp);
		mpz_clear(mod);
		exit(0);
	}else if(k!=3){
		gmp_printf("p-1 is multiples of 8\n");
		mpz_clear(x4);
		mpz_clear(x3);
		mpz_clear(x2);
		mpz_clear(x1);
		mpz_clear(tmp);
		mpz_clear(mod);
		exit(0);
	}

	mpz_clear(x4);
	mpz_clear(x3);
	mpz_clear(x2);
	mpz_clear(x1);
	mpz_clear(tmp);
	mpz_clear(mod);

	//set b

	struct EFp P;
	int legendle;
	struct Fp a,x;
	mpz_t tmp_b;
	Fp_init(&a);
	EFp_init(&P);
	Fp_init(&x);
	mpz_init(tmp_b);

	mpz_set_ui(tmp_b,0);

	for(;;){
		mpz_add_ui(tmp_b,tmp_b,1);
		Fp_set_ui(&x,1);
		legendle=0;
		while(legendle!=1){
			mpz_powm_ui(a.x0,x.x0,3,p);
			mpz_add(a.x0,a.x0,tmp_b);
			if((legendle=mpz_legendre(a.x0,p))==1){
				Fp_sqrt(&P.y,&a);
				Fp_set(&P.x,&x);
				EFp_SCM(ANS,&P,r);
				if(ANS->PaI==TRUE){
					mpz_set(b,tmp_b);
					mpz_clear(tmp_b);
					Fp_clear(&a);
					Fp_clear(&x);
					EFp_clear(&P);
					return;
				}
			}
			Fp_add_ui(&x,&x,1);
		}
	}
	return;
}
//-----------------------------------------------------------------------------------------

void Miller_algo(struct Fp12 *ANS,struct EFp12 *P,struct EFp12 *Q,mpz_t roop){
	struct Fp12 f,l_sum,v_sum;
	Fp12_init(&f);
	Fp12_init(&l_sum);
	Fp12_init(&v_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);
	Fp_set_ui(&v_sum.x0.x0.x0,1);

	struct EFp12 T;
	EFp12_init(&T);
	EFp12_set(&T,P);


	struct Fp12 ltt,ltp,v2t,vtp;
	Fp12_init(&ltt);
	Fp12_init(&ltp);
	Fp12_init(&v2t);
	Fp12_init(&vtp);

	int i;
	struct Fp12 tmp1,lambda;
	Fp12_init(&tmp1);
	Fp12_init(&lambda);

	int r_bit;//bit数
	r_bit= (int)mpz_sizeinbase(roop,2);

	for(i=r_bit-2;i>=0;i--){
		if(mpz_tstbit(roop,i)==1){
			Fp12_mul(&l_sum,&l_sum,&l_sum);
			Fp12_mul(&v_sum,&v_sum,&v_sum);

			ltt_q(&ltt,&T,Q);
			Fp12_mul(&l_sum,&l_sum,&ltt);

			EFp12_ECD(&T,&T);

			v2t_q(&v2t,&T,Q);
			Fp12_mul(&v_sum,&v_sum,&v2t);
			Fp12_mul(&f,&f,&tmp1);

			ltp_q(&ltp,&T,P,Q);
			Fp12_mul(&l_sum,&l_sum,&ltp);

			EFp12_ECA(&T,&T,P);

			vtp_q(&vtp,&T,Q);
			Fp12_mul(&v_sum,&v_sum,&vtp);
		}else{
			Fp12_mul(&l_sum,&l_sum,&l_sum);
			Fp12_mul(&v_sum,&v_sum,&v_sum);

			ltt_q(&ltt,&T,Q);
			Fp12_mul(&l_sum,&l_sum,&ltt);

			EFp12_ECD(&T,&T);

			v2t_q(&v2t,&T,Q);
			Fp12_mul(&v_sum,&v_sum,&v2t);
			Fp12_mul(&f,&f,&tmp1);
		}
	}

	Fp12_div(&f,&l_sum,&v_sum);

	// Fp12_set(ANS,&f);
	Fp12_set(ANS,&l_sum);

	Fp12_clear(&f);
	EFp12_clear(&T);
	Fp12_clear(&tmp1);
	Fp12_clear(&lambda);
}
void Final_Exp(struct Fp12 *ANS,struct Fp12 *A){
	mpz_t p12;
	mpz_init(p12);

	mpz_pow_ui(p12,p,12);
	mpz_sub_ui(p12,p12,1);
	mpz_div(p12,p12,r);

	Fp12_pow(ANS,A,p12);

	mpz_clear(p12);
}
void Tate_Pairing(struct Fp12 *ANS,struct EFp12 *G1,struct EFp12 *G2){
	struct Fp12 t_ans;
	Fp12_init(&t_ans);

	Miller_algo(&t_ans,G2,G1,r);
	Final_Exp(&t_ans,&t_ans);
	Fp12_set(ANS,&t_ans);

	Fp12_clear(&t_ans);
}
void Ate_Pairing(struct Fp12 *ANS,struct EFp12 *G1,struct EFp12 *G2){
	struct Fp12 t_ans;
	Fp12_init(&t_ans);

	mpz_t tm1;
	mpz_init(tm1);
	mpz_sub_ui(tm1,t,1);

	Miller_algo(&t_ans,G2,G1,tm1);
	Final_Exp(&t_ans,&t_ans);
	Fp12_set(ANS,&t_ans);

	Fp12_clear(&t_ans);
}
void Optimal_Ate_Pairing(struct Fp12 *ANS,struct EFp12 *G1,struct EFp12 *G2){
	mpz_t p3;
	mpz_init(p3);
	mpz_mul_ui(p3,p,3);
	struct EFp12 zQ,p3Q;
	EFp12_init(&zQ);
	EFp12_init(&p3Q);

	struct Fp12 ltp;
	Fp12_init(&ltp);

	struct Fp12 Miller_X,Miller_3;
	Fp12_init(&Miller_X);
	Fp12_init(&Miller_3);

	struct Fp12 t_ans;
	Fp12_init(&t_ans);
	Miller_algo(&Miller_X,G2,G1,X);

	mpz_t set_3;
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	Miller_algo(&Miller_3,G2,G1,set_3);
	Fp12_pow(&Miller_3,&Miller_3,p);

	EFp12_SCM(&zQ,G2,X);
	EFp12_SCM(&p3Q,G2,p3);

	ltp_q(&ltp,&zQ,&p3Q,G1);

	Fp12_mul(&t_ans,&Miller_X,&Miller_3);
	Fp12_mul(&t_ans,&t_ans,&ltp);

	Final_Exp(ANS,&t_ans);

}
void ltt_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *Q){
	struct Fp12 tmp1,tmp2,tmp3,lambda,ltt;
	Fp12_init(&tmp1);
	Fp12_init(&tmp2);
	Fp12_init(&tmp3);
	Fp12_init(&lambda);
	Fp12_init(&ltt);

	Fp12_mul(&tmp1,&T->x,&T->x);//xt^2
	Fp12_add(&tmp2,&tmp1,&tmp1);
	Fp12_add(&tmp1,&tmp2,&tmp1);//3xt^3
	Fp12_add(&tmp2,&T->y,&T->y);//2yt

	Fp12_div(&lambda,&tmp1,&tmp2);//lambda=3xt^2/2yt
	Fp12_sub(&tmp3,&Q->x,&T->x);//tmp3=xq-xt
	Fp12_mul(&tmp3,&tmp3,&lambda);//tmp3=lambda(xq-xt)

	Fp12_sub(&ltt,&Q->y,&T->y);//yq-yt
	Fp12_sub(&ltt,&ltt,&tmp3);//ltt=yq-yt-lambda(xq-xt)

	Fp12_set(ANS,&ltt);

	Fp12_clear(&tmp1);
	Fp12_clear(&tmp2);
	Fp12_clear(&tmp3);
	Fp12_clear(&lambda);
	Fp12_clear(&ltt);
}
void v2t_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *Q){
	struct Fp12 v2t;
	Fp12_init(&v2t);

	Fp12_sub(&v2t,&Q->x,&T->x);//v2t=xq-xt
	Fp12_set(ANS,&v2t);

	Fp12_clear(&v2t);
}
void ltp_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *P,struct EFp12 *Q){
	struct Fp12 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp12_init(&tmp1);
	Fp12_init(&tmp2);
	Fp12_init(&tmp3);
	Fp12_init(&tmp4);
	Fp12_init(&lambda);
	Fp12_init(&ltp);

	if((Fp12_cmp(&T->x,&P->x))==0&&(Fp12_cmp(&T->y,&P->y))!=0){//xt==xp&&yt!=yp
		Fp12_sub(&ltp,&Q->x,&T->x);
		Fp12_set(ANS,&ltp);

		return;
	}

	Fp12_sub(&tmp1,&T->x,&P->x);//xt-xp
	Fp12_sub(&tmp2,&T->y,&P->y);//yt-yp
	Fp12_div(&lambda,&tmp2,&tmp1);//lambda=(yt-tp)/(xt-xp)

	Fp12_sub(&tmp3,&Q->x,&T->x);//tmp3=(xq-xt)
	Fp12_mul(&tmp4,&tmp3,&lambda);//tmp4=lambda(xq-xt)

	Fp12_sub(&ltp,&Q->y,&T->y);//ltp=yq-yt
	Fp12_sub(&ltp,&ltp,&tmp4);//ltp=yq-yt-lambda(xq-xt)

	Fp12_set(ANS,&ltp);

	Fp12_clear(&tmp1);
	Fp12_clear(&tmp2);
	Fp12_clear(&tmp3);
	Fp12_clear(&tmp4);
	Fp12_clear(&lambda);
	Fp12_clear(&ltp);
}
void vtp_q(struct Fp12 *ANS,struct EFp12 *T,struct EFp12 *Q){
	struct Fp12 vtp;
	Fp12_init(&vtp);
	if(T->PaI==1){//if T is infity
		Fp12_set_ui(ANS,0);
		Fp_set_ui(&ANS->x0.x0.x0,1);
		return;
	}

	Fp12_sub(&vtp,&Q->x,&T->x);
	Fp12_set(ANS,&vtp);

	Fp12_clear(&vtp);
}
void check_Pairing(void){
	struct EFp tmp;
	EFp_init(&tmp);

	struct EFp12 P,Q,R,S;
	EFp12_init(&P);
	EFp12_init(&Q);
	EFp12_init(&R);
	EFp12_init(&S);

	struct Fp12 ans,tmp1,tmp2,tmp3;
	Fp12_init(&ans);
	Fp12_init(&tmp1);
	Fp12_init(&tmp2);
	Fp12_init(&tmp3);

	mpz_t a,b,ab;
	mpz_init(a);
	mpz_init(b);
	mpz_init(ab);

	mpz_set_ui(a,31);
	mpz_set_ui(b,11);
	mpz_mul(ab,a,b);
	//----------Tate Pairing-----------------------------------------
	// EFp_random_set(&tmp);
	// EFp12_set_EFp(&P,&tmp);
	//
	// EFp12_random_set(&Q);
	// printf("G1=\n");
	//
	// Tate_Pairing(&tmp1,&P,&Q);
	// printf("G1=");
	// EFp12_printf(&P);
	// printf("G2=");
	// EFp12_printf(&Q);
	//
	// Fp12_pow(&tmp1,&tmp1,ab);
	// printf("f^ab=");
	// Fp12_printf(&tmp1);
	//
	// EFp12_SCM(&R,&P,a);
	// EFp12_SCM(&S,&Q,b);
	//
	// Tate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp12_printf(&tmp2);
	//----------Ate Pairing------------------------------------------
	EFp_random_set(&tmp);
	EFp12_set_EFp(&P,&tmp);

	EFp12_random_set_for_Ate(&Q);
	printf("G1=");
	EFp12_printf(&P);
	printf("G2=");
	EFp12_printf(&Q);

	Ate_Pairing(&tmp1,&P,&Q);

	Fp12_pow(&tmp1,&tmp1,ab);
	printf("f^ab=");
	Fp12_printf(&tmp1);

	EFp12_SCM(&R,&P,a);
	EFp12_SCM(&S,&Q,b);

	Ate_Pairing(&tmp2,&R,&S);

	printf("f'  =");
	Fp12_printf(&tmp2);
	//----------Optimal Ate Pairing----------------------------------
	// EFp_random_set(&tmp);
	// EFp12_set_EFp(&P,&tmp);
	//
	// EFp12_random_set_for_Ate(&Q);
	//
	// Optimal_Ate_Pairing(&tmp1,&P,&Q);
	// printf("G1=");
	// EFp12_printf(&P);
	// printf("G2=");
	// EFp12_printf(&Q);
	//
	//
	// Fp12_pow(&tmp1,&tmp1,ab);
	// printf("f^ab=");
	// Fp12_printf(&tmp1);
	//
	// EFp12_SCM(&R,&P,a);
	// EFp12_SCM(&S,&Q,b);
	//
	// Optimal_Ate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp12_printf(&tmp2);
	//----------------------------------------------------


	EFp12_clear(&Q);
	Fp12_clear(&ans);
	Fp12_clear(&tmp1);
	Fp12_clear(&tmp2);
	Fp12_clear(&tmp3);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(ab);
}
