#include"embedding_degree18.h"
// #include"Finite_Field.c"
#include"Finite_Field.c"
#include"Elliptic_Curve.c"

int main(void){
	mpz_init(X);
	mpz_init(prime);
	mpz_init(order);
	mpz_init(order_EFp);
	mpz_init(trace);
	mpz_init(b);

	generate_X();

	EFp_set_EC_parameter();//generate p,r,t,b

	gmp_printf("p=%Zd\n",prime);
	gmp_printf("r=%Zd\n",order);
	gmp_printf("t=%Zd\n",trace);
	gmp_printf("#E(Fp)=%Zd\n",order_EFp);

	printf("p = %dbit\n",(int)mpz_sizeinbase(prime,2));
	printf("r = %dbit\n",(int)mpz_sizeinbase(order,2));
	printf("t = %dbit\n",(int)mpz_sizeinbase(trace,2));
	printf("X = %dbit\n",(int)mpz_sizeinbase(X,2));
	gmp_printf("b = %Zd\n",b);


	// struct EFp18 A,B;
	// EFp18_init(&A);
	// EFp18_init(&B);
	// EFp18_random_set(&A);

	check_Pairing();//pairing check
	// Masure_pairing_time();
	mpz_clear(X);
	mpz_clear(prime);
	mpz_clear(order);
	mpz_clear(order_EFp);
	mpz_clear(trace);
	mpz_clear(b);
	return 0;
}
void generate_X(){
	//set generater X
	//----------------------------------------------------------------------
	//c=3
	// mpz_set_str(X,"17592190238210",10);//348bit
	// X_bit_binary[44]=1;
	// X_bit_binary[22]=1;
	// X_bit_binary[9]=-1;
	// X_bit_binary[1]=1;

	//c=2
	// mpz_set_str(X,"661934",10);
	// X_bit_binary[19]=1;
	// X_bit_binary[17]=1;
	// X_bit_binary[12]=1;
	// X_bit_binary[11]=1;
	// X_bit_binary[8]=1;
	// X_bit_binary[7]=1;
	// X_bit_binary[5]=1;
	// X_bit_binary[3]=1;
	// X_bit_binary[2]=1;
	// X_bit_binary[1]=1;
	// mpz_set_str(X,"-17652315455488",10);//348bit
	// X_bit_binary[44]=-1;
	// X_bit_binary[36]=-1;
	// X_bit_binary[33]=1;
	// X_bit_binary[17]=1;
	// mpz_set_str(X,"-18448925504779055104",10);//508bit
	X_bit_binary[64]=-1;
	X_bit_binary[51]=-1;
	X_bit_binary[46]=1;
	X_bit_binary[12]=1;

	mpz_t tmp,set_2;
	mpz_init(tmp);
	mpz_init(set_2);
	mpz_set_ui(set_2,2);

	int i;
	for(i=x_bit;i>=0;i--){
		printf("%d",X_bit_binary[i]);
		if(X_bit_binary[i]==1){
			mpz_pow_ui(tmp,set_2,i);
			mpz_add(X,X,tmp);
		}else if(X_bit_binary[i]==-1){
			mpz_pow_ui(tmp,set_2,i);
			mpz_sub(X,X,tmp);
		}
	}
	printf("\n");
	mpz_out_str(stdout,2,X);
	printf("\n");
	return;
}
//-----------------------------------------------------------------------------------------

void EFp_set_EC_parameter(void){
	//generate p,r,t
	mpz_t p_tmp,r_tmp,t_tmp;
	mpz_t xpow2,xpow3,xpow4,xpow5,xpow6,xpow7,xpow8;
	mpz_t tmp1,tmp2;

	mpz_init(p_tmp);
	mpz_init(r_tmp);
	mpz_init(t_tmp);

	mpz_init(xpow2);
	mpz_init(xpow3);
	mpz_init(xpow4);
	mpz_init(xpow5);
	mpz_init(xpow6);
	mpz_init(xpow7);
	mpz_init(xpow8);

	mpz_init(tmp1);
	mpz_init(tmp2);

	mpz_mul(xpow2,X,X);
	mpz_mul(xpow3,xpow2,X);
	mpz_mul(xpow4,xpow2,xpow2);
	mpz_mul(xpow5,xpow4,X);
	mpz_mul(xpow6,xpow3,xpow3);
	mpz_mul(xpow7,xpow6,X);
	mpz_mul(xpow8,xpow4,xpow4);

	//t=1/7(x^4+16x+7)
	//-----------------------------------------------------------------------------------------
	mpz_mul_ui(tmp1,X,16);
	mpz_add_ui(tmp2,xpow4,7);
	mpz_add(t_tmp,tmp1,tmp2);

	mpz_div_ui(trace,t_tmp,7);
	//-----------------------------------------------------------------------------------------
	//r=x^6+37x^3+343
	mpz_mul_ui(tmp1,xpow3,37);
	mpz_add_ui(tmp2,xpow6,343);
	mpz_add(r_tmp,tmp1,tmp2);

	mpz_div_ui(order,r_tmp,343);
	//-----------------------------------------------------------------------------------------
	//p=1/21(x^8+5x^7+7x^6+37x^5+188x^4+259x^3+343x^2+1763x+2401)
	mpz_mul_ui(tmp1,xpow7,5);
	mpz_add(p_tmp,tmp1,xpow8);
	mpz_mul_ui(tmp1,xpow6,7);
	mpz_add(p_tmp,tmp1,p_tmp);
	mpz_mul_ui(tmp1,xpow5,37);
	mpz_add(p_tmp,tmp1,p_tmp);
	mpz_mul_ui(tmp1,xpow4,188);
	mpz_add(p_tmp,tmp1,p_tmp);
	mpz_mul_ui(tmp1,xpow3,259);
	mpz_add(p_tmp,tmp1,p_tmp);
	mpz_mul_ui(tmp1,xpow2,343);
	mpz_add(p_tmp,tmp1,p_tmp);
	mpz_mul_ui(tmp1,X,1763);
	mpz_add(p_tmp,tmp1,p_tmp);
	mpz_add_ui(p_tmp,p_tmp,2401);

	mpz_div_ui(prime,p_tmp,21);

	mpz_add_ui(order_EFp,prime,1);
	mpz_sub(order_EFp,order_EFp,trace);
	//-----------------------------------------------------------------------------------------

	if(mpz_probab_prime_p(prime,25)==0){
		gmp_printf("p:%Zd\n",prime);
		printf("not  prime number!\n");
		exit(0);
	}


	// set curve paramater y^2=x^3+b
	struct EFp P,ANS;
	int legendle;
	struct Fp rational_point,x;
	mpz_t tmp_b;
	Fp_init(&rational_point);
	EFp_init(&P);
	EFp_init(&ANS);
	Fp_init(&x);
	mpz_init(tmp_b);
	mpz_set_si(tmp_b,0);

	for(;;){
		mpz_add_ui(tmp_b,tmp_b,1);
		Fp_set_ui(&x,1);
		legendle=0;
		while(legendle!=1){
			mpz_powm_ui(rational_point.x0,x.x0,3,prime);
			mpz_add(rational_point.x0,rational_point.x0,tmp_b);
			if((legendle=mpz_legendre(rational_point.x0,prime))==1){
				Fp_sqrt(&P.y,&rational_point);
				Fp_set(&P.x,&x);
				EFp_SCM(&ANS,&P,order_EFp);
				if(ANS.infity==TRUE){
					mpz_set(b,tmp_b);
					mpz_clear(tmp_b);
					Fp_clear(&rational_point);
					Fp_clear(&x);
					EFp_clear(&P);
					EFp_clear(&ANS);
					return;
				}
			}
			Fp_add_ui(&x,&x,1);
		}
	}
	return;
}

//-----------------------------------------------------------------------------------------
void Miller_algo(struct Fp18 *ANS,struct EFp18 *P,struct EFp18 *Q,mpz_t loop){
	struct Fp18 l_sum,v_sum;
	Fp18_init(&l_sum);
	Fp18_init(&v_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);
	Fp_set_ui(&v_sum.x0.x0.x0,1);

	struct EFp18 T;
	EFp18_init(&T);
	EFp18_set(&T,P);


	struct Fp18 ltt,ltp,v2t,vtp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);
	Fp18_init(&v2t);
	Fp18_init(&vtp);

	int i;
	struct Fp18 tmp1,lambda;
	Fp18_init(&tmp1);
	Fp18_init(&lambda);
	int r_bit;//bit数

	r_bit= (int)mpz_sizeinbase(loop,2);

	for(i=r_bit-2;i>=0;i--){
		if(mpz_tstbit(loop,i)==1){
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Fp18_mul(&v_sum,&v_sum,&v_sum);

			ltt_q(&ltt,&T,Q);
			Fp18_mul(&l_sum,&l_sum,&ltt);

			EFp18_ECD(&T,&T);
			v2t_q(&v2t,&T,Q);
			Fp18_mul(&v_sum,&v_sum,&v2t);

			ltp_q(&ltp,&T,P,Q);
			Fp18_mul(&l_sum,&l_sum,&ltp);

			EFp18_ECA(&T,&T,P);
			vtp_q(&vtp,&T,Q);
			Fp18_mul(&v_sum,&v_sum,&vtp);
		}else{
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Fp18_mul(&v_sum,&v_sum,&v_sum);

			ltt_q(&ltt,&T,Q);
			Fp18_mul(&l_sum,&l_sum,&ltt);

			EFp18_ECD(&T,&T);
			v2t_q(&v2t,&T,Q);
			Fp18_mul(&v_sum,&v_sum,&v2t);
		}
	}


	// EFp18_printf(&T);

	Fp18_div(ANS,&l_sum,&v_sum);
	// Fp18_set(ANS,&l_sum);

	Fp18_clear(&l_sum);
	Fp18_clear(&v_sum);
	EFp18_clear(&T);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
	Fp18_clear(&v2t);
	Fp18_clear(&vtp);
	Fp18_clear(&tmp1);
	Fp18_clear(&lambda);
}
void Optimal_Miller(struct Fp18 *ANS,struct EFp18 *P,struct EFp18 *Q,mpz_t loop){
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp18 T,EFp_tmp;
	EFp18_init(&T);
	EFp18_init(&EFp_tmp);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	struct Fp18 Px_neg;
	Fp18_init(&Px_neg);

	Fp18_neg(&Px_neg,&P->x);

	int i;
	struct Fp18 tmp1,lambda;
	Fp18_init(&tmp1);
	Fp18_init(&lambda);

	struct EFp18 Q_neg;
	EFp18_init(&Q_neg);
	Fp18_neg(&Q_neg.y,&Q->y);
	Fp18_set(&Q_neg.x,&Q->x);
	if(X_bit_binary[x_bit]==-1){
		EFp18_set(&T,&Q_neg);
	}else{
		EFp18_set(&T,Q);
	}

	for(i=x_bit-1;i>=0;i--){
		switch (X_bit_binary[i]){
			case 0:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			DBL_LINE(&ltt,&T,&T,P,&Px_neg);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			break;
			case 1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			Fp18_mul(&l_sum,&l_sum,&ltp);
			break;
			case -1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			Fp18_mul(&l_sum,&l_sum,&ltp);
			break;
		}
	}

	mpz_mul_ui(p3,prime,3);
	EFp18_SCM(&EFp_tmp,Q,p3);

	ltp_q(&ltp,&T,&EFp_tmp,P);
	Fp18_mul(&l_sum,&l_sum,&ltp);

	// Fp18_printf(&ltp);
	Fp18_set(ANS,&l_sum);

	// Fp18_random(&l_sum);
	// struct timeval opt_1,opt_2;
	//
	// gettimeofday(&opt_1, NULL);
	// ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
	// gettimeofday(&opt_2, NULL);
	// opt_add+=((double)(opt_2.tv_sec - opt_1.tv_sec)+ (double)(opt_2.tv_usec-opt_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&opt_1, NULL);
	// DBL_LINE(&ltp,&T,&T,P,&Px_neg);
	// gettimeofday(&opt_2, NULL);
	// opt_dbl+=((double)(opt_2.tv_sec - opt_1.tv_sec)+ (double)(opt_2.tv_usec-opt_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&opt_1, NULL);
	// Fp18_mul(&l_sum,&l_sum,&ltp);
	// gettimeofday(&opt_2, NULL);
	// opt_mul+=((double)(opt_2.tv_sec - opt_1.tv_sec)+ (double)(opt_2.tv_usec-opt_1.tv_usec)*1.0E-6);

	// EFp18_clear(&Q_neg);
	Fp18_clear(&l_sum);
	EFp18_clear(&T);
	EFp18_clear(&EFp_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
	Fp18_clear(&tmp1);
	Fp18_clear(&lambda);
}
void Tate_Pairing(struct Fp18 *ANS,struct EFp18 *G1,struct EFp18 *G2){
	struct Fp18 t_ans;
	Fp18_init(&t_ans);

	Miller_algo(&t_ans,G1,G2,order_EFp);
	Final_Exp(&t_ans,&t_ans);
	Fp18_set(ANS,&t_ans);

	Fp18_clear(&t_ans);
}
void Ate_Pairing(struct Fp18 *ANS,struct EFp18 *G1,struct EFp18 *G2){
	struct Fp18 t_ans;
	Fp18_init(&t_ans);

	mpz_t tm1;
	mpz_init(tm1);
	mpz_sub_ui(tm1,trace,1);

	Miller_algo(&t_ans,G2,G1,tm1);
	Final_Exp(&t_ans,&t_ans);
	Fp18_set(ANS,&t_ans);

	Fp18_clear(&t_ans);
}
void Optimal_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2){
	mpz_t p3;
	mpz_init(p3);
	struct EFp18 zQ,EFp18_G1;
	EFp18_init(&zQ);
	EFp18_init(&EFp18_G1);

	struct Fp18 ltp;
	Fp18_init(&ltp);

	struct Fp18 Miller_X,Miller_3,t_ans;
	Fp18_init(&Miller_X);
	Fp18_init(&Miller_3);
	Fp18_init(&t_ans);

	mpz_t set_3;
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	EFp18_set_EFp(&EFp18_G1,G1);

	Optimal_Miller(&Miller_X,&EFp18_G1,G2,X);

	Miller_algo(&Miller_3,G2,&EFp18_G1,set_3);
	Fp18_pow(&Miller_3,&Miller_3,prime);

	Fp18_mul(&t_ans,&Miller_X,&Miller_3);

	Final_Exp(ANS,&t_ans);

	mpz_clear(p3);
	EFp18_clear(&zQ);
	Fp18_clear(&ltp);
	Fp18_clear(&Miller_X);
	Fp18_clear(&Miller_3);
	Fp18_clear(&t_ans);
	mpz_clear(set_3);
}
void ltt_q(struct Fp18 *ANS,struct EFp18 *T,struct EFp18 *Q){//T:G1,Q:G2
	struct Fp18 tmp1,tmp2,tmp3,lambda,ltt;
	Fp18_init(&tmp1);
	Fp18_init(&tmp2);
	Fp18_init(&tmp3);
	Fp18_init(&lambda);
	Fp18_init(&ltt);

	Fp18_mul(&tmp1,&T->x,&T->x);//xt^2
	Fp18_add(&tmp2,&tmp1,&tmp1);
	Fp18_add(&tmp1,&tmp2,&tmp1);//3xt^3
	Fp18_add(&tmp2,&T->y,&T->y);//2yt
	Fp18_div(&lambda,&tmp1,&tmp2);//lambda=3xt^2/2yt

	Fp18_sub(&tmp3,&Q->x,&T->x);//tmp3=xq-xt
	Fp18_mul(&tmp3,&tmp3,&lambda);//tmp3=lambda(xq-xt)

	Fp18_sub(&ltt,&Q->y,&T->y);//yq-yt
	Fp18_sub(&ltt,&ltt,&tmp3);//ltt=yq-yt-lambda(xq-xt)

	Fp18_set(ANS,&ltt);

	Fp18_clear(&tmp1);
	Fp18_clear(&tmp2);
	Fp18_clear(&tmp3);
	Fp18_clear(&lambda);
	Fp18_clear(&ltt);
}
void v2t_q(struct Fp18 *ANS,struct EFp18 *T,struct EFp18 *Q){//T:G1,Q:G2
	struct Fp18 v2t;
	Fp18_init(&v2t);

	Fp18_sub(&v2t,&Q->x,&T->x);//v2t=xq-xt
	Fp18_set(ANS,&v2t);

	Fp18_clear(&v2t);
}
void ltp_q(struct Fp18 *ANS,struct EFp18 *T,struct EFp18 *P,struct EFp18 *Q){//T:G1,P:G1,Q:G2
	struct Fp18 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp18_init(&tmp1);
	Fp18_init(&tmp2);
	Fp18_init(&tmp3);
	Fp18_init(&tmp4);
	Fp18_init(&lambda);
	Fp18_init(&ltp);

	if((Fp18_cmp(&T->x,&P->x))==0&&(Fp18_cmp(&T->y,&P->y))!=0){//xt==xp&&yt!=yp
		Fp18_sub(&ltp,&Q->x,&T->x);
		Fp18_set(ANS,&ltp);

		return;
	}

	Fp18_sub(&tmp1,&T->x,&P->x);//xt-xp
	Fp18_sub(&tmp2,&T->y,&P->y);//yt-yp
	Fp18_div(&lambda,&tmp2,&tmp1);//lambda=(yt-tp)/(xt-xp)

	Fp18_sub(&tmp3,&Q->x,&T->x);//tmp3=(xq-xt)
	Fp18_mul(&tmp4,&tmp3,&lambda);//tmp4=lambda(xq-xt)

	Fp18_sub(&ltp,&Q->y,&T->y);//ltp=yq-yt
	Fp18_sub(&ltp,&ltp,&tmp4);//ltp=yq-yt-lambda(xq-xt)
	// Fp18_printf(&ltp);
	Fp18_set(ANS,&ltp);

	Fp18_clear(&tmp1);
	Fp18_clear(&tmp2);
	Fp18_clear(&tmp3);
	Fp18_clear(&tmp4);
	Fp18_clear(&lambda);
	Fp18_clear(&ltp);
}
void vtp_q(struct Fp18 *ANS,struct EFp18 *T,struct EFp18 *Q){//T:G1,Q:G2
	struct Fp18 vtp;
	Fp18_init(&vtp);
	if(T->infity==1){//if T is infity
		Fp18_set_ui(ANS,0);
		Fp_set_ui(&ANS->x0.x0.x0,1);
		return;
	}

	Fp18_sub(&vtp,&Q->x,&T->x);
	Fp18_set(ANS,&vtp);

	Fp18_clear(&vtp);
}
void ADD_LINE(struct Fp18 *l_ANS,struct EFp18 *T_ANS,struct EFp18 *T,struct EFp18 *P,struct EFp18 *Q,struct Fp18 *Qx_neg){
	struct Fp18 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp18_init(&tmp1);
	Fp18_init(&tmp2);
	Fp18_init(&tmp3);
	Fp18_init(&tmp4);
	Fp18_init(&lambda);
	Fp18_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp18 x,y,tmp;
	Fp18_init(&x);
	Fp18_init(&y);
	Fp18_init(&tmp);

	struct EFp18 x3_tmp;
	EFp18_init(&x3_tmp);
	struct Fp18 A,B,C,D,E,F;
	Fp18_init(&A);
	Fp18_init(&B);
	Fp18_init(&C);
	Fp18_init(&D);
	Fp18_init(&E);
	Fp18_init(&F);


	Fp18_sub(&A,&P->x,&T->x);//xt-xp
	Fp18_sub(&B,&P->y,&T->y);//yt-yp
	Fp18_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp18_add(&D,&T->x,&P->x);
	Fp18_mul(&tmp1,&C,&C);
	Fp18_sub(&x3_tmp.x,&tmp1,&D);

	Fp18_mul(&tmp2,&C,&T->x);
	Fp18_sub(&E,&tmp2,&T->y);

	Fp18_mul(&tmp3,&C,&x3_tmp.x);
	Fp18_sub(&x3_tmp.y,&E,&tmp3);

	Fp18_set(&l_tmp,&Q->y);

	Fp18_add(&l_tmp,&l_tmp,&E);

	Fp18_mul(&F,&C,Qx_neg);
	Fp18_add(&l_tmp,&l_tmp,&F);

	Fp18_set(l_ANS,&l_tmp);
	EFp18_set(T_ANS,&x3_tmp);

	Fp18_clear(&tmp1);
	Fp18_clear(&tmp2);
	Fp18_clear(&tmp3);
	Fp18_clear(&tmp4);
	Fp18_clear(&lambda);
	Fp18_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp18_clear(&x);
	Fp18_clear(&y);
	Fp18_clear(&tmp);
	EFp18_clear(&x3_tmp);
	Fp18_clear(&A);
	Fp18_clear(&B);
	Fp18_clear(&C);
	Fp18_clear(&D);
	Fp18_clear(&E);
	Fp18_clear(&F);
}
void DBL_LINE(struct Fp18 *l_ANS,struct EFp18 *T_ANS,struct EFp18 *T,struct EFp18 *Q,struct Fp18 *Qx_neg){
	struct Fp18 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp18_init(&tmp1);
	Fp18_init(&tmp2);
	Fp18_init(&tmp3);
	Fp18_init(&tmp4);
	Fp18_init(&lambda);
	Fp18_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp18 x,y,tmp;
	Fp18_init(&x);
	Fp18_init(&y);
	Fp18_init(&tmp);

	struct EFp18 x3_tmp;
	EFp18_init(&x3_tmp);
	struct Fp18 A,B,C,D,E,F;
	Fp18_init(&A);
	Fp18_init(&B);
	Fp18_init(&C);
	Fp18_init(&D);
	Fp18_init(&E);
	Fp18_init(&F);


	Fp18_add(&A,&T->y,&T->y);//xt-xp
	Fp18_mul(&B,&T->x,&T->x);
	Fp18_mul_ui(&B,&B,3);
	Fp18_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp18_add(&D,&T->x,&T->x);
	Fp18_mul(&tmp1,&C,&C);
	Fp18_sub(&x3_tmp.x,&tmp1,&D);

	Fp18_mul(&tmp2,&C,&T->x);
	Fp18_sub(&E,&tmp2,&T->y);

	Fp18_mul(&tmp3,&C,&x3_tmp.x);
	Fp18_sub(&x3_tmp.y,&E,&tmp3);

	Fp18_set(&l_tmp,&Q->y);
	// Fp_set_ui(&l_tmp.x0.x0.x0,1);

	Fp18_add(&l_tmp,&l_tmp,&E);

	Fp18_mul(&F,&C,Qx_neg);
	Fp18_add(&l_tmp,&l_tmp,&F);

	Fp18_set(l_ANS,&l_tmp);
	EFp18_set(T_ANS,&x3_tmp);

	if(T->infity==TRUE){
		EFp18_set(T_ANS,T);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp18_cmp_mpz(&T->y,cmp)==0){//P.y==0
		EFp18_set_infity(T_ANS);
		return;
	}
	Fp18_clear(&tmp1);
	Fp18_clear(&tmp2);
	Fp18_clear(&tmp3);
	Fp18_clear(&tmp4);
	Fp18_clear(&lambda);
	Fp18_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp18_clear(&x);
	Fp18_clear(&y);
	Fp18_clear(&tmp);
	EFp18_clear(&x3_tmp);
	Fp18_clear(&A);
	Fp18_clear(&B);
	Fp18_clear(&C);
	Fp18_clear(&D);
	Fp18_clear(&E);
	Fp18_clear(&F);
	mpz_clear(cmp);
}
//-----------------------------------------------------------------------------------------
void Sparse_Ate_Pairing(struct Fp18 *ANS,struct EFp3 *G1,struct EFp3 *G2){
	struct Fp18 t_ans;
	Fp18_init(&t_ans);

	mpz_t tm1;
	mpz_init(tm1);
	mpz_sub_ui(tm1,trace,1);

	Sparse_type1_Miller(&t_ans,G1,G2,tm1);
	Final_Exp(&t_ans,&t_ans);
	Fp18_set(ANS,&t_ans);

	Fp18_clear(&t_ans);
	mpz_clear(tm1);
}
void Sparse_Optimal_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2){
	// mpz_t p3;
	// mpz_init(p3);
	// struct EFp18 zQ;
	// EFp18_init(&zQ);
	//
	// struct Fp18 ltp;
	// Fp18_init(&ltp);
	//
	// struct Fp18 Miller_X,Miller_3;
	// Fp18_init(&Miller_X);
	// Fp18_init(&Miller_3);
	//
	// struct Fp18 t_ans;
	// Fp18_init(&t_ans);
	//
	// mpz_t set_3;
	// mpz_init(set_3);
	// mpz_set_ui(set_3,3);
	//
	// Sparse_type1_Optimal_Miller(&Miller_X,G1,G2,X);
	//
	// Sparse_type1_Miller(&Miller_3,G1,G2,set_3);
	// Fp18_pow(&Miller_3,&Miller_3,prime);
	//
	// Fp18_mul(&t_ans,&Miller_X,&Miller_3);
	//
	// Final_Exp(ANS,&t_ans);
	//
	// mpz_clear(p3);
	// EFp18_clear(&zQ);
	// Fp18_clear(&ltp);
	// Fp18_clear(&Miller_X);
	// Fp18_clear(&Miller_3);
	// Fp18_clear(&t_ans);
	// mpz_clear(set_3);
	mpz_t p3,set_3;
	mpz_init(p3);
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	struct EFp3 EFp3_G1,EFp3_G2;
	EFp3_init(&EFp3_G1);
	EFp3_init(&EFp3_G2);


	struct Fp18 ltp,Miller_X,Miller_3,t_ans;
	Fp18_init(&ltp);
	Fp18_init(&Miller_X);
	Fp18_init(&Miller_3);
	Fp18_init(&t_ans);
	int i;

	EFp3_set_EFp(&EFp3_G1,G1);
	i=EFp3_set_EFp18_Sparse(&EFp3_G2,G2);
	if(i==1){
		Sparse_type1_Optimal_Miller(&Miller_X,&EFp3_G1,&EFp3_G2,X);
		Sparse_type1_Miller(&Miller_3,&EFp3_G1,&EFp3_G2,set_3);
	}else if(i==2){
		Sparse_type2_Optimal_Miller(&Miller_X,&EFp3_G1,&EFp3_G2,X);
		Sparse_type2_Miller(&Miller_3,&EFp3_G1,&EFp3_G2,set_3);
	}else{
		printf("G2 rational point error\n");
	}

	Fp18_pow(&Miller_3,&Miller_3,prime);
	Fp18_mul(&t_ans,&Miller_X,&Miller_3);

	Final_Exp(ANS,&t_ans);

	mpz_clear(p3);
	Fp18_clear(&ltp);
	Fp18_clear(&Miller_X);
	Fp18_clear(&Miller_3);
	Fp18_clear(&t_ans);
	mpz_clear(set_3);
}
//Sparse type 1 (G2.x.x2.x0,G2.y.x0.x1)
void Sparse_type1_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,EFp_tmp;
	EFp3_init(&T);
	EFp3_init(&EFp_tmp);

	struct Fp3 Px_neg;
	Fp3_init(&Px_neg);

	Fp3_neg(&Px_neg,&P->x);

	mpz_t p3;
	mpz_init(p3);

	EFp3_set(&T,Q);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;
	// EFp3_printf(Q);

	int r_bit;//bit数
	r_bit= (int)mpz_sizeinbase(loop,2);

	for(i=r_bit-2;i>=0;i--){
		if(mpz_tstbit(loop,i)==1){
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			// EFp3_printf(&T);
			// Fp18_printf(&ltt);
			Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
			// EFp3_printf(&T);
			// Fp18_printf(&ltp);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			Fp18_mul(&l_sum,&l_sum,&ltp);
		}else{
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			Fp18_mul(&l_sum,&l_sum,&ltt);
		}
	}
	// EFp3_printf(&T);
	Fp18_set(ANS,&l_sum);

	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&EFp_tmp);
	Fp3_clear(&Px_neg);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
}
void Sparse_type1_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,EFp_tmp;
	EFp3_init(&T);
	EFp3_init(&EFp_tmp);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;

	struct Fp3 Px_neg;
	Fp3_init(&Px_neg);

	Fp3_neg(&Px_neg,&P->x);

	struct EFp3 Q_neg;
	EFp3_init(&Q_neg);
	Fp3_neg(&Q_neg.y,&Q->y);
	Fp3_set(&Q_neg.x,&Q->x);
	if(X_bit_binary[x_bit]==-1){
		EFp3_set(&T,&Q_neg);
	}else{
		EFp3_set(&T,Q);
	}
	for(i=x_bit-1;i>=0;i--){
		switch (X_bit_binary[i]){
			case 0:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);

			Sparse_type1_mul(&l_sum,&l_sum,&ltt);
			break;
			case 1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);

			Sparse_type1_mul(&l_sum,&l_sum,&ltt);
			Sparse_type1_mul(&l_sum,&l_sum,&ltp);
			break;
			case -1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			// EFp3_printf(&T);
			Sparse_type1_ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
			// EFp3_printf(&T);

			Sparse_type1_mul(&l_sum,&l_sum,&ltt);
			Sparse_type1_mul(&l_sum,&l_sum,&ltp);
			break;
		}
	}

	mpz_mul_ui(p3,prime,3);
	EFp3_SCM(&EFp_tmp,Q,p3);
	// EFp3_printf(&EFp_tmp);

	// ltp_q_Sparse_type1(&ltp,&T,&EFp_tmp,P);
	Sparse_type1_ADD_LINE(&ltp,&T,&T,&EFp_tmp,P,&Px_neg);
	Fp18_mul(&l_sum,&l_sum,&ltp);

	Fp18_set(ANS,&l_sum);

	// Fp18_random(&l_sum);
	// struct timeval sps_1,sps_2;
	//
	// gettimeofday(&sps_1, NULL);
	// Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
	// gettimeofday(&sps_2, NULL);
	// sps_add+=((double)(sps_2.tv_sec - sps_1.tv_sec)+ (double)(sps_2.tv_usec-sps_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&sps_1, NULL);
	// Sparse_type1_DBL_LINE(&ltp,&T,&T,P,&Px_neg);
	// gettimeofday(&sps_2, NULL);
	// sps_dbl+=((double)(sps_2.tv_sec - sps_1.tv_sec)+ (double)(sps_2.tv_usec-sps_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&sps_1, NULL);
	// Sparse_type1_mul(&l_sum,&l_sum,&ltp);
	// gettimeofday(&sps_2, NULL);
	// sps_mul+=((double)(sps_2.tv_sec - sps_1.tv_sec)+ (double)(sps_2.tv_usec-sps_1.tv_usec)*1.0E-6);



	EFp3_clear(&Q_neg);
	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&EFp_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
	Fp3_clear(&Px_neg);
}
void Sparse_type1_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *Qx_neg){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E,F;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);
	Fp3_init(&F);


	Fp3_sub(&A,&P->x,&T->x);//xt-xp
	Fp3_sub(&B,&P->y,&T->y);//yt-yp
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&P->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp3_set(&l_tmp.x0.x0,&Q->y);

	Fp3_set(&l_tmp.x0.x1,&E);

	Fp3_mul(&F,&C,Qx_neg);
	Fp3_set(&l_tmp.x1.x0,&F);

	Fp18_set(l_ANS,&l_tmp);
	EFp3_set(T_ANS,&x3_tmp);
	// if((Fp3_cmp(&T->x,&Q->x))==0&&(Fp3_cmp(&T->y,&Q->y))!=0){//xt==xp&&yt!=yp
	// 	Fp3_sub(&ltp,&P->x,&T->x);
	// 	Fp3_set(&l_ANS->x0.x0,&ltp);
	// }
	// if(T->infity==TRUE){//if P2==inf
	// 	EFp3_set(T_ANS,P);
	// 	return;
	// }
	// else if(P->infity==TRUE){//if P1==inf
	// 	EFp3_set(T_ANS,T);
	// 	return;
	// }
	// else if(Fp3_cmp(&T->x,&P->x)==0&&Fp3_cmp(&T->y,&P->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
	// 	EFp3_set_infity(T_ANS);
	// 	return;
	// }
	// else if(EFp3_cmp(T,P)==0){ // P=P
	// 	EFp3_ECD(T_ANS,T);
	// 	return;
	// }

	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
	Fp3_clear(&F);
}
void Sparse_type1_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *Qx_neg){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E,F;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);
	Fp3_init(&F);


	Fp3_add(&A,&T->y,&T->y);//xt-xp
	Fp3_mul(&B,&T->x,&T->x);
	Fp3_mul_ui(&B,&B,3);
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&T->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp3_set(&l_tmp.x0.x0,&Q->y);
	// Fp_set_ui(&l_tmp.x0.x0.x0,1);

	Fp3_set(&l_tmp.x0.x1,&E);

	Fp3_mul(&F,&C,Qx_neg);
	Fp3_set(&l_tmp.x1.x0,&F);

	Fp18_set(l_ANS,&l_tmp);
	EFp3_set(T_ANS,&x3_tmp);

	if(T->infity==TRUE){
		EFp3_set(T_ANS,T);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp3_cmp_mpz(&T->y,cmp)==0){//P.y==0
		EFp3_set_infity(T_ANS);
		return;
	}
	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
	Fp3_clear(&F);
	mpz_clear(cmp);
}
void Sparse_type1_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp6 tmp02,tmp13,tmp45;
	Fp6_init(&tmp02);
	Fp6_init(&tmp13);
	Fp6_init(&tmp45);

	struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&tmp5);
	Fp3_init(&tmp6);

	struct Fp3 theta0,theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8;
	Fp3_init(&theta0);
	Fp3_init(&theta1);
	Fp3_init(&theta2);
	Fp3_init(&theta3);
	Fp3_init(&theta4);
	Fp3_init(&theta5);
	Fp3_init(&theta6);
	Fp3_init(&theta7);
	Fp3_init(&theta8);

	struct Fp18 t_ans,tmp;;
	Fp18_init(&t_ans);
	Fp18_init(&tmp);

	struct Fp6 B_tmp;
	Fp6_init(&B_tmp);

	Fp3_set(&tmp02.x0,&A->x0.x0);
	Fp3_set(&tmp02.x1,&A->x2.x0);
	Fp3_set(&tmp13.x0,&A->x1.x0);
	Fp3_set(&tmp13.x1,&A->x0.x1);
	Fp3_set(&tmp45.x0,&A->x1.x1);
	Fp3_set(&tmp45.x1,&A->x2.x1);

	Fp3_set(&B_tmp.x1,&B->x0.x1);
	Fp3_set(&B_tmp.x0,&B->x1.x0);

	Fp3_mul(&theta1,&tmp02.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta5,&tmp02.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp02.x0,&tmp02.x1);//a+b
	Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta1);//(a+b)(c+d)-ac-bd
	Fp3_sub(&theta3,&tmp3,&theta5);

	Fp3_mul(&theta2,&tmp13.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta6,&tmp13.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp13.x0,&tmp13.x1);//a+b
	// Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta2);
	Fp3_sub(&theta4,&tmp3,&theta6);

	Fp3_mul(&tmp1,&tmp45.x0,&B_tmp.x0);
	Fp3_add(&theta5,&tmp1,&theta5);
	Fp3_mul(&tmp1,&tmp45.x1,&B_tmp.x0);
	Fp3_add(&theta6,&tmp1,&theta6);
	Fp3_mul(&theta7,&tmp45.x0,&B_tmp.x1);
	Fp3_mul(&theta8,&tmp45.x1,&B_tmp.x1);

	Fp3_mul_xi(&theta0,&theta6);
	Fp3_mul_xi(&tmp2,&theta7);
	Fp3_add(&theta1,&theta1,&tmp2);
	Fp3_mul_xi(&tmp3,&theta8);
	Fp3_add(&theta2,&theta2,&tmp3);

	Fp3_set(&t_ans.x0.x0,&theta0);
	Fp3_set(&t_ans.x1.x0,&theta1);
	Fp3_set(&t_ans.x2.x0,&theta2);
	Fp3_set(&t_ans.x0.x1,&theta3);
	Fp3_set(&t_ans.x1.x1,&theta4);
	Fp3_set(&t_ans.x2.x1,&theta5);

	Fp18_mul_Fp(&tmp,A,&B->x0.x0.x0);
	Fp18_add(ANS,&t_ans,&tmp);

	Fp6_clear(&tmp02);
	Fp6_clear(&tmp13);
	Fp6_clear(&tmp45);
	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&tmp5);
	Fp3_clear(&tmp6);
	Fp3_clear(&theta0);
	Fp3_clear(&theta1);
	Fp3_clear(&theta2);
	Fp3_clear(&theta3);
	Fp3_clear(&theta4);
	Fp3_clear(&theta5);
	Fp3_clear(&theta6);
	Fp3_clear(&theta7);
	Fp3_clear(&theta8);
	Fp18_clear(&t_ans);
	Fp6_clear(&B_tmp);
}
//Sparse type 2 (G2.x.x1.x1,G2.y.x0.x1)
void Sparse_type2_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,EFp_tmp;
	EFp3_init(&T);
	EFp3_init(&EFp_tmp);

	struct Fp3 Px_neg;
	Fp3_init(&Px_neg);

	Fp3_neg(&Px_neg,&P->x);

	mpz_t p3;
	mpz_init(p3);

	EFp3_set(&T,Q);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;
	// EFp3_printf(Q);

	int r_bit;//bit数
	r_bit= (int)mpz_sizeinbase(loop,2);

	for(i=r_bit-2;i>=0;i--){
		if(mpz_tstbit(loop,i)==1){
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Sparse_type2_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			// EFp3_printf(&T);
			// Fp18_printf(&ltt);
			Sparse_type2_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
			// EFp3_printf(&T);
			// Fp18_printf(&ltp);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			Fp18_mul(&l_sum,&l_sum,&ltp);
		}else{
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Sparse_type2_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			Fp18_mul(&l_sum,&l_sum,&ltt);
		}
	}
	// EFp3_printf(&T);
	Fp18_set(ANS,&l_sum);

	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&EFp_tmp);
	Fp3_clear(&Px_neg);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
}
void Sparse_type2_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,EFp_tmp;
	EFp3_init(&T);
	EFp3_init(&EFp_tmp);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;

	struct Fp3 Px_neg;
	Fp3_init(&Px_neg);

	Fp3_neg(&Px_neg,&P->x);

	struct EFp3 Q_neg;
	EFp3_init(&Q_neg);
	Fp3_neg(&Q_neg.y,&Q->y);
	Fp3_set(&Q_neg.x,&Q->x);
	if(X_bit_binary[x_bit]==-1){
		EFp3_set(&T,&Q_neg);
	}else{
		EFp3_set(&T,Q);
	}

	for(i=x_bit-1;i>=0;i--){
		switch (X_bit_binary[i]){
			case 0:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Sparse_type2_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			Fp18_mul(&l_sum,&l_sum,&ltt);

			break;
			case 1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Sparse_type2_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			Sparse_type2_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			Fp18_mul(&l_sum,&l_sum,&ltp);
			break;
			case -1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Sparse_type2_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
			Sparse_type2_ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);

			Fp18_mul(&l_sum,&l_sum,&ltt);
			Fp18_mul(&l_sum,&l_sum,&ltp);
			break;
		}
	}

	mpz_mul_ui(p3,prime,3);
	EFp3_type2_SCM(&EFp_tmp,Q,p3);

	Sparse_type2_ADD_LINE(&ltp,&T,&T,&EFp_tmp,P,&Px_neg);
	Fp18_mul(&l_sum,&l_sum,&ltp);

	Fp18_set(ANS,&l_sum);

	// Fp18_random(&l_sum);
	// struct timeval sps_1,sps_2;
	//
	// gettimeofday(&sps_1, NULL);
	// Sparse_type2_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
	// gettimeofday(&sps_2, NULL);
	// sps_add+=((double)(sps_2.tv_sec - sps_1.tv_sec)+ (double)(sps_2.tv_usec-sps_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&sps_1, NULL);
	// Sparse_type2_DBL_LINE(&ltp,&T,&T,P,&Px_neg);
	// gettimeofday(&sps_2, NULL);
	// sps_dbl+=((double)(sps_2.tv_sec - sps_1.tv_sec)+ (double)(sps_2.tv_usec-sps_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&sps_1, NULL);
	// Sparse_type2_mul(&l_sum,&l_sum,&ltp);
	// gettimeofday(&sps_2, NULL);
	// sps_mul+=((double)(sps_2.tv_sec - sps_1.tv_sec)+ (double)(sps_2.tv_usec-sps_1.tv_usec)*1.0E-6);



	EFp3_clear(&Q_neg);
	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&EFp_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
	Fp3_clear(&Px_neg);
}
void Sparse_type2_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *Qx_neg){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E,F;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);
	Fp3_init(&F);


	Fp3_sub(&A,&P->x,&T->x);//xt-xp
	Fp3_sub(&B,&P->y,&T->y);//yt-yp
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&P->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_mul_xi_inv(&tmp1,&tmp1);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp3_set(&l_tmp.x0.x0,&Q->y);

	Fp3_set(&l_tmp.x0.x1,&E);

	Fp3_mul_xi_inv(&F,Qx_neg);
	Fp3_mul(&F,&C,&F);
	Fp3_set(&l_tmp.x2.x1,&F);

	Fp18_set(l_ANS,&l_tmp);
	EFp3_set(T_ANS,&x3_tmp);

	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
	Fp3_clear(&F);
}
void Sparse_type2_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *Qx_neg){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E,F;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);
	Fp3_init(&F);


	Fp3_add(&A,&T->y,&T->y);//xt-xp
	Fp3_mul(&B,&T->x,&T->x);
	Fp3_mul_ui(&B,&B,3);
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&T->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_mul_xi(&tmp1,&tmp1);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_mul_xi(&tmp2,&tmp2);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_mul_xi(&tmp3,&tmp3);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp3_set(&l_tmp.x0.x0,&Q->y);

	Fp3_set(&l_tmp.x0.x1,&E);

	Fp3_mul(&F,&C,Qx_neg);
	Fp3_set(&l_tmp.x2.x1,&F);

	Fp18_set(l_ANS,&l_tmp);
	EFp3_set(T_ANS,&x3_tmp);

	if(T->infity==TRUE){
		EFp3_set(T_ANS,T);
		return;
	}
	mpz_t cmp;
	mpz_init(cmp);
	mpz_set_ui(cmp,0);
	if(Fp3_cmp_mpz(&T->y,cmp)==0){//P.y==0
		EFp3_set_infity(T_ANS);
		return;
	}
	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
	Fp3_clear(&F);
	mpz_clear(cmp);
}
void Sparse_type2_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp6 tmp02,tmp13,tmp45;
	Fp6_init(&tmp02);
	Fp6_init(&tmp13);
	Fp6_init(&tmp45);

	struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&tmp5);
	Fp3_init(&tmp6);

	struct Fp3 theta0,theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8;
	Fp3_init(&theta0);
	Fp3_init(&theta1);
	Fp3_init(&theta2);
	Fp3_init(&theta3);
	Fp3_init(&theta4);
	Fp3_init(&theta5);
	Fp3_init(&theta6);
	Fp3_init(&theta7);
	Fp3_init(&theta8);

	struct Fp18 t_ans,tmp;;
	Fp18_init(&t_ans);
	Fp18_init(&tmp);

	struct Fp6 B_tmp;
	Fp6_init(&B_tmp);

	Fp3_set(&tmp02.x0,&A->x0.x0);
	Fp3_set(&tmp02.x1,&A->x2.x0);
	Fp3_set(&tmp13.x0,&A->x1.x0);
	Fp3_set(&tmp13.x1,&A->x0.x1);
	Fp3_set(&tmp45.x0,&A->x1.x1);
	Fp3_set(&tmp45.x1,&A->x2.x1);

	Fp3_set(&B_tmp.x1,&B->x0.x1);
	Fp3_set(&B_tmp.x0,&B->x1.x0);

	Fp3_mul(&theta1,&tmp02.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta5,&tmp02.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp02.x0,&tmp02.x1);//a+b
	Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta1);//(a+b)(c+d)-ac-bd
	Fp3_sub(&theta3,&tmp3,&theta5);

	Fp3_mul(&theta2,&tmp13.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta6,&tmp13.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp13.x0,&tmp13.x1);//a+b
	// Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta2);
	Fp3_sub(&theta4,&tmp3,&theta6);

	Fp3_mul(&tmp1,&tmp45.x0,&B_tmp.x0);
	Fp3_add(&theta5,&tmp1,&theta5);
	Fp3_mul(&tmp1,&tmp45.x1,&B_tmp.x0);
	Fp3_add(&theta6,&tmp1,&theta6);
	Fp3_mul(&theta7,&tmp45.x0,&B_tmp.x1);
	Fp3_mul(&theta8,&tmp45.x1,&B_tmp.x1);

	Fp3_mul_xi(&theta0,&theta6);
	Fp3_mul_xi(&tmp2,&theta7);
	Fp3_add(&theta1,&theta1,&tmp2);
	Fp3_mul_xi(&tmp3,&theta8);
	Fp3_add(&theta2,&theta2,&tmp3);

	Fp3_set(&t_ans.x0.x0,&theta0);
	Fp3_set(&t_ans.x1.x0,&theta1);
	Fp3_set(&t_ans.x2.x0,&theta2);
	Fp3_set(&t_ans.x0.x1,&theta3);
	Fp3_set(&t_ans.x1.x1,&theta4);
	Fp3_set(&t_ans.x2.x1,&theta5);

	Fp18_mul_Fp(&tmp,A,&B->x0.x0.x0);
	Fp18_add(ANS,&t_ans,&tmp);

	Fp6_clear(&tmp02);
	Fp6_clear(&tmp13);
	Fp6_clear(&tmp45);
	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&tmp5);
	Fp3_clear(&tmp6);
	Fp3_clear(&theta0);
	Fp3_clear(&theta1);
	Fp3_clear(&theta2);
	Fp3_clear(&theta3);
	Fp3_clear(&theta4);
	Fp3_clear(&theta5);
	Fp3_clear(&theta6);
	Fp3_clear(&theta7);
	Fp3_clear(&theta8);
	Fp18_clear(&t_ans);
	Fp6_clear(&B_tmp);
}

//-----------------------------------------------------------------------------------------

void Pseudo_Sparse_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2){
	struct Fp18 t_ans;
	Fp18_init(&t_ans);
	struct EFp3 EFp3_G1,EFp3_G2;
	EFp3_init(&EFp3_G1);
	EFp3_init(&EFp3_G2);
	mpz_t tm1;
	mpz_init(tm1);
	mpz_sub_ui(tm1,trace,1);
	int i;

	EFp3_set_EFp(&EFp3_G1,G1);
	i=EFp3_set_EFp18_Sparse(&EFp3_G2,G2);
	if(i==1){
		Pseudo_type1_Miller(&t_ans,&EFp3_G1,&EFp3_G2,tm1);
	}else if(i==2){
		Pseudo_type2_Miller(&t_ans,&EFp3_G1,&EFp3_G2,tm1);
	}else{
		printf("G2 rational point error\n");
	}


	Final_Exp(&t_ans,&t_ans);
	Fp18_set(ANS,&t_ans);

	Fp18_clear(&t_ans);
	mpz_clear(tm1);
}
void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp18 *ANS,struct EFp *G1,struct EFp18 *G2){
	mpz_t p3,set_3;
	mpz_init(p3);
	mpz_init(set_3);
	mpz_set_ui(set_3,3);

	struct EFp3 EFp3_G1,EFp3_G2;
	EFp3_init(&EFp3_G1);
	EFp3_init(&EFp3_G2);


	struct Fp18 ltp,Miller_X,Miller_3,t_ans;
	Fp18_init(&ltp);
	Fp18_init(&Miller_X);
	Fp18_init(&Miller_3);
	Fp18_init(&t_ans);
	int i;

	EFp3_set_EFp(&EFp3_G1,G1);
	i=EFp3_set_EFp18_Sparse(&EFp3_G2,G2);
	if(i==1){
		Pseudo_type1_Optimal_Miller(&Miller_X,&EFp3_G1,&EFp3_G2,X);
		Pseudo_type1_Miller(&Miller_3,&EFp3_G1,&EFp3_G2,set_3);
	}else if(i==2){
		Pseudo_type2_Optimal_Miller(&Miller_X,&EFp3_G1,&EFp3_G2,X);
		Pseudo_type2_Miller(&Miller_3,&EFp3_G1,&EFp3_G2,set_3);
	}else{
		printf("G2 rational point error\n");
	}

	Fp18_pow(&Miller_3,&Miller_3,prime);
	Fp18_mul(&t_ans,&Miller_X,&Miller_3);

	Final_Exp(ANS,&t_ans);

	mpz_clear(p3);
	Fp18_clear(&ltp);
	Fp18_clear(&Miller_X);
	Fp18_clear(&Miller_3);
	Fp18_clear(&t_ans);
	mpz_clear(set_3);
}
//Sparse type 1 (G2.x.x2.x0,G2.y.x0.x1)
void Pseudo_type1_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){//Q:G2,P:G1
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,P_map,Q_map,EFp3_tmp;
	EFp3_init(&T);
	EFp3_init(&P_map);
	EFp3_init(&Q_map);
	EFp3_init(&EFp3_tmp);

	struct Fp3 L,xy,xy_2,y_inv,tmp,y_tmp;
	Fp3_init(&L);
	Fp3_init(&xy);
	Fp3_init(&xy_2);
	Fp3_init(&y_inv);
	Fp3_init(&tmp);
	Fp3_init(&y_tmp);

	Fp3_invert(&y_inv,&P->y);
	Fp3_mul(&xy,&P->x,&y_inv);

	Fp3_mul(&xy_2,&xy,&xy);
	Fp3_mul(&P_map.x,&xy_2,&P->x);
	Fp3_set(&P_map.y,&P_map.x);

	Fp3_mul(&y_tmp,&xy_2,&xy);
	Fp3_mul(&Q_map.y,&y_tmp,&Q->y);
	Fp3_mul(&Q_map.x,&xy_2,&Q->x);

	EFp3_set(&T,&Q_map);
	Fp3_invert(&L,&P_map.y);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;
	// EFp3_printf(Q);

	int r_bit;//bit数
	r_bit= (int)mpz_sizeinbase(loop,2);

	for(i=r_bit-2;i>=0;i--){
		if(mpz_tstbit(loop,i)==1){
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			// EFp3_printf(&T);
			Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
			// EFp3_printf(&T);

			Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
			Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
		}else{
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);

			Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
		}
	}
	// EFp3_printf(&T);
	Fp18_set(ANS,&l_sum);

	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&P_map);
	EFp3_clear(&Q_map);
	EFp3_clear(&EFp3_tmp);
	Fp3_clear(&L);
	Fp3_clear(&xy);
	Fp3_clear(&xy_2);
	Fp3_clear(&y_inv);
	Fp3_clear(&tmp);
	Fp3_clear(&y_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
}
void Pseudo_type1_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){//Q:G2,P:G1
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,P_map,Q_map,EFp3_tmp;
	EFp3_init(&T);
	EFp3_init(&P_map);
	EFp3_init(&Q_map);
	EFp3_init(&EFp3_tmp);

	struct Fp3 L,xy,xy_2,y_inv,tmp,y_tmp;
	Fp3_init(&L);
	Fp3_init(&xy);
	Fp3_init(&xy_2);
	Fp3_init(&y_inv);
	Fp3_init(&tmp);
	Fp3_init(&y_tmp);

	Fp3_invert(&y_inv,&P->y);
	Fp3_mul(&xy,&P->x,&y_inv);

	Fp3_mul(&xy_2,&xy,&xy);
	Fp3_mul(&P_map.x,&xy_2,&P->x);
	Fp3_set(&P_map.y,&P_map.x);

	Fp3_mul(&y_tmp,&xy_2,&xy);
	Fp3_mul(&Q_map.y,&y_tmp,&Q->y);
	Fp3_mul(&Q_map.x,&xy_2,&Q->x);

	Fp3_invert(&L,&P_map.y);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;

	struct EFp3 Q_neg;
	EFp3_init(&Q_neg);
	Fp3_neg(&Q_neg.y,&Q_map.y);
	Fp3_set(&Q_neg.x,&Q_map.x);
	if(X_bit_binary[x_bit]==-1){
		EFp3_set(&T,&Q_neg);
	}else{
		EFp3_set(&T,&Q_map);
	}
	for(i=x_bit-1;i>=0;i--){
		switch (X_bit_binary[i]){
			case 0:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			// Fp18_printf(&l_sum);
			Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
			// Fp18_printf(&l_sum);
			break;
			case 1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			// Fp18_printf(&l_sum);

			Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);

			Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
			// Fp18_printf(&l_sum);
			Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
			// Fp18_printf(&l_sum);
			break;
			case -1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			// Fp18_printf(&l_sum);

			Pseudo_type1_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_neg,&P_map,&L);

			Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
			// Fp18_printf(&l_sum);
			Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
			// Fp18_printf(&l_sum);
			break;
		}
	}

	mpz_mul_ui(p3,prime,3);
	EFp3_SCM(&EFp3_tmp,&Q_map,p3);

	Pseudo_type1_ADD_LINE(&ltp,&T,&T,&EFp3_tmp,&P_map,&L);
	Fp18_mul(&l_sum,&l_sum,&ltp);

	Fp18_set(ANS,&l_sum);



	// Fp18_random(&l_sum);
	// struct timeval pse_1,pse_2;
	//
	// gettimeofday(&pse_1, NULL);
	// Pseudo_type1_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
	// gettimeofday(&pse_2, NULL);
	// pse_add+=((double)(pse_2.tv_sec - pse_1.tv_sec)+ (double)(pse_2.tv_usec-pse_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&pse_1, NULL);
	// Sparse_type1_DBL_LINE(&ltp,&T,&T,&P_map,&L);
	// gettimeofday(&pse_2, NULL);
	// pse_dbl+=((double)(pse_2.tv_sec - pse_1.tv_sec)+ (double)(pse_2.tv_usec-pse_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&pse_1, NULL);
	// Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
	// gettimeofday(&pse_2, NULL);
	// pse_mul+=((double)(pse_2.tv_sec - pse_1.tv_sec)+ (double)(pse_2.tv_usec-pse_1.tv_usec)*1.0E-6);



	EFp3_clear(&Q_neg);
	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&P_map);
	EFp3_clear(&Q_map);
	EFp3_clear(&EFp3_tmp);
	Fp3_clear(&L);
	Fp3_clear(&xy);
	Fp3_clear(&xy_2);
	Fp3_clear(&y_inv);
	Fp3_clear(&tmp);
	Fp3_clear(&y_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
}
void Pseudo_type1_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *L){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);


	Fp3_sub(&A,&P->x,&T->x);//xt-xp
	Fp3_sub(&B,&P->y,&T->y);//yt-yp
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&P->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp_set_ui(&l_tmp.x0.x0.x0,1);

	Fp3_mul(&l_tmp.x0.x1,&E,L);

	Fp3_neg(&l_tmp.x1.x0,&C);

	Fp18_set(l_ANS,&l_tmp);
	EFp3_set(T_ANS,&x3_tmp);
	// if((Fp3_cmp(&T->x,&Q->x))==0&&(Fp3_cmp(&T->y,&Q->y))!=0){//xt==xp&&yt!=yp
	// 	Fp3_sub(&ltp,&P->x,&T->x);
	// 	Fp3_set(&l_ANS->x0.x0,&ltp);
	// }
	// if(T->infity==TRUE){//if P2==inf
	// 	EFp3_set(T_ANS,P);
	// 	return;
	// }
	// else if(P->infity==TRUE){//if P1==inf
	// 	EFp3_set(T_ANS,T);
	// 	return;
	// }
	// else if(Fp3_cmp(&T->x,&P->x)==0&&Fp3_cmp(&T->y,&P->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
	// 	EFp3_set_infity(T_ANS);
	// 	return;
	// }
	// else if(EFp3_cmp(T,P)==0){ // P=P
	// 	EFp3_ECD(T_ANS,T);
	// 	return;
	// }


	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
}
void Pseudo_type1_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *L){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltt;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltt);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);

	Fp3_add(&A,&T->y,&T->y);//xt-xp
	Fp3_mul(&B,&T->x,&T->x);
	Fp3_mul_ui(&B,&B,3);
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&T->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp_set_ui(&l_tmp.x0.x0.x0,1);

	Fp3_mul(&l_tmp.x0.x1,&E,L);

	Fp3_neg(&l_tmp.x1.x0,&C);
	Fp18_set(l_ANS,&l_tmp);

	EFp3_set(T_ANS,&x3_tmp);

	// if(T->infity==TRUE){
	// 	EFp3_set(T_ANS,T);
	// 	return;
	// }
	// mpz_t cmp;
	// mpz_init(cmp);
	// mpz_set_ui(cmp,0);
	// if(Fp3_cmp_mpz(&T->y,cmp)==0){//P.y==0
	// 	EFp3_set_infity(T_ANS);
	// 	return;
	// }

	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltt);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
}
void Pseudo_type1_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp6 tmp02,tmp13,tmp45;
	Fp6_init(&tmp02);
	Fp6_init(&tmp13);
	Fp6_init(&tmp45);

	Fp3_set(&tmp02.x0,&A->x0.x0);
	Fp3_set(&tmp02.x1,&A->x2.x0);
	Fp3_set(&tmp13.x0,&A->x1.x0);
	Fp3_set(&tmp13.x1,&A->x0.x1);
	Fp3_set(&tmp45.x0,&A->x1.x1);
	Fp3_set(&tmp45.x1,&A->x2.x1);

	struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&tmp5);
	Fp3_init(&tmp6);

	struct Fp3 theta0,theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8;
	Fp3_init(&theta0);
	Fp3_init(&theta1);
	Fp3_init(&theta2);
	Fp3_init(&theta3);
	Fp3_init(&theta4);
	Fp3_init(&theta5);
	Fp3_init(&theta6);
	Fp3_init(&theta7);
	Fp3_init(&theta8);

	struct Fp18 t_ans;
	Fp18_init(&t_ans);

	struct Fp6 B_tmp;
	Fp6_init(&B_tmp);
	Fp3_set(&B_tmp.x1,&B->x0.x1);
	Fp3_set(&B_tmp.x0,&B->x1.x0);

	Fp3_mul(&theta1,&tmp02.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta5,&tmp02.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp02.x0,&tmp02.x1);//a+b
	Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta1);//(a+b)(c+d)-ac-bd
	Fp3_sub(&theta3,&tmp3,&theta5);

	Fp3_mul(&theta2,&tmp13.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta6,&tmp13.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp13.x0,&tmp13.x1);//a+b
	// Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta2);
	Fp3_sub(&theta4,&tmp3,&theta6);

	Fp3_mul(&tmp1,&tmp45.x0,&B_tmp.x0);
	Fp3_add(&theta5,&tmp1,&theta5);
	Fp3_mul(&tmp1,&tmp45.x1,&B_tmp.x0);
	Fp3_add(&theta6,&tmp1,&theta6);
	Fp3_mul(&theta7,&tmp45.x0,&B_tmp.x1);
	Fp3_mul(&theta8,&tmp45.x1,&B_tmp.x1);

	Fp3_mul_xi(&theta0,&theta6);
	Fp3_mul_xi(&tmp2,&theta7);
	Fp3_add(&theta1,&theta1,&tmp2);
	Fp3_mul_xi(&tmp3,&theta8);
	Fp3_add(&theta2,&theta2,&tmp3);

	Fp3_set(&t_ans.x0.x0,&theta0);
	Fp3_set(&t_ans.x1.x0,&theta1);
	Fp3_set(&t_ans.x2.x0,&theta2);
	Fp3_set(&t_ans.x0.x1,&theta3);
	Fp3_set(&t_ans.x1.x1,&theta4);
	Fp3_set(&t_ans.x2.x1,&theta5);

	Fp18_add(ANS,&t_ans,A);

	Fp6_clear(&tmp02);
	Fp6_clear(&tmp13);
	Fp6_clear(&tmp45);
	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&tmp5);
	Fp3_clear(&tmp6);
	Fp3_clear(&theta0);
	Fp3_clear(&theta1);
	Fp3_clear(&theta2);
	Fp3_clear(&theta3);
	Fp3_clear(&theta4);
	Fp3_clear(&theta5);
	Fp3_clear(&theta6);
	Fp3_clear(&theta7);
	Fp3_clear(&theta8);
	Fp18_clear(&t_ans);
	Fp6_clear(&B_tmp);
}
//Sparse type 2 (G2.x.x1.x1,G2.y.x0.x1)
void Pseudo_type2_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){//Q:G2,P:G1
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,P_map,Q_map,EFp3_tmp;
	EFp3_init(&T);
	EFp3_init(&P_map);
	EFp3_init(&Q_map);
	EFp3_init(&EFp3_tmp);

	struct Fp3 L,xy,xy_2,y_inv,tmp,y_tmp;
	Fp3_init(&L);
	Fp3_init(&xy);
	Fp3_init(&xy_2);
	Fp3_init(&y_inv);
	Fp3_init(&tmp);
	Fp3_init(&y_tmp);

	Fp3_invert(&y_inv,&P->y);
	Fp3_mul(&xy,&P->x,&y_inv);

	Fp3_mul(&xy_2,&xy,&xy);
	Fp3_mul(&P_map.x,&xy_2,&P->x);
	Fp3_set(&P_map.y,&P_map.x);

	Fp3_mul(&y_tmp,&xy_2,&xy);
	Fp3_mul(&Q_map.y,&y_tmp,&Q->y);
	Fp3_mul(&Q_map.x,&xy_2,&Q->x);

	EFp3_set(&T,&Q_map);
	Fp3_invert(&L,&P_map.y);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;
	// EFp3_printf(Q);

	int r_bit;//bit数
	r_bit= (int)mpz_sizeinbase(loop,2);

	for(i=r_bit-2;i>=0;i--){
		if(mpz_tstbit(loop,i)==1){
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Pseudo_type2_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			// EFp3_printf(&T);
			Pseudo_type2_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
			// EFp3_printf(&T);

			Pseudo_type2_mul(&l_sum,&l_sum,&ltt);
			Pseudo_type2_mul(&l_sum,&l_sum,&ltp);
		}else{
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Pseudo_type2_DBL_LINE(&ltt,&T,&T,&P_map,&L);

			Pseudo_type2_mul(&l_sum,&l_sum,&ltt);
		}
	}
	// EFp3_printf(&T);
	Fp18_set(ANS,&l_sum);

	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&P_map);
	EFp3_clear(&Q_map);
	EFp3_clear(&EFp3_tmp);
	Fp3_clear(&L);
	Fp3_clear(&xy);
	Fp3_clear(&xy_2);
	Fp3_clear(&y_inv);
	Fp3_clear(&tmp);
	Fp3_clear(&y_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
}
void Pseudo_type2_Optimal_Miller(struct Fp18 *ANS,struct EFp3 *P,struct EFp3 *Q,mpz_t loop){//Q:G2,P:G1
	struct Fp18 l_sum;
	Fp18_init(&l_sum);
	Fp_set_ui(&l_sum.x0.x0.x0,1);

	struct EFp3 T,P_map,Q_map,EFp3_tmp;
	EFp3_init(&T);
	EFp3_init(&P_map);
	EFp3_init(&Q_map);
	EFp3_init(&EFp3_tmp);

	struct Fp3 L,xy,xy_2,y_inv,tmp,y_tmp;
	Fp3_init(&L);
	Fp3_init(&xy);
	Fp3_init(&xy_2);
	Fp3_init(&y_inv);
	Fp3_init(&tmp);
	Fp3_init(&y_tmp);

	Fp3_invert(&y_inv,&P->y);
	Fp3_mul(&xy,&P->x,&y_inv);

	Fp3_mul(&xy_2,&xy,&xy);
	Fp3_mul(&P_map.x,&xy_2,&P->x);
	Fp3_set(&P_map.y,&P_map.x);

	Fp3_mul(&y_tmp,&xy_2,&xy);
	Fp3_mul(&Q_map.y,&y_tmp,&Q->y);
	Fp3_mul(&Q_map.x,&xy_2,&Q->x);

	Fp3_invert(&L,&P_map.y);

	mpz_t p3;
	mpz_init(p3);

	struct Fp18 ltt,ltp;
	Fp18_init(&ltt);
	Fp18_init(&ltp);

	int i;

	struct EFp3 Q_neg;
	EFp3_init(&Q_neg);
	Fp3_neg(&Q_neg.y,&Q_map.y);
	Fp3_set(&Q_neg.x,&Q_map.x);
	if(X_bit_binary[x_bit]==-1){
		EFp3_set(&T,&Q_neg);
	}else{
		EFp3_set(&T,&Q_map);
	}
	for(i=x_bit-1;i>=0;i--){
		switch (X_bit_binary[i]){
			case 0:
			Fp18_mul(&l_sum,&l_sum,&l_sum);
			Pseudo_type2_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			Pseudo_type2_mul(&l_sum,&l_sum,&ltt);
			break;
			case 1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Pseudo_type2_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			Pseudo_type2_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);

			Pseudo_type2_mul(&l_sum,&l_sum,&ltt);
			Pseudo_type2_mul(&l_sum,&l_sum,&ltp);
			break;
			case -1:
			Fp18_mul(&l_sum,&l_sum,&l_sum);

			Pseudo_type2_DBL_LINE(&ltt,&T,&T,&P_map,&L);
			Pseudo_type2_ADD_LINE(&ltp,&T,&T,&Q_neg,&P_map,&L);

			Pseudo_type2_mul(&l_sum,&l_sum,&ltt);
			Pseudo_type2_mul(&l_sum,&l_sum,&ltp);
			break;
		}
	}

	mpz_mul_ui(p3,prime,3);
	EFp3_type2_SCM(&EFp3_tmp,&Q_map,p3);

	Pseudo_type2_ADD_LINE(&ltp,&T,&T,&EFp3_tmp,&P_map,&L);
	Fp18_mul(&l_sum,&l_sum,&ltp);

	Fp18_set(ANS,&l_sum);

	// Fp18_random(&l_sum);
	// struct timeval pse_1,pse_2;
	//
	// gettimeofday(&pse_1, NULL);
	// Pseudo_type2_ADD_LINE(&ltp,&T,&T,&Q_map,&P_map,&L);
	// gettimeofday(&pse_2, NULL);
	// pse_add+=((double)(pse_2.tv_sec - pse_1.tv_sec)+ (double)(pse_2.tv_usec-pse_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&pse_1, NULL);
	// Sparse_DBL_LINE(&ltp,&T,&T,&P_map,&L);
	// gettimeofday(&pse_2, NULL);
	// pse_dbl+=((double)(pse_2.tv_sec - pse_1.tv_sec)+ (double)(pse_2.tv_usec-pse_1.tv_usec)*1.0E-6);
	//
	// gettimeofday(&pse_1, NULL);
	// Pseudo_type2_mul(&l_sum,&l_sum,&ltp);
	// gettimeofday(&pse_2, NULL);
	// pse_mul+=((double)(pse_2.tv_sec - pse_1.tv_sec)+ (double)(pse_2.tv_usec-pse_1.tv_usec)*1.0E-6);

	EFp3_clear(&Q_neg);
	Fp18_clear(&l_sum);
	EFp3_clear(&T);
	EFp3_clear(&P_map);
	EFp3_clear(&Q_map);
	EFp3_clear(&EFp3_tmp);
	Fp3_clear(&L);
	Fp3_clear(&xy);
	Fp3_clear(&xy_2);
	Fp3_clear(&y_inv);
	Fp3_clear(&tmp);
	Fp3_clear(&y_tmp);
	mpz_clear(p3);
	Fp18_clear(&ltt);
	Fp18_clear(&ltp);
}
void Pseudo_type2_ADD_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *P,struct EFp3 *Q,struct Fp3 *L){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltp);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);


	Fp3_sub(&A,&P->x,&T->x);//xt-xp
	Fp3_sub(&B,&P->y,&T->y);//yt-yp
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&P->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_mul_xi_inv(&tmp1,&tmp1);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp_set_ui(&l_tmp.x0.x0.x0,1);

	Fp3_mul(&l_tmp.x0.x1,&E,L);

	Fp3_mul_xi_inv(&tmp4,&C);
	Fp3_neg(&l_tmp.x2.x1,&tmp4);

	Fp18_set(l_ANS,&l_tmp);
	EFp3_set(T_ANS,&x3_tmp);

	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltp);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
}
void Pseudo_type2_DBL_LINE(struct Fp18 *l_ANS,struct EFp3 *T_ANS,struct EFp3 *T,struct EFp3 *Q,struct Fp3 *L){
	struct Fp3 tmp1,tmp2,tmp3,tmp4,lambda,ltt;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&lambda);
	Fp3_init(&ltt);

	struct Fp18 l_tmp;
	Fp18_init(&l_tmp);

	struct Fp3 x,y,tmp;
	Fp3_init(&x);
	Fp3_init(&y);
	Fp3_init(&tmp);

	struct EFp3 x3_tmp;
	EFp3_init(&x3_tmp);
	struct Fp3 A,B,C,D,E;
	Fp3_init(&A);
	Fp3_init(&B);
	Fp3_init(&C);
	Fp3_init(&D);
	Fp3_init(&E);

	Fp3_add(&A,&T->y,&T->y);//xt-xp
	Fp3_mul(&B,&T->x,&T->x);
	Fp3_mul_ui(&B,&B,3);
	Fp3_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)

	Fp3_add(&D,&T->x,&T->x);
	Fp3_mul(&tmp1,&C,&C);
	Fp3_mul_xi(&tmp1,&tmp1);
	Fp3_sub(&x3_tmp.x,&tmp1,&D);

	Fp3_mul(&tmp2,&C,&T->x);
	Fp3_mul_xi(&tmp2,&tmp2);
	Fp3_sub(&E,&tmp2,&T->y);

	Fp3_mul(&tmp3,&C,&x3_tmp.x);
	Fp3_mul_xi(&tmp3,&tmp3);
	Fp3_sub(&x3_tmp.y,&E,&tmp3);

	Fp_set_ui(&l_tmp.x0.x0.x0,1);

	Fp3_mul(&l_tmp.x0.x1,&E,L);

	Fp3_neg(&l_tmp.x2.x1,&C);
	Fp18_set(l_ANS,&l_tmp);

	EFp3_set(T_ANS,&x3_tmp);

	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&lambda);
	Fp3_clear(&ltt);
	Fp18_clear(&l_tmp);
	Fp3_clear(&x);
	Fp3_clear(&y);
	Fp3_clear(&tmp);
	EFp3_clear(&x3_tmp);
	Fp3_clear(&A);
	Fp3_clear(&B);
	Fp3_clear(&C);
	Fp3_clear(&D);
	Fp3_clear(&E);
}
void Pseudo_type2_mul(struct Fp18 *ANS,struct Fp18 *A,struct Fp18 *B){
	struct Fp6 tmp02,tmp13,tmp45;
	Fp6_init(&tmp02);
	Fp6_init(&tmp13);
	Fp6_init(&tmp45);
	Fp3_set(&tmp02.x0,&A->x0.x0);
	Fp3_set(&tmp02.x1,&A->x2.x0);
	Fp3_set(&tmp13.x0,&A->x1.x0);
	Fp3_set(&tmp13.x1,&A->x0.x1);
	Fp3_set(&tmp45.x0,&A->x1.x1);
	Fp3_set(&tmp45.x1,&A->x2.x1);

	struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	Fp3_init(&tmp1);
	Fp3_init(&tmp2);
	Fp3_init(&tmp3);
	Fp3_init(&tmp4);
	Fp3_init(&tmp5);
	Fp3_init(&tmp6);

	struct Fp3 theta0,theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8,theta9,theta10;
	Fp3_init(&theta0);
	Fp3_init(&theta1);
	Fp3_init(&theta2);
	Fp3_init(&theta3);
	Fp3_init(&theta4);
	Fp3_init(&theta5);
	Fp3_init(&theta6);
	Fp3_init(&theta7);
	Fp3_init(&theta8);
	Fp3_init(&theta9);
	Fp3_init(&theta10);

	struct Fp18 t_ans;
	Fp18_init(&t_ans);

	struct Fp6 B_tmp;
	Fp6_init(&B_tmp);
	Fp3_set(&B_tmp.x0,&B->x0.x1);
	Fp3_set(&B_tmp.x1,&B->x2.x1);

	Fp3_mul(&theta3,&tmp02.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta7,&tmp02.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp02.x0,&tmp02.x1);//a+b
	Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta3);//(a+b)(c+d)-ac-bd
	Fp3_sub(&theta5,&tmp3,&theta7);

	Fp3_mul(&theta4,&tmp13.x0,&B_tmp.x0);//a*c
	Fp3_mul(&theta8,&tmp13.x1,&B_tmp.x1);//b*d
	Fp3_add(&tmp1,&tmp13.x0,&tmp13.x1);//a+b
	// Fp3_add(&tmp2,&B_tmp.x0,&B_tmp.x1);//c+d
	Fp3_mul(&tmp6,&tmp1,&tmp2);//(a+b)(c+d)
	Fp3_sub(&tmp3,&tmp6,&theta4);
	Fp3_sub(&theta6,&tmp3,&theta8);

	Fp3_mul(&tmp1,&tmp45.x0,&B_tmp.x0);
	Fp3_add(&theta7,&tmp1,&theta7);
	Fp3_mul(&tmp1,&tmp45.x1,&B_tmp.x0);
	Fp3_add(&theta8,&tmp1,&theta8);
	Fp3_mul(&theta9,&tmp45.x0,&B_tmp.x1);
	Fp3_mul(&theta10,&tmp45.x1,&B_tmp.x1);

	Fp3_mul_xi(&theta0,&theta6);
	Fp3_mul_xi(&theta1,&theta7);
	Fp3_mul_xi(&theta2,&theta8);
	Fp3_mul_xi(&tmp4,&theta9);
	Fp3_add(&theta3,&theta3,&tmp4);
	Fp3_mul_xi(&tmp5,&theta10);
	Fp3_add(&theta4,&theta4,&tmp5);

	Fp3_set(&t_ans.x0.x0,&theta0);
	Fp3_set(&t_ans.x1.x0,&theta1);
	Fp3_set(&t_ans.x2.x0,&theta2);
	Fp3_set(&t_ans.x0.x1,&theta3);
	Fp3_set(&t_ans.x1.x1,&theta4);
	Fp3_set(&t_ans.x2.x1,&theta5);

	Fp18_add(ANS,&t_ans,A);

	Fp6_clear(&tmp02);
	Fp6_clear(&tmp13);
	Fp6_clear(&tmp45);
	Fp3_clear(&tmp1);
	Fp3_clear(&tmp2);
	Fp3_clear(&tmp3);
	Fp3_clear(&tmp4);
	Fp3_clear(&tmp5);
	Fp3_clear(&tmp6);
	Fp3_clear(&theta0);
	Fp3_clear(&theta1);
	Fp3_clear(&theta2);
	Fp3_clear(&theta3);
	Fp3_clear(&theta4);
	Fp3_clear(&theta5);
	Fp3_clear(&theta6);
	Fp3_clear(&theta7);
	Fp3_clear(&theta8);
	Fp3_clear(&theta9);
	Fp3_clear(&theta10);
	Fp18_clear(&t_ans);
	Fp6_clear(&B_tmp);
}

//-----------------------------------------------------------------------------------------
void Final_Exp(struct Fp18 *ANS,struct Fp18 *A){
	int pattern;
	int bit_prime;
	bit_prime=mpz_sizeinbase(prime,2);
	int pow_bit_separate[6][bit_prime];
	int i=0,j=0;

	mpz_t t_mpz;
	mpz_t pow_q,pow_r;
	struct Fp18 A_inv,t_ans,tmp;
	struct Fp18 t[64];

	mpz_init(t_mpz);
	mpz_init(pow_q);
	mpz_init(pow_r);
	Fp18_init(&A_inv);
	Fp18_init(&t_ans);
	Fp18_init(&tmp);
	for(i=0;i<64;i++){
		Fp18_init(&t[i]);
	}
	//A^(p^9-1)^p^3+1
	Fp18_invert(&A_inv,A);

	Fp18_frobenius_map(&t_ans,A,9);
	Fp18_mul(&A_inv,&t_ans,&A_inv);

	Fp18_frobenius_map(&t_ans,&A_inv,3);
	Fp18_mul(&t_ans,&t_ans,&A_inv);

	Fp18_set(&tmp,&t_ans);

	//(p^6-p^3+1)/r
	mpz_pow_ui(t_mpz,prime,3);
	mpz_mul(pow_q,t_mpz,t_mpz);
	mpz_sub(pow_q,pow_q,t_mpz);
	mpz_add_ui(pow_q,pow_q,1);
	mpz_tdiv_q(pow_q,pow_q,order);

	//transrate prime adic representation
	i=0;
	while(mpz_cmp_ui(pow_q,0)!=0){
		mpz_tdiv_qr(pow_q,pow_r,pow_q,prime);
		for(j=0;j<bit_prime;j++){
			pow_bit_separate[i][j]=mpz_tstbit(pow_r,j);
		}
		i++;
	}

	//precomputing
	Fp18_set(&t[1],&tmp);
	Fp18_frobenius_map(&t[2],&t_ans,1);
	Fp18_frobenius_map(&t[4],&t[2],1);
	Fp18_frobenius_map(&t[8],&t[4],1);
	Fp18_frobenius_map(&t[16],&t[8],1);
	Fp18_frobenius_map(&t[32],&t[16],1);

	for(i=1;i<64;i=i*2){
		for(j=1;j<i;j++){
			Fp18_mul(&t[i+j],&t[i],&t[j]);
		}
	}

	//main computation
	Fp18_set_ui(&t_ans,0);
	Fp_set_ui(&t_ans.x0.x0.x0,1);
	for(i=bit_prime-1;i>=0;i--){
		pattern=pow_bit_separate[0][i];
		pattern+=pow_bit_separate[1][i]*2;
		pattern+=pow_bit_separate[2][i]*4;
		pattern+=pow_bit_separate[3][i]*8;
		pattern+=pow_bit_separate[4][i]*16;
		pattern+=pow_bit_separate[5][i]*32;
		Fp18_mul(&t_ans,&t_ans,&t_ans);
		if(pattern!=0){
			Fp18_mul(&t_ans,&t_ans,&t[pattern]);
		}
	}

	Fp18_set(ANS,&t_ans);

	mpz_clear(t_mpz);
	mpz_clear(pow_q);
	mpz_clear(pow_r);
	Fp18_clear(&A_inv);
	Fp18_clear(&t_ans);
	Fp18_clear(&tmp);
	for(i=0;i<64;i++){
		Fp18_clear(&t[i]);
	}
}

void check_Pairing(void){
	struct EFp P,R;
	EFp_init(&P);
	EFp_init(&R);

	struct EFp18 Q,S;
	EFp18_init(&Q);
	EFp18_init(&S);

	struct Fp18 ans,tmp1,tmp2,tmp3;
	Fp18_init(&ans);
	Fp18_init(&tmp1);
	Fp18_init(&tmp2);
	Fp18_init(&tmp3);

	mpz_t a,b,ab;
	mpz_init(a);
	mpz_init(b);
	mpz_init(ab);

	mpz_set_ui(a,31);
	mpz_set_ui(b,11);
	mpz_mul(ab,a,b);

	EFp_random_set(&P);
	EFp18_random_set_G2(&Q);

	// mpz_set_str(tmp.x.x0,"1646176220195929079712050957626783481753869858",10);
	// mpz_set_str(tmp.y.x0,"428008951585620392582406991174654065346629094",10);
	//
	// mpz_set_str(Q.x.x2.x0.x0.x0,"1670519066210195597576851690522033200402307978",10);
	// mpz_set_str(Q.x.x2.x0.x1.x0,"1740081911946076277121569226204477504834988801",10);
	// mpz_set_str(Q.x.x2.x0.x2.x0,"130641026103545360394319616885244749166339863",10);
	//
	// mpz_set_str(Q.y.x0.x1.x0.x0,"1256093386145159259760780362140192333294638484",10);
	// mpz_set_str(Q.y.x0.x1.x1.x0,"888307084556604336179130159561062687061426574",10);
	// mpz_set_str(Q.y.x0.x1.x2.x0,"1479453669659142908289832956903177318716071880",10);


	//---------------------------------------------------
	// printf("Tate Pairing\n");
	// // EFp_random_set(&tmp);
	// EFp18_set_EFp(&P,&tmp);
	//
	// EFp18_random_set(&S);
	//
	// Tate_Pairing(&tmp1,&P,&S);
	// printf("G1=");
	// EFp18_printf(&P);
	// printf("G2=");
	// EFp18_printf(&S);
	//
	// Fp18_pow(&tmp1,&tmp1,ab);
	// printf("f^ab=");
	// Fp18_printf(&tmp1);
	//
	// EFp18_SCM(&R,&P,a);
	// EFp18_SCM(&S,&S,b);
	//
	// Tate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);
	//----------------------------------------------------
	// printf("Ate Pairing\n");
	// printf("G1=");
	// EFp_printf(&P);
	// printf("G2=");
	// EFp18_printf(&Q);
	// EFp18_set_EFp(&S,&P);
	//
	// Ate_Pairing(&tmp1,&S,&Q);
	//
	// Fp18_pow(&tmp1,&tmp1,ab);
	// printf("f^ab=");
	// Fp18_printf(&tmp1);
	//
	// EFp18_SCM(&S,&S,a);
	// EFp18_SCM(&Q,&Q,b);
	//
	// Ate_Pairing(&tmp2,&S,&Q);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);
	//----------------------------------------------------
	printf("Optimal Ate Pairing\n");
	// printf("G1=");
	// EFp_printf(&P);
	// printf("G2=");
	// EFp18_printf(&Q);

	Optimal_Ate_Pairing(&tmp1,&P,&Q);

	// Fp18_pow(&tmp1,&tmp1,ab);
	printf("f^ab=");
	Fp18_printf(&tmp1);

	// EFp_SCM(&R,&P,a);
	// EFp18_SCM(&S,&Q,b);
	//
	// Optimal_Ate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);
	//----------------------------------------------------
	// printf("Ate Pairing Sparse\n");
	//
	// EFp3_set_EFp(&P_EFp3,&tmp);
	// EFp3_set_EFp18_Sparse(&Q_EFp3,&Q);
	// EFp18_printf(&Q);
	//
	// // printf("G1=");
	// // EFp3_printf(&P_EFp3);
	// // printf("G2=");
	// // EFp3_printf(&Q_EFp3);
	//
	// Sparse_Ate_Pairing(&tmp1,&P_EFp3,&Q_EFp3);
	//
	// Fp18_pow(&tmp1,&tmp1,ab);
	// printf("\nf^ab=");
	// Fp18_printf(&tmp1);
	//
	// EFp3_SCM(&R_EFp3,&P_EFp3,a);
	// EFp3_SCM(&S_EFp3,&Q_EFp3,b);
	//
	// Sparse_Ate_Pairing(&tmp2,&R_EFp3,&S_EFp3);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);
	//----------------------------------------------------
	printf("Optimal Ate Pairing Sparse\n");

	// EFp3_set_EFp(&P_EFp3,&tmp);
	// EFp3_set_EFp18_Sparse(&Q_EFp3,&Q);
	// EFp18_printf(&Q);
	//
	// printf("G1=");
	// EFp3_printf(&P_EFp3);
	// printf("G2=");
	// EFp3_printf(&Q_EFp3);
	//
	Sparse_Optimal_Ate_Pairing(&tmp1,&P,&Q);

	// Fp18_pow(&tmp1,&tmp1,ab);
	printf("\nf^ab=");
	Fp18_printf(&tmp1);

	// EFp_SCM(&R,&P,a);
	// EFp18_SCM(&S,&Q,b);
	//
	// Sparse_Optimal_Ate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);
	//----------------------------------------------------
	// printf("Ate Pairing Pseudo Sparse\n");
	//
	// printf("G1=");
	// EFp_printf(&P);
	// printf("G2=");
	// EFp18_printf(&Q);
	// //
	// Sparse_Optimal_Ate_Pairing(&tmp1,&P,&Q);
	//
	// Fp18_pow(&tmp1,&tmp1,ab);
	// printf("\nf^ab=");
	// Fp18_printf(&tmp1);
	//
	// EFp_SCM(&R,&P,a);
	// EFp18_SCM(&S,&Q,b);
	//
	// Sparse_Optimal_Ate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);

	//----------------------------------------------------
	printf("Optimal Ate Pairing Pseudo Sparse\n");
	//
	// printf("G1=");
	// EFp_printf(&P);
	// printf("G2=");
	// EFp18_printf(&Q);
	// //
	Pseudo_Sparse_Optimal_Ate_Pairing(&tmp1,&P,&Q);
	//
	// Fp18_pow(&tmp1,&tmp1,ab);
	printf("\nf^ab=");
	Fp18_printf(&tmp1);

	// EFp_SCM(&R,&P,a);
	// EFp18_SCM(&S,&Q,b);
	//
	// Pseudo_Sparse_Optimal_Ate_Pairing(&tmp2,&R,&S);
	//
	// printf("f'  =");
	// Fp18_printf(&tmp2);

	// ----------------------------------------------------

	EFp18_clear(&Q);
	Fp18_clear(&ans);
	Fp18_clear(&tmp1);
	Fp18_clear(&tmp2);
	Fp18_clear(&tmp3);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(ab);
}
void Masure_pairing_time(void){
	struct EFp P;
	EFp_init(&P);

	struct EFp18 Q;
	EFp18_init(&Q);

	struct Fp18 ans,tmp1;
	Fp18_init(&ans);
	Fp18_init(&tmp1);


	int i;
	int loop=100;
	// struct timeval opt_1,opt_2;
	// struct timeval sparse_1,sparse_2;
	// struct timeval pseudo_1,pseudo_2;
	// double opt_sum=0,sparse_sum=0,pseudo_sum=0;

	mpz_t set_2;
	mpz_init(set_2);
	mpz_set_ui(set_2,2);

	EFp_random_set(&P);
	EFp18_random_set_G2(&Q);

	for(i=1;i<=loop;i++){
		// Fp18_random(&tmp1);
		// gettimeofday(&opt_1, NULL);
		// Fp18_frobenius_map(&ans,&tmp1,1);
		// gettimeofday(&opt_2, NULL);
		// opt_sum+=((double)(opt_2.tv_sec - opt_1.tv_sec)+ (double)(opt_2.tv_usec-opt_1.tv_usec)*1.0E-6);
		// gettimeofday(&pseudo_1, NULL);
		// Fp18_frobenius_map_old(&ans,&tmp1);
		// gettimeofday(&pseudo_2, NULL);
		// pseudo_sum+=((double)(pseudo_2.tv_sec - pseudo_1.tv_sec)+ (double)(pseudo_2.tv_usec-pseudo_1.tv_usec)*1.0E-6);


		EFp_SCM(&P,&P,set_2);
		EFp18_SCM(&Q,&Q,set_2);
		//----------------------------------------------------
		// // printf("Optimal Ate Pairing\n");
		//
		// EFp18_set_EFp(&P,&tmp);
		//
		// gettimeofday(&opt_1, NULL);
		// Optimal_Ate_Pairing(&tmp1,&P,&Q);
		// gettimeofday(&opt_2, NULL);
		// opt_sum+=((double)(opt_2.tv_sec - opt_1.tv_sec)+ (double)(opt_2.tv_usec-opt_1.tv_usec)*1.0E-6);
		// // //----------------------------------------------------
		// // printf("Optimal Ate Pairing Pseudo Sparse\n");
		//
		// EFp3_set_EFp(&P_EFp3,&tmp);
		// EFp3_set_EFp18_Sparse(&Q_EFp3,&Q);
		//
		// gettimeofday(&pseudo_1, NULL);
		Pseudo_Sparse_Optimal_Ate_Pairing(&tmp1,&P,&Q);
		// gettimeofday(&pseudo_2, NULL);
		// pseudo_sum+=((double)(pseudo_2.tv_sec - pseudo_1.tv_sec)+ (double)(pseudo_2.tv_usec-pseudo_1.tv_usec)*1.0E-6);
		//----------------------------------------------------
		// printf("Optimal Ate Pairing Sparse\n");

		// gettimeofday(&sparse_1, NULL);
		Sparse_Optimal_Ate_Pairing(&tmp1,&P,&Q);
		// gettimeofday(&sparse_2, NULL);
		// sparse_sum+=((double)(sparse_2.tv_sec - sparse_1.tv_sec)+ (double)(sparse_2.tv_usec-sparse_1.tv_usec)*1.0E-6);
		//----------------------------------------------------
		printf("%d\n",i);

	}

	// for(i=0;i<loop;i++){
	// 	Fp18_random(&tmp1);
	// 	gettimeofday(&sparse_1, NULL);
	// 	Final_Exp(&ans,&tmp1);
	// 	gettimeofday(&sparse_2, NULL);
	// 	opt_sum+=((double)(sparse_2.tv_sec - sparse_1.tv_sec)+ (double)(sparse_2.tv_usec-sparse_1.tv_usec)*1.0E-6);
	// 	gettimeofday(&sparse_1, NULL);
	// 	Final_Exp_old1(&ans,&tmp1);
	// 	gettimeofday(&sparse_2, NULL);
	// 	sparse_sum+=((double)(sparse_2.tv_sec - sparse_1.tv_sec)+ (double)(sparse_2.tv_usec-sparse_1.tv_usec)*1.0E-6);
	// 	// gettimeofday(&sparse_1, NULL);
	// 	// Final_Exp_old2(&ans,&tmp1);
	// 	// gettimeofday(&sparse_2, NULL);
	// 	// pseudo_sum+=((double)(sparse_2.tv_sec - sparse_1.tv_sec)+ (double)(sparse_2.tv_usec-sparse_1.tv_usec)*1.0E-6);
	// 	printf("%d\n",i);
	// }



	// printf("opt_miller = %lf\n",(opt_sum/loop));
	// printf("sps_miller = %lf\n",(sparse_sum/loop));
	// printf("pse_miller = %lf\n",(pseudo_sum/loop));



	printf("opt_add = %lf\n",(opt_add/loop));
	printf("opt_dbl = %lf\n",(opt_dbl/loop));
	printf("opt_mul = %lf\n",(opt_mul/loop));
	printf("sps_add = %lf\n",(sps_add/loop));
	printf("sps_dbl = %lf\n",(sps_dbl/loop));
	printf("sps_mul = %lf\n",(sps_mul/loop));
	printf("pse_add = %lf\n",(pse_add/loop));
	printf("pse_dbl = %lf\n",(pse_dbl/loop));
	printf("pse_mul = %lf\n",(pse_mul/loop));

	// printf("time = %lf\n",(double)(t2-t1));

	EFp_clear(&P);
	EFp18_clear(&Q);
	Fp18_clear(&ans);
	Fp18_clear(&tmp1);
}
