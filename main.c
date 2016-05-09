//
//  main.c
//  Fp18_Arith
//
//  Created by Khandaker Md. Al-Amin on 4/15/16.
//  Copyright © 2016 Khandaker Md. Al-Amin. All rights reserved.
//

#include"f18.h"
#include "parameters.h"
#include <assert.h>

#define C1 2
#define C2 -1
#define C1_sq 4
int fp_mul,fp_add;
int main(void)
{
    mpz_init(X);
    mpz_set_str(X,"18446893747415302274",10);//64bit22818362
//    mpz_set_str(X,"22818362",10);
    
    mpz_init(prime);
//    mpz_set_str(prime,"7",10);
mpz_set_str(prime,"638508488753389092195439720208817843316062995515124574544010783086949264599518591426924070465759719017041310253064889352571109595156203400346010025740431",10);
    mpz_init(r_order);
    mpz_init(r_order_EFp);
    mpz_init(t_trace);
    mpz_init(b);
    mpz_init(c1_leg);
    mpz_init(c1_leg_bar);
    mpz_init(c1_omega);
    mpz_init(c1_omega_bar);
    
    generate_parameters();
    gmp_printf("p=%Zd\n",prime);
    gmp_printf("r=%Zd\n",r_order);
    gmp_printf("t=%Zd\n",t_trace);
    gmp_printf("#E(Fp)=%Zd\n",r_order_EFp);
    printf("p = %d bit\n",(int)mpz_sizeinbase(prime,2));
    printf("r = %d bit\n",(int)mpz_sizeinbase(r_order,2));
    printf("t = %d bit\n",(int)mpz_sizeinbase(t_trace,2));
    
    fp_mul = 0;
    fp_add = 0;
    
        struct Fp3 P,RES,RES1,Q;
        Fp3_init(&P);
        Fp3_init(&Q);
        Fp3_init(&RES);
        Fp3_init(&RES1);
        Fp3_take_input(&P);
        Fp3_invert(&RES, &P);
        printf("new\n");
        Fp3_printf(&RES);
        printf("FP Mul =%d , FP ADD =%d\n",fp_mul,fp_add);
        Fp3_invert_prev(&RES1, &P);
        printf("prev\n");
        Fp3_printf(&RES1);
    
    //        struct Fp6 P,RES,RES1,Q;
    //        Fp6_init(&P);
    //        Fp6_init(&Q);
    //        Fp6_init(&RES);
    //        Fp6_init(&RES1);
    //        Fp6_take_input(&P);
    //        Fp6_invert(&RES, &P);
    //        printf("new\n");
    //        Fp6_printf(&RES);
    //    printf("FP Mul =%d , FP ADD =%d\n",fp_mul,fp_add);
    //        Fp6_invert_prev(&RES1, &P);
    //        printf("prev\n");
    //        Fp6_printf(&RES1);
    
//    struct Fp18 P,RES,RES1,Q;
//    Fp18_init(&P);
//    Fp18_init(&Q);
//    Fp18_init(&RES);
//    Fp18_init(&RES1);
//    Fp18_take_input(&P);
//    //        Fp18_take_input(&Q);
//    //        Fp18_mul(&RES, &P, &Q);
//    Fp18_invert(&RES, &P);
//    printf("new\n");
//    Fp18_printf(&RES);
//    printf("FP Mul =%d , FP ADD =%d\n",fp_mul,fp_add);
    //        Fp18_invert_prev(&RES1, &P);
    //        printf("prev\n");
    //        Fp18_printf(&RES1);
    
    
    
    
//    struct EFp18 P,Q;
//    EFp18_init(&P);
//    EFp18_init(&Q);
//    
//    EFp18_printf(&Q);
//    EFp18_SCM(&P,&P,prime);
//    EFp18_printf(&P);
//    printf("%d\n",EFp18_cmp(&P,&Q));
    
    mpz_clear(prime);
    mpz_clear(r_order);
    mpz_clear(b);
    mpz_clear(c1_leg);
    mpz_clear(c1_leg_bar);
    return 0;
}

#pragma mark Input from Terminals methods
void Fp_take_input(struct Fp *RES)
{
    struct Fp tmp;
    Fp_init(&tmp);
    
    char inputStr[1025];
    int flag;
    printf ("Enter Fp element:\n");
    scanf("%1024s",inputStr);
    mpz_set_ui(tmp.x_0,0);
    
    flag = mpz_set_str(tmp.x_0,inputStr,10);
    assert(flag == 0);
    
    Fp_set(RES,&tmp);
    //    Fp_printf(RES);
    Fp_clear(&tmp);
}

void Fp3_take_input(struct Fp3 *RES)
{
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_take_input(&tmp.a0);
    Fp_take_input(&tmp.a1);
    Fp_take_input(&tmp.a2);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}

void Fp6_take_input(struct Fp6 *RES)
{
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_take_input(&tmp.a0);
    Fp3_take_input(&tmp.a1);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}

void Fp18_take_input(struct Fp18 *RES)
{
    struct Fp18 tmp;
    Fp18_init(&tmp);
    
    Fp6_take_input(&tmp.m0);
    Fp6_take_input(&tmp.m1);
    Fp6_take_input(&tmp.m2);
    
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}

#pragma mark Fp method implementations
void Fp_init(struct Fp *A){
    mpz_init(A->x_0);
}
void Fp_set(struct Fp *RES,struct Fp *A){
    mpz_set(RES->x_0,A->x_0);
}
void Fp_set_mpz(struct Fp *RES,mpz_t a){
    mpz_set(RES->x_0,a);
}
void Fp_set_ui(struct Fp *A,signed long int B){
    mpz_set_ui(A->x_0,B);
}
void Fp_random(struct Fp *A){
    mpz_random(A->x_0,10);
    mpz_mod(A->x_0,A->x_0,prime);
}
void Fp_clear(struct Fp *A){
    mpz_clear(A->x_0);
}
void Fp_printf(struct Fp *A){
    gmp_printf("%Zd\n",A->x_0);
}
void Fp_add(struct Fp *RES,struct Fp *a,struct Fp *b){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_add(tmp.x_0,a->x_0,b->x_0);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    fp_add++;
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_add_ui(struct Fp *RES,struct Fp *a,unsigned long int b){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_add_ui(tmp.x_0,a->x_0,b);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_sub(struct Fp *RES,struct Fp *a,struct Fp *b){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_sub(tmp.x_0,a->x_0,b->x_0);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_sub_ui(struct Fp *RES,struct Fp *a,unsigned long int b){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_sub_ui(tmp.x_0,a->x_0,b);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    
    Fp_set(RES,&tmp);
    
    Fp_clear(&tmp);
}
void Fp_mul(struct Fp *RES,struct Fp *a,struct Fp *b){
    struct Fp tmp;
    Fp_init(&tmp);
    mpz_mod(b->x_0,b->x_0,prime);
    mpz_mul(tmp.x_0,a->x_0,b->x_0);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    fp_mul++;
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_mul_ui(struct Fp *RES,struct Fp *a,unsigned long int b){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_mul_ui(tmp.x_0,a->x_0,b);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_mul_c1(struct Fp *RES,struct Fp *a){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_mul_2exp(tmp.x_0,a->x_0,1);//C1=2, so 1 bit shift
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_inv(struct Fp *RES,struct Fp *a){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_invert(tmp.x_0,a->x_0,prime);
    
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_div(struct Fp *RES,struct Fp *a,struct Fp *b){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_invert(tmp.x_0,b->x_0,prime);
    mpz_mul(tmp.x_0,a->x_0,tmp.x_0);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
void Fp_pow(struct Fp *RES,struct Fp *A,mpz_t B){
    unsigned long long i,length;
    length= (unsigned long long)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    //    printf("Length: %lld\n",length);
    //    printf("p: %s\n",&B_binary[0]);
    
    struct Fp tmp,dumy;
    Fp_init(&tmp);
    Fp_init(&dumy);
    Fp_set(&tmp,A);
    Fp_set(&dumy,A);
    
    for(i=1;B_binary[i]!='\0';i++){
        Fp_mul(&tmp,&tmp,&tmp);
        Fp_mul(&dumy,&tmp,A);
        if(B_binary[i]=='1'){
            Fp_mul(&tmp,&tmp,A);
            Fp_mul(&dumy,&dumy,&dumy);
        }
    }
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
    Fp_clear(&dumy);
}

void Fp_sqrt(struct Fp *RES,struct Fp *A){
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
    
    while(mpz_legendre(n_tmp.x_0,prime)!=-1){
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
        // gmp_printf("%Zd,%Zd\n",x.a0.a0,x.a1.a0);
        mpz_sub_ui(tmp_mpz,r_tmp,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,prime);
        Fp_pow(&t_tmp,&y_tmp,tmp_mpz);
        Fp_pow(&y_tmp,&t_tmp,set_2);
        mpz_set_ui(r_tmp,m);
        Fp_mul(&x_tmp,&x_tmp,&t_tmp);
        Fp_mul(&b_tmp,&b_tmp,&y_tmp);
    }
    
    Fp_set(RES,&x_tmp);
    
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
}
void Fp_neg(struct Fp *RES,struct Fp *A){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_sub(tmp.x_0,prime,A->x_0);
    mpz_mod(tmp.x_0,tmp.x_0,prime);
    
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}
int Fp_cmp(struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x_0,B->x_0)==0){
        return 0;
    }
    return 1;
}
int Fp_cmp_mpz(struct Fp *A,mpz_t B){
    if(mpz_cmp(A->x_0,B)==0){
        return 0;
    }
    return 1;
}


#pragma mark ==Fp3== method implementations
void Fp3_init(struct Fp3 *A){
    Fp_init(&A->a0);
    Fp_init(&A->a1);
    Fp_init(&A->a2);
}
void Fp3_set(struct Fp3 *RES,struct Fp3 *A){
    Fp_set(&RES->a0,&A->a0);
    Fp_set(&RES->a1,&A->a1);
    Fp_set(&RES->a2,&A->a2);
}
void Fp3_set_ui(struct Fp3 *A,signed long int B){
    Fp_set_ui(&A->a0,B);
    Fp_set_ui(&A->a1,B);
    Fp_set_ui(&A->a2,B);
}
void Fp3_random(struct Fp3 *A){
    Fp_random(&A->a0);
    Fp_random(&A->a1);
    Fp_random(&A->a2);
}
void Fp3_clear(struct Fp3 *A){
    Fp_clear(&A->a0);
    Fp_clear(&A->a1);
    Fp_clear(&A->a2);
}
void Fp3_printf(struct Fp3 *A){
    gmp_printf("%Zd,%Zd,%Zd\n",A->a0.x_0,A->a1.x_0,A->a2.x_0);
}
void Fp3_add(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_add(&tmp.a0,&A->a0,&B->a0);
    Fp_add(&tmp.a1,&A->a1,&B->a1);
    Fp_add(&tmp.a2,&A->a2,&B->a2);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}
void Fp3_add_ui(struct Fp3 *RES,struct Fp3 *A,unsigned long int B){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_add_ui(&tmp.a0,&A->a0,B);
    Fp_add_ui(&tmp.a1,&A->a1,B);
    Fp_add_ui(&tmp.a2,&A->a2,B);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}
void Fp3_sub(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_sub(&tmp.a0,&A->a0,&B->a0);
    Fp_sub(&tmp.a1,&A->a1,&B->a1);
    Fp_sub(&tmp.a2,&A->a2,&B->a2);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}
void Fp3_mul(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B){
    //(a0,a1,a2)*(b0,b1,b2)=(x0y0+xi((a1+a2)(b1+b2)-a1b1-x2y2),xix2y2+(a0+a1)(b0+b1)-x0y0-a1b1,a1b1+(a0+a2)(b0+b2)-x0y0-x2y2)
    struct Fp A0,A1,A2,tmp_a01,tmp_a12,tmp_a20,tmp_b01,tmp_b12,tmp_b20,A3,A4,A5,t0,t1,t2,t3,t4,tmp;
    
    struct Fp3 t_RES;
    Fp_init(&A0);
    Fp_init(&A1);
    Fp_init(&A2);
    Fp_init(&tmp_a01);
    Fp_init(&tmp_a12);
    Fp_init(&tmp_a20);
    Fp_init(&tmp_b01);
    Fp_init(&tmp_b12);
    Fp_init(&tmp_b20);
    Fp_init(&A3);
    Fp_init(&A4);
    Fp_init(&A5);
    Fp_init(&t0);
    Fp_init(&t1);
    Fp_init(&t2);
    Fp_init(&t3);
    Fp_init(&t4);
    Fp_init(&tmp);
    Fp3_init(&t_RES);
    
    Fp_mul(&A0,&A->a0,&B->a0);//A0=a0*b0
    Fp_mul(&A1,&A->a1,&B->a1);//A1=a1*b1
    Fp_mul(&A2,&A->a2,&B->a2);//A2=a2*b2
    
    Fp_add(&tmp_a01,&A->a0,&A->a1);//a0+a1
    Fp_add(&tmp_a12,&A->a1,&A->a2);//a1+a2
    Fp_add(&tmp_a20,&A->a0,&A->a2);//a2+a0
    Fp_add(&tmp_b01,&B->a0,&B->a1);//b0+b1
    Fp_add(&tmp_b12,&B->a1,&B->a2);//b1+b2
    Fp_add(&tmp_b20,&B->a0,&B->a2);//b2+b0
    
    Fp_mul(&A3,&tmp_a01,&tmp_b01);//A3=(a0+a1)(b0+b1)
    Fp_mul(&A4,&tmp_a20,&tmp_b20);//A4=(a0+a2)(b0+b2)
    Fp_mul(&A5,&tmp_a12,&tmp_b12);//A5=(a1+a2)(b1+b2)
    
    Fp_set(&t0,&A0);//t0=A0
    Fp_set(&t4, &A2);//t4= A2
    Fp_sub(&A3,&A3,&A0);//A3-A0
    Fp_sub(&t1,&A3,&A1);//t1 = A3-A1-A0
    Fp_sub(&A4,&A4,&A0);//A4-A0
    Fp_sub(&A4,&A4,&A2);//A4-A2-A0
    Fp_add(&t2,&A1,&A4);//t2 = A4−A2−A0+A1
    Fp_sub(&A5,&A5,&A1);//A5-A1
    Fp_sub(&t3,&A5,&A2);//t3 = A5 −A1 −A2
    
    //ab = (t0 + c1t3) + (t1 + c1t4)ω + t2ω2.
    Fp_mul_c1(&tmp, &t3);
    //    Fp_mul_ui(&tmp,&t3,C1);
    Fp_add(&t_RES.a0, &tmp, &t0);
    Fp_mul_c1(&tmp, &t4);
    //    Fp_mul_ui(&tmp,&t4,C1);
    Fp_add(&t_RES.a1, &tmp, &t1);
    Fp_set(&t_RES.a2, &t2);
    
    Fp3_set(RES,&t_RES);
    
    Fp_clear(&A0);
    Fp_clear(&A1);
    Fp_clear(&A2);
    Fp_clear(&tmp_a01);
    Fp_clear(&tmp_a12);
    Fp_clear(&tmp_a20);
    Fp_clear(&tmp_b01);
    Fp_clear(&tmp_b12);
    Fp_clear(&tmp_b20);
    Fp_clear(&A3);
    Fp_clear(&A4);
    Fp_clear(&A5);
    Fp_clear(&t0);
    Fp_clear(&t1);
    Fp_clear(&t2);
    Fp_clear(&t3);
    Fp_clear(&t4);
    Fp_clear(&tmp);
}

void Fp3_mul_omega(struct Fp3 *RES,struct Fp3 *A){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_mul_c1(&tmp.a0,&A->a2);//omega^3=C1
    Fp_set(&tmp.a1,&A->a0);
    Fp_set(&tmp.a2,&A->a1);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}

void Fp3_mul_ui(struct Fp3 *RES,struct Fp3 *A,unsigned long int B){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_mul_ui(&tmp.a0,&A->a0,B);
    Fp_mul_ui(&tmp.a1,&A->a1,B);
    Fp_mul_ui(&tmp.a2,&A->a2,B);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}
void Fp3_mul_Fp(struct Fp3 *RES,struct Fp3 *A,struct Fp *B){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    
    Fp_mul(&tmp.a0,&A->a0,B);
    Fp_mul(&tmp.a1,&A->a1,B);
    Fp_mul(&tmp.a2,&A->a2,B);
    
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}
void Fp3_frobenius_map(struct Fp3 *RES,struct Fp3 *A){
    //a^−1 = n(a)^−1(a^q*a^q^2 ),
    struct Fp3 t_RES;
    Fp3_init(&t_RES);
    
    struct Fp tmp_c1;
    Fp_init(&tmp_c1);
    
    mpz_set(tmp_c1.x_0,c1_leg);
    Fp_set(&t_RES.a0,&A->a0);
    Fp_mul(&t_RES.a1,&A->a1,&tmp_c1);
    Fp_set_mpz(&tmp_c1, c1_leg_bar);
    Fp_mul(&t_RES.a2,&A->a2,&tmp_c1);
    
    Fp3_set(RES,&t_RES);
    Fp3_clear(&t_RES);
    Fp_clear(&tmp_c1);
}


void Fp3_invert(struct Fp3 *RES, struct Fp3 *A){
    struct Fp3 t_RES,Aq,Aqsq,t_Aq,AqAqsq,t_norm;
    Fp3_init(&t_RES);
    Fp3_init(&Aq);
    Fp3_init(&Aqsq);
    Fp3_init(&t_Aq);
    Fp3_init(&AqAqsq);
    Fp3_init(&t_norm);
    
    Fp3_frobenius_map(&Aq, A);
    Fp3_set(&t_Aq, &Aq);
    Fp3_frobenius_map(&Aqsq,&t_Aq);
    Fp3_mul(&AqAqsq, &Aq, &Aqsq);
    Fp3_mul(&t_norm, &AqAqsq, A);
    
    struct Fp norm,T0;
    Fp_init(&norm);
    Fp_init(&T0);
    printf("Fp3 NORM=");
    Fp3_printf(&t_norm);
    if (mpz_cmp_ui(t_norm.a1.x_0,0) != 0 || mpz_cmp_ui(t_norm.a2.x_0,0) != 0)
    {
        printf("Fp3 inverse calculation goes wrong\n");
        exit(0);
    }
    
    mpz_set(T0.x_0,t_norm.a0.x_0);
    Fp_inv(&norm, &T0);
    Fp_mul(&t_RES.a0, &AqAqsq.a0, &norm);
    Fp_mul(&t_RES.a1, &AqAqsq.a1, &norm);
    Fp_mul(&t_RES.a2, &AqAqsq.a2, &norm);
    
    Fp3_set(RES,&t_RES);
    
    Fp3_clear(&t_RES);
    Fp3_clear(&Aq);
    Fp3_clear(&Aqsq);
    Fp3_clear(&AqAqsq);
    Fp3_clear(&t_norm);
    Fp_clear(&T0);
    Fp_clear(&norm);
}
void Fp3_div(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B){
    struct Fp3 tmp,t_RES;
    Fp3_init(&tmp);
    Fp3_init(&t_RES);
    
    Fp3_invert(&tmp,B);
    Fp3_mul(&t_RES,A,&tmp);
    
    Fp3_set(RES,&t_RES);
    
    Fp3_clear(&tmp);
    Fp3_clear(&t_RES);
}
void Fp3_pow(struct Fp3 *RES,struct Fp3 *A,mpz_t B){
    int i;
    int r;//bit数
    r= (int)mpz_sizeinbase(B,2);
    //    printf("r= %d\n",r);
    
    struct Fp3 res_tmp;
    Fp3_init(&res_tmp);
    Fp3_set(&res_tmp,A);
    
    struct Fp3 in_tmp;
    Fp3_init(&in_tmp);
    Fp3_set(&in_tmp,A);
    
    for(i=r-2;i>=0;i--){
        if(mpz_tstbit(B,i)==1){
            Fp3_mul(&res_tmp,&res_tmp,&res_tmp);//a*2
            Fp3_mul(&res_tmp,&res_tmp,&in_tmp);//*a
        }else{
            Fp3_mul(&res_tmp,&res_tmp,&res_tmp);//a*2
        }
    }
    
    Fp3_set(RES,&res_tmp);
    
    Fp3_clear(&res_tmp);
    Fp3_clear(&in_tmp);
}
void Fp3_sqrt(struct Fp3 *RES,struct Fp3 *A){
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
    }
    
    Fp3_set(RES,&x);
    
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
}
int Fp3_cmp(struct Fp3 *A,struct Fp3 *B){
    if(Fp_cmp(&A->a0,&B->a0)==0 && Fp_cmp(&A->a1,&B->a1)==0 && Fp_cmp(&A->a2,&B->a2)==0){
        return 0;
    }
    return 1;
}
int Fp3_cmp_mpz(struct Fp3 *A,mpz_t B){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    if(Fp_cmp_mpz(&A->a0,B)==0 && Fp_cmp(&A->a1,&tmp.a1)==0 && Fp_cmp(&A->a2,&tmp.a2)==0){
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
        mpz_clear(i);
        return 1;
    }else{
        Fp3_clear(&tmp);
        mpz_clear(i);
        return -1;
    }
}
void Fp3_neg(struct Fp3 *RES,struct Fp3 *a){
    struct Fp3 tmp;
    Fp3_init(&tmp);
    Fp_neg(&tmp.a0,&a->a0);
    Fp_neg(&tmp.a1,&a->a1);
    Fp_neg(&tmp.a2,&a->a2);
    Fp3_set(RES,&tmp);
    Fp3_clear(&tmp);
}

#pragma mark ==Fp6== method implementations
void Fp6_init(struct Fp6 *A){
    Fp3_init(&A->a0);
    Fp3_init(&A->a1);
}
void Fp6_set(struct Fp6 *RES,struct Fp6 *A){
    Fp3_set(&RES->a0,&A->a0);
    Fp3_set(&RES->a1,&A->a1);
}
void Fp6_set_ui(struct Fp6 *A,signed long int B){
    Fp3_set_ui(&A->a0,B);
    Fp3_set_ui(&A->a1,B);
}
void Fp6_random(struct Fp6 *A){
    Fp3_random(&A->a0);
    Fp3_random(&A->a1);
}
void Fp6_clear(struct Fp6 *A){
    Fp3_clear(&A->a0);
    Fp3_clear(&A->a1);
}
void Fp6_printf(struct Fp6 *A){
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->a0.a0.x_0,A->a0.a1.x_0,A->a0.a2.x_0,A->a1.a0.x_0,A->a1.a1.x_0,A->a1.a2.x_0);
}
void Fp6_add(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_add(&tmp.a0,&A->a0,&B->a0);
    Fp3_add(&tmp.a1,&A->a1,&B->a1);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_add_ui(struct Fp6 *RES,struct Fp6 *A,unsigned long int B){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_add_ui(&tmp.a0,&A->a0,B);
    Fp3_add_ui(&tmp.a1,&A->a1,B);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_sub(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_sub(&tmp.a0,&A->a0,&B->a0);
    Fp3_sub(&tmp.a1,&A->a1,&B->a1);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_mul(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B){
    //x^2-c2=0;c2=-1
    //mn= (a0b0 + c2a1b1) + (a0 + a1)(b0 + b1)τ−(a0b0)τ − (a1b1)τ.
    struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5;
    Fp3_init(&tmp1);
    Fp3_init(&tmp2);
    Fp3_init(&tmp3);
    Fp3_init(&tmp4);
    Fp3_init(&tmp5);
    
    struct Fp6 t_RES;
    Fp6_init(&t_RES);
    
    Fp3_mul(&tmp1,&A->a0,&B->a0);//a0*b0
    Fp3_mul(&tmp2,&A->a1,&B->a1);//a1*b1
    Fp3_sub(&t_RES.a0, &tmp1, &tmp2);//c2*a1*b1; a0b0-a1b1
    Fp3_add(&tmp3,&A->a0,&A->a1);//a0+a1
    Fp3_add(&tmp4,&B->a0,&B->a1);//b0+b1
    Fp3_mul(&tmp5,&tmp3,&tmp4);//(a0+a1)(b0+b1)
    Fp3_sub(&t_RES.a1,&tmp5,&tmp1);//(a0+a1)(b0+b1)-a0b0
    Fp3_sub(&t_RES.a1,&t_RES.a1,&tmp2);//(a0+a1)(b0+b1)-a0b0-a1b1
    
    Fp6_set(RES,&t_RES);
    
    Fp3_clear(&tmp1);
    Fp3_clear(&tmp2);
    Fp3_clear(&tmp3);
    Fp3_clear(&tmp4);
    Fp3_clear(&tmp5);
    Fp6_clear(&t_RES);
}
void Fp6_mul_v_prev(struct Fp6 *RES,struct Fp6 *A){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_mul_omega(&tmp.a0,&A->a1);
    Fp3_set(&tmp.a1,&A->a0);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_mul_tau(struct Fp6 *RES,struct Fp6 *A){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_neg(&tmp.a0,&A->a1);
    Fp3_set(&tmp.a1,&A->a0);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_mul_omega(struct Fp6 *RES,struct Fp6 *A){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_mul_omega(&tmp.a0,&A->a0);
    Fp3_mul_omega(&tmp.a1,&A->a1);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_mul_ui(struct Fp6 *RES,struct Fp6 *A,unsigned long int B){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_mul_ui(&tmp.a0,&A->a0,B);
    Fp3_mul_ui(&tmp.a1,&A->a1,B);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_mul_Fp(struct Fp6 *RES,struct Fp6 *A,struct Fp *B){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    Fp3_mul_Fp(&tmp.a0,&A->a0,B);
    Fp3_mul_Fp(&tmp.a1,&A->a1,B);
    
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_invert(struct Fp6 *RES, struct Fp6 *A){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    //(m^q^3 )= a0−a1τ
    Fp3_set(&tmp.a0,&A->a0);
    Fp3_neg(&tmp.a1,&A->a1);
    
    struct Fp3 norm,a,b;
    Fp3_init(&norm);
    Fp3_init(&a);
    Fp3_init(&b);
    
    //n(m)=a0^2+a1^2
    Fp3_mul(&a,&A->a0,&A->a0); // a=a0^2
    Fp3_mul(&b,&A->a1,&A->a1); // b=a1^2
    Fp3_add(&norm,&a,&b);
    Fp_printf(&norm.a0);
    Fp_printf(&norm.a1);
    Fp_printf(&norm.a2);
    Fp3_invert(&norm,&norm);
    
    //m^−1 = n(m)^−1(m^q^3 )
    Fp3_mul(&tmp.a0,&tmp.a0,&norm);
    Fp3_mul(&tmp.a1,&tmp.a1,&norm);
    
    Fp6_set(RES,&tmp);
    
    Fp3_clear(&norm);
    Fp3_clear(&a);
    Fp3_clear(&b);
    Fp6_clear(&tmp);
}
void Fp6_div(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B){
    struct Fp6 tmp,t_RES;
    Fp6_init(&tmp);
    Fp6_init(&t_RES);
    
    Fp6_invert(&tmp,B);
    Fp6_mul(&t_RES,A,&tmp);
    
    Fp6_set(RES,&t_RES);
    
    Fp6_clear(&tmp);
    Fp6_clear(&t_RES);
}
void Fp6_pow(struct Fp6 *RES,struct Fp6 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    
    struct Fp6 res_tmp;
    Fp6_init(&res_tmp);
    Fp6_set(&res_tmp,A);
    
    struct Fp6 in_tmp;
    Fp6_init(&in_tmp);
    Fp6_set(&in_tmp,A);
    
    for(i=length-2;i>=0;i--){
        if(mpz_tstbit(B,i)==1){
            Fp6_mul(&res_tmp,&res_tmp,&res_tmp);//a*2
            Fp6_mul(&res_tmp,&res_tmp,&in_tmp);//*a
        }else{
            Fp6_mul(&res_tmp,&res_tmp,&res_tmp);//a*2
        }
    }
    
    Fp6_set(RES,&res_tmp);
    Fp6_clear(&res_tmp);
    Fp6_clear(&in_tmp);
}
void Fp6_sqrt(struct Fp6 *RES,struct Fp6 *A){
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
        // gmp_printf("%Zd,%Zd,\n",y.a0.a0.a0,y.a0.a1.a0);
        mpz_set_ui(r,m);
        Fp6_mul(&x,&x,&t);
        Fp6_mul(&b,&b,&y);
    }
    
    Fp6_set(RES,&x);
    
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
int Fp6_legendre(struct Fp6 *a){
    mpz_t i,cmp;
    struct Fp6 tmp;
    Fp6_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,prime,6);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp6_pow(&tmp,a,i);
    
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
    if(Fp3_cmp(&A->a0,&B->a0)==0 && Fp3_cmp(&A->a1,&B->a1)==0){
        return 0;
    }
    return 1;
}
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t B){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    if(Fp3_cmp_mpz(&A->a0,B)==0 && Fp3_cmp(&A->a1,&tmp.a1)==0){
        Fp6_clear(&tmp);
        return 0;
    }
    Fp6_clear(&tmp);
    return 1;
}
void Fp6_neg(struct Fp6 *RES,struct Fp6 *a){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    Fp3_neg(&tmp.a0,&a->a0);
    Fp3_neg(&tmp.a1,&a->a1);
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}


#pragma mark ==Fp18== method implementations

void Fp18_init(struct Fp18 *A){
    Fp6_init(&A->m0);
    Fp6_init(&A->m1);
    Fp6_init(&A->m2);
}
void Fp18_set(struct Fp18 *RES,struct Fp18 *A){
    Fp6_set(&RES->m0,&A->m0);
    Fp6_set(&RES->m1,&A->m1);
    Fp6_set(&RES->m2,&A->m2);
}
void Fp18_set_ui(struct Fp18 *A,signed long int B){
    Fp6_set_ui(&A->m0,B);
    Fp6_set_ui(&A->m1,B);
    Fp6_set_ui(&A->m2,B);
}
void Fp18_random(struct Fp18 *A){
    Fp6_random(&A->m0);
    Fp6_random(&A->m1);
    Fp6_random(&A->m2);
}
void Fp18_clear(struct Fp18 *A){
    Fp6_clear(&A->m0);
    Fp6_clear(&A->m1);
    Fp6_clear(&A->m2);
}
void Fp18_printf(struct Fp18 *A){
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->m0.a0.a0.x_0,A->m0.a0.a1.x_0,A->m0.a0.a2.x_0,A->m0.a1.a0.x_0,A->m0.a1.a1.x_0,A->m0.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->m1.a0.a0.x_0,A->m1.a0.a1.x_0,A->m1.a0.a2.x_0,A->m1.a1.a0.x_0,A->m1.a1.a1.x_0,A->m1.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n\n",A->m2.a0.a0.x_0,A->m2.a0.a1.x_0,A->m2.a0.a2.x_0,A->m2.a1.a0.x_0,A->m2.a1.a1.x_0,A->m2.a1.a2.x_0);
}
void Fp18_add(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    
    Fp6_add(&tmp.m0,&A->m0,&B->m0);
    Fp6_add(&tmp.m1,&A->m1,&B->m1);
    Fp6_add(&tmp.m2,&A->m2,&B->m2);
    
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}
void Fp18_add_ui(struct Fp18 *RES,struct Fp18 *A,unsigned long int B){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    
    Fp6_add_ui(&tmp.m0,&A->m0,B);
    Fp6_add_ui(&tmp.m1,&A->m1,B);
    Fp6_add_ui(&tmp.m2,&A->m2,B);
    
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}
void Fp18_sub(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    
    Fp6_sub(&tmp.m0,&A->m0,&B->m0);
    Fp6_sub(&tmp.m1,&A->m1,&B->m1);
    Fp6_sub(&tmp.m2,&A->m2,&B->m2);
    
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}
void Fp18_mul(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B){
    struct Fp6 A0,A1,A2,tmp_A01,tmp_A12,tmp_A20,tmp_B01,tmp_B12,tmp_B20,A3,A5,A4,tmp;
    struct Fp18 t_RES;
    Fp6_init(&A0);
    Fp6_init(&A1);
    Fp6_init(&A2);
    Fp6_init(&tmp_A01);
    Fp6_init(&tmp_A12);
    Fp6_init(&tmp_A20);
    Fp6_init(&tmp_B01);
    Fp6_init(&tmp_B12);
    Fp6_init(&tmp_B20);
    Fp6_init(&A3);
    Fp6_init(&A5);
    Fp6_init(&A4);
    Fp6_init(&tmp);
    Fp18_init(&t_RES);
    
    Fp6_mul(&A0,&A->m0,&B->m0);//a0*b0
    Fp6_mul(&A1,&A->m1,&B->m1);//a1*b1
    Fp6_mul(&A2,&A->m2,&B->m2);//a2*b2
    
    Fp6_add(&tmp_A01,&A->m0,&A->m1);//a0+a1
    Fp6_add(&tmp_A12,&A->m1,&A->m2);//a1+a2
    Fp6_add(&tmp_A20,&A->m0,&A->m2);//a2+a0
    Fp6_add(&tmp_B01,&B->m0,&B->m1);//b0+b1
    Fp6_add(&tmp_B12,&B->m1,&B->m2);//b1+b2
    Fp6_add(&tmp_B20,&B->m0,&B->m2);//b2+b0
    
    Fp6_mul(&A3,&tmp_A01,&tmp_B01);//(a0+a1)(b0+b1)
    Fp6_mul(&A4,&tmp_A20,&tmp_B20);//(a2+a0)(b2+b0)
    Fp6_mul(&A5,&tmp_A12,&tmp_B12);//(a1+a2)(b1+b2)
    
    Fp6_sub(&A5,&A5,&A1);
    Fp6_sub(&A5,&A5,&A2);
    Fp6_mul_omega(&tmp,&A5);//wt3
    Fp6_add(&t_RES.m0,&A0,&tmp);
    
    Fp6_sub(&A3,&A3,&A0);
    Fp6_sub(&A3,&A3,&A1);//t1
    Fp6_mul_omega(&tmp,&A2);
    Fp6_add(&t_RES.m1,&tmp,&A3);
    
    Fp6_sub(&A4,&A4,&A0);
    Fp6_sub(&A4,&A4,&A2);
    Fp6_add(&t_RES.m2,&A1,&A4);
    
    Fp18_set(RES,&t_RES);
    
    Fp6_clear(&A0);
    Fp6_clear(&A1);
    Fp6_clear(&A2);
    Fp6_clear(&tmp_A01);
    Fp6_clear(&tmp_A12);
    Fp6_clear(&tmp_A20);
    Fp6_clear(&tmp_B01);
    Fp6_clear(&tmp_B12);
    Fp6_clear(&tmp_B20);
    Fp6_clear(&A3);
    Fp6_clear(&A5);
    Fp6_clear(&A4);
    Fp6_clear(&tmp);
}
void Fp18_mul_ui(struct Fp18 *RES,struct Fp18 *A,unsigned long int B){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    
    Fp6_mul_ui(&tmp.m0,&A->m0,B);
    Fp6_mul_ui(&tmp.m1,&A->m1,B);
    Fp6_mul_ui(&tmp.m2,&A->m2,B);
    
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}
void Fp18_mul_Fp(struct Fp18 *RES,struct Fp18 *A,struct Fp *B){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    
    Fp6_mul_Fp(&tmp.m0,&A->m0,B);
    Fp6_mul_Fp(&tmp.m1,&A->m1,B);
    Fp6_mul_Fp(&tmp.m2,&A->m2,B);
    
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}
void Fp18_frobenius_map(struct Fp18 *RES,struct Fp18 *A){
    struct Fp18 t_RES;
    Fp18_init(&t_RES);
    struct Fp tmp;
    Fp_init(&tmp);
    Fp_set_mpz(&tmp,c1_omega);
    
    Fp6_set(&t_RES.m0, &A->m0);
    Fp6_mul_Fp(&t_RES.m1, &A->m1, &tmp);
    Fp_set_mpz(&tmp, c1_omega_bar);
    Fp6_mul_Fp(&t_RES.m2, &A->m2, &tmp);
    
    Fp18_set(RES,&t_RES);
    Fp18_clear(&t_RES);
    Fp_clear(&tmp);
}
void Fp18_invert(struct Fp18 *RES, struct Fp18 *A){
    struct Fp18 t_RES,Aq,Aqsq,t_Aq,AqAqsq,t_norm;
    Fp18_init(&t_RES);
    Fp18_init(&Aq);
    Fp18_init(&Aqsq);
    Fp18_init(&t_Aq);
    Fp18_init(&AqAqsq);
    Fp18_init(&t_norm);
    
    Fp18_frobenius_map(&Aq, A);
    Fp18_set(&t_Aq, &Aq);
    Fp18_frobenius_map(&Aqsq,&t_Aq);
    Fp18_mul(&AqAqsq, &Aq, &Aqsq);
    Fp18_mul(&t_norm, &AqAqsq, A);
    
    printf("M2 \n");
    Fp18_printf(&t_norm);
//    printf("M2 \n");
//    //    printf("SIZE %d\n",(int)mpz_sizeinbase(t_norm.m1.a0.a0.x_0,2));
//    Fp6_printf(&t_norm.m2);
    
    struct Fp6 norm,T0;
    Fp6_init(&norm);
    Fp6_init(&T0);
    mpz_t zero;
    mpz_init(zero);
    mpz_set_ui(zero, 0);
    
    if (Fp6_cmp_mpz(&t_norm.m1, zero)!= 0 || Fp6_cmp_mpz(&t_norm.m2, zero)!= 0)
    {
        printf("Fp18 Norm func calculation goes wrong\n");
        exit(0);
    }
    
    Fp6_set(&T0, &t_norm.m0);
    Fp6_invert(&norm, &T0);
    Fp6_mul(&t_RES.m0, &AqAqsq.m0, &norm);
    Fp6_mul(&t_RES.m1, &AqAqsq.m1, &norm);
    Fp6_mul(&t_RES.m2, &AqAqsq.m2, &norm);
    
    Fp18_set(RES,&t_RES);
    
    Fp18_clear(&t_RES);
    Fp18_clear(&Aq);
    Fp18_clear(&Aqsq);
    Fp18_clear(&AqAqsq);
    Fp18_clear(&t_norm);
    Fp6_clear(&T0);
    Fp6_clear(&norm);
}

void Fp18_div(struct Fp18 *RES,struct Fp18 *A,struct Fp18 *B){
    struct Fp18 tmp,t_RES;
    Fp18_init(&tmp);
    Fp18_init(&t_RES);
    
    Fp18_invert(&tmp,B);
    Fp18_mul(&t_RES,A,&tmp);
    
    Fp18_set(RES,&t_RES);
    
    Fp18_clear(&tmp);
    Fp18_clear(&t_RES);
}
void Fp18_pow(struct Fp18 *RES,struct Fp18 *A,mpz_t B){
    int i;
    int r;//bit数
    r= (int)mpz_sizeinbase(B,2);
    //printf("r= %d\n",r);
    
    struct Fp18 res_tmp;
    Fp18_init(&res_tmp);
    Fp18_set(&res_tmp,A);
    
    struct Fp18 in_tmp;
    Fp18_init(&in_tmp);
    Fp18_set(&in_tmp,A);
    
    for(i=r-2;i>=0;i--){
        if(mpz_tstbit(B,i)==1){
            Fp18_mul(&res_tmp,&res_tmp,&res_tmp);//a*2
            Fp18_mul(&res_tmp,&res_tmp,&in_tmp);//*a
        }else{
            Fp18_mul(&res_tmp,&res_tmp,&res_tmp);//a*2
        }
    }
    
    Fp18_set(RES,&res_tmp);
    
    Fp18_clear(&res_tmp);
    Fp18_clear(&in_tmp);
}
void Fp18_sqrt(struct Fp18 *RES,struct Fp18 *A){
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
        // printf("%d\n",Fp18_legendre(&n));
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
    }
    
    Fp18_set(RES,&x);
    
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
}
int Fp18_cmp(struct Fp18 *A,struct Fp18 *B){
    if(Fp6_cmp(&A->m0,&B->m0)==0 && Fp6_cmp(&A->m1,&B->m1)==0 && Fp6_cmp(&A->m2,&B->m2)==0){
        return 0;
    }
    return 1;
}
int Fp18_cmp_mpz(struct Fp18 *A,mpz_t B){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    if(Fp6_cmp_mpz(&A->m0,B)==0 && Fp6_cmp(&A->m1,&tmp.m1)==0 && Fp6_cmp(&A->m2,&tmp.m2)==0){
        Fp18_clear(&tmp);
        return 0;
    }
    Fp18_clear(&tmp);
    return 1;
}
int Fp18_legendre(struct Fp18 *a){
    mpz_t i,cmp;
    struct Fp18 tmp;
    mpz_init(i);
    mpz_init(cmp);
    Fp18_init(&tmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,prime,18);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp18_pow(&tmp,a,i);
    
    if((Fp18_cmp_mpz(&tmp,cmp))==0){
        Fp18_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp18_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
void Fp18_neg(struct Fp18 *RES,struct Fp18 *a){
    struct Fp18 tmp;
    Fp18_init(&tmp);
    Fp6_neg(&tmp.m0,&a->m0);
    Fp6_neg(&tmp.m1,&a->m1);
    Fp6_neg(&tmp.m2,&a->m2);
    Fp18_set(RES,&tmp);
    Fp18_clear(&tmp);
}

#pragma mark **EFp** method implementations

void EFp_init(struct EFp *A){
    Fp_init(&A->px);
    Fp_init(&A->py);
    A->isInfinity=FALSE;
}
void EFp_set(struct EFp *A,struct EFp *B){
    Fp_set(&A->px,&B->px);
    Fp_set(&A->py,&B->py);
    A->isInfinity=B->isInfinity;
}
void EFp_set_infity(struct EFp *A){
    Fp_set_ui(&A->px,0);
    Fp_set_ui(&A->py,0);
    A->isInfinity=TRUE;
}
void EFp_clear(struct EFp *A){
    Fp_clear(&A->px);
    Fp_clear(&A->py);
}
void EFp_printf(struct EFp *A){
    gmp_printf("(%Zd,%Zd)\n",A->px.x_0,A->py.x_0);
}
void EFp_SCM(struct EFp *RES, struct EFp *P,mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    // int eca=0,ecd=0;
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
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
    EFp_set(RES,&Q);
    // printf("%d,%d\n",eca,ecd);
    EFp_clear(&Q);
    EFp_clear(&R);
    return;
}
void EFp_ECD(struct EFp *RES, struct EFp *P){
    if(P->isInfinity==TRUE){
        EFp_set(RES,P);
        return;
    }
    if(mpz_sgn(P->py.x_0)==0){//P.y==0
        EFp_set_infity(RES);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_RES;
    Fp_init(&x);
    Fp_init(&lambda);
    Fp_init(&tmp);
    Fp_init(&y);
    EFp_init(&t_RES);
    
    Fp_mul(&x,&P->px,&P->px);
    Fp_add(&tmp,&x,&x);
    Fp_add(&x,&tmp,&x);//3x^2
    Fp_add(&y,&P->py,&P->py);//2y
    
    Fp_div(&lambda,&x,&y);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P->px,&P->px);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P->px,&x);
    
    Fp_set(&t_RES.px,&x);
    
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_RES.py,&tmp,&P->py);
    
    EFp_set(RES,&t_RES);
    
    Fp_clear(&x);
    Fp_clear(&lambda);
    Fp_clear(&y);
    Fp_clear(&tmp);
    EFp_clear(&t_RES);
}
void EFp_ECA(struct EFp *RES, struct EFp *P, struct EFp *Q){
    if(Q->isInfinity==TRUE){//if P2==inf
        EFp_set(RES,P);
        return;
    }
    else if(P->isInfinity==TRUE){//if P1==inf
        EFp_set(RES,Q);
        return;
    }
    else if(Fp_cmp(&P->px,&Q->px)==0&&Fp_cmp(&P->py,&Q->py)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp_set_infity(RES);
        return;
    }
    else if(EFp_cmp(P,Q)==0){ // P=Q
        EFp_ECD(RES,P);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_RES;
    
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&lambda);
    Fp_init(&tmp);
    EFp_init(&t_RES);
    
    Fp_sub(&x,&Q->px,&P->px);
    Fp_sub(&y,&Q->py,&P->py);
    Fp_div(&lambda,&y,&x);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P->px,&Q->px);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P->px,&x);
    Fp_set(&t_RES.px,&x);
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_RES.py,&tmp,&P->py);
    
    EFp_set(RES,&t_RES);
    
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&lambda);
    Fp_clear(&tmp);
    EFp_clear(&t_RES);
}
int EFp_cmp(struct EFp *A,struct EFp *B){
    if(Fp_cmp(&A->px,&B->px)==0 && Fp_cmp(&A->py,&B->py)==0){
        return 0;
    }
    return 1;
}

#pragma mark **EFp3** method implementations
void EFp3_init(struct EFp3 *A){
    Fp3_init(&A->p3x);
    Fp3_init(&A->p3y);
    A->isInfinity=FALSE;
}
void EFp3_set(struct EFp3 *A,struct EFp3 *B){
    Fp3_set(&A->p3x,&B->p3x);
    Fp3_set(&A->p3y,&B->p3y);
    A->isInfinity=B->isInfinity;
}
void EFp3_set_infity(struct EFp3 *A){
    Fp3_set_ui(&A->p3x,0);
    Fp3_set_ui(&A->p3y,0);
    A->isInfinity=TRUE;
}
void EFp3_clear(struct EFp3 *A){
    Fp3_clear(&A->p3x);
    Fp3_clear(&A->p3y);
}
void EFp3_printf(struct EFp3 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",A->p3x.a0.x_0,A->p3x.a1.x_0,A->p3x.a2.x_0,A->p3y.a0.x_0,A->p3y.a1.x_0,A->p3y.a2.x_0);
}
void EFp3_ECD(struct EFp3 *RES, struct EFp3 *P){
    if(P->isInfinity==TRUE){
        EFp3_set(RES,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp3_cmp_mpz(&P->p3y,cmp)==0){//P.y==0
        EFp3_set_infity(RES);
        return;
    }
    
    struct Fp3 x,y,lambda,tmp;
    struct EFp3 t_RES;
    Fp3_init(&x);
    Fp3_init(&lambda);
    Fp3_init(&tmp);
    Fp3_init(&y);
    EFp3_init(&t_RES);
    
    Fp3_mul(&x,&P->p3x,&P->p3x);
    Fp3_add(&tmp,&x,&x);
    Fp3_add(&x,&tmp,&x);
    Fp3_add(&y,&P->p3y,&P->p3y);
    Fp3_div(&lambda,&x,&y);
    Fp3_mul(&tmp,&lambda,&lambda);
    Fp3_add(&x,&P->p3x,&P->p3x);
    Fp3_sub(&x,&tmp,&x);
    Fp3_sub(&tmp,&P->p3x,&x);
    Fp3_set(&t_RES.p3x,&x);
    Fp3_mul(&tmp,&tmp,&lambda);
    Fp3_sub(&t_RES.p3y,&tmp,&P->p3y);
    
    EFp3_set(RES,&t_RES);
    
    Fp3_clear(&x);
    Fp3_clear(&lambda);
    Fp3_clear(&y);
    Fp3_clear(&tmp);
    EFp3_clear(&t_RES);
}
void EFp3_ECA(struct EFp3 *RES, struct EFp3 *P, struct EFp3 *Q){
    if(Q->isInfinity==TRUE){//if P2==inf
        EFp3_set(RES,P);
        return;
    }
    else if(P->isInfinity==TRUE){//if P1==inf
        EFp3_set(RES,Q);
        return;
    }
    else if(Fp3_cmp(&P->p3x,&Q->p3x)==0&&Fp3_cmp(&P->p3y,&Q->p3y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp3_set_infity(RES);
        return;
    }
    else if(EFp3_cmp(P,Q)==0){ // P=Q
        EFp3_ECD(RES,P);
        return;
    }
    
    struct Fp3 x,y,lambda,tmp;
    struct EFp3 t_RES;
    
    Fp3_init(&x);
    Fp3_init(&y);
    Fp3_init(&lambda);
    Fp3_init(&tmp);
    EFp3_init(&t_RES);
    
    Fp3_sub(&x,&Q->p3x,&P->p3x);
    Fp3_sub(&y,&Q->p3y,&P->p3y);
    Fp3_div(&lambda,&y,&x);
    Fp3_mul(&tmp,&lambda,&lambda);
    Fp3_add(&x,&P->p3x,&Q->p3x);
    Fp3_sub(&x,&tmp,&x);
    Fp3_sub(&tmp,&P->p3x,&x);
    Fp3_set(&t_RES.p3x,&x);
    Fp3_mul(&tmp,&tmp,&lambda);
    Fp3_sub(&t_RES.p3y,&tmp,&P->p3y);
    
    EFp3_set(RES,&t_RES);
    
    Fp3_clear(&x);
    Fp3_clear(&y);
    Fp3_clear(&lambda);
    Fp3_clear(&tmp);
    EFp3_clear(&t_RES);
}
void EFp3_SCM(struct EFp3 *RES, struct EFp3 *P, mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    // printf("%d\n",length);
    mpz_get_str(r_binary,2,scalar);
    struct EFp3 Q,R;
    EFp3_init(&Q);
    EFp3_set(&Q,P);
    EFp3_init(&R);
    for(i=1;r_binary[i]!='\0';i++){
        EFp3_ECD(&Q,&Q);
        if(r_binary[i]=='1'){
            EFp3_ECA(&Q,&Q,P);
        }
    }
    EFp3_set(RES,&Q);
    
    EFp3_clear(&Q);
    EFp3_clear(&R);
    return;
}
int EFp3_cmp(struct EFp3 *A,struct EFp3 *B){
    if(Fp3_cmp(&A->p3x,&B->p3x)==0 && Fp3_cmp(&A->p3y,&B->p3y)==0){
        return 0;
    }
    return 1;
}
void EFp3_set_EFp(struct EFp3 *A,struct EFp *B){
    Fp3_set_ui(&A->p3x,0);
    Fp3_set_ui(&A->p3y,0);
    
    Fp_set(&A->p3x.a0,&B->px);
    Fp_set(&A->p3y.a0,&B->py);
    A->isInfinity=B->isInfinity;
}
void EFp3_neg(struct EFp3 *RES, struct EFp3 *A){
    struct EFp3 tmp;
    EFp3_init(&tmp);
    Fp3_neg(&tmp.p3y,&A->p3y);
    Fp3_set(&tmp.p3x,&A->p3x);
    
    EFp3_set(RES,&tmp);
    EFp3_clear(&tmp);
}


#pragma mark **EFp6** method implementations
void EFp6_init(struct EFp6 *A){
    Fp6_init(&A->p6x);
    Fp6_init(&A->p6y);
    A->isInfinity=FALSE;
}
void EFp6_set(struct EFp6 *A,struct EFp6 *B){
    Fp6_set(&A->p6x,&B->p6x);
    Fp6_set(&A->p6y,&B->p6y);
    A->isInfinity=B->isInfinity;
}
void EFp6_set_infity(struct EFp6 *A){
    Fp6_set_ui(&A->p6x,0);
    Fp6_set_ui(&A->p6y,0);
    A->isInfinity=TRUE;
}
void EFp6_clear(struct EFp6 *A){
    Fp6_clear(&A->p6x);
    Fp6_clear(&A->p6y);
}
void EFp6_printf(struct EFp6 *A){
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->p6x.a0.a0.x_0,A->p6x.a0.a1.x_0,A->p6x.a0.a2.x_0,A->p6x.a1.a0.x_0,A->p6x.a1.a1.x_0,A->p6x.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->p6y.a0.a0.x_0,A->p6y.a0.a1.x_0,A->p6y.a0.a2.x_0,A->p6y.a1.a0.x_0,A->p6y.a1.a1.x_0,A->p6y.a1.a2.x_0);
}
void EFp6_ECD(struct EFp6 *RES, struct EFp6 *P){
    if(P->isInfinity==TRUE){
        EFp6_set(RES,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp6_cmp_mpz(&P->p6y,cmp)==0){//P.y==0
        EFp6_set_infity(RES);
        return;
    }
    
    struct Fp6 x,y,lambda,tmp;
    struct EFp6 t_RES;
    Fp6_init(&x);
    Fp6_init(&lambda);
    Fp6_init(&tmp);
    Fp6_init(&y);
    EFp6_init(&t_RES);
    
    Fp6_mul(&x,&P->p6x,&P->p6x);
    Fp6_add(&tmp,&x,&x);
    Fp6_add(&x,&tmp,&x);
    Fp6_add(&y,&P->p6y,&P->p6y);
    Fp6_div(&lambda,&x,&y);
    Fp6_mul(&tmp,&lambda,&lambda);
    Fp6_add(&x,&P->p6x,&P->p6x);
    Fp6_sub(&x,&tmp,&x);
    Fp6_sub(&tmp,&P->p6x,&x);
    Fp6_set(&t_RES.p6x,&x);
    Fp6_mul(&tmp,&tmp,&lambda);
    Fp6_sub(&t_RES.p6y,&tmp,&P->p6y);
    
    EFp6_set(RES,&t_RES);
    
    Fp6_clear(&x);
    Fp6_clear(&lambda);
    Fp6_clear(&y);
    Fp6_clear(&tmp);
    EFp6_clear(&t_RES);
}
void EFp6_ECA(struct EFp6 *RES, struct EFp6 *P1, struct EFp6 *P2){
    if(P2->isInfinity==TRUE){//if P2==inf
        EFp6_set(RES,P1);
        return;
    }
    else if(P1->isInfinity==TRUE){//if P1==inf
        EFp6_set(RES,P2);
        return;
    }
    else if(Fp6_cmp(&P1->p6x,&P2->p6x)==0&&Fp6_cmp(&P1->p6y,&P2->p6y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp6_set_infity(RES);
        return;
    }
    else if(EFp6_cmp(P1,P2)==0){ // P=Q
        EFp6_ECD(RES,P1);
        return;
    }
    
    struct Fp6 x,y,lambda,tmp;
    struct EFp6 t_RES;
    
    Fp6_init(&x);
    Fp6_init(&y);
    Fp6_init(&lambda);
    Fp6_init(&tmp);
    EFp6_init(&t_RES);
    
    Fp6_sub(&x,&P2->p6x,&P1->p6x);
    Fp6_sub(&y,&P2->p6y,&P1->p6y);
    Fp6_div(&lambda,&y,&x);
    Fp6_mul(&tmp,&lambda,&lambda);
    Fp6_add(&x,&P1->p6x,&P2->p6x);
    Fp6_sub(&x,&tmp,&x);
    Fp6_sub(&tmp,&P1->p6x,&x);
    Fp6_set(&t_RES.p6x,&x);
    Fp6_mul(&tmp,&tmp,&lambda);
    Fp6_sub(&t_RES.p6y,&tmp,&P1->p6y);
    
    EFp6_set(RES,&t_RES);
    
    Fp6_clear(&x);
    Fp6_clear(&y);
    Fp6_clear(&lambda);
    Fp6_clear(&tmp);
    EFp6_clear(&t_RES);
}
void EFp6_SCM(struct EFp6 *RES,struct EFp6 *P,mpz_t S){
    int i,length;
    length= (int)mpz_sizeinbase(S,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,S);
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
    EFp6_set(RES,&Q);
    
    EFp6_clear(&Q);
    EFp6_clear(&R);
    return;
}
int EFp6_cmp(struct EFp6 *A,struct EFp6 *B){
    if(Fp6_cmp(&A->p6x,&B->p6x)==0 && Fp6_cmp(&A->p6y,&B->p6y)==0){
        return 0;
    }
    return 1;
}

#pragma mark **EFp18** method implementations
void EFp18_init(struct EFp18 *A){
    Fp18_init(&A->p18x);
    Fp18_init(&A->p18y);
    A->isInfinity=FALSE;
}
void EFp18_set(struct EFp18 *A,struct EFp18 *B){
    Fp18_set(&A->p18x,&B->p18x);
    Fp18_set(&A->p18y,&B->p18y);
    A->isInfinity=B->isInfinity;
}
void EFp18_set_infity(struct EFp18 *A){
    Fp18_set_ui(&A->p18x,0);
    Fp18_set_ui(&A->p18y,0);
    A->isInfinity=TRUE;
}
void EFp18_clear(struct EFp18 *A){
    Fp18_clear(&A->p18x);
    Fp18_clear(&A->p18y);
}
void EFp18_printf(struct EFp18 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->p18x.m0.a0.a0.x_0,A->p18x.m0.a0.a1.x_0,A->p18x.m0.a0.a2.x_0,A->p18x.m0.a1.a0.x_0,A->p18x.m0.a1.a1.x_0,A->p18x.m0.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->p18x.m1.a0.a0.x_0,A->p18x.m1.a0.a1.x_0,A->p18x.m1.a0.a2.x_0,A->p18x.m1.a1.a0.x_0,A->p18x.m1.a1.a1.x_0,A->p18x.m1.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n",A->p18x.m2.a0.a0.x_0,A->p18x.m2.a0.a1.x_0,A->p18x.m2.a0.a2.x_0,A->p18x.m2.a1.a0.x_0,A->p18x.m2.a1.a1.x_0,A->p18x.m2.a1.a2.x_0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->p18y.m0.a0.a0.x_0,A->p18y.m0.a0.a1.x_0,A->p18y.m0.a0.a2.x_0,A->p18y.m0.a1.a0.x_0,A->p18y.m0.a1.a1.x_0,A->p18y.m0.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",A->p18y.m1.a0.a0.x_0,A->p18y.m1.a0.a1.x_0,A->p18y.m1.a0.a2.x_0,A->p18y.m1.a1.a0.x_0,A->p18y.m1.a1.a1.x_0,A->p18y.m1.a1.a2.x_0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,%Zd,%Zd)\n\n",A->p18y.m2.a0.a0.x_0,A->p18y.m2.a0.a1.x_0,A->p18y.m2.a0.a2.x_0,A->p18y.m2.a1.a0.x_0,A->p18y.m2.a1.a1.x_0,A->p18y.m2.a1.a2.x_0);
}
void EFp18_ECD(struct EFp18 *RES, struct EFp18 *P){
    if(P->isInfinity==TRUE){
        EFp18_set(RES,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp18_cmp_mpz(&P->p18y,cmp)==0){//P.y==0
        EFp18_set_infity(RES);
        return;
    }
    
    struct Fp18 x,y,lambda,tmp;
    struct EFp18 t_RES;
    Fp18_init(&x);
    Fp18_init(&lambda);
    Fp18_init(&tmp);
    Fp18_init(&y);
    EFp18_init(&t_RES);
    
    Fp18_mul(&x,&P->p18x,&P->p18x);
    Fp18_add(&tmp,&x,&x);
    Fp18_add(&x,&tmp,&x);
    Fp18_add(&y,&P->p18y,&P->p18y);
    Fp18_div(&lambda,&x,&y);
    Fp18_mul(&tmp,&lambda,&lambda);
    Fp18_add(&x,&P->p18x,&P->p18x);
    Fp18_sub(&x,&tmp,&x);
    Fp18_sub(&tmp,&P->p18x,&x);
    Fp18_set(&t_RES.p18x,&x);
    Fp18_mul(&tmp,&tmp,&lambda);
    Fp18_sub(&t_RES.p18y,&tmp,&P->p18y);
    
    EFp18_set(RES,&t_RES);
    
    Fp18_clear(&x);
    Fp18_clear(&lambda);
    Fp18_clear(&y);
    Fp18_clear(&tmp);
    EFp18_clear(&t_RES);
}
void EFp18_ECA(struct EFp18 *RES, struct EFp18 *P1, struct EFp18 *P2){
    if(P2->isInfinity==TRUE){//if P2==inf
        EFp18_set(RES,P1);
        return;
    }
    else if(P1->isInfinity==TRUE){//if P1==inf
        EFp18_set(RES,P2);
        return;
    }
    else if(Fp18_cmp(&P1->p18x,&P2->p18x)==0&&Fp18_cmp(&P1->p18y,&P2->p18y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp18_set_infity(RES);
        return;
    }
    else if(EFp18_cmp(P1,P2)==0){ // P=Q
        EFp18_ECD(RES,P1);
        return;
    }
    
    struct Fp18 x,y,lambda,tmp;
    struct EFp18 t_RES;
    
    Fp18_init(&x);
    Fp18_init(&y);
    Fp18_init(&lambda);
    Fp18_init(&tmp);
    EFp18_init(&t_RES);
    
    Fp18_sub(&x,&P2->p18x,&P1->p18x);
    Fp18_sub(&y,&P2->p18y,&P1->p18y);
    Fp18_div(&lambda,&y,&x);
    Fp18_mul(&tmp,&lambda,&lambda);
    Fp18_add(&x,&P1->p18x,&P2->p18x);
    Fp18_sub(&x,&tmp,&x);
    Fp18_sub(&tmp,&P1->p18x,&x);
    Fp18_set(&t_RES.p18x,&x);
    Fp18_mul(&tmp,&tmp,&lambda);
    Fp18_sub(&t_RES.p18y,&tmp,&P1->p18y);
    
    EFp18_set(RES,&t_RES);
    
    Fp18_clear(&x);
    Fp18_clear(&y);
    Fp18_clear(&lambda);
    Fp18_clear(&tmp);
    EFp18_clear(&t_RES);
}
void EFp18_SCM(struct EFp18 *RES,struct EFp18 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp18 Q,R;
    EFp18_init(&Q);
    EFp18_set(&Q,P);
    EFp18_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp18_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp18_ECA(&Q,&Q,P);
        }
    }
    EFp18_set(RES,&Q);
    
    EFp18_clear(&Q);
    EFp18_clear(&R);
    return;
}
int EFp18_cmp(struct EFp18 *A,struct EFp18 *B){
    if(Fp18_cmp(&A->p18x,&B->p18x)==0 && Fp18_cmp(&A->p18y,&B->p18y)==0){
        return 0;
    }
    return 1;
}
void EFp18_set_EFp(struct EFp18 *A,struct EFp *B){
    Fp18_set_ui(&A->p18x,0);
    Fp18_set_ui(&A->p18y,0);
    
    Fp_set(&A->p18x.m0.a0.a0,&B->px);
    Fp_set(&A->p18y.m0.a0.a0,&B->py);
    A->isInfinity=B->isInfinity;
}


#pragma mark Parameters method implementations
void generate_parameters(void){
    mpz_t p_tmp,r_tmp,t_tmp,mod4,p_m1;
    mpz_t xpow2,xpow3,xpow4,xpow5,xpow6,xpow7,xpow8;
    mpz_t tmp1,tmp2;
    
    mpz_init(p_tmp);
    mpz_init(r_tmp);
    mpz_init(t_tmp);
    mpz_init(mod4);
    mpz_init(p_m1);
    
    mpz_init(xpow2);
    mpz_init(xpow3);
    mpz_init(xpow4);
    mpz_init(xpow5);
    mpz_init(xpow6);
    mpz_init(xpow7);
    mpz_init(xpow8);
    
    mpz_init(tmp1);
    mpz_init(tmp2);
    
    mpz_tdiv_r_ui(tmp1,X,42); //tmp1 = x/42 | x = 14(mod 42)
    while(mpz_cmp_ui(tmp1,14)!=0){
        mpz_add_ui(X,X,1);
        mpz_tdiv_r_ui(tmp1,X,42);
    }
    
    mpz_mul(xpow2,X,X);
    mpz_mul(xpow3,xpow2,X);
    mpz_mul(xpow4,xpow2,xpow2);
    mpz_mul(xpow5,xpow4,X);
    mpz_mul(xpow6,xpow3,xpow3);
    mpz_mul(xpow7,xpow6,X);
    mpz_mul(xpow8,xpow4,xpow4);
    
    //trace t(x)=1/7(x^4+16x+7)
    mpz_mul_ui(tmp1,X,16);
    mpz_add_ui(tmp2,xpow4,7);
    mpz_add(t_tmp,tmp1,tmp2);
    mpz_div_ui(t_trace,t_tmp,7);
    
    //order r(x)=x^6+37x^3+343
    mpz_mul_ui(tmp1,xpow3,37);
    mpz_add_ui(tmp2,xpow6,343);
    mpz_add(r_tmp,tmp1,tmp2);
    mpz_set(r_order,r_tmp);
    
CHECK_PRIME:
    //r = p+1-t
    mpz_add_ui(r_order_EFp,prime,1);
    mpz_sub(r_order_EFp,r_order_EFp,t_trace);
    
    mpz_mod_ui(mod4,prime,4);
    unsigned long m4 = mpz_get_ui(mod4);
    
    mpz_sub_ui(p_m1,prime,1);
    mpz_cdiv_r_ui(tmp1,p_m1,3);
    unsigned long r = mpz_get_ui(tmp1);
    
    if(mpz_probab_prime_p(prime,25) >=0 && m4 == 3 && r == 0){
        
        mpz_set(p_tmp,prime);
        mpz_t prime_m1,p_m1_d3,c2_tmp;
        mpz_init(prime_m1);
        mpz_init(c2_tmp);
        mpz_init(p_m1_d3);
        
        mpz_set_ui(c1_leg,C1);
        mpz_set_ui(c2_tmp,C2);
        
        mpz_sub_ui(prime_m1,p_tmp,1);
        mpz_cdiv_q_ui(p_m1_d3,prime_m1,3);
        mpz_powm(c1_leg,c1_leg,p_m1_d3,p_tmp);
        mpz_powm_ui(c1_leg_bar,c1_leg,2,p_tmp);
        printf("c1_leg = %dbit \n",(int)mpz_sizeinbase(c1_leg,2));
        gmp_printf("C1 = %Zd,\nC1^2 = %Zd\n",c1_leg,c1_leg_bar);
        
        mpz_t p6_m1,q_p_m1_d9;
        mpz_init(p6_m1);
        mpz_init(q_p_m1_d9);
        
        mpz_set_ui(c1_omega,C1);
        
        mpz_pow_ui(p6_m1,p_tmp,6);
        mpz_sub_ui(p6_m1,p6_m1,1);
//        mpz_mod(p6_m1,p6_m1,prime);
//         mpz_cdiv_r_ui(q_p_m1_d9,p6_m1,9);
//        gmp_printf("r= %Zd\n",q_p_m1_d9);
        
        mpz_cdiv_q_ui(q_p_m1_d9,p6_m1,9);
        mpz_powm(c1_omega,c1_omega,q_p_m1_d9,p_tmp);
        mpz_powm_ui(c1_omega_bar,c1_omega,2,p_tmp);
        gmp_printf("C1_omg = %Zd,\nC1_omg^2 = %Zd\n",c1_omega,c1_omega_bar);
        
        int leg = mpz_legendre(c2_tmp,p_tmp);
        
        mpz_clear(prime_m1);
        mpz_clear(p_m1_d3);
        mpz_clear(p6_m1);
        mpz_clear(q_p_m1_d9);
        
        if ((mpz_cmp_ui(c1_leg,0) != 0 &&  mpz_cmp_ui(c1_leg,1) != 0) && leg == -1)
        {
            gmp_printf("%Zd\n",X);
            mpz_set(prime,p_tmp);
            gmp_printf("p:%Zd\n",p_tmp);
            printf("p = %dbit\n\n\n",(int)mpz_sizeinbase(p_tmp,2));
        }
    }
    else
    {
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
        
        if(mpz_probab_prime_p(prime,25)==0){
            gmp_printf("p:%Zd\n",prime);
            printf("not  prime number!\n");
            exit(0);
        }
        else
        {
            goto CHECK_PRIME;
        }
    }
    
    mpz_clear(p_tmp);
    mpz_clear(r_tmp);
    mpz_clear(t_tmp);
    mpz_clear(mod4);
    mpz_clear(p_m1);
    mpz_clear(xpow2);
    mpz_clear(xpow3);
    mpz_clear(xpow4);
    mpz_clear(xpow5);
    mpz_clear(xpow6);
    mpz_clear(xpow7);
    mpz_clear(xpow8);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    return;
}

void get_C1_C1bar(){
    mpz_t prime_m1,p_m1_d3;
    mpz_init(prime_m1);
    mpz_init(p_m1_d3);
    
    mpz_set_ui(c1_leg,C1);
    
    mpz_sub_ui(prime_m1,prime,1);
    mpz_cdiv_q_ui(p_m1_d3,prime_m1,3);
    mpz_powm(c1_leg,c1_leg,p_m1_d3,prime);
    mpz_powm_ui(c1_leg_bar,c1_leg,2,prime);
    
    gmp_printf("C1=%Zd,\n C1^2=%Zd\n",c1_leg,c1_leg_bar);
    
    mpz_clear(prime_m1);
    mpz_clear(p_m1_d3);
}

void get_C1omega_C1omegabar(){
    mpz_t p6_m1,q_p_m1_d9;
    mpz_init(p6_m1);
    mpz_init(q_p_m1_d9);
    
    mpz_set_ui(c1_omega,C1);
    
    mpz_pow_ui(p6_m1,prime,6);
    //    gmp_printf("p6_m1=%Zd,\n c1_omega=%Zd\n",p6_m1,c1_omega);
    mpz_sub_ui(p6_m1,p6_m1,1);
    mpz_mod(p6_m1,p6_m1,prime);
    mpz_cdiv_q_ui(q_p_m1_d9,p6_m1,9);
    mpz_powm(c1_omega,c1_omega,q_p_m1_d9,prime);
    mpz_powm_ui(c1_omega_bar,c1_omega,2,prime);
    //    gmp_printf("C1_omg=%Zd,\n C1_omg^2=%Zd\n",c1_omega,c1_omega_bar);
    
    mpz_clear(p6_m1);
    mpz_clear(q_p_m1_d9);
}

#pragma mark Previous Methods
void Fp_pow_bin(struct Fp *RES,struct Fp *A,mpz_t B){//RES = A^B
    unsigned long long i,length;
    length= (unsigned long long)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    printf("Length: %lld\n",length);
    //    printf("p: %s\n",&B_binary[0]);
    
    struct Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,A);
    
    for(i=1;B_binary[i]!='\0';i++){
        Fp_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp_mul(&tmp,&tmp,A);
        }
    }
    Fp_set(RES,&tmp);
    Fp_clear(&tmp);
}

void Fp3_mul_prev(struct Fp3 *RES,struct Fp3 *A,struct Fp3 *B){
    //(a0,a1,a2)*(b0,b1,b2)=(x0y0+xi((a1+a2)(b1+b2)-a1b1-x2y2),xix2y2+(a0+a1)(b0+b1)-x0y0-a1b1,a1b1+(a0+a2)(b0+b2)-x0y0-x2y2)
    struct Fp A0,A1,A2,tmp_a01,tmp_a12,tmp_a20,tmp_b01,tmp_b12,tmp_b20,A3,A4,A5,tmp;
    struct Fp3 t_RES;
    Fp_init(&A0);
    Fp_init(&A1);
    Fp_init(&A2);
    Fp_init(&tmp_a01);
    Fp_init(&tmp_a12);
    Fp_init(&tmp_a20);
    Fp_init(&tmp_b01);
    Fp_init(&tmp_b12);
    Fp_init(&tmp_b20);
    Fp_init(&A3);
    Fp_init(&A4);
    Fp_init(&A5);
    Fp_init(&tmp);
    Fp3_init(&t_RES);
    
    Fp_mul(&A0,&A->a0,&B->a0);//A0=a0*b0
    Fp_mul(&A1,&A->a1,&B->a1);//A1=a1*b1
    Fp_mul(&A2,&A->a2,&B->a2);//A2=a2*b2
    
    Fp_add(&tmp_a01,&A->a0,&A->a1);//a0+a1
    Fp_add(&tmp_a12,&A->a1,&A->a2);//a1+a2
    Fp_add(&tmp_a20,&A->a0,&A->a2);//a2+a0
    Fp_add(&tmp_b01,&B->a0,&B->a1);//b0+b1
    Fp_add(&tmp_b12,&B->a1,&B->a2);//b1+b2
    Fp_add(&tmp_b20,&B->a0,&B->a2);//b2+b0
    
    Fp_mul(&A3,&tmp_a01,&tmp_b01);//A3=(a0+a1)(b0+b1)
    Fp_mul(&A4,&tmp_a12,&tmp_b12);//A4=(a1+a2)(b1+b2)
    Fp_mul(&A5,&tmp_a20,&tmp_b20);//A5=(a2+a0)(b2+b0)
    
    //(a0,a1,a2)*(b0,b1,b2)=(a0b0 + ai((a1+a2)(b1+b2)-a1b1-a2b2), aia2b2+(a0+a1)(b0+b1)-a0b0-a1b1, a1b1+(a0+a2)(b0+b2)-a0b0-a2b2)
    
    //t1 = A3 −A1 −A0
    //t2 = A4−A2−A0+A1
    //t3 = A5 −A1 −A2
    //ab = (t0 + c1t3) + (t1 + c1t4)ω + t2ω2.
    Fp_sub(&A4,&A4,&A1);//A4-A1
    Fp_sub(&A4,&A4,&A2);//A4-A2-A1        /(a1+a2)(b1+b2)-a1b1-a2b2
    Fp_mul_ui(&tmp,&A4,C1);
    Fp_add(&t_RES.a0,&A0,&tmp);
    
    Fp_sub(&A3,&A3,&A0);//A3-A0
    Fp_sub(&A3,&A3,&A1);//t1 = A3-A1-A0
    Fp_mul_ui(&tmp,&A2,C1);//A2*C1
    Fp_add(&t_RES.a1,&tmp,&A3);
    
    Fp_sub(&A5,&A5,&A0);//A5-A0
    Fp_sub(&A5,&A5,&A2);//A5-A2-A0
    Fp_add(&t_RES.a2,&A1,&A5);//A5-A2-A0 + A1;;t2 = A4−A2−A0+A1
    
    Fp3_set(RES,&t_RES);
    
    Fp_clear(&A0);
    Fp_clear(&A1);
    Fp_clear(&A2);
    Fp_clear(&tmp_a01);
    Fp_clear(&tmp_a12);
    Fp_clear(&tmp_a20);
    Fp_clear(&tmp_b01);
    Fp_clear(&tmp_b12);
    Fp_clear(&tmp_b20);
    Fp_clear(&A3);
    Fp_clear(&A4);
    Fp_clear(&A5);
    Fp_clear(&tmp);
}
void Fp3_invert_prev(struct Fp3 *ANS, struct Fp3 *A){
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
    Fp_mul(&T0,&A->a0,&A->a0);
    Fp_mul_ui(&t0,&A->a1,C1);
    
    Fp_mul(&T1,&t0,&A->a2);
    Fp_sub(&t1,&T0,&T1); // t1=(a^2-bci) mod q
    
    Fp_mul(&T0,&A->a2,&A->a2);
    Fp_mul_ui(&T0,&T0,C1);
    Fp_mul(&T1,&A->a0,&A->a1);
    Fp_sub(&t2,&T0,&T1); // t2=(c^2i-ab) mod q
    
    Fp_mul(&T0,&A->a1,&A->a1);
    Fp_mul(&T1,&A->a0,&A->a2);
    Fp_sub(&t3,&T0,&T1); // t3=(b^2-ac) mod q
    
    Fp_mul(&T0,&t0,&t3);
    Fp_mul(&T1,&A->a0,&t1);
    Fp_add(&T0,&T0,&T1); // T0={bi(b^2-ac)+a(a^2-bci)} mod q
    
    Fp_mul_ui(&t0,&A->a2,C1);
    Fp_mul(&T1,&t0,&t2);
    Fp_add(&t0,&T0,&T1); // t0={ci(c^2i-ab)+{bi(b^2-ac)+a(a^2-bci)}} mod q .0
    
    mpz_invert(t0.x_0,t0.x_0,prime);
    
    Fp_mul(&t_ans.a0,&t1,&t0);
    Fp_mul(&t_ans.a1,&t2,&t0);
    Fp_mul(&t_ans.a2,&t3,&t0);
    
    Fp3_set(ANS,&t_ans);
    
    Fp3_clear(&t_ans);
    Fp_clear(&T0);
    Fp_clear(&T1);
    Fp_clear(&t0);
    Fp_clear(&t1);
    Fp_clear(&t2);
    Fp_clear(&t3);
}
void Fp6_mul_prev(struct Fp6 *RES,struct Fp6 *A,struct Fp6 *B){
    //x^2-c2=0;c2=-1
    //= (a0b0 + c2a1b1) + (a0 + a1)(b0 + b1)τ−(a0b0)τ − (a1b1)τ.
    struct Fp3 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    Fp3_init(&tmp1);
    Fp3_init(&tmp2);
    Fp3_init(&tmp3);
    Fp3_init(&tmp4);
    Fp3_init(&tmp5);
    Fp3_init(&tmp6);
    
    struct Fp6 t_RES;
    Fp6_init(&t_RES);
    
    Fp3_mul(&tmp1,&A->a0,&B->a0);//a0*b0
    Fp3_mul(&tmp2,&A->a1,&B->a1);//a1*b1
    Fp3_mul_omega(&tmp3,&tmp2);//c2*a1*b1
    //    Fp3_mul_omega(&tmp3,&tmp2);
    Fp3_add(&t_RES.a0,&tmp1,&tmp3);//a*c+b*d*v
    Fp3_add(&tmp4,&A->a0,&A->a1);//a+b
    Fp3_add(&tmp5,&B->a0,&B->a1);//c+d
    Fp3_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
    Fp3_sub(&t_RES.a1,&tmp6,&tmp1);
    Fp3_sub(&t_RES.a1,&t_RES.a1,&tmp2);
    
    Fp6_set(RES,&t_RES);
    
    Fp3_clear(&tmp1);
    Fp3_clear(&tmp2);
    Fp3_clear(&tmp3);
    Fp3_clear(&tmp4);
    Fp3_clear(&tmp5);
    Fp3_clear(&tmp6);
    Fp6_clear(&t_RES);
}
void Fp6_pow_prev(struct Fp6 *RES,struct Fp6 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    for(i=1;B_binary[i]!='\0';i++){
        Fp6_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp6_mul(&tmp,&tmp,A);
        }
    }
    Fp6_set(RES,&tmp);
    Fp6_clear(&tmp);
}
void Fp6_invert_prev(struct Fp6 *ANS, struct Fp6 *A){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    
    // tmp=A^(q^6)=(x0,-x1)
    Fp3_set(&tmp.a0,&A->a0);
    Fp3_neg(&tmp.a1,&A->a1);
    
    struct Fp3 c,a,b;
    Fp3_init(&c);
    Fp3_init(&a);
    Fp3_init(&b);
    
    Fp3_mul(&a,&A->a0,&A->a0); // a=x0^2
    Fp3_mul(&b,&A->a1,&A->a1); // b=x1^2
    Fp3_mul_omega(&b,&b); // b=x1^2*v
    Fp3_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp3_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp6-element and tmp is a vector A Fp6
    Fp3_mul(&tmp.a0,&tmp.a0,&c);
    Fp3_mul(&tmp.a1,&tmp.a1,&c);
    
    Fp6_set(ANS,&tmp);
    
    Fp3_clear(&c);
    Fp3_clear(&a);
    Fp3_clear(&b);
    Fp6_clear(&tmp);
}
void Fp18_invert_prev(struct Fp18 *RES, struct Fp18 *A){
    struct Fp18 t_RES;
    Fp18_init(&t_RES);
    
    struct Fp6 T0,T1,t0,t1,t2,t3;
    Fp6_init(&T0);
    Fp6_init(&T1);
    Fp6_init(&t0);
    Fp6_init(&t1);
    Fp6_init(&t2);
    Fp6_init(&t3);
    
    // An optimized version of Grewal's Algo. 3   (a,b,c)
    Fp6_mul(&T0,&A->m0,&A->m0);
    Fp6_mul_tau(&t0,&A->m1);//TODO
    
    Fp6_mul(&T1,&t0,&A->m2);
    Fp6_sub(&t1,&T0,&T1); // t1=(a^2-bci) mod q
    
    Fp6_mul(&T0,&A->m2,&A->m2);
    Fp6_mul_tau(&T0,&T0);
    Fp6_mul(&T1,&A->m0,&A->m1);
    Fp6_sub(&t2,&T0,&T1); // t2=(c^2i-ab) mod q
    
    Fp6_mul(&T0,&A->m1,&A->m1);
    Fp6_mul(&T1,&A->m0,&A->m2);
    Fp6_sub(&t3,&T0,&T1); // t3=(b^2-ac) mod q
    
    Fp6_mul(&T0,&t0,&t3);
    Fp6_mul(&T1,&A->m0,&t1);
    Fp6_add(&T0,&T0,&T1); // T0={bi(b^2-ac)+a(a^2-bci)} mod q
    
    Fp6_mul_tau(&t0,&A->m2);
    Fp6_mul(&T1,&t0,&t2);
    Fp6_add(&t0,&T0,&T1); // t0={ci(c^2i-ab)+{bi(b^2-ac)+a(a^2-bci)}} mod q .0
    
    Fp6_invert(&t0,&t0);
    Fp6_printf(&t0);
    Fp6_mul(&t_RES.m0,&t1,&t0);
    Fp6_mul(&t_RES.m1,&t2,&t0);
    Fp6_mul(&t_RES.m2,&t3,&t0);
    
    Fp18_set(RES,&t_RES);
    
    Fp18_clear(&t_RES);
    Fp6_clear(&T0);
    Fp6_clear(&T1);
    Fp6_clear(&t0);
    Fp6_clear(&t1);
    Fp6_clear(&t2);
    Fp6_clear(&t3);
}