//gcc -I/usr/local/include -L/usr/local/lib -Wall -O3 -o degree18.out degree18.c -lgmp
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
int main(int argc, char *argv[]){

	if (argc != 2) {
		printf("Usage: %s num\n", argv[0]);
		exit(1);
	}
	mpz_t X,X2;
	mpz_init(X);
	mpz_init(X2);

	// mpz_set_str(X,argv[1], 10);
	// mpz_set_str(X,"14",10);
	int X_int;
	X_int=atoi(argv[1]);
	X_int=X_int/8;
	printf("%d\n",X_int);
	mpz_setbit(X,X_int);
	mpz_div_ui(X2,X,4);
	mpz_add(X,X2,X);


	mpz_t p_tmp,r_tmp;
	mpz_t xpow2,xpow3,xpow4,xpow5,xpow6,xpow7,xpow8;
	mpz_t tmp1,tmp2;

	mpz_init(p_tmp);
	mpz_init(r_tmp);

	mpz_init(xpow2);
	mpz_init(xpow3);
	mpz_init(xpow4);
	mpz_init(xpow5);
	mpz_init(xpow6);
	mpz_init(xpow7);
	mpz_init(xpow8);

	mpz_init(tmp1);
	mpz_init(tmp2);

	mpz_tdiv_r_ui(tmp1,X,42);
	while(mpz_cmp_ui(tmp1,14)!=0){
		mpz_add_ui(X,X,1);
		mpz_tdiv_r_ui(tmp1,X,42);
	}

	mpz_t tmp,mod;
	mpz_init(tmp);
	mpz_init(mod);
	int i,j,k,l;
	mpz_set_ui(mod,3);

	for(;;){
		mpz_mul(xpow2,X,X);
		mpz_mul(xpow3,xpow2,X);
		mpz_mul(xpow4,xpow2,xpow2);
		mpz_mul(xpow5,xpow4,X);
		mpz_mul(xpow6,xpow3,xpow3);
		mpz_mul(xpow7,xpow6,X);
		mpz_mul(xpow8,xpow4,xpow4);

		mpz_mul_ui(tmp1,xpow3,37);
		mpz_add_ui(tmp2,xpow6,343);
		mpz_add(r_tmp,tmp1,tmp2);


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

		mpz_tdiv_r_ui(tmp1,p_tmp,21);


		if(mpz_cmp_ui(tmp1,0)==0){
			mpz_div_ui(p_tmp,p_tmp,21);
			// gmp_printf("%Zd\n",p_tmp);

			mpz_pow_ui(tmp,p_tmp,18);
			mpz_sub_ui(tmp,tmp,1);
			mpz_tdiv_r(tmp2,tmp,r_tmp);

			if(mpz_cmp_ui(tmp2,0)==0){
				// printf("afh\n");
				mpz_sub_ui(tmp,p_tmp,1);
				mpz_mod(tmp,tmp,mod);
				k=mpz_get_ui(tmp);

				j=mpz_probab_prime_p(p_tmp,25);
				mpz_sub_ui(tmp,p_tmp,1);
				mpz_tdiv_r_ui(tmp,tmp,4);
				l=mpz_get_ui(tmp);
				// i=mpz_probab_prime_p(r_tmp,25);
				if(k==0&&j!=0&&l!=0){
					gmp_printf("%Zd\n",X);
					gmp_printf("p:%Zd\n",p_tmp);
					printf("p = %dbit\n",(int)mpz_sizeinbase(p_tmp,2));
				}
			}
		}
		mpz_add_ui(X,X,42);
	}
	return 0;
}
