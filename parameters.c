#define DEFINE_GLOBAL_VARIABLES
#include "parameters.h"
#include "embedding_degree18.h"

// Define global variables for the project
unsigned long int c1 = 2;
char X_bit_binary[x_bit+1];

// Global counters for optimizations
double opt_sum = 0, sparse_sum = 0, pseudo_sum = 0;
double opt_dbl = 0, opt_add = 0, opt_mul = 0;
double sps_dbl = 0, sps_add = 0, sps_mul = 0;
double pse_dbl = 0, pse_add = 0, pse_mul = 0;

void EFp_set_EC_parameter(void) {
    mpz_init(X);
    mpz_init(prime);
    mpz_init(order);
    mpz_init(trace);
    mpz_init(order_EFp);
    mpz_init(b);

    // Parameter initialization
    mpz_set_str(prime, "1209574531123134088399535937915047554087933139044681453126413342342151962271", 10);
    mpz_set_str(order, "1209574531123134088399535922471231530268062181839235806606083532738352871035", 10);
    mpz_set_str(order_EFp, "1209574531123134088399535922471231530268062181839235806606083532738352871035", 10);
    mpz_set_str(b, "4", 10);

    // Generate X parameter
    generate_X();
}

void generate_X(void) {
    int i;
    mpz_set_str(X, "57896044618658097711785492504343953926634992332820282019728792003956564819968", 10);

    // Convert X to binary representation
    for(i = 0; i < x_bit; i++) {
        if(mpz_tstbit(X, i)) {
            X_bit_binary[x_bit-1-i] = '1';
        } else {
            X_bit_binary[x_bit-1-i] = '0';
        }
    }
    X_bit_binary[x_bit] = '\0';
}
