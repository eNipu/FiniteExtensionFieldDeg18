//
//  main.c
//  Fp18_Arith
//
//  Created by Khandaker Md. Al-Amin on 4/15/16.
//  Copyright © 2016 Khandaker Md. Al-Amin. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "embedding_degree18.h" // Include the main header file
#include "parameters.h" // Include parameter definitions/functions

// Global counters (consider moving to a dedicated stats module if complexity grows)
int fp_mul = 0;
int fp_add = 0;

int main(void)
{
    struct Fp A_fp, B_fp, C_fp, EXPECTED_fp;
    mpz_t temp_sum;

    printf("Finite Extension Field Degree 18 Library Test\n");
    printf("============================================\n\n");
    
    // Initialize parameters
    EFp_set_EC_parameter();
    
    // Initialize GMP variables
    mpz_init(temp_sum);
    
    // Initialize field elements
    Fp_init(&A_fp);
    Fp_init(&B_fp);
    Fp_init(&C_fp);
    Fp_init(&EXPECTED_fp);
    
    // Set test values
    Fp_set_ui(&A_fp, 123);
    Fp_set_ui(&B_fp, 456);
    
    // Test addition
    Fp_add(&C_fp, &A_fp, &B_fp);
    
    // Calculate expected result
    mpz_set_ui(temp_sum, 123 + 456);
    mpz_mod(temp_sum, temp_sum, prime);
    Fp_set_mpz(&EXPECTED_fp, temp_sum);
    
    // Display results
    printf("Testing Fp addition: (123 + 456) mod p\n");
    printf("Result:   "); Fp_printf(&C_fp);
    printf("Expected: "); Fp_printf(&EXPECTED_fp);
    
    // Check if the result is correct
    if (Fp_cmp(&C_fp, &EXPECTED_fp) == 0) {
        printf("Addition test: SUCCESS!\n");
    } else {
        printf("Addition test: FAILED!\n");
    }
    
    // Clean up
    Fp_clear(&A_fp);
    Fp_clear(&B_fp);
    Fp_clear(&C_fp);
    Fp_clear(&EXPECTED_fp);
    mpz_clear(temp_sum);
    
    return 0;
}