//
//  main.c
//  Fp18_Arith
//
//  Created by Khandaker Md. Al-Amin on 4/15/16.
//  Copyright © 2016 Khandaker Md. Al-Amin. All rights reserved.
//

#include "embedding_degree18.h"
#include "fp.h"
#include "fp3.h"
#include "fp6.h"
#include "fp18.h"
#include "ec.h"
#include "pairing.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

int main(void) {
    // Initialize random seed
    srand(time(NULL));

    // Initialize KSS18 parameters
    init_kss18_params();

    // Example: create and print a random Fp element
    Fp a;
    fp_init(&a);
    fp_random(&a);
    printf("Random Fp element: ");
    fp_print(&a);
    printf("\n");
    fp_clear(&a);

    // Example: pairing stub usage
    EcFp P; EcFp3 Q; Fp18 gt;
    ecfp_init(&P); ecfp3_init(&Q); fp18_init(&gt);
    optimal_ate_pairing(&gt, &P, &Q);
    printf("Pairing output (stub): ");
    fp18_print(&gt);
    printf("\n");
    fp18_clear(&gt); ecfp_clear(&P); ecfp3_clear(&Q);

    // Clear KSS18 parameters
    clear_kss18_params();
    return 0;
}