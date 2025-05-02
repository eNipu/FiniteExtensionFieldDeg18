//
//  main.c
//  Fp18_Arith
//
//  Created by Khandaker Md. Al-Amin on 4/15/16.
//  Copyright © 2016 Khandaker Md. Al-Amin. All rights reserved.
//

#include "embedding_degree18.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

int main(void) {
    // Initialize random seed
    srand(time(NULL));
    
    // Initialize parameters for KSS degree 18 curve
    init_parameters();
    
    // Generate parameters for the curve
    generate_parameters();
    
    // Run pairing checks
    check_Pairing();
    
    // Measure pairing performance
    Masure_pairing_time();
    
    return 0;
}