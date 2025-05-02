//
//  main.c
//  Fp18_Arith
//
//  Created by Khandaker Md. Al-Amin on 4/15/16.
//  Copyright © 2016 Khandaker Md. Al-Amin. All rights reserved.
//

#include "fp18_arith.h" // Include the consolidated header
#include "parameters.h" // Include parameter definitions/functions

// Global counters (consider moving to a dedicated stats module if complexity grows)
int fp_mul = 0;
int fp_add = 0;

int main(void)
{
    // Initialize parameters
    init_parameters();

    // Set the generator value X (example)
    // This could also be read from input or config
    mpz_set_str(X,"18446893747415302274",10); // Default X from original main.c

    // Generate curve parameters based on X
    generate_parameters();

    // Print generated parameters for verification
    gmp_printf("p = %Zd\n", prime);
    gmp_printf("r = %Zd\n", r_order);
    gmp_printf("t = %Zd\n", t_trace);
    gmp_printf("#E(Fp) = %Zd\n", r_order_EFp);
    gmp_printf("b = %Zd\n", b);
    printf("p = %d bits\n", (int)mpz_sizeinbase(prime, 2));
    printf("r = %d bits\n", (int)mpz_sizeinbase(r_order, 2));
    printf("t = %d bits\n", (int)mpz_sizeinbase(t_trace, 2));
    printf("\n");

    // --- Example Usage / Testing --- 
    // (Keep relevant examples, move extensive tests to test suite)

    printf("--- Fp3 Inversion Test ---\n");
    struct Fp3 P3, RES3;
    Fp3_init(&P3);
    Fp3_init(&RES3);

    // Example input (replace Fp3_take_input for non-interactive use)
    // Fp3_take_input(&P3);
    // Using a sample Fp3 element for demonstration:
    Fp_set_ui(&P3.a0, 5);
    Fp_set_ui(&P3.a1, 8);
    Fp_set_ui(&P3.a2, 2);
    printf("Input P3: "); Fp3_printf(&P3);

    fp_mul = 0; // Reset counters for this operation
    fp_add = 0;
    Fp3_invert(&RES3, &P3);
    printf("Inverted RES3: "); Fp3_printf(&RES3);
    printf("Fp3 Invert: FP Mul = %d, FP Add = %d\n", fp_mul, fp_add);

    // Verification: P3 * RES3 should be 1
    struct Fp3 VERIFY3;
    Fp3_init(&VERIFY3);
    Fp3_mul(&VERIFY3, &P3, &RES3);
    printf("Verification (P3 * RES3): "); Fp3_printf(&VERIFY3);
    printf("\n");

    Fp3_clear(&P3);
    Fp3_clear(&RES3);
    Fp3_clear(&VERIFY3);

    // --- Add other tests as needed (Fp6, Fp18, EC, Pairing) ---

    printf("\n--- Fp Addition Test ---\n");
    struct Fp A_fp, B_fp, R_fp, EXPECTED_fp;
    Fp_init(&A_fp);
    Fp_init(&B_fp);
    Fp_init(&R_fp);
    Fp_init(&EXPECTED_fp);

    // Set values for A and B
    Fp_set_ui(&A_fp, 12345);
    Fp_set_ui(&B_fp, 67890);
    printf("Input A_fp: "); Fp_printf(&A_fp);
    printf("Input B_fp: "); Fp_printf(&B_fp);

    // Calculate expected result (manually or using mpz for large numbers if needed)
    // Expected = (12345 + 67890) mod p
    mpz_t temp_sum;
    mpz_init(temp_sum);
    mpz_add_ui(temp_sum, A_fp.x0, 67890);
    mpz_mod(temp_sum, temp_sum, prime);
    Fp_set_mpz(&EXPECTED_fp, temp_sum);
    printf("Expected R_fp: "); Fp_printf(&EXPECTED_fp);

    // Perform Fp_add
    Fp_add(&R_fp, &A_fp, &B_fp);
    printf("Result R_fp = A_fp + B_fp: "); Fp_printf(&R_fp);

    // Verification
    if (Fp_cmp(&R_fp, &EXPECTED_fp) == 0) {
        printf("Fp_add Verification: PASSED\n");
    } else {
        printf("Fp_add Verification: FAILED\n");
    }
    printf("\n");

    Fp_clear(&A_fp);
    Fp_clear(&B_fp);
    Fp_clear(&R_fp);
    Fp_clear(&EXPECTED_fp);
    mpz_clear(temp_sum);

    // Example: EFp Addition Test
    printf("--- EFp Addition Test ---\n");
    struct EFp P_fp, Q_fp, R_fp;
    EFp_init(&P_fp);
    EFp_init(&Q_fp);
    EFp_init(&R_fp);

    // Find two random points (ensure they are valid and on the curve)
    // EFp_random_set might be complex, using placeholder values for now
    // TODO: Implement robust point generation or use known points
    // Placeholder: Set P_fp and Q_fp manually (assuming they are valid)
    // These values are arbitrary and likely NOT on the curve y^2 = x^3 + b
    // Replace with actual points from EFp_random_set or known test vectors.
    Fp_set_ui(&P_fp.px, 10);
    Fp_set_ui(&P_fp.py, 20); // Placeholder
    P_fp.isInfinity = FALSE;
    Fp_set_ui(&Q_fp.px, 15);
    Fp_set_ui(&Q_fp.py, 25); // Placeholder
    Q_fp.isInfinity = FALSE;

    printf("Input P_fp: "); EFp_printf(&P_fp);
    printf("Input Q_fp: "); EFp_printf(&Q_fp);

    EFp_ECA(&R_fp, &P_fp, &Q_fp);
    printf("Result R_fp = P_fp + Q_fp: "); EFp_printf(&R_fp);
    printf("\n");

    EFp_clear(&P_fp);
    EFp_clear(&Q_fp);
    EFp_clear(&R_fp);

    // Cleanup parameters
    clear_parameters();

    return 0;
}

// Remove all function definitions (Fp_*, Fp3_*, etc.) as they are now
// expected to be in separate implementation files (Finite_Field.c, Elliptic_Curve.c, etc.)
// and their prototypes are in fp18_arith.h.

// Remove Input functions (Fp_take_input, etc.) or move them to a dedicated
// test utility file if interactive input is desired for testing.