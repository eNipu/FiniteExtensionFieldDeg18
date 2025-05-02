#include "../fp18_arith.h"
#include "../parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// Number of iterations for benchmarks
#define ITERATIONS 1000

// Function to get current time in microseconds
long get_microseconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000000 + tv.tv_usec;
}

// Benchmark Fp arithmetic
void benchmark_fp_arithmetic() {
    printf("\nBenchmarking Fp arithmetic operations...\n");
    
    struct Fp a, b, result;
    Fp_init(&a);
    Fp_init(&b);
    Fp_init(&result);
    
    // Generate random values
    Fp_random(&a);
    Fp_random(&b);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp_add(&result, &a, &b);
    }
    long end_time = get_microseconds();
    printf("Fp_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp_mul(&result, &a, &b);
    }
    end_time = get_microseconds();
    printf("Fp_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        Fp_invert(&result, &a);
    }
    end_time = get_microseconds();
    printf("Fp_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    Fp_clear(&a);
    Fp_clear(&b);
    Fp_clear(&result);
}

// Benchmark Fp3 arithmetic
void benchmark_fp3_arithmetic() {
    printf("\nBenchmarking Fp3 arithmetic operations...\n");
    
    struct Fp3 a3, b3, result3;
    Fp3_init(&a3);
    Fp3_init(&b3);
    Fp3_init(&result3);
    
    // Generate random values
    Fp3_random(&a3);
    Fp3_random(&b3);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp3_add(&result3, &a3, &b3);
    }
    long end_time = get_microseconds();
    printf("Fp3_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp3_mul(&result3, &a3, &b3);
    }
    end_time = get_microseconds();
    printf("Fp3_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        Fp3_invert(&result3, &a3);
    }
    end_time = get_microseconds();
    printf("Fp3_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    Fp3_clear(&a3);
    Fp3_clear(&b3);
    Fp3_clear(&result3);
}

// Benchmark Fp6 arithmetic
void benchmark_fp6_arithmetic() {
    printf("\nBenchmarking Fp6 arithmetic operations...\n");
    
    struct Fp6 a6, b6, result6;
    Fp6_init(&a6);
    Fp6_init(&b6);
    Fp6_init(&result6);
    
    // Generate random values
    Fp6_random(&a6);
    Fp6_random(&b6);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp6_add(&result6, &a6, &b6);
    }
    long end_time = get_microseconds();
    printf("Fp6_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp6_mul(&result6, &a6, &b6);
    }
    end_time = get_microseconds();
    printf("Fp6_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        Fp6_invert(&result6, &a6);
    }
    end_time = get_microseconds();
    printf("Fp6_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    Fp6_clear(&a6);
    Fp6_clear(&b6);
    Fp6_clear(&result6);
}

// Benchmark Fp18 arithmetic
void benchmark_fp18_arithmetic() {
    printf("\nBenchmarking Fp18 arithmetic operations...\n");
    
    struct Fp18 a18, b18, result18;
    Fp18_init(&a18);
    Fp18_init(&b18);
    Fp18_init(&result18);
    
    // Generate random values
    Fp18_random(&a18);
    Fp18_random(&b18);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp18_add(&result18, &a18, &b18);
    }
    long end_time = get_microseconds();
    printf("Fp18_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        Fp18_mul(&result18, &a18, &b18);
    }
    end_time = get_microseconds();
    printf("Fp18_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        Fp18_invert(&result18, &a18);
    }
    end_time = get_microseconds();
    printf("Fp18_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    Fp18_clear(&a18);
    Fp18_clear(&b18);
    Fp18_clear(&result18);
}

// Benchmark elliptic curve operations
void benchmark_elliptic_operations() {
    printf("\nBenchmarking elliptic curve operations...\n");
    
    struct EFp P, Q, R;
    EFp_init(&P);
    EFp_init(&Q);
    EFp_init(&R);
    
    // Generate random points
    EFp_random_set(&P);
    EFp_random_set(&Q);
    
    // Benchmark point addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        EFp_ECA(&R, &P, &Q);
    }
    long end_time = get_microseconds();
    printf("EFp_ECA (point addition): %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark point doubling
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        EFp_ECD(&R, &P);
    }
    end_time = get_microseconds();
    printf("EFp_ECD (point doubling): %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark scalar multiplication
    mpz_t scalar;
    mpz_init(scalar);
    mpz_set_ui(scalar, 0xABCDEF); // Large enough scalar
    
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) { // Less iterations due to cost
        EFp_SCM(&R, &P, scalar);
    }
    end_time = get_microseconds();
    printf("EFp_SCM (scalar multiplication): %.2f μs per operation\n", 
           (double)(end_time - start_time) / (ITERATIONS / 100));
    
    EFp_clear(&P);
    EFp_clear(&Q);
    EFp_clear(&R);
    mpz_clear(scalar);
}

// Benchmark pairing operations
void benchmark_pairing() {
    printf("\nBenchmarking pairing operations...\n");
    
    struct EFp P;
    struct EFp18 Q;
    struct Fp18 result;
    
    EFp_init(&P);
    EFp18_init(&Q);
    Fp18_init(&result);
    
    // Generate random points
    EFp_random_set(&P);
    EFp18_random_set_G2(&Q);
    
    // Benchmark optimal ate pairing
    long start_time = get_microseconds();
    // Due to complexity, run fewer iterations
    for (int i = 0; i < 10; i++) {
        Optimal_Ate_Pairing(&result, &P, &Q);
    }
    long end_time = get_microseconds();
    printf("Optimal_Ate_Pairing: %.2f ms per operation\n", 
           (double)(end_time - start_time) / (10 * 1000)); // Convert to milliseconds
    
    EFp_clear(&P);
    EFp18_clear(&Q);
    Fp18_clear(&result);
}

int main() {
    printf("Running benchmarks for Finite Extension Field Degree 18 Library\n");
    printf("Each operation is repeated multiple times to get reliable averages\n");
    
    // Initialize parameters
    init_parameters();
    
    // Set the generator value X
    mpz_set_str(X, "18446893747415302274", 10);
    
    // Generate curve parameters based on X
    generate_parameters();
    
    // Run all benchmark suites
    benchmark_fp_arithmetic();
    benchmark_fp3_arithmetic();
    benchmark_fp6_arithmetic();
    benchmark_fp18_arithmetic();
    benchmark_elliptic_operations();
    benchmark_pairing();
    
    // Report completion
    printf("\nBenchmark complete\n");
    
    // Cleanup
    clear_parameters();
    
    return 0;
}