#include "fp.h"
#include "fp3.h"
#include "fp6.h"
#include "fp18.h"
#include "ec.h"
#include "pairing.h"
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
    
    Fp a, b, result;
    fp_init(&a);
    fp_init(&b);
    fp_init(&result);
    
    // Generate random values
    fp_random(&a);
    fp_random(&b);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp_add(&result, &a, &b);
    }
    long end_time = get_microseconds();
    printf("fp_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp_mul(&result, &a, &b);
    }
    end_time = get_microseconds();
    printf("fp_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        fp_invert(&result, &a);
    }
    end_time = get_microseconds();
    printf("fp_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    fp_clear(&a);
    fp_clear(&b);
    fp_clear(&result);
}

// Benchmark Fp3 arithmetic
void benchmark_fp3_arithmetic() {
    printf("\nBenchmarking Fp3 arithmetic operations...\n");
    
    Fp3 a3, b3, result3;
    fp3_init(&a3);
    fp3_init(&b3);
    fp3_init(&result3);
    
    // Generate random values
    fp3_random(&a3);
    fp3_random(&b3);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp3_add(&result3, &a3, &b3);
    }
    long end_time = get_microseconds();
    printf("fp3_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp3_mul(&result3, &a3, &b3);
    }
    end_time = get_microseconds();
    printf("fp3_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        fp3_invert(&result3, &a3);
    }
    end_time = get_microseconds();
    printf("fp3_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    fp3_clear(&a3);
    fp3_clear(&b3);
    fp3_clear(&result3);
}

// Benchmark Fp6 arithmetic
void benchmark_fp6_arithmetic() {
    printf("\nBenchmarking Fp6 arithmetic operations...\n");
    
    Fp6 a6, b6, result6;
    fp6_init(&a6);
    fp6_init(&b6);
    fp6_init(&result6);
    
    // Generate random values
    fp6_random(&a6);
    fp6_random(&b6);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp6_add(&result6, &a6, &b6);
    }
    long end_time = get_microseconds();
    printf("fp6_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp6_mul(&result6, &a6, &b6);
    }
    end_time = get_microseconds();
    printf("fp6_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        fp6_invert(&result6, &a6);
    }
    end_time = get_microseconds();
    printf("fp6_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    fp6_clear(&a6);
    fp6_clear(&b6);
    fp6_clear(&result6);
}

// Benchmark Fp18 arithmetic
void benchmark_fp18_arithmetic() {
    printf("\nBenchmarking Fp18 arithmetic operations...\n");
    
    Fp18 a18, b18, result18;
    fp18_init(&a18);
    fp18_init(&b18);
    fp18_init(&result18);
    
    // Generate random values
    fp18_random(&a18);
    fp18_random(&b18);
    
    // Benchmark addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp18_add(&result18, &a18, &b18);
    }
    long end_time = get_microseconds();
    printf("fp18_add: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark multiplication
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        fp18_mul(&result18, &a18, &b18);
    }
    end_time = get_microseconds();
    printf("fp18_mul: %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark inversion (less iterations due to cost)
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) {
        fp18_invert(&result18, &a18);
    }
    end_time = get_microseconds();
    printf("fp18_invert: %.2f μs per operation\n", (double)(end_time - start_time) / (ITERATIONS / 100));
    
    fp18_clear(&a18);
    fp18_clear(&b18);
    fp18_clear(&result18);
}

// Benchmark elliptic curve operations
void benchmark_elliptic_operations() {
    printf("\nBenchmarking elliptic curve operations...\n");
    
    ECPoint P, Q, R;
    ec_init(&P);
    ec_init(&Q);
    ec_init(&R);
    
    // Generate random points
    ec_random_set(&P);
    ec_random_set(&Q);
    
    // Benchmark point addition
    long start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        ec_add(&R, &P, &Q);
    }
    long end_time = get_microseconds();
    printf("ec_add (point addition): %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark point doubling
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS; i++) {
        ec_double(&R, &P);
    }
    end_time = get_microseconds();
    printf("ec_double (point doubling): %.2f μs per operation\n", (double)(end_time - start_time) / ITERATIONS);
    
    // Benchmark scalar multiplication
    mpz_t scalar;
    mpz_init(scalar);
    mpz_set_ui(scalar, 0xABCDEF); // Large enough scalar
    
    start_time = get_microseconds();
    for (int i = 0; i < ITERATIONS / 100; i++) { // Less iterations due to cost
        ec_scalar_mul(&R, &P, scalar);
    }
    end_time = get_microseconds();
    printf("ec_scalar_mul (scalar multiplication): %.2f μs per operation\n", 
           (double)(end_time - start_time) / (ITERATIONS / 100));
    
    ec_clear(&P);
    ec_clear(&Q);
    ec_clear(&R);
    mpz_clear(scalar);
}

// Benchmark pairing operations
void benchmark_pairing() {
    printf("\nBenchmarking pairing operations...\n");
    
    ECPoint P;
    ECPoint18 Q;
    Fp18 result;
    
    ec_init(&P);
    ec18_init(&Q);
    fp18_init(&result);
    
    // Generate random points
    ec_random_set(&P);
    ec18_random_set_G2(&Q);
    
    // Benchmark optimal ate pairing
    long start_time = get_microseconds();
    // Due to complexity, run fewer iterations
    for (int i = 0; i < 10; i++) {
        optimal_ate_pairing(&result, &P, &Q);
    }
    long end_time = get_microseconds();
    printf("optimal_ate_pairing: %.2f ms per operation\n", 
           (double)(end_time - start_time) / (10 * 1000)); // Convert to milliseconds
    
    ec_clear(&P);
    ec18_clear(&Q);
    fp18_clear(&result);
}

int main() {
    printf("Running benchmarks for Finite Extension Field Degree 18 Library\n");
    printf("Each operation is repeated multiple times to get reliable averages\n");
    
    // Initialize parameters
    init_kss18_params();
    
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
    clear_kss18_params();
    
    return 0;
}