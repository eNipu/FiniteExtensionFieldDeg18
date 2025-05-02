CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Werror -I. -I./tests -I./bench
LDFLAGS = -lgmp

SRC_DIR = .
BUILD_DIR = ./build
TEST_DIR = ./tests
BENCHMARK_DIR = ./bench

# Modular source files
SRCS = fp.c fp3.c fp6.c fp18.c ec.c pairing.c parameters.c main.c
OBJS = $(SRCS:.c=.o)

TEST_SRCS = tests/test_fp.c tests/test_main.c
TEST_OBJS = $(TEST_SRCS:.c=.o)

BENCH_SRCS = bench/benchmark.c
BENCH_OBJS = $(BENCH_SRCS:.c=.o)

.PHONY: all clean test bench dirs

all: finitefield

dirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(TEST_DIR)
	@mkdir -p $(BENCHMARK_DIR)

finitefield: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: finitefield $(TEST_OBJS)
	$(CC) $(CFLAGS) -o test_fp tests/test_fp.c fp.o fp3.o fp6.o fp18.o ec.o pairing.o parameters.o $(LDFLAGS)
	$(CC) $(CFLAGS) -o test_main tests/test_main.c fp.o fp3.o fp6.o fp18.o ec.o pairing.o parameters.o $(LDFLAGS)

bench: finitefield $(BENCH_OBJS)
	$(CC) $(CFLAGS) -o bench/benchmark bench/benchmark.c fp.o fp3.o fp6.o fp18.o ec.o pairing.o parameters.o $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o finitefield test_fp test_main bench/benchmark
