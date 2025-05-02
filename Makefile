CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11
LDFLAGS = -lgmp -lm

SRC_DIR = .
BUILD_DIR = ./build
TEST_DIR = ./tests
BENCHMARK_DIR = ./bench

SRCS = $(wildcard $(SRC_DIR)/*.c)
SRCS := $(filter-out $(SRC_DIR)/main.c, $(SRCS))
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_OBJS = $(TEST_SRCS:$(TEST_DIR)/%.c=$(BUILD_DIR)/%.o)

BENCH_SRCS = $(wildcard $(BENCHMARK_DIR)/*.c)
BENCH_OBJS = $(BENCH_SRCS:$(BENCHMARK_DIR)/%.c=$(BUILD_DIR)/%.o)

MAIN_OBJ = $(BUILD_DIR)/main.o
EXEC = fp18_arith
TEST_EXEC = run_tests
BENCH_EXEC = run_benchmarks

.PHONY: all clean test bench docs dirs

all: dirs $(EXEC)

dirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(TEST_DIR)
	@mkdir -p $(BENCHMARK_DIR)

$(EXEC): $(MAIN_OBJ) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(TEST_EXEC): $(TEST_OBJS) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BENCH_EXEC): $(BENCH_OBJS) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(BENCHMARK_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

test: dirs $(TEST_EXEC)
	./$(TEST_EXEC)

bench: dirs $(BENCH_EXEC)
	./$(BENCH_EXEC)

docs:
	doxygen Doxyfile

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(EXEC) $(TEST_EXEC) $(BENCH_EXEC)
