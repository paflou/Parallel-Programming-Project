CC=gcc

all: statistics

# TODO 1: add flag to build OpenMP code
statistics: statistics.c
		$(CC) -O3 -Wall -fopenmp statistics.c -o statistics

clean:
		rm -f statistics

.PHONY: all clean
