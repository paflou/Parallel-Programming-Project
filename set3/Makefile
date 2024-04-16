CC=gcc
CFLAGS=-O3 -fopenmp -Wall -lm

.PHONY: all
all: mc_serial mc_parallel

mc_serial: mc_serial.c
	$(CC) $(CFLAGS) -o mc_serial mc_serial.c

mc_parallel: mc_parallel.c
	$(CC) $(CFLAGS) -o mc_parallel mc_parallel.c

.PHONY: clean
clean:
	rm -f mc_serial mc_parallel
