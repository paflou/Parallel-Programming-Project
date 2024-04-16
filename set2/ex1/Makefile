CC=gcc
CFLAGS=-Wall -O3 -pthread -lm

all: integral_seq integral_threads

clean:
	rm -f integral_seq 
	rm -f integral_threads

%: %.c
	$(CC) $(CFLAGS) -o $@ $< timer.c
