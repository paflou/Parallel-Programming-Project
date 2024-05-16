CC=gcc
CFLAGS=-Wall -O3 -fopenmp -lmpi
CFLAGS+=
LDLIBS=-lm

#TODO: add the following cases: multistart_mds_omp multistart_mds_omp_tasks multistart_mds_mpi

all: multistart_mds_seq multistart_mds_omp multistart_mds_omp_tasks multistart_mds_omp_tasks_old_mds multistart_mds_mpi

multistart_mds_seq: multistart_mds_seq.c torczon.c Makefile
	$(CC) $(CFLAGS) -o multistart_mds_seq multistart_mds_seq.c torczon.c $(LDLIBS)

multistart_mds_omp: multistart_mds_omp.c torczon.c Makefile
	$(CC) $(CFLAGS) -o multistart_mds_omp multistart_mds_omp.c torczon.c $(LDLIBS)

multistart_mds_omp_tasks: multistart_mds_omp_tasks.c torczon_tasks.c Makefile
	$(CC) $(CFLAGS) -o multistart_mds_omp_tasks multistart_mds_omp_tasks.c torczon_tasks.c $(LDLIBS)

multistart_mds_omp_tasks_old_mds: multistart_mds_omp_tasks.c torczon.c Makefile
	$(CC) $(CFLAGS) -o multistart_mds_omp_tasks_old_mds multistart_mds_omp_tasks.c torczon.c $(LDLIBS)

multistart_mds_mpi: multistart_mds_mpi.c torczon.c Makefile
	$(CC) $(CFLAGS) -o multistart_mds_mpi multistart_mds_mpi.c torczon.c $(LDLIBS)

clean:
	rm -f multistart_mds_seq multistart_mds_omp multistart_mds_omp_tasks multistart_mds_omp_tasks_old_mds multistart_mds_mpi
