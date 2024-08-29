### Parallel Programming Project 2024

This project focuses on parallelizing a sequential multi-start multidimensional scaling (MDS) algorithm (`multistart_mds_seq.c`) using different parallel programming paradigms: OpenMP, OpenMP tasks, and MPI. The aim is to explore and compare the performance and scalability of these parallelization techniques when applied to a computationally intensive algorithm.

#### Files and Parallelization Techniques

1. **`multistart_mds_omp.c`**:
   - **Parallelization Method**: OpenMP

2. **`multistart_mds_omp_tasks.c`**:
   - **Parallelization Method**: OpenMP Tasks

3. **`multistart_mds_mpi.c`**:
   - **Parallelization Method**: MPI (Message Passing Interface)

#### Key Objectives

- **Performance Analysis**: Measure and compare the performance of each parallelization technique, focusing on speedup, efficiency, and scalability.
- **Scalability Testing**: Evaluate how each method scales with increasing numbers of processors/threads.
- **Optimization**: Investigate potential optimizations in each parallelization method to maximize performance.

#### How to Run

Run make on the directory.

#### Dependencies

- **OpenMP**: Requires a compiler with OpenMP support (e.g., GCC, Clang).
- **MPI**: Requires an MPI implementation (e.g., MPICH, OpenMPI).
