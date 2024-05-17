//DO NOT COMPILE - JUST NOTETAKING

----------  PROCESSES VS POSIX Threads  ------------

ID:       pid_t, getpid() -- pthread_t, pthread_self()

Creation:      fork()     -- pthread_create

start:         exec()

wait:          wait()     -- pthread_join()

exit:          exit()     -- pthread_exit()


---------------- Mutex Management --------------------

Declaration and static initialization of a mutex:
pthread_mutex_t mymutex = PTHREAD_MUTEX_INITIALIZER;

Declaration and dynamic initialization of a mutex:
pthread_mutex_t mymutex;
pthread_mutex_init(&mymutex, NULL);

Dynamic Allocation of Mutex:
pthread_mutex_t *mymutex;
mymutex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));

Locking (acquiring) the mutex:
pthread_mutex_lock(&mymutex);

Unlocking (releasing) the mutex after the critical section:
pthread_mutex_unlock(&mymutex);

Destroying mutex:
pthread_mutex_destroy(mymutex);


--------------- Check and Lock ------------------------

int pthread_mutex_trylock(pthread_mutex_t *m);
• Allows a thread to try to lock a mutex
• If the mutex is available then the thread locks the mutex
• If the mutex is locked then the function informs the user by
returning a special value (EBUSY)
• This approach allows for implementations of spinlocks


------------------- Barriers --------------------------

• Barrier: synchronization mechanism
• No thread can cross the barrier until all the threads have reached it

int pthread_barrier_init(pthread_barrier_t * bar,
                         const pthread_barrierattr_t*attr,
                         usigned int count); // Count: No. of threads that must call barrier
int pthread_barrier_wait(pthread_barrier_t *bar);
int pthread_barrier_destroy(pthread_barrier_t *bar);


//////////////////////////////////////////////////////////
----------------------- OpenMP ---------------------------
//////////////////////////////////////////////////////////

----------------------- Parallel -------------------------

double A[20];
#pragma omp parallel [clause ...]
{
    //code
}

•Each threads runs the same code
•All threads share A
•Execution continues when all threads have finished their work

Clauses for omp parallel:

if (scalar_expression)              Only parallelize if the expression is true. Can be
                                    used to stop parallelization if the work is too little

num_threads (integer-expression)    Set the number of threads

private (list)                      The specified variables are thread-private

shared (list)                       The specified variables are shared among all threads

firstprivate (list)                 The specified variables are thread-private and
                                    initialized from the master thread

reduction (operator: list)          Perform a reduction on the thread-local variables
                                    and assign it to the master thread

default (shared | none)             Unspecified variables are shared or not


----------------- Static / Dynamic Mode ------------------

• Dynamic mode (default):
    - The number of threads can differ between parallel regions of the same program
    - The specified number of threads actually defines the maximum number,
                                    the actual number of threads can be smaller
• Static mode:
    - The number of threads is fixed and exactly equal to the number specified by the programmer
• OpenMP supports nested parallel regions but…
    - The compiler is allowed to serialize all the inner levels
    - This means that it uses a single OpenMP thread for those parallel regions


---------------------- Sections construct ----------------

• The sections construct gives a different structured block to
each thread


------------------- Critical Section ---------------------

Prevents multiple threads from accessing a section of code 
at the same time.

#pragma omp critical [(name)]
{
    code_block
}


----------------------- Barriers -------------------------
-------------- Overridable with nowait -------------------
Explicit:
#pragma omp barrier

Implicit:
- Parallel
- Loop
- Single


------------------ Amdahl's Law ----------------------------

How much a computation can be sped up by running part of a program in parallel

T = Total time of serial execution

B = Total time of non-parallelizable part

T - B = Total time of parallelizable part (when executed serially)

N = The number of threads

- T(N) = B + (T - B) / N


- Speedup = Original Execution Time / Execution time after Enhancement


Example:
T = 1
Non-parallelizable part = 40% of T
Parallelization factor = 2

T(2) = 0.4 + (1 - 0.4) / 2
     = 0.4 + 0.3 = 0.7

Speedup = 1 / 0.7



SERIAL CODE FINAL RESULTS:
Elapsed time = 42.432 s
Total number of trials = 64
Total number of function evaluations = 640308
Best result at trial 40 used 1250 iterations, 10005 function calls and returned
x[  0] =   1.0074912e+00 
x[  1] =   1.0150409e+00 
x[  2] =   1.0303266e+00 
x[  3] =   1.0617248e+00 
f(x) =   1.2043904e-03


TIME:


Sequential:             42.096s

Omp4t (5 run average):  10.357s (4.064 times faster)
[10.373s, 10.330s, 10.392s, 10.368s, 10.323s]

mpi4c (5 run average):
[10.701s, 10.700s, 11.098s, 10.704s, 10.693s]