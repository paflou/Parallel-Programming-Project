// numerical integration / sequential code 
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include "timer.h"

#define THREAD_NUM 8

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
double a = 1.0;
double b = 4.0;
unsigned long const n = 1e9;
double S = 0;

typedef struct {
    int start;
    int end;
} ThreadArgs;

double f(double x)
{
	return log(x)*sqrt(x);
}

void *func(void * arg)
{
    ThreadArgs *args = (ThreadArgs*) arg;

    const double dx = (b-a)/n;
    double S_local = 0;
	for (int i = args->start; i < args->end; i++) {
		double xi = a + (i + 0.5)*dx;
        S_local += f(xi);
        //printf("%f\n",S_local);
	}
    S_local *=dx;
    pthread_mutex_lock(&mutex);
	S += S_local;
    pthread_mutex_unlock(&mutex);
    printf("thread %ld finished\n",pthread_self());
    return NULL;
}

// WolframAlpha: integral_1^4 log(x) sqrt(x) dx = 4/9 (4 log(64)-7) ~~ 4.28245881486164
int main(int argc, char *argv[])
{
    pthread_t id[THREAD_NUM];
    int iterations_per_thread = n/THREAD_NUM;
    ThreadArgs thread_args[THREAD_NUM];
    
	double t0 = get_wtime();
    for (int i = 0; i < THREAD_NUM; i++) {
        thread_args[i].start = i * iterations_per_thread;
        thread_args[i].end = (i + 1) * iterations_per_thread;

        pthread_create(&id[i], NULL, func, (void *)&thread_args[i]);
    }
    for (long i = 0; i < THREAD_NUM; i++) {
        pthread_join(id[i], NULL);
    }

	double t1 = get_wtime();

	printf("Time=%lf seconds, Result=%.8f\n", t1-t0, S);
    pthread_mutex_destroy(&mutex);
	return 0;
}
