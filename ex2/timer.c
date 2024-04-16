// Timer for the HPCSE I course
//

#include <sys/time.h>
#include <stddef.h>

double get_wtime(void)
{
	struct timeval t;

	gettimeofday(&t, NULL);

	return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}
