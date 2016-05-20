#include "time.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_statistics_double.h"
#include "math.h"
#include "stdio.h"
#include "stdarg.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_rng.h>

//One-fold derivation function for integrated marginal objective function


int main(int argc, char *argv[])
{


	
	gsl_vector * v = gsl_vector_alloc(11);
	for (int i = 0; i < 10; i++)
	{
		gsl_vector_set(v, i, 1.23 + i);
	}
	for (int i = 0; i < 10; i++)
	{
	//	printf("v_%d = %g\n", i, gsl_vector_get(v, i));
	//	printf("gsl_vector_max_index %d\n", gsl_vector_max_index(v));
	}
	
	//printf("gsl_vector_max_index %d\n", v->size);
	//gsl_vector_add_constant(v, 1);
	//printf("v_%d = %g\n", gsl_vector_get(v, 9));
	
	double data[5] = { 17.2, 18.1, 16.5, 18.3, 12.6 };
	double t = gsl_stats_variance(data, 1, 5);
	printf(" variance = %f", t);
	printf(" variance = %d", v->owner);

	gsl_vector_free(v);
	getchar();
	


	double dur;
	clock_t start, end;
	int times = 300, i = 0;
	size_t len = 200;
	gsl_vector *vec1 = gsl_vector_alloc(len), *vec2 = gsl_vector_alloc(len), *vec3 = NULL;

	start = clock();
	for (i = 0; i < times; ++i) {
		vec3 = gsl_vector_alloc(vec1->size);
		gsl_vector_memcpy(vec3, vec1);
		gsl_vector_add(vec3, vec2);
		gsl_vector_free(vec3);
	}
	end = clock();
	dur = (double)(end - start);
	printf("using vec3: %f\n", dur / CLOCKS_PER_SEC);

	start = clock();
	for (i = 0; i < times; ++i) {
		gsl_vector_add(vec1, vec2);
		gsl_vector_sub(vec1, vec2);
	}
	end = clock();
	dur = (double)(end - start);
	printf("using vec1: %f\n", dur / CLOCKS_PER_SEC);

	start = clock();
	for (i = 0; i < 100; ++i) {
		gsl_cdf_gaussian_Pinv(0.5, 1);
	}
	end = clock();
	dur = (double)(end - start);
	printf("time for ppf: %f\n", dur / CLOCKS_PER_SEC);

	printf("sizeof double: %ld\n", sizeof(double));

	return 0;
}
