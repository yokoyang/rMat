#ifndef TRANS_MYFUNC
#define TRANS_MYFUNC


#include <stdarg.h>
#include "gsl/gsl_vector.h"
#include "../include/type.h"

double myfunc_individual(const double x[], va_list argv);

void myfunc_individual_der(const double x[], double res[], va_list argv);

double myfunc_marginal_integrated(const double x[], va_list argv);

void myfunc_marginal_integrated_der(const double x[], double res[], va_list argv);

int vec2psi(gsl_vector* psi, gsl_vector *inc, gsl_vector *skp,
            int inclu_len, int skip_len);

void mle_result_set(mle_result * mle, double sum, gsl_vector * psi, double beta0, double beta1, double var);

int MLE_marginal_iteration(gsl_vector* inc, gsl_vector* skp, gsl_vector* gene_exp, 
	const int inclu_len, const int skip_len, mle_result* mle);
 
int MLE_marginal_iteration_constrain(gsl_vector* inc, gsl_vector* skp, gsl_vector* gene_exp, 
    	const int inclu_len, const int skip_len, mle_result * mle, double base_line);

double* likelihood_test(gsl_vector* inc, gsl_vector* skp, gsl_vector* gene_exp, 
        int inclu_len, int skip_len, double delta, int flag, char* id);

void* thread_wrapper_for_LT(void* arg);

void* batch_wrapper_for_LT(void* arg);

#endif
