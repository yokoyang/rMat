#ifndef TRANS_MYFUNC
#define TRANS_MYFUNC


#include "stdarg.h"
#include "gsl/gsl_vector.h"
#include "type.h"
#include "gsl/gsl_math.h"
#include <gsl/gsl_statistics.h>  
#include <gsl/gsl_fit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>


void mle_result_set(mle_result* mle, double sum, gsl_vector* psi1, gsl_vector* psi2,
	double beta0, double beta1, double var1, double var2);

void log_vector(gsl_vector *vec);

double myfunc_multivar(const double x[], va_list argv);

void myfunc_multivar_der(const double x[], double res[], va_list argv);

double myfunc_1_2(const double x[], va_list argv);

void myfunc_der_1_2(const double x[], double res[], va_list argv);

double myfunc_marginal_1_2(const double x[], va_list argv);

void myfunc_marginal_1_2_der(const double x[], double res[], va_list argv);

double myfunc_individual(const double x[], va_list argv);

void myfunc_individual_der(const double x[], double res[], va_list argv);

double myfunc_marginal(const double x[], va_list argv);

void myfunc_marginal_der(const double x[], double res[], va_list argv);

int MLE_marginal_iteration(gsl_vector* i1, gsl_vector* i2,
	gsl_vector* s1, gsl_vector* s2,
	const int inclu_len, const int skip_len,
	mle_result* mle);

int MLE_marginal_iteration_constrain(gsl_vector* i1, gsl_vector* i2,
	gsl_vector* s1, gsl_vector* s2,
	const int inclu_len, const int skip_len,
	mle_result* mle);

void* thread_wrapper_for_LT(void* arg);

void* batch_wrapper_for_LT(void* arg);

double likelihood_test(gsl_vector *i1, gsl_vector *i2, gsl_vector *s1, gsl_vector *s2,
	int inclu_len, int skip_len, int flag, char* id);

gsl_vector* vec2psi(gsl_vector* psi, gsl_vector *inc, gsl_vector *skp,
	int inclu_len, int skip_len);

///////////////
double logit();



inline double logit(double in)
{

	return in;
}


void mle_result_set(mle_result * mle, double sum, gsl_vector * psi1, gsl_vector * psi2, double beta0, double beta1, double var1, double var2)
{
	mle->sum = sum;
	mle->params.beta0 = beta0;
	mle->params.beta1 = beta1;
	mle->params.psi1 = psi1;
	mle->params.psi2 = psi2;
	mle->params.var1 = var1;
	mle->params.var2 = var2;

}

int MLE_marginal_iteration(gsl_vector * i1, gsl_vector * i2, gsl_vector * s1, gsl_vector * s2, const int inclu_len, const int skip_len, mle_result * mle)
{
	double slope, intercept;
	double cov00 = 0.0;
	double cov01 = 0.0;
	double cov11 = 0.0;
	double sumsq = 0.0;
	gsl_vector* psi = vec2psi(i1, i2, s2, inclu_len, skip_len);
	gsl_vector *ltemp = gsl_vector_calloc(psi->size);
	for (int k = 0; k != psi->size; k++) {
		double temp = gsl_vector_get(psi, k);
		temp = logit(temp);
		gsl_vector_set(ltemp, k, temp);
	}
	double var=gsl_stats_variance(ltemp->data,1,ltemp->size);
	if (var <= 0.01) { var = 0.01; }
		//int gsl_fit_linear(const double * x, const size_t xstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq)
		//x是自变量
		// xstrside是步长
		// y是因变量
		//ystride是y的步长
		//n是数据组数
		//c0 储存计算出来的a （f(x) = a * x + b）
		// c1 储存计算出来的b （f(x) = a * x + b）
		// cov00 cov01 cov11是点的协方差 我的博客中有介绍协方差
		// 没有cov10是因为它和cov01是相等的
		//sumsq是残差

	gsl_fit_linear(s1->data, 1, ltemp->data, 1, ltemp->size, &slope, &intercept, &cov00, &cov01, &cov11, &sumsq);
	double beta[2] = { slope ,intercept };
	int iter_cutoff = 1, iter_maxrun = 100, count = 0, previous_sum = 0;
	while ((iter_cutoff>0.01)&(count <= iter_maxrun))
	{
		count += 1;


	}



	return 0;
}




double likelihood_test(gsl_vector * i1, gsl_vector * i2, gsl_vector * s1, gsl_vector * s2, int inclu_len, int skip_len, int flag, char * id)
{

	//double result= MLE_marginal_iteration(i1,i2,s1,s2, inclu_len, skip_len, )
	//

	//return 0.0;
}

#endif
