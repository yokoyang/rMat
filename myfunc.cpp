#include "myfunc.h"
int vec2psi(gsl_vector* psi, gsl_vector *inc, gsl_vector *skp,
	int inclu_len, int skip_len) {
	size_t idx;
	for (idx = 0; idx < inc->size; ++idx) {
		if (gsl_vector_get(inc, idx) == 0 && gsl_vector_get(skp, idx) == 0) {
			gsl_vector_set(psi, idx, 0.5);
			//ignore the difference between int and float
		}
		else if (gsl_vector_get(inc, idx) == 0) {
			gsl_vector_set(psi, idx, 0.01);
		}
		else if (gsl_vector_get(skp, idx) == 0) {
			gsl_vector_set(psi, idx, 0.99);
		}
		else {
			gsl_vector_set(psi, idx, gsl_vector_get(inc, idx) / inclu_len /
				(gsl_vector_get(inc, idx) / inclu_len + gsl_vector_get(skp, idx) / skip_len));
		}
	}
	return 0;
}

inline double logit(double input_x)
{
	double value;
	if (input_x < 0.01)
	{
		input_x = 0.01;
	}
	else if (input_x > 0.99) {
		input_x = 0.99;
	}
	value = gsl_log1p(1 - (input_x / (1 - input_x)));
	return input_x;
}


void mle_result_set(mle_result * mle, double sum, gsl_vector * psi, double beta0, double beta1, double var)
{
	mle->sum = sum;
	mle->params.psi = psi;
	mle->params.beta0 = beta0;
	mle->params.beta1 = beta1;
	mle->params.var = var;

}



double myfunc_individual(const double x[], va_list argv)
{
	double I = va_arg(argv, double), S = va_arg(argv, double);
	double beta0 = va_arg(argv, double), beta1 = va_arg(argv, double);
	double var = va_arg(argv, double);
	double gene_exp = va_arg(argv, double);  //Don't know the type of gene_exp
	int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int);
	double new_psi = inclu_len * x[0] / (inclu_len * x[0] + skip_len * (1 - x[0]));

	// TODO This change the result.
	return 1.0 * (pow((logit(x[0]) - (beta0 + beta1 * gene_exp)), 2) / (2 * var) - (I * log(new_psi) + S * log(1 - new_psi) + log(x[0]) + log(1 - x[0]) + log(sqrt(var))));
}


void myfunc_individual_der(const double x[], double res[], va_list argv) {
	double I = va_arg(argv, double), S = va_arg(argv, double);
	double beta0 = va_arg(argv, double), beta1 = va_arg(argv, double);
	double var = va_arg(argv, double), gene_exp = va_arg(argv, double);
	int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int);
	double part1 = (beta0 + beta1 * gene_exp - logit(x[0])) / (var *(double)x[0] * (1 - x[0]) - 1.0 / x[0] + 1.0 / (1 - x[0]));
	//为了和part1的类型统一，把float(x*(1-x))改成double(x*(1-x))
	double part2 = I * 1.0 / x[0] + S * 1.0 / (x[0] - 1) - (I + S) *((inclu_len - skip_len) / (inclu_len * x[0] + skip_len * (1 - x[0])));
	//因为part2和part3有很多重复的部分所以整合了一下
	res[0] = -1.0 * (part1 + part2);
	return;
}

int MLE_marginal_iteration(gsl_vector * i1, gsl_vector * i2, gsl_vector * s1, gsl_vector * s2, const int inclu_len, const int skip_len, mle_result * mle)
{
	double slope, intercept;
	double cov00 = 0.0;
	double cov01 = 0.0;
	double cov11 = 0.0;
	double sumsq = 0.0;
	int i;
	gsl_vector* psi = gsl_vector_calloc(i1->size);
	if (!vec2psi(psi, i1, s1, inclu_len, skip_len))
	{
		return 0;//failed
	}
	gsl_vector *ltemp = gsl_vector_calloc(psi->size);
	for (i = 0; i != psi->size; i++) {
		double temp = gsl_vector_get(psi, i);
		temp = logit(temp);
		gsl_vector_set(ltemp, i, temp);
	}
	double var = gsl_stats_variance(ltemp->data, 1, ltemp->size) / 2;
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
	double iter_cutoff = 1.0, iter_maxrun = 100.0, count = 0.0, previous_sum = 0.0, current_sum = 0.0;
	while ((iter_cutoff>0.01)&(count <= iter_maxrun))
	{
		count += 1;


		/*xopt = fmin_l_bfgs_b(myfunc_marginal_integrated, beta, myfunc_marginal_integrated_der,
		args = [inc, skc, psi, var, gene_expression, effective_inclusion_length, effective_skipping_length],
		bounds = [[-99999999.0, 99999999.0], [-99999999.0, 99999999.0]], iprint = -1);


		beta = xopt[0];*/


		gsl_vector *new_psi = gsl_vector_calloc(psi->size);
		current_sum = 0;
		for (int i = 0; i != psi->size; i++) {
			/*xopt = fmin_l_bfgs_b(myfunc_individual, [psi[i]], myfunc_individual_der, args = [inc[i], skc[i], beta[0], beta[1], var,
			gene_expression[i], effective_inclusion_length, effective_skipping_length], bounds = [[0.01, 0.99]], iprint = -1);
			new_psi.append(float(xopt[0])); current_sum+=float(xopt[1]);*/
		}
		psi = new_psi;
		iter_cutoff = fabs(previous_sum - current_sum);
		previous_sum = current_sum;
	}
	if (count > iter_maxrun) {
		mle_result_set(mle, current_sum, psi, 0, 0, var);
	}
	else
	{
		mle_result_set(mle, current_sum, psi, beta[0], beta[1], var);
	}


	return 1;
}

int MLE_marginal_iteration_constrain(gsl_vector * i1, gsl_vector * i2, gsl_vector * s1, gsl_vector * s2, const int inclu_len, const int skip_len, mle_result * mle,double base_line)
{
	double slope, intercept;
	double cov00 = 0.0;
	double cov01 = 0.0;
	double cov11 = 0.0;
	double sumsq = 0.0;
	int i;
	gsl_vector* psi = gsl_vector_calloc(i1->size);
	if (!vec2psi(psi, i1, s1, inclu_len, skip_len))
	{
		return 0;//failed
	}
	gsl_vector *ltemp = gsl_vector_calloc(psi->size);
	for (i = 0; i != psi->size; i++) {
		double temp = gsl_vector_get(psi, i);
		temp = logit(temp);
		gsl_vector_set(ltemp, i, temp);
	}
	double var = gsl_stats_variance(ltemp->data, 1, ltemp->size) / 2;
	if (var <= 0.01) { var = 0.01; }
	gsl_fit_linear(s1->data, 1, ltemp->data, 1, ltemp->size, &slope, &intercept, &cov00, &cov01, &cov11, &sumsq);
	double beta[2] = { slope ,intercept };
	double iter_cutoff = 1.0, iter_maxrun = 100.0, count = 0.0, previous_sum = 0.0, current_sum = 0.0;
	while ((iter_cutoff > 0.01)&(count <= iter_maxrun)) 
	{
		count++;
		/*xopt = fmin_l_bfgs_b(myfunc_marginal_integrated, beta, myfunc_marginal_integrated_der,
			args = [inc, skc, psi, var, gene_expression, effective_inclusion_length, effective_skipping_length],
			bounds = [[-9999999.0, 9999999.0], [base_line, base_line]], iprint = -1);
			beta=xopt[0];*/
		gsl_vector *new_psi = gsl_vector_calloc(psi->size);
		current_sum = 0;
		for (int i = 0; i != psi->size; i++) {
			
			//注意里面需要用到base_line


			/*xopt = fmin_l_bfgs_b(myfunc_individual, [psi[i]], myfunc_individual_der, args = [inc[i], skc[i], beta[0], beta[1], var,
			gene_expression[i], effective_inclusion_length, effective_skipping_length], bounds = [[0.01, 0.99]], iprint = -1);
			new_psi.append(float(xopt[0])); current_sum+=float(xopt[1]);*/

		}
		psi = new_psi;
		iter_cutoff = fabs(previous_sum - current_sum);
		previous_sum = current_sum;
	}



	return 0;
}




double* likelihood_test(gsl_vector * i1, gsl_vector * i2, gsl_vector * s1, gsl_vector * s2, int inclu_len, int skip_len, double Delta)
{
	mle_result *result;
	result = (mle_result*)malloc(sizeof(mle_result));
	mle_result *result_constrain;
	result_constrain = (mle_result*)malloc(sizeof(mle_result));
	double *likelireasult;
	likelireasult = (double*)malloc(sizeof(double) * 3);
	MLE_marginal_iteration(i1, i2, s1, s2, inclu_len, skip_len, result);
	double Base = fabs(Delta);
	if (Base < 0.01) 
	{
		if (fabs(result->params.beta1 - 0.0) < 0.0001)
		{
			likelireasult[0] = 1;
			likelireasult[1] = result->params.beta0;
			likelireasult[2] = result->params.beta1;
			return likelireasult;
		}
		else
		{
			MLE_marginal_iteration_constrain(i1, i2, s1, s2, inclu_len, skip_len, result_constrain,0.0);
			likelireasult[0]=gsl_cdf_chisq_P(2 * (fabs(result_constrain->sum - result->sum)), 1);
			likelireasult[1] = result->params.beta0;
			likelireasult[2] = result->params.beta1;
			return likelireasult;
		}
	}
	else {
		if (fabs(result->params.beta1) < Base)
		{
			likelireasult[0] = 1;
			likelireasult[1] = result->params.beta0;
			likelireasult[2] = result->params.beta1;
			return likelireasult;
		}
		else {
			if (result->params.beta1 < 0.0) {
				Base = -Base;
			}
			MLE_marginal_iteration_constrain(i1, i2, s1, s2, inclu_len, skip_len, result_constrain, Base);
			likelireasult[0] = gsl_cdf_chisq_P(2 * (fabs(result_constrain->sum - result->sum)), 1);
			likelireasult[1] = result->params.beta0;
			likelireasult[2] = result->params.beta1;
			return likelireasult;
		}
	}
}