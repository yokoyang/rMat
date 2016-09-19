#ifndef TRANS_UTIL
#define TRANS_UTIL


#define prefetch(x) __builtin_prefetch(x)
#define prefetch_for_r(x) __builtin_prefetch(x, 0, 3)
#define prefetch_for_w(x) __builtin_prefetch(x, 1, 3)


#include <stdarg.h>
#include "gsl/gsl_vector.h"
#include "../include/type.h"

double cuscumsum(gsl_vector *vec, double (*fun) (double, va_list), int argc, ...);

double* cuscumsum_der(gsl_vector *vec, void (*fun) (double, double* , va_list), int argc, ...) ;

double logit(double i);

double sum_for_marginal_integrated(const double i, va_list argv);

void sum_for_marginal_integrated_der(const double i, double* tmp, va_list argv);

double myfunc_marginal_2der(const double x, const double I, const double S,
                            const double beta0, const double beta1,
                            const double var, const double gene_exp,
                            const int inclu_len, const int skip_len);

/*
double sum_for_marginal_integrated_1_der(const double i, va_list argv);

double sum_for_marginal_integrated_2_der(const double i, va_list argv);
*/

double l_bfgs_b_wrapper(integer n, integer m, doublereal x[], doublereal l[],
                        doublereal u[], integer nbd[],
                        double (*fp) (const double x[], va_list argv),
                        void (*gp) (const double x[], double res[], va_list argv),
                        doublereal factr, doublereal pgtol, integer iprint,
                        int maxfun, int maxiter, int argc, ...);

void mp_threadpool(int nthread, int ntask, void* (*func)(void *), void** datum, void **ret);

//int parse_file(const char* file_exon, const char* file_gene, diff_list_node* list, char** title_list, int gene_size);
int parse_file(FILE* fin_exon, FILE* fin_gene, diff_list_node* list, char** title_list, int gene_size);

void parse_line_exon(char* str, char* id, gsl_vector** inc, gsl_vector** skp, int* inclu_len, int* skip_len);

void parse_line_gene(char* str, gsl_vector** gene, int gene_size);

int str_to_vector(char* input, gsl_vector** vec);

odiff* diff_alloc(gsl_vector* inc, gsl_vector* skp, gsl_vector* gene_exp, int inclu_len, int skip_len, int flag, char* id);

int diff_append(diff_list_node* header, odiff* data);

int diff_insert(diff_list_node* header, odiff* data, int idx);

int diff_get_next(diff_list_node* header, odiff* data);

int diff_get_at(diff_list_node* header, odiff* data, int idx);

//void rm_NA_Samples(diff_list_node* node);

#endif
