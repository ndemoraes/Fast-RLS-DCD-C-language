//
// Created by Naomi Sutcliffe de Moraes
//


#ifndef ECHO_CANCELING_H
#include "matrix_operations/Standard_operations.h"

double *filt(double *y, const double *h, const double *x, long num_time_steps, long M);

void RLS(double *w, double *e, double *x, const double *s, long M, double lambda, double delta, long num_time_steps);


void test_DCBF(const gsl_rng *r);
void test_fRLSDCD(long num_time_steps, const gsl_rng *r);
void test_RLS(long num_time_steps, const gsl_rng *r);
void test_RLSDCD(long num_time_steps, const gsl_rng *r);

#define ECHO_CANCELING_H

#endif //ECHO_CANCELING_H
