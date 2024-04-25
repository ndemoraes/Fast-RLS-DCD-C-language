//
// Created by Naomi Nascimento on 3/28/24.
//

#ifndef FASTDCD_STANDARD_RLSDCD_H
#define FASTDCD_STANDARD_RLSDCD_H
#include "../matrix_operations/Standard_operations.h"

void DCD(double H, long B, long num_dcd_iterations, REGULAR_MATRIX *R_phi, double *res, double *delta_w);
void RLSDCD(double *w, double *e, double*x, const double *s, double lambda, double delta, long num_dcd_iterations, long num_time_steps);
#endif //FASTDCD_STANDARD_RLSDCD_H
