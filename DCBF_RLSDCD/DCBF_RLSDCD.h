//
// Created by Naomi J Sutcliffe de Moraes
//

#ifndef FASTDCD_DCBF_RLSDCD_H
#define FASTDCD_DCBF_RLSDCD_H
#import "../matrix_operations/DCBF_operations.h"

void fDCD(double H, long B, long num_dcd_iterations, struct DCBF_Matrix *R_phi, double *res, double *delta_w);

void fRLSDCD(double *w, double *e, double*x, const double *s, double lambda, double delta, long num_dcd_iterations, long num_time_steps);
#endif //FASTDCD_DCBF_RLSDCD_H
