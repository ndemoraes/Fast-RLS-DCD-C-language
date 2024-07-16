//
// Created by Naomi Sutcliffe de Moraes
//
#ifndef DCBF_OPERATIONS_H
#include <stdlib.h>
#include <stdio.h>
#include "Standard_operations.h"

// see notes in standard_operations.h, as the functions below are similar

// structure to store dmtrx in DIAGONAL, CIRCULAR-BUFFER FORMAT (DCBF)
struct DCBF_Matrix
{
    double dmtrx[ARRAY_SIZE]; // 1-D array representation of DCBF matrix
    double *top[MATRIX_SIZE]; // array of pointers to first element of each diagonal at beginning, these should not change
    long cb_head[MATRIX_SIZE]; // relative location with respect to the top element (head) of the circular buffer for each diagonal, 0..MATRIX_SIZE - 1, vary as circular buffer filled
    long diag_num_elements[MATRIX_SIZE]; // Number of elements in each diagonal
};

long julia_mod1(long dividend, long divisor); // to imitate the mod1 function in Julia
double *pointer_to_head(long diagonal, const struct DCBF_Matrix *mtrx);

struct DCBF_Matrix *initialize_DCBF_Matrix(double delta);

// constructor for DCBF matrix, used in unit tests
void initialize_DCBF_Matrix_with_RegularMatrix(REGULAR_MATRIX *regular_matrix, struct DCBF_Matrix *initialized );

void print_DCBF_Matrix(struct DCBF_Matrix *mtrx);
void print_vector(double vtr[]);

void update_DCBF_Matrix(struct DCBF_Matrix *mtrx, const double *vtr, double scalar);
void update_DCBF_matrix_w_backward_vector(struct DCBF_Matrix *mtrx, const double *vtr, long start_slice, double scalar);

void update_DCBF_Vector(double *vtr, const struct DCBF_Matrix *mtrx, double scalar, long col_numb);

#define DCBF_OPERATIONS_H
#endif // DCBF_OPERATIONS_H