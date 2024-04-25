#ifndef STANDARD_OPERATIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define MATRIX_SIZE 1000
#define ARRAY_SIZE 500500  // number of elements in the upper triangle of the matrix
// for MATRIX_SIZE = 10,  ARRAY_SIZE = 10 + 9 + 8 + 7 + 6 + 5+4+3+2+1 (or 11 * 5) = 55
// for MATRIX_SIZE = 100, ARRAY_SIZE = (100 + 1) + (99 + 2) + ... or 101 * 50 = 5050
// for MATRIX_SIZE = 1000, ARRAY_SIZE = (1000 + 1) + ... or 1001 * 500 = 500,500


typedef double REGULAR_MATRIX[MATRIX_SIZE][MATRIX_SIZE]; // [row][col], C stores matrices by rows
// to create a variable of this type, use
//     REGULAR_MATRIX my_matrix;
// to allocate memory for a variable of this type, use:
//     REGULAR_MATRIX *regular_matrix_pt = malloc(sizeof(REGULAR_MATRIX));
// to set all matrix elements to zero, call the zero_regular_matrix function.
// alternatively, the initialize_matrix function sets all diagonal elements to a
// scalar value, and all other elements to zero.
// to pass the pointer to a function (instead of copying a large matrix), use:
//     void function_name((REGULAR_MATRIX *mtrx) { ... }
// inside the function, you will have to use (*mtrx)[row][col] to access elements


REGULAR_MATRIX *initialize_5x5();
REGULAR_MATRIX *initialize_matrix(double scalar);
void zero_regular_matrix(REGULAR_MATRIX * mtrx);

// functions for algorithm
void update_matrix_w_backward_vector(REGULAR_MATRIX *mtrx, const double *vtr, long start_slice, const double scalar);
void update_matrix(REGULAR_MATRIX *mtrx, const double *vtr, double scalar);
void update_vector(double *vtr, const REGULAR_MATRIX *mtrx, double scalar, long col_num);

// extra linear algebra and math functions
void print_regular_matrix(REGULAR_MATRIX *regular_matrix);

void print_array(double *array, long size);
void zero_vector(double *vector, long size);
//void matrix_x_vector(double **matrix, double *vector, double *result);

// backwards vector functions (to avoid copying a slice of x to the vector u
void matrix_x_backward_vector(const REGULAR_MATRIX *matrix, const double *vector, long start_index, double *result);
double backward_dot_product(const double *backward, const double *forward, long start_index);
double get_backward_element(const double *vector, long start_slice, long element);


void outer_product(const double *column_vector, const double *row_vector, REGULAR_MATRIX *matrix);
double Euclidean_distance(double *vector1, double *vector2, long length);
double mean_of_vector(const double *vector, long length);
double minimum_of_vector(const double *vector, long length);
long argmax(double *vector, long length);
double sign(double value);

double gaussian_rand(const gsl_rng *r, double sigma);
double *gaussian_rand_vector(const gsl_rng *r, double sigma, long size);
long min(long a, long b);

#define STANDARD_OPERATIONS_H

#endif // STANDARD_OPERATIONS_H