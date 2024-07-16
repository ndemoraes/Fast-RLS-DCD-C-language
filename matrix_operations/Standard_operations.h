//
// Created by Naomi Sutcliffe de Moraes
//
#ifndef STANDARD_OPERATIONS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

// The matrix size for all simulations is hard-coded here (MATRIX_SIZE)
// In a typical application, the matrix size would be constant and known at compile time.
// For the DCBF matrix, in order to save space, we store the upper triangle of the matrix in a one-dimensional array
// The size of the one-dimensional array must also be hard-coded (ARRAY_SIZE)
#define MATRIX_SIZE 1000
#define ARRAY_SIZE 500500  // number of elements in the upper triangle of the matrix (including the diagonal)

// for MATRIX_SIZE = 10,  ARRAY_SIZE = 10 + 9 + 8 + 7 + 6 + 5 + 4 + 3 + 2 + 1 (or 11 * 10/2) = 55
// for MATRIX_SIZE = 100, ARRAY_SIZE = (100 + 1) + (99 + 2) + ... or 101 * (100/2) = 5050
// for MATRIX_SIZE = 1000, ARRAY_SIZE = (1000 + 1) + ... or 1001 * (1000/2) = 500,500

// This is the type for a regular matrix. The DCBF matrix is defined in DCBF_operations.h
typedef double REGULAR_MATRIX[MATRIX_SIZE][MATRIX_SIZE]; // [row][col], C stores matrices by rows

// to create a variable of this type, use:
//     REGULAR_MATRIX *regular_matrix_pt = malloc(sizeof(REGULAR_MATRIX));
// to set all matrix elements to zero, call the zero_regular_matrix function.
// alternatively, the initialize_matrix function sets all diagonal elements to a
// scalar value, and all other elements to zero.
// to pass the pointer to a function (instead of copying a large matrix), use:
//     void function_name((REGULAR_MATRIX *mtrx)) { ... }
// inside the function, you will have to use (*mtrx)[row][col] to access elements

// The functions below are specific to the regular matrix type
// or are generic math functions that can be used with arrays or numbers

//REGULAR_MATRIX *initialize_5x5();  // used for unit tests
REGULAR_MATRIX *initialize_matrix(double scalar);
void zero_regular_matrix(REGULAR_MATRIX * mtrx);

// functions for algorithm
void update_matrix_w_backward_vector(REGULAR_MATRIX *mtrx, const double *vtr, long start_slice, const double scalar);
//void update_matrix(REGULAR_MATRIX *mtrx, const double *vtr, double scalar); // used in unit tests
void update_vector(double *vtr, const REGULAR_MATRIX *mtrx, double scalar, long col_num);

// Print functions for debugging
void print_regular_matrix(REGULAR_MATRIX *regular_matrix);  // used for debugging
void print_array(double *array, long size); // used for debugging

// extra linear algebra and math functions
void zero_vector(double *vector, long size);
//void matrix_x_vector(double **matrix, double *vector, double *result);  // used in unit tests

// There are two versions of some functions. Some have the word BACKWARDS in the name, the corresponding similar
// function does not.
// In the Julia code, the variable u is used as an alias to a backwards slice of the vector x. (see the Julia code)

// In an earlier version of the C code, we used a separate vector u to store the reverse of a section of x.
// These functions do not have the word "backward" in them.

// In the current C version, in order to avoid the cost
// of creating a new vector (u in Julia) that is the reverse of a section of x,
// we continue using x and access it backwards using the get_backward_element function. All functions
// with "backward" in the name are used to access a section of x backwards.

// The original "non-backward" functions have not been eliminated because they are helpful for debugging.

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