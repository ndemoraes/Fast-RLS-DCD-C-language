//
// Created by Naomi Sutcliffe de Moraes
//
#include "Standard_operations.h"

// create regular dmtrx for 5x5 test
REGULAR_MATRIX *initialize_5x5()
{
    REGULAR_MATRIX *regular_matrix_pt = (REGULAR_MATRIX *)malloc(sizeof(REGULAR_MATRIX));
    (*regular_matrix_pt)[0][0] = 11.0;
    (*regular_matrix_pt)[0][1] = 12.0;
    (*regular_matrix_pt)[1][0] = 12.0;
    (*regular_matrix_pt)[0][2] = 13.0;
    (*regular_matrix_pt)[2][0] = 13.0;
    (*regular_matrix_pt)[0][3] = 14.0;
    (*regular_matrix_pt)[3][0] = 14.0;
    (*regular_matrix_pt)[0][4] = 15.0;
    (*regular_matrix_pt)[4][0] = 15.0;

    (*regular_matrix_pt)[1][1] = 22.0;
    (*regular_matrix_pt)[1][2] = 23.0;
    (*regular_matrix_pt)[2][1] = 23.0;
    (*regular_matrix_pt)[1][3] = 24.0;
    (*regular_matrix_pt)[3][1] = 24.0;
    (*regular_matrix_pt)[1][4] = 25.0;
    (*regular_matrix_pt)[4][1] = 25.0;

    (*regular_matrix_pt)[2][2] = 33.0;
    (*regular_matrix_pt)[2][3] = 34.0;
    (*regular_matrix_pt)[3][2] = 34.0;
    (*regular_matrix_pt)[2][4] = 35.0;
    (*regular_matrix_pt)[4][2] = 35.0;

    (*regular_matrix_pt)[3][3] = 44.0;
    (*regular_matrix_pt)[3][4] = 45.0;
    (*regular_matrix_pt)[4][3] = 45.0;

    (*regular_matrix_pt)[4][4] = 55.0;

    return regular_matrix_pt;
}

REGULAR_MATRIX *initialize_matrix(double scalar)
{
    REGULAR_MATRIX *regular_matrix_pt = malloc(sizeof(REGULAR_MATRIX));
    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        for (long col = 0; col < MATRIX_SIZE; ++col)
        {

            if (row == col) {
                (*regular_matrix_pt)[row][col] = scalar;
            }
            else
            {
                (*regular_matrix_pt)[row][col] = 0.0;
            }
        }
    }

    return regular_matrix_pt;
}

void zero_regular_matrix(REGULAR_MATRIX * regular_matrix_pt)
{
    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        for (long col = 0; col < MATRIX_SIZE; ++col)
        {
                (*regular_matrix_pt)[row][col] = 0.0;
        }
    }
}

void print_regular_matrix(REGULAR_MATRIX *regular_matrix)
{

    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        for (long col = 0; col < MATRIX_SIZE; ++col)
        {
            printf("%f  ", (*regular_matrix)[row][col]);
        }
        printf("\n"); // return at the end of each row
    }
}

void print_array(double *array, long size)
{
    for (long element = 0; element < size; ++element)
    {
        printf("%f  ", array[element]);
    }
    printf("\n");
}
void zero_vector(double *vector, long size)
{
    for (long i = 0; i < size; ++i) {vector[i] = 0.0;}
}

void matrix_x_vector(double **matrix, const double *vector, double *result)
{
    long i, j;
    for (i = 0; i < MATRIX_SIZE; ++i)
    {
        result[i] = 0.0;
        for (j = 0; j < MATRIX_SIZE; ++j)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    } // end for
} // end matrix_x_vector

void matrix_x_backward_vector(const REGULAR_MATRIX *matrix, const double *vector, long start_index, double *result)
{
    // there is a matrix of size MATRIX_SIZE X MATRIX_SIZE, a vector that is size >= MATRIX_SIZE. You pick a starting point inside the
    // vector and step through it backwards, as if it were in reverse order.
    long row, col;
    for (row = 0; row < MATRIX_SIZE; ++row)
    {
        result[row] = 0.0;
        for (col = 0; col < MATRIX_SIZE; ++col)
        {
            result[row] += (*matrix)[row][col] * vector[start_index - col];
        }
    }
} // end matrix_x_backward_vector

double backward_dot_product(const double *backward, const double *forward, long start_index)
{
    // dot product between backward and forward vectors
    double y = 0.0;
    for (long elem = 0; elem < MATRIX_SIZE; ++elem)
    {
        y = y + forward[elem] * backward[start_index - elem];
    }
    return y;
}

double get_backward_element(const double *vector, long start_slice, long element)
{
    // you have a vector of length N, and you want to obtain an element that is offset from
    // element start_slice, going backwards.
    return vector[start_slice - element];
}


void outer_product(const double * column_vector, const double *row_vector, REGULAR_MATRIX *matrix)
{
    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        for (long col = 0; col < MATRIX_SIZE; ++col)
            (*matrix)[row][col] = column_vector[row] * row_vector[col];
    }
}

double Euclidean_distance(double *vector1, double *vector2, long length)
{
    double distance = 0;
    for (long i = 0; i < length; i++)
    {
        distance = distance + pow(  (vector1[i] - vector2[i]), 2  );
    }
    distance = sqrt(distance);
    return distance;
}

double mean_of_vector(const double *vector, long length)
{
    double mean = 0.0;
    for (long i = 0; i < length; ++i)
    {
        mean = mean + vector[i];
    }
    mean = mean / (double)length;
    return mean;
}

double minimum_of_vector(const double *vector, long length)
{
    double minimum = vector[0];
    for (long i = 1; i < length; ++i)
    {
        if (vector[i] < minimum)
        {
            minimum = vector[i];
        }
    }
    return minimum;
}

long argmax(double *vector, long length)
{
    double max_value = fabs(vector[0]);
    long max_index = 0;
    for (long j = 1; j < length; ++j)
    {
        if (  fabs(vector[j]) > max_value  )
        {
            max_value = fabs(vector[j]);
            max_index = j;
        }
    }
    return max_index;
}

double sign(double value)
{
    if (value == 0.0) {return 0.0;}
    else if (value > 0.0) {return 1.0;}
    else {return -1.0;}
}
// *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
void update_matrix_w_backward_vector(REGULAR_MATRIX *mtrx, const double *vtr, long start_slice, const double scalar)
{
    // format is mtrx[row][col]

    // store first column (index 0) for later use
    double first_col[MATRIX_SIZE];
    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        first_col[row] = (*mtrx)[row][0];
    }

    // loop backwards from MATRIX_SIZE - 2 to 0 in increments of -1
    for (long row = MATRIX_SIZE - 2; row >= 0; --row)
    {
        for (long col = MATRIX_SIZE - 2; col >= 0; --col)
        {
            // shift elements down 1 row and one col to the right
            (*mtrx)[row+1][col+1] = (*mtrx)[row][col];
        }
    }
    double vtr_zero = get_backward_element(vtr, start_slice, 0);
    // define first column and row
    for (long row = MATRIX_SIZE - 2; row >= 0; --row)
    {
        (*mtrx)[row+1][0] = scalar * first_col[row+1] + vtr_zero * get_backward_element(vtr,start_slice,row+1);
        // make matrix symmetric
        (*mtrx)[0][row+1] = (*mtrx)[row+1][0];
    }

    (*mtrx)[0][0] = scalar * first_col[0] + vtr_zero * vtr_zero;
}


void update_matrix(REGULAR_MATRIX *mtrx, const double *vtr, const double scalar)
{
    //format is  mtrx[row][col]

    // store first column (index 0) for later use
    double first_col[MATRIX_SIZE];
    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        first_col[row] = (*mtrx)[row][0];
    }

    // loop backwards from MATRIX_SIZE - 2 to 0 in increments of -1
    for (long row = MATRIX_SIZE - 2; row >= 0; --row)
    {
        for (long col = MATRIX_SIZE - 2; col >= 0; --col)
        {
        // shift elements down 1 row and one col to the right
            (*mtrx)[row+1][col+1] = (*mtrx)[row][col];
        }

    }
    // define first column and row
    for (long row = MATRIX_SIZE - 2; row >= 0; --row)
    {
        (*mtrx)[row+1][0] = scalar * first_col[row+1] + vtr[0] * vtr[row+1];
        // make matrix symmetric
        (*mtrx)[0][row+1] = (*mtrx)[row+1][0];
    }

    (*mtrx)[0][0] = scalar * first_col[0] + vtr[0] * vtr[0];
}

void update_vector(double *vtr, const REGULAR_MATRIX *mtrx, const double scalar, const long col_num)
{
    // format is mtrx[row][col]
    for (long row = 0; row < MATRIX_SIZE; ++row)
    {
        vtr[row] = vtr[row] - scalar * (*mtrx)[row][col_num];
    }
}

double gaussian_rand(const gsl_rng *r, double sigma)
{
    // generate gaussian random number with mean zero and standard deviation sigma
    // r is a pointer to the random number generator
    // calling routine must set the seed with srand()

    double rand_num = gsl_ran_gaussian(r, sigma);
    return rand_num;
}

double *gaussian_rand_vector(const gsl_rng *r, double sigma, long size)
{

    // create a vector of size "size" with gaussian random numbers
    double *vector_pt = (double *)malloc(sizeof(double) * size);
    for (long element = 0; element < size; ++element)
    {
        vector_pt[element] = gsl_ran_gaussian(r, sigma);
    }
    return vector_pt;
    // don't forget to release memory in calling routine
}

long min(long a, long b)
{
    return (a<b ? a : b);
}