//
// Created by Naomi Sutcliffe de Moraes
//
#include "DCBF_operations.h"

long julia_mod1(long dividend, long divisor) {
    long result = (dividend % divisor);
    if (result < 0) // If the result is negative, add divisor to get into the range [0, divisor - 1]
        result += divisor;
    return result;
}
double * pointer_to_head(const long diagonal, const struct DCBF_Matrix *mtrx)
{
    // given the DCBF dmtrx and a diagonal index, return a pointer to the head of the circular buffer of that diagonal
    // pointer to head = pointer to top + offset to head

    double * pt_to_head = mtrx->top[diagonal] + mtrx->cb_head[diagonal];
    return pt_to_head;
}

struct DCBF_Matrix *initialize_DCBF_Matrix(double delta)
{
    // Allocate memory for the structure
    struct DCBF_Matrix *initialized = (struct DCBF_Matrix *)malloc(sizeof(struct DCBF_Matrix));

    // Initialize longest diagonal to all ones, rest of dmtrx to zeroes
    for (long element = 0; element < MATRIX_SIZE; ++element)  // first diagonal
    {
        initialized->dmtrx[element] = 1.0 * delta;
    }
    for (long element = MATRIX_SIZE; element < ARRAY_SIZE; ++element)  // remaining diagonals
    {
        initialized->dmtrx[element] = 0.0;
    }
    // initialize cb_head vector to be all zeroes (we start with circular buffers aligned at the top of the diagonal)
    for (long diagonal = 0; diagonal < MATRIX_SIZE; ++diagonal)
    {
        initialized->cb_head[diagonal] = 0;
    }


    long diagonal_size = MATRIX_SIZE;
    long index = 0;
    for (long diagonal = 0; diagonal < MATRIX_SIZE; ++diagonal)
    {
        // initialize first diagonal to have MATRIX_SIZE elements
        // the second diagonal will have MATRIX_SIZE - 1 elements
        // and so on until the last diagonal has just 1 element
        initialized->diag_num_elements[diagonal] = diagonal_size;

        // Initialize the indices to the head elements of the circular buffer for each diagonal, 0..MATRIX_SIZE - 1
        //initialized.cb_head[diagonal] = index;

        // find addresses corresponding to first element of each diagonal
        initialized->top[diagonal] = &initialized->dmtrx[index];

        index = index + diagonal_size;
        diagonal_size-- ;

    }
    return initialized;
    // don't forget to free the memory in the calling program
}



void initialize_DCBF_Matrix_with_RegularMatrix(REGULAR_MATRIX *regular_matrix, struct DCBF_Matrix *initialized )
{
    long reg_row = 0;
    long reg_col = 0;
    long array_element = 0;
    long diagonal_size = MATRIX_SIZE; // decreases each loop

    while (array_element < ARRAY_SIZE) {
        for (long element = 0; element < diagonal_size; ++element)  // loop for each diagonal
        {
            initialized->dmtrx[array_element] = (*regular_matrix)[reg_row][reg_col];
            reg_col++;
            reg_row++;
            array_element++;
        }
        // reset counters to start at next diagonal
        diagonal_size--;
        reg_col = reg_col - diagonal_size;
        reg_row = 0;
    }
}



void print_DCBF_Matrix(struct DCBF_Matrix *mtrx)
{
    printf("\nthe mtrx *******************************\n");
    long matrix_size = mtrx->diag_num_elements[0]; // Number of elements in largest [0] diagonal
    long diagonal_size = matrix_size; //Number of elements in diagonal decreases each loop
    long element = 0; // element steps through the entire 1-D array representing the DCBF dmtrx
    for (long diagonal = 0; diagonal < matrix_size; ++diagonal)
    {
        printf("diagonal %ld: ", diagonal);
        for (long counter = 0; counter < diagonal_size; ++counter) // counter resets each diagonal
        {
            printf("%f  ", mtrx->dmtrx[element]);
            element++;
        }
        diagonal_size--;
        printf("\n"); // return at the end of each diagonal
    }
    printf("\nHead of each diagonal (index with respect to TOP pointer), mtrx->cb_head: ");
    printf("\nshould be zero or a positive number less than the diagonal size\n");
    for (long diagonal = 0;  diagonal < matrix_size; ++diagonal)
    {
        printf("%ld ", mtrx->cb_head[diagonal]);
    }
    printf("\n\n");

    printf("Top of each diagonal (value in 1-D array): ");
    for (long diagonal = 0;  diagonal < matrix_size; ++diagonal)
    {
        printf("%f ", *mtrx->top[diagonal]);
    }
    printf("\n");
    printf("Number of elements in each diagonal: ");
    for (long diagonal = 0;  diagonal < matrix_size; ++diagonal)
    {
        printf("%ld ", mtrx->diag_num_elements[diagonal]);
    }
    printf("\n");

    for (long diagonal = 0; diagonal < matrix_size; ++diagonal) {
        double *head_diag = pointer_to_head(diagonal, mtrx);
        printf("value of head of diagonal %ld = %f \n", diagonal, *head_diag);
    }
    printf("\nEnd of the mtrx data ***********************\n");
}  // end print_DCBF_Matrix

void print_vector(double vtr[])
{
    printf("Vector: \n");
    for (long element = 0; element < MATRIX_SIZE; ++element)
    {
        printf("%.2f ", vtr[element]);
    }
    printf("\n");
}

// ****************************************************************
// update DIAGONAL CIRCULAR-BUFFER FORMAT MATRIX by shifting down
// and to the right by one column and one row
// (used in fRLSDCD function)
// ****************************************************************
// input: DCBF matrix, vector and scalar for update equation
// output: updated DCBF matrix
// ****************************************************************


void update_DCBF_matrix_w_backward_vector(struct DCBF_Matrix *mtrx, const double *vtr, long start_slice, const double scalar)
{
    for(long i = 0; i < MATRIX_SIZE; i++)  // i is diagonal
    {
        // value of element at head of circular buffer before updating
        double old_head = *pointer_to_head(i, mtrx);


        // shift diagonal up by moving pointer within diagonal (mod is to make it wrap around)
        // Formula is different than in Julia code because indices start at zero in C, but 1 in Julia
        mtrx->cb_head[i] = (mtrx->cb_head[i] - 1) % mtrx->diag_num_elements[i];
        if (mtrx->cb_head[i] < 0)
        {
            mtrx->cb_head[i] = mtrx->diag_num_elements[i] + mtrx->cb_head[i];
        }

        // update a single element (store new value in element representing the head of the diagonal)
        double *pt_to_new_head = pointer_to_head(i, mtrx);
        *pt_to_new_head = scalar * old_head + get_backward_element(vtr, start_slice, 0) * get_backward_element(vtr, start_slice, i);

    }
}

void update_DCBF_Matrix(struct DCBF_Matrix *mtrx, const double *vtr, const double scalar)
{
    for(long i = 0; i < MATRIX_SIZE; i++)  // i is diagonal
    {
        // value of element at head of circular buffer before updating
        double old_head = *pointer_to_head(i, mtrx);


        // shift diagonal up by moving pointer within diagonal (mod is to make it wrap around)
        // Formula is different than in Julia code because indices start at zero in C, but 1 in Julia
        mtrx->cb_head[i] = (mtrx->cb_head[i] - 1) % mtrx->diag_num_elements[i];
        if (mtrx->cb_head[i] < 0)
        {
            mtrx->cb_head[i] = mtrx->diag_num_elements[i] + mtrx->cb_head[i];
        }

        // update a single element (store new value in element representing the head of the diagonal)
        double *pt_to_new_head = pointer_to_head(i, mtrx);
        *pt_to_new_head = scalar * old_head + vtr[0] * vtr[i];

    }
}

// *********************************************************************
// multiply column of the DIAGONAL CIRCULAR-BUFFER FORMAT MATRIX by a
// factor and add result to a vector (used in fDCD function)
// *********************************************************************
// input: vector, matrix, factor, column number
// output: vector <- vector + factor * matrix[:,column number]
// ********************************************************************
void update_DCBF_Vector(double *vtr, const struct DCBF_Matrix *mtrx, const double scalar, const long col_numb)
{
    long diagonal = col_numb;  // diagonal starts out as column in the original standard matrix format
    long diagonal_step = -1;   // the algorithm steps through the columns until it reaches the first ([0],
                               // then changes direction
    long cb_offset = 0;        // cb_offset is the offset from the head element
                               // inside the circular buffer for that diagonal
    long kstep = 1;

    // iterate over elements in the vector
    for(long i = 0; i < MATRIX_SIZE; i++)  // i is element in vector
    {
        // given the diagonal, select the element in the circular buffer (row)
        long cb_element = julia_mod1((mtrx->cb_head[diagonal] + cb_offset), mtrx->diag_num_elements[diagonal]);

        // update the vector
        // vtr[i] = vtr[i] - scalar * mtrx->dmtrx[diagonal][cb_element]; is equivalent in regular matrix

        vtr[i] = vtr[i] - scalar * (*(mtrx->top[diagonal] + cb_element));

        // select the next diagonal
        diagonal = diagonal + diagonal_step;
        cb_offset = cb_offset + kstep;
        if (diagonal == -1)  // changed from 0 to -1 because Julia arrays start at 1 and C arrays start at 0
        {
            diagonal = 1;  // changed from 2 to 1 because Julia arrays start at 1 and C arrays start at 0
            diagonal_step = 1;
            cb_offset = col_numb;  // changed from 2 to 1 because Julia arrays start at 1 and C arrays start at 0
            kstep = 0;
        }  // end if

    } // end for


}