#include <stdio.h>
//#include <stdlib.h>
//#include "matrix_operations/Standard_operations.h"
//#include "matrix_operations/DCBF_operations.h"
#include "operations_unit_tests.h"
#include "echo_canceling.h"


int main() {
    printf("Start of program!\n");
    // create 5x5 dmtrx for test
/*    REGULAR_MATRIX *regular_matrix_pt = initialize_5x5();
    printf("\n printing 5x5 regular dmtrx for 5x5 test\n");
    print_regular_matrix(regular_matrix_pt);
    test_DCBF_matrix(regular_matrix_pt);

    //test_initialization();
    test_matrix(regular_matrix_pt);

    // free up variables
    free(regular_matrix_pt);

    // create random vector
    double *vect = gaussian_rand_vector(MATRIX_SIZE);
    printf("\n random vector test \n");
    print_array(vect, MATRIX_SIZE);
    free(vect);*/

    //test_DCBF();
    // test_get_backward_element();


    // RANDOM NUMBER GENERATOR PREPARATION
    const gsl_rng_type * T;
    gsl_rng * r;
    // Initialize the random number generator
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Set the seed for the random number generator
    gsl_rng_set(r, 123654);

    // END OF RANDOM NUMBER GENERATOR PREPARATION

    //test_RLS(40000, r);
    //test_fRLSDCD(40000, r);
    //test_RLSDCD(40000, r);
    test_DCBF(r);


    gsl_rng_free(r);
    return 0;
}
