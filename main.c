//
// Created by Naomi Sutcliffe de Moraes
//
#include <stdio.h>
#include "echo_canceling.h"


int main() {

    // INSTRUCTIONS
    // The matrix size is hard-coded in the file standard_operations.h
    // Set the value there first, before running the test(s)

    printf("Tests running\n");

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
    // to test the RLS algorithm with 40,000 time steps, uncomment the following line
    //test_RLS(40000, r);

    // to test the fast RLS-DCD algorithm (using the DCBF) with 40,000 time steps, uncomment the following line
    //test_fRLSDCD(40000, r);

    // to test the RLS-DCD algorithms (using the standard matrix format) with 40,000 time steps, uncomment the following line
    //test_RLSDCD(40000, r);

    // to test all the algorithms with the same parameters
    // and print mean and minimum times, uncomment the following line
    test_DCBF(r);


    gsl_rng_free(r);
    return 0;
}
