//
// Created by Naomi Sutcliffe de Moraes
//


#include "matrix_operations/Standard_operations.h"
#include "matrix_operations/DCBF_operations.h"
#include "standard_RLSDCD/standard_RLSDCD.h"
#include "DCBF_RLSDCD/DCBF_RLSDCD.h"

void filt(double *y, const double *h, const double *xx, long num_time_steps, long M)
{
    // h is an array of length M
    // x and y are  arrays of length num_time_steps
    for (long n = 0; n < num_time_steps; ++n)
    {
        y[n] = h[0] * xx[M - 1 + n ];
        for (long k = 1; k < M; k++)
        {
            y[n] = y[n] + h[k] * xx[M - 1 + n - k];
        }

    }
}


void RLS(double *w, double *e,  double *x, const double *s, long M, double lambda, double delta, long num_time_steps)
{
    // this function modifies two vectors, w, e
    // initialize P to delta * Identity matrix
    REGULAR_MATRIX *out_product = malloc(sizeof(REGULAR_MATRIX));
    //REGULAR_MATRIX  out_product;
    zero_regular_matrix(out_product);


    REGULAR_MATRIX *P = initialize_matrix(delta);


    double lambda_i = 1.0 / lambda;
    for (long step = 0; step < num_time_steps; ++step)
    {

        // start and end of slice are inclusive in Julia
        long start_slice = step + M - 1;
        double y = backward_dot_product(x, w, start_slice);

        e[step] = s[step] - y;

        double *g = malloc(sizeof(double) * MATRIX_SIZE);

        matrix_x_backward_vector(P, x, step + M - 1, g);


        double gamma = lambda + backward_dot_product(x, g, step + M - 1);

        double gamma_inv = 1.0 / gamma;

        for (long j = 0; j < M; ++j)
        {
            w[j] = w[j] + (e[step] * gamma_inv) * g[j];
        }


        for (long k = 0; k < M; ++k)
        {
            g[k] = g[k] * sqrt(fabs(gamma_inv));
        }

        outer_product(g, g, out_product);
        //printf(" \n outer product\n");
        //print_regular_matrix(out_product);


        for (long row = 0; row < MATRIX_SIZE; ++row)
        {
            for (long col = 0; col < MATRIX_SIZE; ++col)
            {
                (*P)[row][col] = lambda_i * (   (*P)[row][col] - (*out_product)[row][col]  );

            }
        } // end for row

    } // end for step


    free(P);
    free(out_product);

} // end RLS

void test_DCBF(const gsl_rng *r)
{
    printf("inside test_DCBF \n");
    FILE *file = fopen("C_RLSDCD.txt", "w");

    double sigma_v = 0.1;
    double sigma_one = 1.0;
    double lambda = 0.99;
    double delta = 1.0;

    long num_simulations = 100; // number of simulations to run, then average results - normally 100
    long num_time_steps = 5000;  // normally 5000
    long num_dcd_iterations;

    // create arrays to store running times for different simulations
    // array[row][col] = array[different M values][number of simulations run, average of them is taken later
    double *TRLS = malloc(sizeof (double) * num_simulations );
    double *TRLSDCD = malloc(sizeof (double) * num_simulations );
    double *TRLSDCD4 = malloc(sizeof (double) * num_simulations );
    double *TfRLSDCD = malloc(sizeof (double) * num_simulations );
    double *TfRLSDCD4 = malloc(sizeof (double) * num_simulations );


    double *hi = gaussian_rand_vector(r, sigma_one, MATRIX_SIZE);
    double *s = malloc(sizeof(double) * num_time_steps);

    for (long i = 0; i < num_time_steps; ++i)
    {
        s[i] += gaussian_rand(r, sigma_v);
    }

    long longer_len = MATRIX_SIZE + num_time_steps - 1;  // length of longer_x vector
    double *longer_x = malloc(sizeof(double) * longer_len);

    // the vector longer_x has MATRIX_SIZE - 1 zeros, followed by (num_time_steps) random numbers
    for (long elem = 0; elem < MATRIX_SIZE - 1; ++elem)
    {
        longer_x[elem] = 0.0;
    }

    for (long elem = MATRIX_SIZE-1; elem < longer_len; ++elem)
    {
        longer_x[elem] = gaussian_rand(r, sigma_one);
    }
    filt(s, hi, longer_x, num_time_steps, MATRIX_SIZE);                  // s is the return value

    double *e = malloc(sizeof(double) * num_time_steps);
    printf( " Matrix size = %d \n", MATRIX_SIZE);
    // loop over number of simulations to average for a given value of MATRIX_SIZE

    for (long i = 0; i < num_simulations; ++i)
    {

        double w[MATRIX_SIZE] = {0};  // need to zero each simulation because it is updated inside RLS


        clock_t start, end;

        // RLS time test
        start = clock();
        long M = (long)MATRIX_SIZE;

        RLS(w, e, longer_x, s, M, lambda, delta, num_time_steps);
        end = clock();
        TRLS[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
        //printf(  "\n i: %ld, TRLS %lf ", i, TRLS[i]  );

       // RLSDCD time test
        num_dcd_iterations = 1;
        start = clock();
        RLSDCD(w, e, longer_x, s, lambda, delta, num_dcd_iterations, num_time_steps);
        end = clock();
        TRLSDCD[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
        //printf(  "\n i: %ld, TRLSDCD %lf ", i, TRLSDCD[i]  );


        // RLSDCD 4 time test
        num_dcd_iterations = 4;
        start = clock();
        RLSDCD(w, e, longer_x, s, lambda, delta, num_dcd_iterations, num_time_steps);
        end = clock();
        TRLSDCD4[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
        //printf(  "\n i: %ld, TRLSDCD4 %lf ", i, TRLSDCD4[i]  );


        // fast RLSDCD time test
        num_dcd_iterations = 1;
        start = clock();
        fRLSDCD(w, e, longer_x, s, lambda, delta, num_dcd_iterations, num_time_steps);
        end = clock();
        TfRLSDCD[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
        //printf(  "\n i: %ld, TfRLSDCD4 %lf ", i, TfRLSDCD4[i]  );

        // fast RLSDCD4 time test
        num_dcd_iterations = 4;
        start = clock();
        fRLSDCD(w, e, longer_x, s, lambda, delta, num_dcd_iterations, num_time_steps);
        end = clock();
        TfRLSDCD4[i] = ((double) (end - start)) / CLOCKS_PER_SEC;
        //printf(  "\n i: %ld, TfRLSDCD4 %lf ", i, TfRLSDCD4[i]  );
    } // end for num_simulations

    double MeRLS = mean_of_vector(TRLS, num_simulations);
    double MinRLS = minimum_of_vector(TRLS, num_simulations);
    printf("\n Mean RLS time: %lf ", MeRLS);
    printf("\n Minimum RLS time: %lf", MinRLS);

    double MeRLSDCD = mean_of_vector(TRLSDCD, num_simulations);
    double MinRLSDCD = minimum_of_vector(TRLSDCD, num_simulations);
    printf("\n Mean RLSDCD time: %lf ", MeRLSDCD);
    printf("\n Minimum RLSDCD time: %lf", MinRLSDCD);

    double MeRLSDCD4 = mean_of_vector(TRLSDCD4, num_simulations);
    double MinRLSDCD4 = minimum_of_vector(TRLSDCD4, num_simulations);
    printf("\n Mean RLSDCD4 time: %lf ", MeRLSDCD4);
    printf("\n Minimum RLSDCD4 time: %lf", MinRLSDCD4);

    double MefRLSDCD = mean_of_vector(TfRLSDCD, num_simulations);
    double MinfRLSDCD = minimum_of_vector(TfRLSDCD, num_simulations);
    printf("\n Mean fRLSDCD time: %lf ", MefRLSDCD);
    printf("\n Minimum RfLSDCD time: %lf", MinfRLSDCD);

    double MefRLSDCD4 = mean_of_vector(TfRLSDCD4, num_simulations);
    double MinfRLSDCD4 = minimum_of_vector(TfRLSDCD4, num_simulations);
    printf("\n Mean fRLSDCD4 time: %lf ", MefRLSDCD4);
    printf("\n Minimum fRLSDCD4 time: %lf", MinfRLSDCD4);


    if (file != NULL)
    {
        fprintf(file, "%lf %lf %ld %ld\n", MeRLS, MinRLS, (long)MATRIX_SIZE, num_time_steps);
    }


    free(s);
    free(longer_x);
    free(e);


    free(TRLS);
    free(TRLSDCD);
    free(TRLSDCD4);
    free(TfRLSDCD);
    free(TfRLSDCD4);


}  // end test_DCBF


void test_fRLSDCD(long num_time_steps, const gsl_rng *r)
{
    double sigma_v = 0.01;
    double sigma_one = 1.0;
    double lambda = 0.999;
    double delta = 1.0;
    long num_dcd_iterations = 10;

    double *hi = gaussian_rand_vector(r, sigma_one, MATRIX_SIZE);
    double w[MATRIX_SIZE] = {0};
    double *s = malloc(sizeof(double) * num_time_steps);
    double *e = malloc(sizeof(double) * num_time_steps);

    long longer_len = MATRIX_SIZE + num_time_steps - 1;  // length of longer_x vector
    double *longer_x = malloc(sizeof(double) * longer_len);
    // the vector longer_x has MATRIX_SIZE-1 zeros, followed by (num_time_steps) random numbers
    for (long elem = 0; elem < MATRIX_SIZE - 1; ++elem)
    {
        longer_x[elem] = 0.0;
    }

    for (long elem = MATRIX_SIZE-1; elem < longer_len; ++elem)
    {
        longer_x[elem] = gaussian_rand(r, sigma_one);
    }
    filt(s, hi, longer_x, num_time_steps, MATRIX_SIZE);        // s is the return value
    for (long i=0; i<num_time_steps; i++)
    {
        s[i] += gaussian_rand(r, sigma_v);
    }
    fRLSDCD(w, e, longer_x, s, lambda, delta, num_dcd_iterations, num_time_steps);
/*    for (long i=0; i<MATRIX_SIZE; i++)
    {
        printf("hi[%ld] = %lf, w[%ld] = %lf \n", i, hi[i], i, w[i]);
    }*/
    double MSD = pow(Euclidean_distance(w, hi, MATRIX_SIZE), 2);
    printf("\n fRLSDCD MSD = %lg ", MSD);


    free(hi);
    free(longer_x);
    free(s);
    free(e);

} // end test_fRLSDCD

void test_RLSDCD(long num_time_steps, const gsl_rng *r)
{
    double sigma_v = 0.01;
    double sigma_one = 1.0;
    double lambda = 0.999;
    double delta = 1.0;
    long num_dcd_iterations = 10;

    double *hi = gaussian_rand_vector(r, sigma_one, MATRIX_SIZE);
    double w[MATRIX_SIZE] = {0};
    double *s = malloc(sizeof(double) * num_time_steps);
    double *e = malloc(sizeof(double) * num_time_steps);

    long longer_len = MATRIX_SIZE + num_time_steps - 1;  // length of longer_x vector
    double *longer_x = malloc(sizeof(double) * longer_len);
    // the vector longer_x has MATRIX_SIZE-1 zeros, followed by (num_time_steps) random numbers
    for (long elem = 0; elem < MATRIX_SIZE - 1; ++elem)
    {
        longer_x[elem] = 0.0;
    }

    for (long elem = MATRIX_SIZE-1; elem < longer_len; ++elem)
    {
        longer_x[elem] = gaussian_rand(r, sigma_one);
    }
    filt(s, hi, longer_x, num_time_steps, MATRIX_SIZE);        // s is the return value
    for (long i=0; i<num_time_steps; i++)
    {
        s[i] += gaussian_rand(r, sigma_v);
    }
    RLSDCD(w, e, longer_x, s, lambda, delta, num_dcd_iterations, num_time_steps);
/*    for (long i=0; i<MATRIX_SIZE; i++)
    {
        printf("hi[%ld] = %lf, w[%ld] = %lf \n", i, hi[i], i, w[i]);
    }*/
    double MSD = pow(Euclidean_distance(w, hi, MATRIX_SIZE), 2);
    printf("\n RLSDCD MSD = %lg ", MSD);


    free(hi);
    free(longer_x);
    free(s);
    free(e);

} // end test_RLSDCD

void test_RLS(long num_time_steps, const gsl_rng *r)
{
    double sigma_v = 0.01;
    double sigma_one = 1.0;
    double lambda = 0.999;
    double delta = 1.0;

    double *hi = gaussian_rand_vector(r, 1.0, MATRIX_SIZE);
    double w[MATRIX_SIZE] = {0};
    double *s = malloc(sizeof(double) * num_time_steps);
    double *e = malloc(sizeof(double) * num_time_steps);

    long longer_len = MATRIX_SIZE + num_time_steps - 1;  // length of longer_x vector
    double *longer_x = malloc(sizeof(double) * longer_len);
    // the vector longer_x has MATRIX_SIZE-1 zeros, followed by (num_time_steps) random numbers
    for (long elem = 0; elem < MATRIX_SIZE - 1; ++elem)
    {
        longer_x[elem] = 0.0;
    }

    for (long elem = MATRIX_SIZE-1; elem < longer_len; ++elem)
    {
        longer_x[elem] = gaussian_rand(r, sigma_one);
    }
    filt(s, hi, longer_x, num_time_steps, MATRIX_SIZE);        // s is the return value
    for (long i=0; i < num_time_steps; i++)
    {
        s[i] += gaussian_rand(r, sigma_v);
    }
    RLS(w, e, longer_x, s, MATRIX_SIZE, lambda, delta, num_time_steps);

/*    for (long i=0; i<MATRIX_SIZE; i++)
    {
        printf("hi[%ld] = %lf, w[%ld] = %lf \n", i, hi[i], i, w[i]);
    }*/
    double MSD = pow(Euclidean_distance(w, hi, MATRIX_SIZE), 2);
    printf("\n RLS MSD = %lg ", MSD);


    free(hi);
    free(longer_x);
    free(s);
    free(e);

} // end test_RLS