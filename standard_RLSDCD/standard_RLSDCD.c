//
// Created by Naomi J S de Moraes
//

#include "standard_RLSDCD.h"
void DCD(double H, long B, long num_dcd_iterations, REGULAR_MATRIX *R_phi, double *res, double *delta_w)
{
    // function modifies res and delta_w
    double h = H / 2.0;
    long b = 1;
    zero_vector(delta_w, MATRIX_SIZE);
    for (long k = 0; k < num_dcd_iterations; ++k)
    {
        long p = argmax(res, MATRIX_SIZE);  // a value from 0 to M-1
        while (   fabs(res[p]) <= (  (h/2.0) * (*R_phi)[ p][p ]  )  ) {
            b++;
            h = h / 2.0;
            if (b > B) { break; }
        } // end while

        double delta_wp =  ( (sign(res[p])) * h );
        delta_w[p] = delta_w[p] + delta_wp;

        update_vector(res, R_phi, delta_wp, p);

    } // end for
}  // end fn DCD


void RLSDCD(double *w, double *e, double*x, const double *s, double lambda, double delta, long num_dcd_iterations, long num_time_steps)
{
    zero_vector(w, MATRIX_SIZE);
    double beta[MATRIX_SIZE] = {0};
    double delta_w[MATRIX_SIZE] = {0};

    REGULAR_MATRIX *R_phi = initialize_matrix(delta);
    zero_vector(e, num_time_steps);

    for (long step = 0; step < num_time_steps; ++step)
    {
        // note: start and end of slice are inclusive in Julia, so must take this into account
        long start_slice = step + MATRIX_SIZE - 1;
        double y = backward_dot_product(x, w, start_slice);
        e[step] = s[step] - y;

        update_matrix_w_backward_vector(R_phi, x, step + MATRIX_SIZE - 1, lambda);

        for (long element = 0; element < MATRIX_SIZE; ++element)
        {
            beta[element] = lambda * beta[element] + e[step] * get_backward_element(x,step + MATRIX_SIZE - 1, element);
        }
        DCD(4.0, 16, num_dcd_iterations, R_phi, beta, delta_w);

        for (long i = 0; i < MATRIX_SIZE; ++i)
        {
            w[i] = w[i] + delta_w[i];
        }


    } // end for step
    free(R_phi);

} // end RLSDCD