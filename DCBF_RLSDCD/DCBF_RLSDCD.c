//
// Created by Naomi Nascimento on 3/31/24.
//

#include "DCBF_RLSDCD.h"

void fDCD(double H, long B, long num_dcd_iterations, struct DCBF_Matrix *R_phi, double *res, double *delta_w)
{
    // function modifies res and delta_w
    double h = H / 2.0;
    long b = 1;
    zero_vector(delta_w, MATRIX_SIZE);

    for (long k = 0; k < num_dcd_iterations; ++k)
    {
        // p is the index corresponding to the maximum residual. I want to update the pth column of Rphi
        long p = argmax(res, MATRIX_SIZE);  // a value from 0 to M-1
        long p_cb_position = julia_mod1( (R_phi->cb_head[0] + p), R_phi->diag_num_elements[0]);

        // R_phi->dmtrx[p_cb_position] is the cb position for the main diagonal
        while (    fabs(res[p]) <= (   (h / 2.0) * R_phi->dmtrx[p_cb_position]    ))
        {
            b++;
            h = h / 2.0;
            if (b > B) { break; }
        } // end while
        double delta_wp =  ( sign(res[p])) * h ;
        delta_w[p] = delta_w[p] + delta_wp;

        update_DCBF_Vector(res, R_phi, delta_wp, p);

    }  // end loop over num_dcd_iterations


}  // end fDCD

void fRLSDCD(double *w, double *e, double*x, const double *s, double lambda, double delta, long num_dcd_iterations, long num_time_steps)
{
    zero_vector(w, MATRIX_SIZE);
    double beta[MATRIX_SIZE] = {0};
    double delta_w[MATRIX_SIZE] = {0};

    struct DCBF_Matrix *R_phi = initialize_DCBF_Matrix(delta);
    zero_vector(e, num_time_steps);

    for (long step = 0; step < num_time_steps; ++step)
    {
        // start and end of slice are inclusive in Julia
        long start_slice = step + MATRIX_SIZE - 1;
        double y = backward_dot_product(x, w, start_slice);
        e[step] = s[step] - y;

        update_DCBF_matrix_w_backward_vector(R_phi, x, step + MATRIX_SIZE - 1, lambda);

        // TO DO - check if indices are correct
        for (long element = 0; element < MATRIX_SIZE; ++element)
        {
            beta[element] = lambda * beta[element] + e[step] * get_backward_element(x, step + MATRIX_SIZE - 1, element);
        }
        fDCD(4.0, 16, num_dcd_iterations, R_phi, beta, delta_w);

        for (long i = 0; i < MATRIX_SIZE; ++i)
        {
            w[i] = w[i] + delta_w[i];
        }
    } // end for
    free(R_phi);

}  // end RLSDCD