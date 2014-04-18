/*
 * gravsim - simulates masses acting under their gravity
 *  Copyright (C) 2014 Michael J. Beer <michael@ubeer.org>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*---------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/*---------------------------------------------------------------------------*/
/* Solving the equations occuring by implicit GAUSS-LEGENDRE integration */
#include <gsl/gsl_multiroots.h>
gsl_multiroot_fsolver *multiRootSolver;
size_t iteration;
/*---------------------------------------------------------------------------*/
/**
 *
 */
struct exponential_params {
    double dt;
    const gsl_vector * y_last;
    const gsl_vector * a;
};
/*---------------------------------------------------------------------------*/
int exponential(struct exponential_params * params, gsl_vector * new_value, const gsl_vector * k)
{
   fprintf(stderr, "\n\n\nexponential: iteration %d \nk : ", iteration);
   iteration++;
   gsl_vector_fprintf(stderr, k, "%f");
   fprintf(stderr, "exponential: y_last : ");
   gsl_vector_fprintf(stderr, params->y_last, "%f");
   gsl_vector_memcpy(new_value, k);
   gsl_vector_scale(new_value,  0.5);
   gsl_vector_add(new_value, params->y_last);
   gsl_vector_scale(new_value, params->dt);
   gsl_vector_mul(new_value, params->a);
   fprintf(stderr, "exponential: new_value : ");
   gsl_vector_fprintf(stderr, new_value, "%f");
   return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------*/
int exponential_zero_function(const gsl_vector * x, void * params, gsl_vector * f)
{
    struct exponential_params *exp_params = (struct exponential_params *) params;
    exponential(params, f, x);
   fprintf(stderr, "exponential_zero_function: f : ");
   gsl_vector_fprintf(stderr, f, "%f");
    gsl_vector_sub(f, x);
   fprintf(stderr, "exponential_zero_function: f : ");
   gsl_vector_fprintf(stderr, f, "%f");
    return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------*/
int integrate_exponential(double start_t, double dt, double end_t,const gsl_vector * y_0, const gsl_vector * a)
{
    printf("%f\n", start_t);
    /* For Gauss-Legendre 2nd order */
    double time = start_t;
    double max_err = 0.0000000000000001;
    size_t size_dim = (size_t) y_0->size;
    gsl_vector * init_val = gsl_vector_alloc(size_dim);
    size_t n;
    for(n = 0; n < size_dim; n++) 
    { 
        gsl_vector_set(init_val, n, 1.0);
    }
    gsl_vector * y_implicit = gsl_vector_alloc(size_dim);
    gsl_vector_memcpy(y_implicit, y_0);
    struct exponential_params params;
    params.dt = dt;
    params.y_last = y_implicit;
    params.a = a;
    multiRootSolver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 2);
    gsl_multiroot_function zero_func;
    zero_func.f = &exponential_zero_function;
    zero_func.n = 2;
    zero_func.params = &params;
    gsl_multiroot_fsolver_set(multiRootSolver, &zero_func, init_val);
    /* For common Euler method */
    gsl_vector * y = gsl_vector_alloc(y_0->size);
    gsl_vector_memcpy(y, y_0);
    gsl_vector * temp = gsl_vector_alloc(y_0->size);
    /* And for Gauss-Legendre 'by-hand' */
    gsl_vector * y_by_hand = gsl_vector_alloc(y_0->size);
    gsl_vector * factor_by_hand = gsl_vector_alloc(y_0->size);
    gsl_vector_memcpy(y_by_hand, y_0);
    /*---*/
    printf("%f\n", time);
    while(time < end_t)
    {
        iteration = 0;
        gsl_multiroot_fsolver_set(multiRootSolver, &zero_func, init_val);
        gsl_multiroot_fsolver_iterate (multiRootSolver);
        int test_residual_result;
        while((test_residual_result = gsl_multiroot_test_residual(gsl_multiroot_fsolver_f(multiRootSolver), max_err)) == GSL_CONTINUE)
        {
            int result = gsl_multiroot_fsolver_iterate (multiRootSolver);
            if(result != GSL_SUCCESS)
            {
                fprintf(stderr, "Iteration aborted due to %i\n", gsl_strerror(result));
                exit(1);
            }
            fprintf(stderr, "%d %s %d %s\n", result, gsl_strerror(result), test_residual_result, gsl_strerror(test_residual_result));
        }
        if(test_residual_result != GSL_SUCCESS)
        {
            fprintf(stderr, "Iteration aborted due to %i\n", gsl_strerror(test_residual_result));
            exit(1);
        }
        gsl_vector * k = gsl_multiroot_fsolver_root(multiRootSolver);
        /* Use found k to calculate actual new value */
        gsl_vector_add(y_implicit, k);
        gsl_vector_memcpy(init_val, k);
        /* Now integrate using default Euler method */
        gsl_vector_memcpy(temp, y);
        gsl_vector_mul(temp, a);
        gsl_vector_scale(temp, dt);
        gsl_vector_add(y, temp);
        /* And Gauss-Legendre 'by hand */
        for(n = 0; n < y_by_hand->size; n++)
        {
            double a_value = gsl_vector_get(a, n);
            double factor  = 2.0 * dt * a_value / (2.0  - dt * a_value) + 1.0;
            gsl_vector_set(factor_by_hand, n, factor);
        }
        gsl_vector_mul(y_by_hand, factor_by_hand);
        /* --- */
        time += dt;
        printf("%f ", time);
        for( n = 0; n < size_dim; n++)
        { 
            
            printf("%f ", gsl_vector_get(y_implicit, n));
            printf("%f ", gsl_vector_get(y, n));
            printf("%f ", gsl_vector_get(y_by_hand, n));
        }
        printf("\n");
    }
    gsl_vector_free(init_val);
    gsl_vector_free(y_implicit);
    gsl_vector_free(y);
    gsl_vector_free(temp);
    gsl_vector_free(y_by_hand);
    gsl_vector_free(factor_by_hand);
    gsl_multiroot_fsolver_free(multiRootSolver);
    return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------
 * MAIN - do parameter parsing etc...
 *---------------------------------------------------------------------------*/
int main(int argc, char** argv) 
{
    gsl_vector * y_0 = gsl_vector_alloc(2);
    gsl_vector * a   = gsl_vector_alloc(2);
    gsl_vector_set(y_0, 0, 1.0);
    gsl_vector_set(y_0, 1, 1.0);
    gsl_vector_set(a  , 0, -0.1);
    gsl_vector_set(a  , 1, -2.0);
    integrate_exponential(0.0, 0.1, 50.0, y_0, a);
    gsl_vector_free(y_0);
    gsl_vector_free(a);
}

