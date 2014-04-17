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
/*---------------------------------------------------------------------------*/
/* All masses in 10^6 kg */
/* All distances in 10^6 m */
#define GRAVITATIONAL_CONSTANT_ORIGIN ((double)6.0e2)
#define GRAVITATIONAL_CONSTANT ((double)(GRAVITATIONAL_CONSTANT_ORIGIN))
/*---------------------------------------------------------------------------*/
int      number_particles;
int      number_mass_particles;
double   * masses;
double dt;
double dt_min, dt_max;
/* Should never exceeded */
double v_max;
/* previous states */
double * p_x_prev;
double * p_y_prev;
double * v_x_prev;
double * v_y_prev;
/* temporarily store new p */
double * p_x_temp;
double * p_y_temp;
/* to store result */
double * p_x_current;
double * p_y_current;
double * v_x_current;
double * v_y_current;
/* acceleration */
double * a_x;
double * a_y;
/*---------------------------------------------------------------------------*/
void (* next)(void);
/*---------------------------------------------------------------------------
 * Helpers 
 *---------------------------------------------------------------------------*/
void print_array(char * name, int n, double * array) 
{
    fprintf(stderr, "%s ", name);
    int i;
    for(i = 0; i < n; i++)
    {
        fprintf(stderr, "%f;", array[i]);
    }
    fprintf(stderr, "\n");
}
/*---------------------------------------------------------------------------*/
void exchange(double **a, double **b) 
{
    double *temp = *a;
    *a = *b;
    *b = temp;
}
/*---------------------------------------------------------------------------*/
double dabs(double a)
{
    if(a < 0.0) return -a;
    return a;
}
/*---------------------------------------------------------------------------*/
void recalculate_v() 
{
    int i;
    double dt_invers = 1.0 / dt;
    for( i = 0; i < number_particles; i++)
    {
        v_x_current[i]  = p_x_current[i] - p_x_prev[i];
        v_y_current[i]  = p_y_current[i] - p_y_prev[i];
        v_x_current[i] *= dt_invers;
        v_y_current[i] *= dt_invers;
    }
}
/*---------------------------------------------------------------------------
 * Function typedefs
 *---------------------------------------------------------------------------*/
typedef void(* InitFunction)(int * no_particles, int * no_mass_particles, double ** masses,
        double ** p_x, double ** p_y, double ** v_x, double ** v_y) ;
typedef void(* OutputFunction)(int no_particles, int no_mass_particles, double * masses,
        double time_current, double * p_x, double * p_y, double *a_x, double *a_y) ;
/*---------------------------------------------------------------------------*/
/**
 * Initialize an velocity vector thus that the particle will circle the central
 * body at (0,0) in a circle
 */
void calculate_satelite_orbital_v(double central_mass, double x, double y, double *v_x, double *v_y) 
{
    double r_2 = x * x + y * y;
    double r   = sqrt(r_2);
    double v   = sqrt(central_mass * GRAVITATIONAL_CONSTANT / r);
    if (x == 0) 
    {
        *v_x = v; *v_y = 0;
        return;
    }
    if (y == 0) 
    {
        *v_y = v; *v_x = 0;
        return;
    }
    *v_x =    v   * y / r;
    *v_y = - *v_x * x / y;
}
/*---------------------------------------------------------------------------*/
/**
 * Custom init method that takes care to initialize the number of masses, 
 * the location and velocity vectors of these masses etc.
 */
void initialize(int * no_particles, int * no_mass_particles, double ** masses,
        double ** p_x, double ** p_y, double ** v_x, double ** v_y) 
{
#define CENTRAL_MASS 1e12
#define MEDIUM_MASS  1e7
#define X_0 0.0
#define Y_0 0.0
#define D_X_0 300000.0
#define D_Y_0 D_X_0
#define V_X_0 0.0
#define V_Y_0 0.0    
    int n_p   = *no_particles      =  10; 
    int n_m_p = *no_mass_particles =  1; 
    *masses = (double *)malloc(n_m_p * sizeof(double));
    int n;
    for(n = 0; n < n_m_p; n++) 
    {
        (*masses)[n] = 100.0;

    }
    *p_x = (double *)malloc(4 * n_p * sizeof(double));
    *p_y = *p_x + n_p;
    *v_x = *p_y + n_p;
    *v_y = *v_x + n_p;
    (*p_x)[1] = 0;
    (*p_y)[1] = 300000.0;
    (*v_x)[1] = 1000;
    (*v_y)[1] = .0;
    /* calculate_satelite_orbital_v(CENTRAL_MASS, (*p_x)[1], (*p_y)[1], &(*v_x)[1], &(*v_y)[1]);  */
    /*fprintf(stderr, "x: %f y: %f v_x : %f v_y: %f\n", (*p_x)[1], (*p_y)[1], (*v_x)[1], (*v_y)[1]); */
     for(n = 0; n < n_p; n++)  
     { 
         double dx = (double) D_X_0 * random() / (double) RAND_MAX; 
         double dy = (double) D_Y_0 * random() / (double) RAND_MAX; 
         dx -= D_X_0 * 0.5; 
         dy -= D_X_0 * 0.5; 
         (*p_x)[n] = X_0 + dx; 
         (*p_y)[n] = Y_0 + dy; 
         (*v_x)[n] = V_X_0; 
         (*v_y)[n] = V_Y_0; 
         calculate_satelite_orbital_v(CENTRAL_MASS, (*p_x)[n], (*p_y)[n], &((*v_x)[n]), &((*v_y)[n]));
     } 
    /* First mass is central mass */
    (*masses)[0] = CENTRAL_MASS;
    (*p_x)[0]    = 0;
    (*p_y)[0]    = 0;
    (*v_x)[0]    = 0;
    (*v_y)[0]    = 0;
#undef X_0
#undef Y_0
#undef D_X_0
#undef D_Y_0 
#undef V_X_0
#undef V_Y_0
}
/*---------------------------------------------------------------------------*/
/** 
 * Custom output function that will print out the current time as well as,
 * successively for each particle, its location and acceleration vector
 */
void output_stdout(int n_p, int n_m_p, double * masses, double now, 
        double * p_x, double *p_y, double *a_x, double *a_y) 
{
    printf("%f ", now);
    int i;
    for(i = 0; i < n_p; i++)
    {
        printf("%f %f %f %f ", p_x[i], p_y[i], a_x[i], a_y[i]);
    }
    printf("\n");
}
/*---------------------------------------------------------------------------
 * General initialization methods
 *---------------------------------------------------------------------------*/
void alloc_helper_arrays() 
{

    p_x_prev = (double *)malloc(4 * number_particles * sizeof(double));
    p_y_prev = p_x_prev + number_particles;
    v_x_prev = p_y_prev + number_particles;
    v_y_prev = v_x_prev + number_particles;
    p_x_temp = (double *)malloc(2 * number_particles * sizeof(double));
    p_y_temp = p_x_temp + number_particles;
    a_x      = (double *)malloc(2 * number_particles * sizeof(double));
    a_y      = a_x + number_particles;
}
/*---------------------------------------------------------------------------*/
void initialize_adaptive_solver()
{
    int i;
    double v_max = 0;
    for(i = 0; i < number_particles; i++) 
    {
        if(v_max < dabs(v_x_current[i])) 
        {
            v_max = dabs(v_x_current[i]);
        }
        if(v_max < dabs(v_y_current[i])) 
        {
            v_max = dabs(v_y_current[i]);
        }
    }
    dt_min = 0.000001 * dt;
    dt_max = 100.0 * dt;
    if(dt_max > 1.0)
    {
        dt_max = 0.99;
    }
    fprintf(stderr, "v_max : %f  dt_min: %.10f  dt_max: %.10f\n", v_max, dt_min, dt_max);
}
/*---------------------------------------------------------------------------*/
void finish_initialization() 
{
    exchange(&p_x_prev, &p_x_current);
    exchange(&p_y_prev, &p_y_current);
    int n;
    /* Approximate current location by
     * x_n+1 = x_n + dx / dt * Dt */
    for(n = 0; n < number_particles; n++) 
    {
         p_x_current[n] = p_x_prev[n] + dt * v_x_current[n];
         p_y_current[n] = p_y_prev[n] + dt * v_y_current[n]; 
    }
}
/*---------------------------------------------------------------------------*/
void initialize_solver() 
{
    initialize_adaptive_solver();
    multiRootSolver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 
            4 * number_particles);
    alloc_helper_arrays();
    finish_initialization();
}
/*---------------------------------------------------------------------------
 * Clean up functions 
 *---------------------------------------------------------------------------*/
void free_arrays()
{
   free(p_x_current); 
   free(p_x_prev); 
   free(p_x_temp); 
   free(a_x);
   free(masses);
}
/*---------------------------------------------------------------------------*/
void clean_up()
{
    gsl_multiroot_fsolver_free(multiRootSolver);
    free_arrays();
}
/*---------------------------------------------------------------------------*/
/**
 *
 */
struct exponential_params {
    double dt;
    const gsl_vector * y_last;
    const gsl_vector * a;
};
int exponential(struct exponential_params * params, gsl_vector * new_value, const gsl_vector * k)
{
   gsl_vector_memcpy(new_value, k);
   gsl_vector_scale(new_value,  0.5);
   gsl_vector_add(new_value, params->y_last);
   gsl_vector_scale(new_value, params->dt);
   gsl_vector_mul(new_value, params->a);
   return GSL_SUCCESS;
}
int exponential_zero_function(const gsl_vector * x, void * params, gsl_vector * f)
{
    struct exponential_params *exp_params = (struct exponential_params *) params;
    exponential(params, f, x);
    gsl_vector_sub(f, x);
    return GSL_SUCCESS;
}
int integrate_exponential(double start_t, double dt, double end_t,const gsl_vector * y_0, const gsl_vector * a)
{
    printf("%f\n", start_t);
    double time = start_t;
    double max_err = 0.000001;
    size_t size_dim = (size_t) y_0->size;
    gsl_vector * init_val = gsl_vector_alloc(size_dim);
    gsl_vector * k;
    gsl_vector * new_value = gsl_vector_alloc(size_dim);
    gsl_vector * old_value = gsl_vector_alloc(size_dim);
    size_t n;
    for(n = 0; n < size_dim; n++) 
    { 
        gsl_vector_set(init_val, n, 1.0);
    }
    gsl_vector_memcpy(old_value, y_0);
    struct exponential_params params;
    params.dt = dt;
    params.y_last = old_value;
    params.a = a;
    multiRootSolver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 2);
    gsl_multiroot_function zero_func;
    zero_func.f = &exponential_zero_function;
    zero_func.n = 2;
    zero_func.params = &params;
    gsl_multiroot_fsolver_set(multiRootSolver, &zero_func, init_val);
    printf("%f\n", time);
    while(time < end_t)
    {
        gsl_multiroot_fsolver_iterate (multiRootSolver);
        int test_residual_result;
        while((test_residual_result = gsl_multiroot_test_residual(gsl_multiroot_fsolver_f(multiRootSolver), max_err)) == GSL_CONTINUE)
        {
            int result = gsl_multiroot_fsolver_iterate (multiRootSolver);
            if(result != GSL_SUCCESS)
            {
                fprintf(stderr, "Iteration aborted due to %i\n", result);
                exit(1);
            }
        }
        if(test_residual_result != GSL_SUCCESS)
        {
            fprintf(stderr, "Iteration aborted due to %i\n", test_residual_result);
            exit(1);
        }
        k = gsl_multiroot_fsolver_root(multiRootSolver);
        /* Use found k to calculate actual new value */
        exponential(&params, new_value, k);
        gsl_vector_add(new_value, params.y_last);
        time += dt;
        printf("%f ", time);
        for( n = 0; n < size_dim; n++)
        { 
            double value = gsl_vector_get(new_value, n);
            printf("%f ", value);
        }
        printf("\n");
        gsl_vector_memcpy(old_value, new_value);
    }
    gsl_vector_free(init_val);
    gsl_vector_free(new_value);
    gsl_vector_free(old_value);
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
    gsl_vector_set(a  , 0, 2.0);
    gsl_vector_set(a  , 1, 5.0);
    integrate_exponential(0.0, 0.5, 100.0, y_0, a);
    gsl_vector_free(y_0);
    gsl_vector_free(a);
}

