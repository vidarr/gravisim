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
#include "solver_gauss_legendre_2.h"
/*---------------------------------------------------------------------------*/
/* Solving the equations occuring by implicit GAUSS-LEGENDRE integration */
#include <gsl/gsl_multiroots.h>
/*---------------------------------------------------------------------------
 * HELPER
 *---------------------------------------------------------------------------*/
void print_vector(FILE * stream, double time, const gsl_vector * vector)
{
    int i;
    fprintf(stream, "%f ", time);
    for(i = 0; i < vector->size; i++)
    {
        fprintf(stream, "%f ", gsl_vector_get(vector, i));
    }
    fprintf(stream, "\n");
}
/*---------------------------------------------------------------------------
 * Root solver wrapper
 *---------------------------------------------------------------------------*/
RootProblem * root_problem_initialize(
        int(* func)(const gsl_vector *, void *, gsl_vector *), void * params,
        gsl_vector * init_values, double max_abs_error, double max_rel_error)
{
    size_t dimensions = init_values->size;
    RootProblem * problem = (RootProblem *)malloc(sizeof(RootProblem));
    problem->solver =
        gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, dimensions);
    problem->zero_function =
        (gsl_multiroot_function *)malloc(sizeof(gsl_multiroot_function));
    problem->zero_function->f      = func;
    problem->zero_function->params = params;
    problem->zero_function->n      = dimensions;
    problem->init_values           = init_values;
    problem->max_abs_error         = max_abs_error;
    problem->max_rel_error         = max_rel_error;
    return problem;
}
/*---------------------------------------------------------------------------*/
void root_problem_reset (RootProblem * problem)
{
    gsl_multiroot_fsolver_set(problem->solver,
            problem->zero_function, problem->init_values);
}
/*---------------------------------------------------------------------------*/
void root_problem_free(RootProblem * problem)
{
    free(problem->zero_function);
    gsl_multiroot_fsolver_free(problem->solver);
    free(problem);
}
/*---------------------------------------------------------------------------*/
int  root_problem_solve(RootProblem * problem)
{
    root_problem_reset(problem);
    int iteration_result = gsl_multiroot_fsolver_iterate (problem->solver);
    if(iteration_result != GSL_SUCCESS)
    {
        return iteration_result;
    }
    int test_residual_result = gsl_multiroot_test_residual(
            gsl_multiroot_fsolver_f(problem->solver), problem->max_abs_error);
    while(test_residual_result == GSL_CONTINUE)
    {
        iteration_result = gsl_multiroot_fsolver_iterate (problem->solver);
        if(iteration_result != GSL_SUCCESS)
        {
            return iteration_result;
        }
    test_residual_result = gsl_multiroot_test_residual(
            gsl_multiroot_fsolver_f(problem->solver), problem->max_abs_error);
    }
    if(test_residual_result != GSL_SUCCESS)
    {
        return test_residual_result;
    }
    return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------*/
gsl_vector* root_problem_get_solution(RootProblem * problem)
{
    return gsl_multiroot_fsolver_root(problem->solver);
}
/*---------------------------------------------------------------------------
 * Here goes the actual gravitational stuff
 *---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
typedef struct
{
    gsl_vector * masses;
    gsl_vector * old_coords;
    gsl_vector * temp;
    double       dt;
    double       dt_min;
    double       dt_max;
    double       f_max;
} GravitationalParams;
/*---------------------------------------------------------------------------*/
/**
 * This is the right-hand side of the diff equation of motion within
 * a gravitational field.
 * As the diff equations for both x and y coordinates are de-coupled, this
 * function can be applied separately to coords_x and coords_y.
 */
int gravitational_dif_func(const gsl_vector* x, void* params, gsl_vector* f)
{

    GravitationalParams * grav_params = (GravitationalParams *)params;
    gsl_vector * masses  = grav_params->masses;
    size_t num_masses    = masses->size;
    size_t num_particles = x->size / 4;
    size_t i, n;
    double x_val, y_val, x_diff, y_diff;
    double x_force, y_force;
    double r_invers, r_2;
    for(i = 0; i < num_particles; i++)
    {
        x_val = gsl_vector_get(x, V_X(i, num_particles));
        y_val = gsl_vector_get(x, V_Y(i, num_particles));
        gsl_vector_set(f, X(i, num_particles), x_val);
        gsl_vector_set(f, Y(i, num_particles), y_val);
        x_val = gsl_vector_get(x, X(i, num_particles));
        y_val = gsl_vector_get(x, Y(i, num_particles));
        x_force = 0;
        y_force = 0;
        for(n = 0; n < num_masses; n++)
        {
            if(n == i) continue;
            x_diff = gsl_vector_get(x, X(n, num_particles));
            y_diff = gsl_vector_get(x, Y(n, num_particles));
            x_diff -= x_val;
            y_diff -= y_val;
            r_2 = x_diff * x_diff + y_diff * y_diff;
            r_invers = 1.0 / sqrt(r_2);
            x_diff *= r_invers;
            y_diff *= r_invers;
            x_force += gsl_vector_get(masses, n) * x_diff / r_2;
            y_force += gsl_vector_get(masses, n) * y_diff / r_2;
        }
        x_force *= GRAVITATIONAL_CONSTANT;
        y_force *= GRAVITATIONAL_CONSTANT;
        gsl_vector_set(f, V_X(i, num_particles), x_force);
        gsl_vector_set(f, V_Y(i, num_particles), y_force);
    }
    return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------*/
/**
 * This is the equation we have to find the roots for.
 * GAUSS-LEGENDRE defines the k vector to be
 *    k = h * f(old_y + 1/2 * k)
 */
int gravitational_zero_func(const gsl_vector * x, void * params, gsl_vector * f)
{
    double f_max = 0.0;
    double dt_new = 0.0;
    /* x is the vector of Runge-Kutta k's ! */
    GravitationalParams * grav_params = (GravitationalParams *) params;
    gsl_vector_memcpy(grav_params->temp, x);
    gsl_vector_scale(grav_params->temp, 0.5);
    gsl_vector_add(grav_params->temp, grav_params->old_coords);
    gravitational_dif_func(grav_params->temp, params, f);
    f_max = gsl_vector_max(f);
    if(0.0 == f_max) exit(1);
    grav_params->dt = grav_params->f_max / f_max;
    if(grav_params->dt > grav_params->dt_max)
    {
        grav_params->dt = grav_params->dt_max;
    }
    else if(grav_params->dt < grav_params->dt_min)
    {
        grav_params->dt = grav_params->dt_min;
    }
    gsl_vector_scale(f, grav_params->dt);
    gsl_vector_sub(f, x);
    return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------*/
void integrate_system(void (* init_func)(gsl_vector **, gsl_vector **),
        double dt_min, double dt_max, double dt_out, double end_time,
        double f_max_newton,
        double max_abs_error, double max_rel_error)
{
    double time = 0;
    GravitationalParams * params =
        (GravitationalParams *)malloc(sizeof(GravitationalParams));
    params->dt = dt_max;
    params->dt_min = dt_min;
    params->dt_max = dt_max;
    params->f_max  = f_max_newton;
    gsl_vector * coords;
    init_func(&coords, &(params->masses));
    size_t no_particles = coords->size / 4;
    printf("# NO_PARTICLES = %i X(0) = %i   Y(0) = %i   "
           "V_X(0) = %i   V_Y(0) = %i\n",
            no_particles,
            X(0, no_particles), Y(0, no_particles),
            V_X(0, no_particles), V_Y(0, no_particles));
    params->old_coords = gsl_vector_alloc(coords->size);
    params->temp       = gsl_vector_alloc(coords->size);
    gsl_vector_memcpy(params->old_coords, coords);
    gsl_vector * init_values = gsl_vector_alloc(coords->size);
    gsl_vector_memcpy(init_values, coords);
    RootProblem * problem =
        root_problem_initialize(gravitational_zero_func, params,
        init_values, max_abs_error, max_rel_error);
    double output_time = 0;
    print_vector(stdout,time, params->old_coords);
    while(time < end_time)
    {
        root_problem_solve(problem);
        gsl_vector* k = root_problem_get_solution(problem);
        /* GAUSS-LEGENDRE : y_new = y_old + h * f(y_old + 1/2 * k)
         * where k = h * f(y_old + 1/2 * k)
         * thus
         * y_new = y_old + k                                       */
        gsl_vector_add(params->old_coords, k);
        time += params->dt;
        output_time += params->dt;
        if(output_time > dt_out)
        {
            fprintf(stderr, "%lf   %lf\n", time, params->dt);
            print_vector(stdout,time, params->old_coords);
            output_time = 0;
        }
    }
    gsl_vector_free(coords);
    gsl_vector_free(params->masses);
    gsl_vector_free(params->old_coords);
    gsl_vector_free(params->temp);
    gsl_vector_free(init_values);
    root_problem_free(problem);
    free(params);
}
