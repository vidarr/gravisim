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
typedef struct {
    gsl_multiroot_fsolver  * solver;
    gsl_multiroot_function * zero_function;
    double                   max_abs_error;
    double                   max_rel_error;
    gsl_vector             * init_values;
} RootProblem;
/*---------------------------------------------------------------------------*/
RootProblem * RootProblem_initialize(
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
void RootProblem_reset (RootProblem * problem) 
{
    gsl_multiroot_fsolver_set(problem->solver, 
            problem->zero_function, problem->init_values);
}
/*---------------------------------------------------------------------------*/
void RootProblem_free(RootProblem * problem)
{
    free(problem->zero_function);
    gsl_multiroot_fsolver_free(problem->solver);
    free(problem);
}
/*---------------------------------------------------------------------------*/
int  RootProblem_solve(RootProblem * problem) 
{
    RootProblem_reset(problem);
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
gsl_vector * RootProblem_get_solution(RootProblem * problem)
{
    return gsl_multiroot_fsolver_root(problem->solver);
}
/*---------------------------------------------------------------------------
 * Here goes the actual gravitational stuff
 *---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------
 *
 * All units in 10^3 km, kg and s
 * All masses are given as multiples of 1 base mass, the actual value of 
 * 1 uniform mass is given by the constant BASE_MASS
 * There are two different kind of particles:
 * 1. Particles with neglectible gravitational influence and mass "1"
 * 2. Massive particles with masses > 1, whose gravitational influence on 
 *    other particles is NOT neglectible.
 *---------------------------------------------------------------------------*/
/**
 * Minimum mass of one particle is 1 metric ton = 10e3 kg
 */
#define BASE_MASS 10e3
/**
 * This is the gravitational constant - 6.67384 * 10^(-11) m^3 / kg / s^2 =
 * 6.67384 * 10^(-2) km^3 / kg / s^2
 */
#define GRAVITATIONAL_CONSTANT (6.67384e-2 / BASE_MASS)
/**
 * contains the coordinates of all simulated particles - their position AND 
 * velocities.
 * Thus the vector is 4 * number_particles long.
 * coordinates_x[0 ... number_particles - 1] contains the x positions
 * coordinates_x[number_particles ... 2 * number_particles - 1] contains their
 * velocities. 
 *
 * The rest of the vector consists of the y coordinates, again, location and
 * velociy.
 *
 * Particles with masses different from 1 go at the top positions, thus 
 * coordinates[0 ... number_of_massive_particles - 1] massive particles
 * coordinates[number_of_massive_particles ... number_particles] particles 
 *     with neglectible masses
 * The same holds for coordinates_y.
 *
 * Another vector contains the masses of all particles that DO have a mass. 
 * As there are at most as many massive particles as there are particles all
 * in all, this vector must have a maximum size of the vector coords_x.
 */
/** 
 * Accessing the elements of the coords vector
 */
#define X(i, total_no)   (i)
#define Y(i, total_no)   (i + 2 * (total_no))
#define V_X(i, total_no) (i + (total_no))
#define V_Y(i, total_no) (i + 3 * (total_no))
/*---------------------------------------------------------------------------*/
typedef struct 
{
    gsl_vector * masses;
    gsl_vector * old_coords;
    gsl_vector * temp;
    double       dt;
} GravitationalParams;
/*---------------------------------------------------------------------------*/
/** 
 * This is the right-hand side of the diff equation of motion within 
 * a gravitational field.
 * As the diff equations for both x and y coordinates are de-coupled, this
 * function can be applied separately to coords_x and coords_y.
 */
int gravitational_dif_func(const gsl_vector * x, void * params, gsl_vector * f)
{

    GravitationalParams * grav_params = (GravitationalParams *)params;
    gsl_vector * masses  = grav_params->masses;
    size_t num_masses    = masses->size;
    size_t num_particles = x->size / 4;
    size_t i, n;
    double x_val, y_val, x_diff, y_diff;
    double x_force, y_force;
    double r_2;
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
            x_force += gsl_vector_get(masses, n) * x_diff / r_2;
            y_force += gsl_vector_get(masses, n) * y_diff / r_2;
        }
        if(i < num_masses)
        {
            x_force *= gsl_vector_get(masses, i);
            y_force *= gsl_vector_get(masses, i);
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
    /* x is the vector of Runge-Kutta k's ! */
    GravitationalParams * grav_params = (GravitationalParams *) params;
    gsl_vector_memcpy(grav_params->temp, x);
    gsl_vector_scale(grav_params->temp, 0.5);
    gsl_vector_add(grav_params->temp, grav_params->old_coords);
    gravitational_dif_func(grav_params->temp, params, f);
    gsl_vector_scale(f, grav_params->dt);
    gsl_vector_sub(f, x);
    return GSL_SUCCESS;
}
/*---------------------------------------------------------------------------*/
void integrate_system(void (* init_func)(gsl_vector **, gsl_vector **), 
        double dt, double dt_out, double end_time, 
        double max_abs_error, double max_rel_error)
{
    double time = 0;
    GravitationalParams * params = 
        (GravitationalParams *)malloc(sizeof(GravitationalParams));
    params->dt = dt;
    gsl_vector * coords;
    init_func(&coords, &(params->masses));
    params->old_coords = gsl_vector_alloc(coords->size);
    params->temp       = gsl_vector_alloc(coords->size);
    gsl_vector_memcpy(params->old_coords, coords);
    gsl_vector * init_values = gsl_vector_alloc(coords->size);
    gsl_vector_memcpy(init_values, coords);
    RootProblem * problem = RootProblem_initialize(gravitational_zero_func, params,
        init_values, max_abs_error, max_rel_error);
    double output_time = 0;
    print_vector(stdout,time, params->old_coords);
    while(time < end_time)
    {
        RootProblem_solve(problem);
        gsl_vector * k = RootProblem_get_solution(problem);
        /* GAUSS-LEGENDRE : y_new = y_old + h * f(y_old + 1/2 * k) 
         * where k = h * f(y_old + 1/2 * k)
         * thus 
         * y_new = y_old + k                                       */
        gsl_vector_add(params->old_coords, k);
        time += dt;
        output_time += dt;
        if(output_time > dt_out)
        {
            fprintf(stderr, "%f \n", time);
            print_vector(stdout,time, params->old_coords);
            output_time = 0;
        }
    }
    gsl_vector_free(params->old_coords);
    gsl_vector_free(params->temp);
    gsl_vector_free(init_values);
    RootProblem_free(problem);
    free(params);
}
/*---------------------------------------------------------------------------*/
/**
 * Inits vectors to resemble the Earth - Moon system
 */
void init_earth_moon(gsl_vector ** coord, gsl_vector ** masses)
{
    size_t no = 2;
    *coord = gsl_vector_alloc(4 * no);
    gsl_vector_set(*coord, X(0, no), 0);
    gsl_vector_set(*coord, Y(0, no), 0);
    gsl_vector_set(*coord, V_X(0, no), 0);
    gsl_vector_set(*coord, V_Y(0, no), 0);
    gsl_vector_set(*coord, X(1, no), 300);
    gsl_vector_set(*coord, V_X(1, no), 0);
    gsl_vector_set(*coord, V_Y(1, no), 1);
    *masses = gsl_vector_alloc(1);
    gsl_vector_set(*masses, 0, 1.0 / 0.0123);
}
/*---------------------------------------------------------------------------
 * MAIN - do parameter parsing etc...
 *---------------------------------------------------------------------------*/
int main(int argc, char** argv) 
{
    float time_end = 1000;
    float dt       = 0.1;
    float dt_out   = 0.1;
    double abs_error = 0.000000001;
    char *opts = " ";

    if(argc > 1) 
    {
        sscanf(argv[1], "%f", &time_end);
    }
    if(argc > 2) 
    {
        sscanf(argv[2], "%f", &dt);
    }
    if(argc > 3) 
    {
        sscanf(argv[3], "%f", &dt_out);
    }
    if(argc > 4) 
    { 
        if(argv[4][0] == 'a') 
        { 
        } 
    } 
    fprintf(stderr, "time %f dt %f dt_out %f\n", time_end, dt, dt_out);
    integrate_system(init_earth_moon, dt, dt_out, time_end, abs_error, 0.001);
    return 0;
}

