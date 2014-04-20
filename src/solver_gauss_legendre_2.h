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
#ifndef __SOLVER_GAUSS_LEGENDRE_2_H__
#define __SOLVER_GAUSS_LEGENDRE_2_H__
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/*---------------------------------------------------------------------------
 * The GNU Scientific Library stuff
 *---------------------------------------------------------------------------*/
#include <gsl/gsl_multiroots.h>
/*---------------------------------------------------------------------------
 * HELPERS
 *---------------------------------------------------------------------------*/
void print_vector(FILE * stream, double time, const gsl_vector * vector);
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
RootProblem * root_problem_initialize(
        int(* func)(const gsl_vector *, void *, gsl_vector *), void * params,
        gsl_vector * init_values, double max_abs_error, double max_rel_error);
/*---------------------------------------------------------------------------*/
void root_problem_reset (RootProblem * problem);
/*---------------------------------------------------------------------------*/
void root_problem_free(RootProblem * problem);
/*---------------------------------------------------------------------------*/
int  root_problem_solve(RootProblem * problem);
/*---------------------------------------------------------------------------*/
gsl_vector * root_problem_get_solution(RootProblem * problem);
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
 *
 * Each state vector contains the coordinates of all simulated particles - 
 * their position AND velocities.
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
/**
 * Minimum mass of one particle, is 1 metric ton = 10e3 kg
 */
#define BASE_MASS 10e3
/**
 * This is the gravitational constant - 6.67384 * 10^(-11) m^3 / kg / s^2 =
 * 6.67384 * 10^(-2) km^3 / kg / s^2
 */
#define GRAVITATIONAL_CONSTANT (6.67384e-2 / BASE_MASS)
/*---------------------------------------------------------------------------*/
/**
 * Function to integrate a scenario
 */
void integrate_system(void (* init_func)(gsl_vector **, gsl_vector **), 
        double dt, double dt_out, double end_time, 
        double max_abs_error, double max_rel_error);
/*---------------------------------------------------------------------------*/
#endif
