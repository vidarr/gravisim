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
#include "solver_gauss_legendre_2.h"
/*---------------------------------------------------------------------------
 * HELPERS
 *---------------------------------------------------------------------------*/
/**
 * Initialise velocity vector in a way that the particle would orbit in a 
 * plain circle around a central mass placed at (0,0)
 */
void calculate_satelite_orbital_v(double central_mass,
        double x, double y,
        double *v_x, double *v_y) 
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
/*---------------------------------------------------------------------------
 * Methods for initialising specific scenarios
 *---------------------------------------------------------------------------*/
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
    gsl_vector_set(*coord, Y(1, no), 0);
    gsl_vector_set(*coord, V_X(1, no), 0);
    gsl_vector_set(*coord, V_Y(1, no), 0.001);
    *masses = gsl_vector_alloc(1);
    gsl_vector_set(*masses, 0, 1.0 / 0.0123);
}
/**
 * Inits vectors to resemble the Earth - Moon system with add. satellite
 */
void init_earth_moon_sat(gsl_vector ** coord, gsl_vector ** masses)
{
    size_t no = 3;
    *coord = gsl_vector_alloc(4 * no);
    gsl_vector_set(*coord, X(0, no), 0);
    gsl_vector_set(*coord, Y(0, no), 0);
    gsl_vector_set(*coord, V_X(0, no), 0);
    gsl_vector_set(*coord, V_Y(0, no), 0);
    gsl_vector_set(*coord, X(1, no), 300);
    gsl_vector_set(*coord, Y(1, no), 0);
    gsl_vector_set(*coord, V_X(1, no), 0);
    gsl_vector_set(*coord, V_Y(1, no), 0.01);
    gsl_vector_set(*coord, X(2, no), 0);
    gsl_vector_set(*coord, Y(2, no), 100);
    gsl_vector_set(*coord, V_X(2, no), 0.1);
    gsl_vector_set(*coord, V_Y(2, no), 0);
    *masses = gsl_vector_alloc(2);
    gsl_vector_set(*masses, 0, 1000);
    gsl_vector_set(*masses, 1, 100);
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
    printf("# time %f dt %f dt_out %f abs_error %f\n", time_end, dt, dt_out, abs_error);
    integrate_system(init_earth_moon_sat, dt, dt_out, time_end, abs_error, 0.001);
    return 0;
}

