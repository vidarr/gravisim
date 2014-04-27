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
/*---------------------------------------------------------------------------*/
#define MASS_SUN            (1.989 * 10e33)
#define MASS_EARTH          (5.974 * 10e27)
#define MASS_MOON           (7.349 * 10e25)
#define DISTANCE_SUN_EARTH  (149.6 * 10e6)
#define DISTANCE_EARTH_MOON (0.384 * 10e6)
/*---------------------------------------------------------------------------
 * HELPERS
 *---------------------------------------------------------------------------*/
/**
 * Initialise velocity vector in a way that the particle would orbit in a 
 * plain circle around a central mass placed at (0,0)
 */
void calculate_satellite_orbital_v(double central_mass,
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
    gsl_vector_set(*coord, X(1, no), 300000);
    gsl_vector_set(*coord, Y(1, no), 0);
    gsl_vector_set(*coord, V_X(1, no), 0);
    gsl_vector_set(*coord, V_Y(1, no), 1);
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
    gsl_vector_set(*coord, X(1, no), 300000);
    gsl_vector_set(*coord, Y(1, no), 0);
    gsl_vector_set(*coord, V_X(1, no), 0);
    gsl_vector_set(*coord, V_Y(1, no), 1);
    gsl_vector_set(*coord, X(2, no), 0);
    gsl_vector_set(*coord, Y(2, no), 100000);
    gsl_vector_set(*coord, V_X(2, no), 10);
    gsl_vector_set(*coord, V_Y(2, no), 0);
    *masses = gsl_vector_alloc(2);
    gsl_vector_set(*masses, 0, 1000);
    gsl_vector_set(*masses, 1, 100);
}
/**
 * Inits vectors to resemble the Earth - Moon system with add. satellite 
 * Initialises velocities such that both satellites orbit around the central
 * mass
 */
void init_earth_moon_sat_orbit(gsl_vector ** coord, gsl_vector ** masses)
{
    size_t no = 3;
    double x, y, v_x, v_y;
    *masses = gsl_vector_alloc(2);
    gsl_vector_set(*masses, 0, 1000);
    gsl_vector_set(*masses, 1, 40);
    *coord = gsl_vector_alloc(4 * no);
    gsl_vector_set(*coord, X(0, no), 0);
    gsl_vector_set(*coord, Y(0, no), 0);
    gsl_vector_set(*coord, V_X(0, no), 0);
    gsl_vector_set(*coord, V_Y(0, no), 0);
    gsl_vector_set(*coord, X(1, no), 300000);
    gsl_vector_set(*coord, Y(1, no), 0);
    calculate_satellite_orbital_v(1000, 300000, 0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(1, no), v_x);
    gsl_vector_set(*coord, V_Y(1, no), v_y);
    gsl_vector_set(*coord, X(2, no), 0);
    gsl_vector_set(*coord, Y(2, no), 100000);
    calculate_satellite_orbital_v(1000, 0, 100000, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(2, no), v_x);
    gsl_vector_set(*coord, V_Y(2, no), v_y);
}
/**
 * Inits vectors to resemble the Earth with satellites in different dists.
 * Initialises velocities such that all satellites orbit around the central
 * mass
 */
void init_earth_sat_orbit(gsl_vector ** coord, gsl_vector ** masses)
{
    size_t no = 10;
    printf("# x(0) in col %u, y(0) in col %u, v_x(0) in col %u, v_y(0) in col %u\n", 
            X(0, no), Y(0, no), V_X(0, no), V_Y(0, no));
    double x, y, v_x, v_y;
    /* masses */
    *masses = gsl_vector_alloc(1);
    gsl_vector_set(*masses, 0, 1000000000);
    /* Coordinates */
    *coord = gsl_vector_alloc(4 * no);
    gsl_vector_set(*coord, X(0, no), 0);
    gsl_vector_set(*coord, Y(0, no), 0);
    gsl_vector_set(*coord, V_X(0, no), 0);
    gsl_vector_set(*coord, V_Y(0, no), 0);
    size_t no_satellite;
    for(no_satellite = 1; no_satellite < no; no_satellite++)
    {
        double dist = (10000 + no_satellite * 300);
        gsl_vector_set(*coord, X(no_satellite, no), dist);
        gsl_vector_set(*coord, Y(no_satellite, no), 0);
        calculate_satellite_orbital_v(gsl_vector_get(*masses, 0), 
                dist, 0, &v_x, &v_y);
        gsl_vector_set(*coord, V_X(no_satellite, no), v_x);
        gsl_vector_set(*coord, V_Y(no_satellite, no), v_y);
    }
}
/**
 * Just like init_earth_sat_orbit , but takes the furthest out satellite
 * assigns it a mass of 1/50 of the central mass, and puts it 1/4 furher away
 * from the central mass
 */
void init_earth_sat_moon_orbit(gsl_vector ** coord, gsl_vector ** masses)
{
    init_earth_sat_orbit(coord, masses);
    double central_mass = gsl_vector_get(*masses, 0);
    double moon_mass = 1.0/50.0 * central_mass;
    gsl_vector_free(*masses);
    *masses = gsl_vector_alloc(2);
    gsl_vector_set(*masses, 0, central_mass);
    gsl_vector_set(*masses, 1, moon_mass);
    size_t no_satellites = (*coord)->size / 4;
    double x = gsl_vector_get(*coord, X(no_satellites - 1, no_satellites));
    x *= 1.25;
    double y = gsl_vector_get(*coord, Y(no_satellites - 1, no_satellites));
    y *= 1.25;
    gsl_vector_set(*coord, X(no_satellites - 1, no_satellites), x);
    gsl_vector_set(*coord, Y(no_satellites - 1, no_satellites), y);
    double v_x, v_y;
    calculate_satellite_orbital_v(central_mass, x, y, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(no_satellites - 1, no_satellites), v_x);
    gsl_vector_set(*coord, V_Y(no_satellites - 1, no_satellites), v_y);
    /* at last, we need to swap the new massive moon into position 1 of 
     * the coord array such that the index of the moon is the same for 
     * both the mass and coord array                                    */
    double x_sat   = gsl_vector_get(*coord, X(1, no_satellites));
    double y_sat   = gsl_vector_get(*coord, Y(1, no_satellites));
    double v_x_sat = gsl_vector_get(*coord, V_X(1, no_satellites));
    double v_y_sat = gsl_vector_get(*coord, V_Y(1, no_satellites));
    gsl_vector_set(*coord, X(1, no_satellites), gsl_vector_get(*coord, X(no_satellites - 1, no_satellites)));
    gsl_vector_set(*coord, Y(1, no_satellites), gsl_vector_get(*coord, Y(no_satellites - 1, no_satellites)));
    gsl_vector_set(*coord, V_X(1, no_satellites), gsl_vector_get(*coord, V_X(no_satellites - 1, no_satellites)));
    gsl_vector_set(*coord, V_Y(1, no_satellites), gsl_vector_get(*coord, V_Y(no_satellites - 1, no_satellites)));
    gsl_vector_set(*coord, X(no_satellites - 1, no_satellites), x_sat);
    gsl_vector_set(*coord, Y(no_satellites - 1, no_satellites), y_sat);
    gsl_vector_set(*coord, V_X(no_satellites - 1, no_satellites), v_x_sat);
    gsl_vector_set(*coord, V_Y(no_satellites - 1, no_satellites), v_y_sat);
}
/**
 * Create Sun - Earth - Moon system
 */
void init_sun_earth_moon(gsl_vector ** coord, gsl_vector ** masses)
{
    *masses = gsl_vector_alloc(3);
    *coord  = gsl_vector_alloc(3 * 4);
    double v_x, v_y_moon, v_y_earth;
    double omega; 
    int body, index;
    for(body = 0; body < 3; body++)
    {
        for(index = 0; index < 4; index++)
        {
            gsl_vector_set(*coord, body + index, 0.0);
        }
    }
    gsl_vector_set(*masses, 0, MASS_SUN);
    gsl_vector_set(*masses, 1, MASS_EARTH);
    gsl_vector_set(*masses, 2, MASS_MOON);
    gsl_vector_set(*coord, X(1, 3), DISTANCE_SUN_EARTH) ;
    gsl_vector_set(*coord, X(2, 3), DISTANCE_SUN_EARTH + DISTANCE_EARTH_MOON);
    calculate_satellite_orbital_v(MASS_SUN, DISTANCE_SUN_EARTH, 0.0, &v_x, &v_y_earth);
    gsl_vector_set(*coord, V_X(1, 3), v_x);
    gsl_vector_set(*coord, V_Y(1, 3), v_y_earth);
    omega = v_y_earth / gsl_vector_get(*coord, X(1, 3));
    calculate_satellite_orbital_v(MASS_EARTH, DISTANCE_EARTH_MOON, 0.0, &v_x, &v_y_moon);
    gsl_vector_set(*coord, V_X(2, 3), v_x);
    gsl_vector_set(*coord, V_Y(2, 3), omega * gsl_vector_get(*coord, X(2, 3)) + v_y_moon);
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
    integrate_system(init_sun_earth_moon, dt, dt_out, time_end, abs_error, 0.001);
    return 0;
}

