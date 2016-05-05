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
#include <assert.h>
#include "solver_gauss_legendre_2.h"
/*---------------------------------------------------------------------------*/
#define USAGE "\nUSAGE\n"   \
     "\n\n"        \
     "   ./gravisim [di|da|do|t|fm  VALUE] ...\n" \
     "\n"  \
     " Where 'di' indicates the minimum dt to calculate with\n"  \
     "       'da' indicates the maximum dt to calculate with\n"  \
     "       'fm' indicates the maximum power momentum to calculate with\n"   \
     "       't'  indicates the time to calculate up to\n"       \
     "       'do' indicates the dt to write a line after\n"
/*---------------------------------------------------------------------------*/
/* MASSES in g */
#define MASS_SUN            (1.989e33)
#define MASS_MERCURY        (0.330e27)
#define MASS_VENUS          (4.87e27)
#define MASS_EARTH          (5.974e27)
#define MASS_MOON           (7.349e25)
#define MASS_MARS           (0.642e27)
#define MASS_JUPITER        (1898e27)
#define MASS_SATURN         (568e27)
#define MASS_URANUS         (86.8e27)
#define MASS_NEPTUN         (102e27)
/* DISTANCES in km */
#define DISTANCE_SUN_MERCURY (57.9e6)
#define DISTANCE_SUN_VENUS   (108.2e6)
#define DISTANCE_SUN_EARTH  (149.6e6)
#define DISTANCE_SUN_MARS   (227.9e6)
#define DISTANCE_SUN_JUPITER (778.6e6)
#define DISTANCE_SUN_SATURN  (1433.5e6)
#define DISTANCE_SUN_URANUS  (2872.5e6)
#define DISTANCE_SUN_NEPTUN  (4495.1e6)
#define DISTANCE_EARTH_MOON (0.384e6)
/*---------------------------------------------------------------------------
 * HELPERS
 *---------------------------------------------------------------------------*/
/**
 * Return -1 if s < 0, 0 if s == 0, 1 if s > 0
 */
double sign(double s)
{
    if(0 == s) return 0;
    return s / abs(s);
}
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
        *v_x = sign(y) * v; *v_y = 0;
        return;
    }
    if (y == 0)
    {
        *v_y = sign(x) * v; *v_x = 0;
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
/**
 * Create Sun and all 'planets' along with test masses
 */
void init_solar_system(gsl_vector ** coord, gsl_vector ** masses)
{
    size_t no_satellites = 50;
    size_t no_bodies = 1 + 8 + no_satellites;
    *masses = gsl_vector_alloc(no_bodies - no_satellites);
    *coord  = gsl_vector_alloc((no_bodies) * 4);
    double v_x, v_y;
    int body, index;
    for(body = 0; body < no_bodies; body++)
    {
        for(index = 0; index < 4; index++)
        {
            gsl_vector_set(*coord, body * 4 + index, 0.0);
        }
    }
    gsl_vector_set(*masses, 0, MASS_SUN);
    gsl_vector_set(*masses, 1, MASS_MERCURY);
    gsl_vector_set(*masses, 2, MASS_VENUS);
    gsl_vector_set(*masses, 3, MASS_EARTH);
    gsl_vector_set(*masses, 4, MASS_MARS);
    gsl_vector_set(*masses, 5, MASS_JUPITER);
    gsl_vector_set(*masses, 6, MASS_SATURN);
    gsl_vector_set(*masses, 7, MASS_URANUS);
    gsl_vector_set(*masses, 8, MASS_NEPTUN);
    gsl_vector_set(*coord, X(1, no_bodies),  DISTANCE_SUN_MERCURY) ;
    gsl_vector_set(*coord, X(2, no_bodies), -DISTANCE_SUN_VENUS);
    gsl_vector_set(*coord, X(3, no_bodies),  DISTANCE_SUN_EARTH);
    gsl_vector_set(*coord, X(4, no_bodies), -DISTANCE_SUN_MARS);
    gsl_vector_set(*coord, X(5, no_bodies),  DISTANCE_SUN_JUPITER);
    gsl_vector_set(*coord, X(6, no_bodies), -DISTANCE_SUN_SATURN);
    gsl_vector_set(*coord, X(7, no_bodies),  DISTANCE_SUN_URANUS);
    gsl_vector_set(*coord, X(8, no_bodies), -DISTANCE_SUN_NEPTUN);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_MERCURY,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(1, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(1, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_VENUS,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(2, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(2, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_EARTH,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(3, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(3, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_MARS,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(4, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(4, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_JUPITER,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(5, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(5, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_SATURN,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(6, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(6, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_URANUS,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(7, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(7, no_bodies), v_y);
    calculate_satellite_orbital_v(MASS_SUN,  DISTANCE_SUN_NEPTUN,
        0.0, &v_x, &v_y);
    gsl_vector_set(*coord, V_X(8, no_bodies), v_x);
    gsl_vector_set(*coord, V_Y(8, no_bodies), v_y);
    /* Now init 'test' satellites */
    double dist_increment = DISTANCE_SUN_NEPTUN / (no_satellites + 1.0f);
    double dist = 0.0;
    for(index = 0; index < no_satellites; index++)
    {
        dist = dist_increment * (1 + index);
        gsl_vector_set(*coord, X(9 + index, no_bodies), 0.0);
        gsl_vector_set(*coord, Y(9 + index, no_bodies), dist);
        calculate_satellite_orbital_v(gsl_vector_get(*masses, 0),
                0.0, dist, &v_x, &v_y);
        gsl_vector_set(*coord, V_X(9 + index, no_bodies), v_x);
        gsl_vector_set(*coord, V_Y(9 + index, no_bodies), v_y);
    }
}
/*---------------------------------------------------------------------------
 * MAIN - do parameter parsing etc...
 *---------------------------------------------------------------------------*/
int main(int argc, char** argv)
{
    float time_end = 1000;
    float dt_min  = 0.1;
    float dt_max  = 0.1;
    float dt_out   = 0.1;
    double f_max = 1e5;
    double abs_error = 0.000000001;
    char *opts = " ";
    int i;
    enum { NONE, DT, DT_MIN, DT_MAX, DT_OUT, T_END, F_MAX } arg_expected = NONE;

    if(1 == argc)
    {
        fprintf(stderr, USAGE);
        exit(1);
    }

    for(i = 0; i < argc; i++)
    {
        switch(arg_expected)
        {
            case NONE:
                if('d' == argv[i][0])
                {
                    if(0 != argv[i][1])
                    {
                        if('i' == argv[i][1])
                        {
                            arg_expected = DT_MIN;
                        }
                        else if('a' == argv[i][1])
                        {
                            arg_expected = DT_MAX;
                        }
                        else if('t' == argv[i][1])
                        {
                            arg_expected = DT;
                        }
                        else if('o' == argv[i][1])
                        {
                            arg_expected = DT_OUT;
                        }
                    }
                }
                else if('t' == argv[i][0])
                {
                    arg_expected = T_END;
                }
                else if('f' == argv[i][0])
                {
                    arg_expected = F_MAX;
                }
                break;
            case DT:
                sscanf(argv[i], "%f", &dt_min);
                dt_max = dt_min;
                arg_expected = NONE;
                break;
            case DT_MIN:
                sscanf(argv[i], "%f", &dt_min);
                if(dt_max < dt_min)
                {
                    dt_max = dt_min;
                }
                arg_expected = NONE;
                break;
            case DT_MAX:
                sscanf(argv[i], "%f", &dt_max);
                if(dt_max < dt_min)
                {
                    dt_min = dt_max;
                }
                arg_expected = NONE;
                break;
            case DT_OUT:
                sscanf(argv[i], "%f", &dt_out);
                arg_expected = NONE;
                break;
            case T_END:
                sscanf(argv[i], "%f", &time_end);
                arg_expected = NONE;
                break;
            case F_MAX:
                sscanf(argv[i], "%f", &time_end);
                arg_expected = NONE;
                break;
            default:
                assert(! "THIS SHOULD NEVER EVER HAPPEN");
        };
    }
    if(dt_out < dt_max)
    {
        dt_out = dt_max;
    }
    printf("# time %f dt_min %f dt_max %f   f_max %f   dt_out %f abs_error %f\n",
           time_end, dt_min, dt_max, f_max, dt_out, abs_error);
    integrate_system(init_solar_system, dt_min, dt_max, dt_out, time_end,
            f_max,
            abs_error, 0.001);
    return 0;
}
