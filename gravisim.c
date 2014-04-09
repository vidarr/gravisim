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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


/* All masses in 10^6 kg */
/* All distances in 10^6 m */
#define GRAVITATIONAL_CONSTANT_ORIGIN ((double)6.0e2)
#define GRAVITATIONAL_CONSTANT ((double)(GRAVITATIONAL_CONSTANT_ORIGIN))

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

/* temporarily store new v */
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

void (* next)(void);

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

void exchange(double **a, double **b) 
{
    double *temp = *a;
    *a = *b;
    *b = temp;
}

double dabs(double a)
{
    if(a < 0.0) return -a;
    return a;
}

typedef void(* InitFunction)(int * no_particles, int * no_mass_particles, double ** masses,
        double ** p_x, double ** p_y, double ** v_x, double ** v_y) ;
typedef void(* OutputFunction)(int no_particles, int no_mass_particles, double * masses,
        double time_current, double * p_x, double * p_y, double *a_x, double *a_y) ;


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
    *p_x = (double *)malloc(n_p * sizeof(double));
    *p_y = (double *)malloc(n_p * sizeof(double));
    *v_x = (double *)malloc(n_p * sizeof(double));
    *v_y = (double *)malloc(n_p * sizeof(double));
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

void alloc_helper_arrays() 
{

    p_x_prev = (double *)malloc(number_particles * sizeof(double));
    p_y_prev = (double *)malloc(number_particles * sizeof(double));
    p_x_temp = (double *)malloc(number_particles * sizeof(double));
    p_y_temp = (double *)malloc(number_particles * sizeof(double));
    a_x      = (double *)malloc(number_particles * sizeof(double));
    a_y      = (double *)malloc(number_particles * sizeof(double));
}

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

void recalculate_acceleration() 
{
    int n;
    for(n = 0; n < number_particles; n++)
    {
        a_x[n] = a_y[n] = 0.0;
        int i;
        for(i = 0; i < number_mass_particles; i++) 
        {
            if(n == i) continue;
            double r_x = p_x_current[i] - p_x_current[n];
            double r_y = p_y_current[i] - p_y_current[n];
            double r_2 = r_x * r_x + r_y * r_y;
            double r_inverted = 1.0 / sqrt(r_2);
            double a_abs = GRAVITATIONAL_CONSTANT * masses[i] / r_2;
            a_x[n] += a_abs * r_x * r_inverted;
            a_y[n] += a_abs * r_y * r_inverted;
        }
    }
}

void static_next() 
{
    recalculate_acceleration();
    double dt_square_half =0.5f * dt * dt;
    int n;
    for(n = 0; n < number_particles; n++)
    {
        /* effectively: 
         * x_1 = x_0 + Dt * dx / dt + Dt * Dt * a */
        /* Only valid if Dt does not change */
        p_x_temp[n] = 2.0 * p_x_current[n] - p_x_prev[n] + dt_square_half * a_x[n];
        p_y_temp[n] = 2.0 * p_y_current[n] - p_y_prev[n] + dt_square_half * a_y[n];
    }
    exchange(&p_x_prev, &p_x_current);
    exchange(&p_y_prev, &p_y_current);
    exchange(&p_x_temp, &p_x_current);
    exchange(&p_y_temp, &p_y_current);
}

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

void adapt_dt() 
{
    int    n          = 0;
    double a_max      = 0.0;
    double a_max_temp = 0.0;
    double dt_temp    = 0.0;
    for(n = 0; n < number_particles; n++)
    {
        a_max_temp = dabs(a_x[n]);
        if(a_max_temp > a_max) a_max = a_max_temp;
        a_max_temp = dabs(a_y[n]);
        if(a_max_temp > a_max) a_max = a_max_temp;
    }
    if(a_max == 0.0) 
    {
        dt = dt_max;
        return;
    }
    dt_temp = v_max / a_max;
    if(dt_temp < dt_min) dt_temp = dt_min;
    if(dt_temp > dt_max) dt_max  = dt_max;
    dt = dt_temp;
}

void adaptive_next() 
{
    recalculate_acceleration();
    adapt_dt();
    double dt_square_half =0.5f * dt * dt;
    int n;
    for(n = 0; n < number_particles; n++)
    {
        /* p_x_current[n] = p_x_current[n] +  dt * v_x_current[n] + dt_square_half * a_x[n];
           p_y_current[n] = p_y_current[n] +  dt * v_y_current[n] + dt_square_half * a_y[n]; */
        p_x_temp[n] = 2.0 * p_x_current[n] - p_x_prev[n] + dt_square_half * a_x[n];
        p_y_temp[n] = 2.0 * p_y_current[n] - p_y_prev[n] + dt_square_half * a_y[n];
    }
    exchange(&p_x_prev, &p_x_current);
    exchange(&p_y_prev, &p_y_current);
    exchange(&p_x_temp, &p_x_current);
    exchange(&p_y_temp, &p_y_current);
    /* recalculate_v();*/
}

void free_arrays()
{
   free(p_x_current); 
   free(p_y_current); 
   free(p_x_prev); 
   free(p_y_prev); 
   free(p_x_temp); 
   free(p_y_temp); 
   free(a_x);
   free(a_y);
   free(masses);
}

void output_stdout(int n_p, int n_m_p, double * masses, double now, double * p_x, double *p_y, double *a_x, double *a_y) 
{
    printf("%f ", now);
    int i;
    for(i = 0; i < n_p; i++)
    {
        printf("%f %f %f %f ", p_x[i], p_y[i], a_x[i], a_y[i]);
    }
    printf("\n");
}

void gravitate(InitFunction init, OutputFunction output, double time_end, double initial_dt, double output_dt)
{
    dt = initial_dt;
    double time_current = 0;
    double time_since_output = 0;
    init(&number_particles, &number_mass_particles, &masses,
        &p_x_current, &p_y_current, &v_x_current, &v_y_current); 
    initialize_adaptive_solver();
    alloc_helper_arrays();
    finish_initialization();
    output(number_particles, number_mass_particles, masses, 
        time_current, p_x_current, p_y_current, a_x, a_y);
    while(time_current < time_end) 
    {
        next();
        time_since_output += dt;
        if(time_since_output >= output_dt) 
        {
            fprintf(stderr, "time: %.15f    dt: %.15f\n", time_current, dt);
            output(number_particles, number_mass_particles, masses, 
                time_current, p_x_current, p_y_current, a_x, a_y);
            time_since_output = 0;
        }
        time_current += dt;
    }
    free_arrays();
}


int main(int argc, char** argv) 
{
    float time_end = 1000;
    float dt       = 0.1;
    float dt_out   = 0.1;
    char *opts = " ";
    next = static_next;

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
            next = adaptive_next; 
        } 
    } 
    fprintf(stderr, "time %f dt %f dt_out %f\n", time_end, dt, dt_out);
    fprintf(stderr, "argc %i time %s dt %s dt_out %s\n", argc,  argv[1], argv[2], argv[3]);
    gravitate(initialize, output_stdout, time_end, dt, dt_out);
    return 0;
}

