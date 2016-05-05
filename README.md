# GRAVISIM

## Overview

GRAVISIM aims to provide a simple framework to calculate trajectories of objects
subjected to gravitational force.

By default, it will compile to an executable that will calculate a scenario
of our solar system, i.e. the sun including the 8 planets, as well as 25 test
masses distributed in between the planets up to an arbitray amount of time.

## Compile

In order to compile, you will need the
[GNU Scientific Library](https://www.gnu.org/software/gsl/ "GNU Scientific Library")
headers and libraries installed.

Compilation is done via  
`cd gravisim && make`

This will compile the source files underneath `src/` and create a binary
`gravisim` underneath `bin`.
Calling this binary like `bin/gravisim t 1000` will calculate the default
scenario for 1000 days.
Thereby, it will print to stdout a lengthy line of numbers for each `0.1` days.
This line consists of the x coordinates of all objects involved (e.g.
for the sun, the 8 planets and the 25 test masses), followed by the
x coords of their velocities, followed by the y coords of their places,
followed by the y coords of their velocities.
Moreover, it will print to stderr two numbers, the current time, and the
current delta t used for the calculation.

# Howto Use

TDB

## Beyond the scenes

### The numerics

GRAVISIM basically calculates for each object the gravitational force `(f_x, f_y)`
dragging it by using Newtons Law of Gravity.

It then takes for each object the vector `(x, dx/dt, y, dy/dt)`
and starting with given values for `x_0, dx/dt_0, y_0, dy/dt_0`, it will
calculate the new value of this vector after a small amount of time.
The simplest way would be to just calculate
`(x_1, dx/dt_1, y_1, dy_dt_1) = (x_0, dx/dt_0, y_0, dy/dt_0) + dt * (dx/dt_0, f_x_0, dy/dt_0, f_y_0)` .  
This is called the *Euler Method*, one out of a large number of numerical *schemes*
called *Runge-Kutta schemes*.
However, since *Euler* is not particulary good, GRAVISIM instead uses
the *implicit* *Gauss-Legendre scheme of second order*.

Moreover, it facilitates an adaptive solver, fitting the used `delta t` to
keep the `(dx/dt, f_x, dy/dt, f_y) * dt` just underneath a given value `f_max`.
