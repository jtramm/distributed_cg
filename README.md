# distributed_cg
A parallel conjugate gradient linear solver built on MPI for solving
Poisson's Equation in 2D.

## Compilation
Compilation ettings can be configured at the top of the included makefile.
Settings include turning on/off MPI.

$> make

## Setting Problem Size
The physical domain size (n) can be set at runtime using the first
command line argument of the program.

## Running
example:

$> ./cg 75

will run the solver for a 75 x 75 physical domain size.
