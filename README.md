# pp = Phase Portrait

This program computes the solution orbit of [**a hybrid system**](https://en.wikipedia.org/wiki/Hybrid_system)
which includes dynamical systems and a finite state machine.
The implemented method for numerical integration is [RK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods).

For the detection of the event, we have implemented [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method)
and [the bisection method](https://en.wikipedia.org/wiki/Bisection_method).

## How to build

1. Download all files
2. Move to the root directory where this file exists
3. Run `make pp` and get the pp program in `bin` directory

## How to use

1. Prepare an input file, e.g., `01.pt`, including follows:
   1. Parameters
   2. Initial conditions 
   3. Initial mode
   4. Additional informations you need
2. Run `pp [-m <counts_of_maps>] filename` and get the following files:
   1. `FILENAME.pp.orbit` which has the data of the solution orbit
   2. `FILENAME.pp.poin` which has the data of [the Poincare map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) of the orbit

# fix = FIXed point

This program computes the fixed point of a map.

## How to build

1. Download all files
2. Move to the root directory where this file exists
3. Run `make fix` and get the pp program in `bin` directory

## How to use

1. Prepare an input file, e.g., `01.pt`, including follows:
   1. Parameters 
   2. Initial conditions 
   3. Additional infomation you need
2. Run `fix [-p <period>] filename` and get the following files:
   1. `FILENAME.fix.pt` which has the data of the fixed point
   2. `FILENAME.fix.jac` which has the data of the Jacobian matrix with respect to the fixed point

