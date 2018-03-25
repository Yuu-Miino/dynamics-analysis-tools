# pp = Phase Portrait

This program computes the solution orbit of [**a hybrid system**](https://en.wikipedia.org/wiki/Hybrid_system)
which includes dynamical systems and a finite state machine.
The implemented method for numerical integration is [RK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods).

For the detection of the event, we have implemented [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method)
and [the bisection method](https://en.wikipedia.org/wiki/Bisection_method).

# How to build

1. Download all files
2. Move to the root directory
3. Run `make pp`
   - then you get pp program in `bin` directory

# How to use

1. Prepare an input file
   - including parameters and initial conditions, e.g., `01.pt`
2. Run `pp [-m counts_of_maps] filename`
   - then you get the data of the solution orbit in `stdout` and 
     also get `pp.poin` which gives [the Poincare map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) of the orbit

