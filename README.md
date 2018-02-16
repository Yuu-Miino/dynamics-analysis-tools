# pp = Phase Portrait

This program computes the solution orbit of **a hybrid system**
which includes dynamical systems and a finite state machine.
The implemented method for numerical integration is RK4.

For the detection of the event, we have implemented Newton's method and the binary method.

# How to build

1. Download all files
2. Move to the **root** directory of source files
3. Run `make pp`
   - then you get pp program in **bin** directory

# How to use

1. Prepare an input file
   - including parameters and initial conditions
2. Run `pp (counts of maps) (input file name)`
   - then you get the data of the solution orbit in `stdout` and also get *pp.poin* which include the Poincar\'e map of the orbit

