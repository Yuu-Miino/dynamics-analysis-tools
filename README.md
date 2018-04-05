# Dynamics Analysis Tools

This is a C++ program package containing the tools to simulate and analyze the dynamical systems.

- Compatible systems:
	- [Continuous-time dynamical systems](https://en.wikipedia.org/wiki/Dynamical_system_(definition)) (either autonomous or non-autonomous)
	- [*a hybrid system*](https://en.wikipedia.org/wiki/Hybrid_system)

An input file format is as follows:

```
# filename.pt
0.1 0.1 0.1 	# parameters
1.0 1.0     	# initial conditions
1           	# initial mode (only necessary for hybrid systems)
```

All programs have a option `--help` that displays the usage.

## How to build?

```
$ cd dynamics-analysis-tools
$ make
```

## pp = Phase Portrait

This program computes the solution orbit of the system.
The implemented method for numerical integration is [RK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods).

For the hybrid systems, we have implemented [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method)
and [the bisection method](https://en.wikipedia.org/wiki/Bisection_method)
to detect the event.

### Usage

- `$ ./pp [-m <counts_of_maps>] FILENAME`

### Outputs
- `FILENAME.pp.orbit:` The data of the solution orbit.  
- `FILENAME.pp.poin:` The data of [the Poincare map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) of the orbit.

### Options
- `-m <count_of_maps>:` Set the count of maps to calculate (default: 100).

## fix = FIXed point

This program computes [the fixed point](https://en.wikipedia.org/wiki/Fixed_point_(mathematics)) of a map.
This also supports the periodic point.

### Usage

- `$ ./fix [-p <period>] FILENAME`

### Outputs
- `FILENAME.fix.pt:` The data of the fixed point with the parameters.  
- `FILENAME.fix.jac:` The data of [the Jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) with respect to the fixed point.

If using the parameter continuation,

- `FILENAME.bf2.cont:` The data set to plot a 2-dim bifurcation diagram.

### Options
- `-p <period>:` Set the period of the periodic point (default: 1).  
- `-C <parameter_index>:` Use the parameter continuation with the selected parameter.

## bf1 = 1-dim BiFurcation diagram

This program computes a data set to plot a [1-dimensional bifurcation diagram](https://en.wikipedia.org/wiki/Bifurcation_diagram) with a parameter continuation.

### Usage

```
$ ./bf1 [-m <counts_of_maps>] [-i <index_of_parameter>] \
        [-e <end_value>] [-r <resolution>] filename
```

### Outputs

- `FILENAME.bf1:` The data set to plot a 1-dim bifurcation diagram.


### Options
- `-m <counts_of_maps>:` Set the count of maps to calculate (default: 100).
- `-i <index_of_parameter>:` Select the index of the continuation parameter (default: 0).
- `-e <end_value>:` Set the end value of the continuation.  
- `-r <resolution>:` Set the resolution of the continuation (default: 100).

## bf2 = 2-dim BiFurcation diagram

This program computes a bifurcation parameter of a map.

### Usage

```
$ ./bf2 [-p <period>] [-i <parameter_index>] [-G | -I] \
        [-C <parameter_index> [-s <parameter_step>]] filename
```

### Outputs
- `FILENAME.bf2.pt:` The data of the fixed point with the parameters.  

If using the parameter continuation,

- `FILENAME.bf2.cont:` The data set to plot a 2-dim bifurcation diagram.

### Options
- `-p <period>:` Set the period of the fixed point getting bifurcation.
- `-i <parameter_index>:` Set the bifurcation parameter index.
- `-G | -I:` Set the type of the bifurcation. (G: tangent bifurcation, I: period-doubling bifurcation)
- `-C <parameter_index>:` Use the parameter continuation with the selected parameter.
	- `-s <parameter_step>` Set the step of the parameter continuation.
