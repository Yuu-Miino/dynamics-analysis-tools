#ifndef _PENDULUM_
#define _PENDULUM_

#include "Dynamics.hpp"
#include "ODEsolver.hpp"

#define STATE_DIM  2
#define PARA_DIM   1
#define JAC_MAT_DIM 2
#define JAC_PARA_DIM 2

#define PRINT_DIM     2

// Dynamical system
class pendulum: public Dynamics{
public:
  pendulum(){}
  ~pendulum(){}
  void ode(double* dxdt, const State& state, const Parameter& para);
};

#endif
