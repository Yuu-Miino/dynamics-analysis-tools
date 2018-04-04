#include "pendulum.hpp"
#include <math.h>


void pendulum::ode(double* dxdt, const State& state, const Parameter& para){
  double t = state.getT(), x = state.getX(0), y = state.getX(1);
  double B = para.getValue(0);

  dxdt[0] = y;
  dxdt[1] = -sin(x) + B * cos(t);
}
