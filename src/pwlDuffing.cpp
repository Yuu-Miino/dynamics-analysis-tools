#include "pwlDuffing.hpp"
#include <math.h>

// Internal function
double f(double x, double th, unsigned int mode){
  double th1 = -1.5, th2 = th, th3 = -2, th4 = 0.05;

  switch(mode){
  case 0: return (x * 3.0 + th1);
  case 1: return (x / 3.0 + th2);
  case 2: return (x * 3.0 + th3);
  case 3: return (x / 3.0 + th4);
  default:
    fprintf(stderr,"Error: undefined mode = %d in pwlDuffing::f\n",mode);
    exit(1);
  }
}

// Class member functions 
pwlDuffing::pwlDuffing(unsigned int inMode){
  DIM = 2;
  DIM_state = 2;
  mode = inMode;
}

void pwlDuffing::ode(double* dxdt, const State& state, const Parameter& para){
  double t = state.getT(), x = state.getX(0), y = state.getX(1);
  double k, b0, b, theta;
  double fx;
  
  k  = para.getValue(0);
  b0 = para.getValue(1);
  b  = para.getValue(2);
  theta = para.getValue(3);
  fx = f(x, theta, mode);

  dxdt[0] = y;
  dxdt[1] = -k * y - fx + b0 + b*cos(t);
}
