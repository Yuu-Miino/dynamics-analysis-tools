#include "pwlDuffing.hpp"
#include <math.h>

// Internal function
double f(double x, const Parameter& para, unsigned int mode){
  double th0, th1, th2, th3;
  th0 = para.getValue(3);
  th1 = para.getValue(4);
  th2 = para.getValue(5);
  th3 = para.getValue(6);
  
  switch(mode){
  case 0: return (x * 3.0 + th0);
  case 1: return (x / 3.0 + th1);
  case 2: return (x * 3.0 + th2);
  case 3: return (x / 3.0 + th3);
  default:
    fprintf(stderr,"Error: undefined mode = %d in pwlDuffing::f\n",mode);
    exit(1);
  }
}

double df(double in, int mode){
  switch(mode){
  case 0: return (in * 3.0);
  case 1: return (in / 3.0);
  case 2: return (in * 3.0);
  case 3: return (in / 3.0);
  default:
    fprintf(stderr,"Error: undefined mode = %d in pwlDuffing::f\n",mode);
    exit(1);
  }
}

// Class member functions 
pwlDuffing::pwlDuffing(unsigned int inMode){
  mode = inMode;
}

void pwlDuffing::ode(double* dxdt, const State& state, const Parameter& para){
  double k, b0, b;
  double fx, fxx, fxy, fxz, fxk, fxb0, fxb;

  double x = state.getX(0), y = state.getX(1), z = state.getX(2);
  double xx = state.getX(3), yx = state.getX(4), zx = state.getX(5);
  double xy = state.getX(6), yy = state.getX(7), zy = state.getX(8);
  double xz = state.getX(9), yz = state.getX(10), zz = state.getX(11);
  double xk = state.getX(12), yk = state.getX(13), zk = state.getX(14);
  double xb0 = state.getX(15), yb0  = state.getX(16), zb0 = state.getX(17);
  double xb  = state.getX(18), yb = state.getX(19), zb = state.getX(20);

  k  = para.getValue(0);
  b0 = para.getValue(1);
  b  = para.getValue(2);
  
  fx  = f(x, para, mode);
  fxx = df(xx, mode);
  fxy = df(xy, mode);
  fxz = df(xz, mode);
  fxk = df(xk, mode);
  fxb0 = df(xb0, mode);
  fxb  = df(xb, mode);

  // dx/dt
  dxdt[0] = y;
  dxdt[1] = -k * y - fx + b0 + b * cos(z);
  dxdt[2] = 1;

  // dT/dx
  dxdt[3] = yx;
  dxdt[4] = -k * yx - fxx - b * sin(z) * zx;
  dxdt[5] = 0;

  // dT/dy
  dxdt[6] = yy;
  dxdt[7] = -k * yy - fxy - b * sin(z) * zy;
  dxdt[8] = 0;

  // dT/dz
  dxdt[9]  = yz;
  dxdt[10] = -k * yz - fxz - b * sin(z) * zz;
  dxdt[11] = 0;

  // dT/dk
  dxdt[12] = yk;
  dxdt[13] = -y - k * yk - fxk - b * sin(z) * zk;
  dxdt[14] = 0;

  // dT/db0
  dxdt[15] = yb0;
  dxdt[16] = -k * yb0 - fxb0 + 1 - b * sin(z) * zb0;
  dxdt[17] = 0;

  // dT/db
  dxdt[18] = yb;
  dxdt[19] = -k * yb - fxb - b * sin(z) * zb + cos(z);
  dxdt[20] = 0;
}
