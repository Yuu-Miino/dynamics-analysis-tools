#ifndef _PWL_DUFFING_
#define _PWL_DUFFING_

#include "Dynamics.hpp"
#include "HSODEsolver.hpp"

#define MODE_NUM  4

#define STATE_DIM  21
#define PARA_DIM   4
#define JAC_MAT_DIM 3
#define JAC_PARA_DIM 3

#define EVENT_STATE_INDEX 0 // x
#define TIME_STATE_INDEX  2 // t

#define PRINT_DIM     2

// Dynamical system
class pwlDuffing: public Dynamics{
private:
  unsigned int mode;
public:
  pwlDuffing(unsigned int inMode);
  ~pwlDuffing(){};
  void ode(double* dxdt, const State& state, const Parameter& para);
};

// Mode property
class pwlDuffingMode: public ModeProperty{
public:
  pwlDuffingMode(unsigned int namein);
  ~pwlDuffingMode();

  // Overrides
  void eventFunction(double* EF, const State& state, const Parameter& para);
  void dEFdt(double* dEFdt, const State& state, const Parameter& para, double* dxdt);
  unsigned int modeDist(int eventIndex);
  bool inDomain(const State& state, const Parameter& para);
};

#endif
