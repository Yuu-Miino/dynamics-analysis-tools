#ifndef _ODE_SOLVER_
#define _ODE_SOLVER_

#include <stdio.h>
#include <math.h>
#include <string>

#include"Dynamics.hpp"

#define ZERO 1e-16

class ODEsolver
{
private:
  void RK4(Dynamics& dyna, 
	   const State& current, const Parameter& para,
	   State& next);
protected:
  double h, initStep;
  std::string solverName;

public:

  ODEsolver(std::string name,double step){
    solverName = name;
    h = initStep = step;
  }
  virtual ~ODEsolver(){
  }

  void stepODEsolver(Dynamics& dyna, 
		     const State& current, const Parameter& para, 
		     State& next);
  void runODEsolver(Dynamics& dyna, 
		    const State& init, const Parameter& para, double tfinal, 
		    State& next, FILE* printDist=NULL, int printDim=-1);
};


#endif
