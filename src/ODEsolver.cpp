#include "ODEsolver.hpp"

// Class member functions
void ODEsolver::stepODEsolver(Dynamics& dyna, 
			     const State& current, const Parameter& para, 
			     State& next)
{
  if(solverName == "RK4"){ RK4(dyna,current,para,next);}
  else{
    fprintf(stderr,"Error: undefined solverName = %s\n",solverName.c_str());
    exit(1);
  }
}

void ODEsolver::RK4(Dynamics& dyna, 
		    const State& current, const Parameter& para,
		    State& next)
{
  int dim = current.getDIM();
  double k0[dim], k1[dim], k2[dim], k3[dim];
  State tmp(current);

  // Get K0
  dyna.ode(k0, tmp, para);

  // Get K1
  tmp.setT(current.getT()+h/2.0);
  for( unsigned int i = 0; i < dim; i++)
    tmp.setX(i, current.getX(i) + k0[i]*(h/2.0));
  dyna.ode(k1, tmp, para);
  
  // Get K2
  tmp.setT(current.getT()+h/2.0);
  for( unsigned int i = 0; i < dim; i++)
    tmp.setX(i, current.getX(i) + k1[i]*(h/2.0));
  dyna.ode(k2, tmp, para);
  
  // Get K3
  tmp.setT(current.getT()+h);
  for( unsigned int i = 0; i < dim; i++)
    tmp.setX(i, current.getX(i) + k2[i]*h);
  dyna.ode(k3, tmp, para);
  
  // Get next point
  next.setT(current.getT() + h);
  for( unsigned int i = 0; i < dim; i++){
    next.setX(i, current.getX(i) + (h/6.0)*(k0[i] + 2.0*k1[i] + 2.0*k2[i] + k3[i]));
  }
}

void ODEsolver::runODEsolver(Dynamics& dyna,
			     const State& init, const Parameter& para, double tfinal, 
			     State& dst, FILE* printDist,int printDim)
{
  State current(init), next(init);
  bool FINISH = false;

  if(printDim == -1) printDim = init.getDIM();

  // if print the orbit
  if(printDist != NULL){
    init.printT(printDist); init.printX(printDist,printDim);
    fprintf(printDist,"\n");
  }

  // main loop
  while(!FINISH){
    stepODEsolver(dyna,current,para,next);
    
    if((next.getT() - tfinal)*(h/fabs(h)) > ZERO){
      h = tfinal - current.getT();
      stepODEsolver(dyna,current,para,next);
      h = initStep;
      FINISH = true;
    }

    current = next;
    // if print the orbit
    if(printDist != NULL){
      current.printT(printDist);
      current.printX(printDist,printDim);
      fprintf(printDist,"\n");
    }
  }
  dst = next;
}
