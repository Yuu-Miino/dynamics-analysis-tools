#include "HSODEsolver.hpp"
// Class member functions
void HSODEsolver::runHSODEsolver(ModeProperty& mode, 
				 const State& init, const Parameter& para, double tfinal, 
				 StateWithEvent& out, FILE* printDist, int printDim)
{
  int dim   = init.getDIM();
  if(printDim == -1) printDim = dim;
  int efnum = mode.getNumOfEventFunc();

  double ef0[efnum], ef1[efnum], defdt[efnum];;
  double dxdt[dim];
  State current(init), next(init);
  bool FINISH = false, grazeFlag = false, newtonFlag = true, eventFIN = false;
  int index, dir;

  // if print the orbit
  if(printDist != NULL){
    init.printT(printDist); init.printX(printDist,printDim);
    fprintf(printDist,"%d\n",mode.getMode());
  }

  // main loop
  while(!FINISH){
    index = -1; dir = 0;

    for(unsigned int i = 0; i < efnum; i++){
      ef0[i]=0.0; ef1[i]=0.0;
    }
    stepODEsolver(*mode.dyna,current,para,next);
    if(!teventFlag){
      eventDetect(mode, current, ef0, next, ef1, para, &index, &dir);
    }else{
      teventFlag = false;
    }
    if((next.getT() - tfinal)*(h/fabs(h)) > ZERO){
      h = tfinal - current.getT();
      stepODEsolver(*mode.dyna,current,para,next);
      h = initStep;

      *out.state  = next;
      FINISH = true;

      eventDetect(mode, current, ef0, next, ef1, para, &index, &dir);
      if(index >= 0){// if an event is detected
	FINISH = false;
      }
    }

    if(index >= 0){
      if(newtonFlag){// Newton's method
	while(fabs(ef1[index]) > EEPS){
	  mode.dyna->ode(dxdt, current, para);
	  mode.dEFdt(defdt, current, para, dxdt);
	  if(fabs(defdt[index]) < 1e-3){
	    grazeFlag = true;
	    break;
	  }
	  h = -(ef1[index])/defdt[index];
	  
	  stepODEsolver(*mode.dyna,current,para,next);
	  current = next;
	  mode.eventFunction(ef1, current, para);
	}
	if(!grazeFlag) eventFIN = true;

      }else{ // Bisection method
	if(fabs(ef1[index]) > EEPS){
	  h /= 2.0;
	}else{
	  eventFIN = true;
	}
      }
    
      if(grazeFlag){// Newton's method -> Bisection method
	grazeFlag = false;
	newtonFlag = false;
	h = initStep/2.0;
      }

      if(eventFIN){
	eventFIN   = false;
	newtonFlag = true;
	teventFlag = true;

	h = initStep;

	if(mode.getEventFlag(index)){
	  // Set output
	  *out.state  = next;
	  out.eventFlag  = true;
	  out.eventIndex = index;
	  out.eventDir   = dir;
	  FINISH = true;
	}
      }
    }
    
    // if print the orbit
    if((printDist != NULL) && newtonFlag){
      next.printT(printDist); next.printX(printDist,printDim);
      fprintf(printDist,"%d\n",mode.getMode());
    }
    current = next;
  }// for "while(!FINISH)"
}

void HSODEsolver::eventDetect(ModeProperty& mode, 
			       const State& current, double *ef0, 
			       const State& next, double *ef1, const Parameter& para,
			       int* index, int* dir){
  mode.eventFunction(ef0, current, para);
  mode.eventFunction(ef1, next, para);

  for(unsigned int i = 0; i < mode.getNumOfEventFunc(); i++){
    switch (mode.getEventDir(i)) {
    case 0:
      if((ef1[i] > ZERO) && (ZERO > ef0[i])){
	*index = i; *dir = 1;
      }else if(ef0[i] > ZERO && ZERO > ef1[i]){
	*index = i; *dir = -1;
      }
      break;
    case 1:
      if((ef1[i] > ZERO) && (ZERO > ef0[i])) {
	*index = i; *dir =  1;
      }
      break;
    case -1:
      if((ef0[i] > ZERO) && (ZERO > ef1[i])) {
	*index = i; *dir = -1;
      }
      break;
    default:
      fprintf(stderr,"Error: undefined EDIR = %d\n",mode.getEventDir(i));
      exit(1);
    }
  }
}
