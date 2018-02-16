#include "ODE.hpp"

void ODE::stepRK(){
  double k0[DIM_MAX], k1[DIM_MAX], k2[DIM_MAX], k3[DIM_MAX], tmp[DIM_MAX];
  
  //exe K0
  for( unsigned int i=0; i<dyna->DIM; i++){
    tmp[i] = dyna->x[i];
  }
  dyna->func(k0, tmp, *(dyna->t));
  
  //exeK1
  for( unsigned int i=0; i<dyna->DIM; i++){
    tmp[i] = dyna->x[i] + k0[i]*(h/2.0);
  }
  dyna->func(k1, tmp, *(dyna->t)+0.5*h);
  
  //exeK2
  for( unsigned int i=0; i<dyna->DIM; i++){
    tmp[i] = dyna->x[i] + k1[i]*(h/2.0);
  }
  dyna->func(k2, tmp, *(dyna->t)+0.5*h);
  
  //exeK3
  for( unsigned int i=0; i<dyna->DIM; i++){
    tmp[i] = dyna->x[i] + k2[i]*h;
  }
  dyna->func(k3, tmp, *(dyna->t)+h);
  
  //exe next point
  for( unsigned int i=0; i<dyna->DIM; i++){
    dyna->x[i] += (h/6.0)*(k0[i] + 2.0*k1[i] + 2.0*k2[i] + k3[i]);
  }
  *(dyna->t) += h;
}

void ODE::stepODE(int printFlag)
{
  double xold[DIM_MAX]; double told;
  unsigned int FINISH = 0;

  while(!FINISH){
    for(unsigned int i = 0; i < dyna->DIM; i++){xold[i]=dyna->x[i];}
    told = *(dyna->t);
    stepRK();
    if(fabs(*(dyna->t)) - fabs(dyna->tend) > 1e-12){
      for(unsigned int i=0;i<dyna->DIM;i++){dyna->x[i]=xold[i];}
      *(dyna->t) = told;
      h = dyna->tend - told;
      stepRK();
      FINISH = 1;
      h = initStep;
    }
    if(printFlag){
      printf("%+.12lf ", dyna->tfinal + *(dyna->t));
      for(unsigned int i = 0; i < dyna->DIM; i++)
	printf("%+.12lf ",dyna->x[i]);
      printf("\n");
    }
  }
}

int ODEwithEvent::stepODEwithEvent(int printFlag)
{
  double tmp0[EFNUM_MAX], tmp1[EFNUM_MAX], dx[DIM_MAX], df[DIM_MAX], xold[DIM_MAX], told;
  int index, tmp_dir, grazeFlag = 0, newtonFlag = 1;
  int FINISH = 0, ret = 0, eventFIN = 0;

  // initialize
  for(unsigned int i=0;i<dyna->DIM;i++){ dyna->x[i]=dyna->x0[i]; }
  *(dyna->t) = dyna->tstart;

  // if print the orbit in stdout
  if(printFlag){
    printf("%+.12lf ", dyna->tfinal + *(dyna->t));
    for(unsigned int i = 0; i < dyna->DIM_state; i++)
      printf("%+.12lf ",dyna->x[i]);
    printf("%d \n",dyna->mode);
  }

  // main loop
  while(!FINISH){
    // store the previous values
    for(unsigned int i=0;i<dyna->DIM;i++){ xold[i]=dyna->x[i];} told = *(dyna->t);

    // initialize indexes
    index = -1; tmp_dir = 0;
    for(unsigned int i = 0; i < dyna->event->EFNUM; i++){
      tmp0[i]=0.0; tmp1[i]=0.0;
    }

    stepRK(); // take a step of RK method
    
    // if not the first time after the event
    if(teventFlag != 1){
      eventDetect(xold, told, tmp0, tmp1, &index, &tmp_dir);
    }else{
      h = initStep;
      teventFlag = 0;
    }

    // if t passes through tend
    if(fabs(*(dyna->t)) - fabs(dyna->tend) > ZERO){
      // restore the previous values to the current state
      for(unsigned int i=0;i<dyna->DIM;i++){dyna->x[i]=xold[i];}
      *(dyna->t) = told;
      
      // set the rest time as the step value 
      h = dyna->tend - told;

      stepRK();     // take a step of RK method

      FINISH = 1;

      // never events
	index = -1; // initialize
	eventDetect(xold, told, tmp0, tmp1, &index, &tmp_dir);
	
	if(index >= 0){// if an event is detected
	  crossFlag = 1;
	  FINISH = 0;
	}else{
	  h = initStep; // initialize step
	}
    }

    // if an event
    if(index >= 0){
      if(tevent < ZERO) tevent = *(dyna->t);

      // for the event detection
      if(newtonFlag){// Newton's method
	while(fabs(tmp1[index]) > EEPS){
	  dyna->func(dx, dyna->x, *(dyna->t));
	  dyna->event->deventFuncdt(df, dyna->x, dyna->tfinal+*(dyna->t), dx);
	  if(fabs(df[index]) < 1e-3){
	    grazeFlag = 1;
	    break;
	  }
	  h = -(tmp1[index])/df[index];
	  stepRK();
	  dyna->event->eventFunc(tmp1, dyna->x, dyna->tfinal+*(dyna->t));
	}
	if(!grazeFlag) eventFIN = 1;
      }else{ // Binary method
	if(fabs(tmp1[index]) > EEPS){
	  h /= 2.0;
	  for(unsigned int i=0;i<dyna->DIM;i++){dyna->x[i]=xold[i];}
	  *(dyna->t) = told;
	}else{
	  eventFIN = 1;
	}
      }
      
      if(grazeFlag){// Newton's method -> Binary method
	newtonFlag = 0;
	
	// restore the previous values to the current state
	// set the half step as h
	h = initStep/2.0;
	for(unsigned int i=0;i<dyna->DIM;i++){dyna->x[i]=xold[i];}
	*(dyna->t) = told;
		
	grazeFlag = 0; // reset flag
      }
      if(eventFIN){
	if(fabs(*(dyna->t) - dyna->tend) < ZERO){
	  h = initStep;
	  ret = 0;
	}else{
	  teventFlag = 1; //reset flag
	  h = tevent - *(dyna->t);
	  ret = 1;
	}

	eventFIN = 0;   //reset flag
	newtonFlag = 1; //reset flag

	tevent = 0;
	
	// set indexes
	dyna->event->EiNUM = index;
	dyna->event->EiDIR = tmp_dir;

	if(dyna->event->Eflag[index] == 1){
	  FINISH = 1;
	}
      }
    }
    
    if(printFlag && newtonFlag){
      printf("%+.12lf ", dyna->tfinal+*(dyna->t));
      for(unsigned int i = 0; i < dyna->DIM_state; i++)
	printf("%+.12lf ",dyna->x[i]);
      printf("%d \n",dyna->mode);
    }
  }// for "while(!FINISH)"
  
  return ret;
}

void ODE::reverseStep(){
  h = initStep = -h;
}

void ODEwithEvent::eventDetect(double* xold, double told, double* tmp0, double* tmp1, int* index, int* tmp_dir){
  dyna->event->eventFunc(tmp0, xold, dyna->tfinal+told);
  dyna->event->eventFunc(tmp1, dyna->x, dyna->tfinal+*(dyna->t));
  for(unsigned int i=0;i<dyna->event->EFNUM;i++){
    switch (dyna->event->EDIR[i]) {
    case 0:
      if((tmp1[i] > ZERO) && (ZERO > tmp0[i])){
	*index = i;*tmp_dir = 1;
      }else if(tmp0[i] > ZERO && ZERO > tmp1[i]){
	*index = i;*tmp_dir = -1;
      }
      break;
    case 1:
      if((tmp1[i] > ZERO) && (ZERO > tmp0[i])) {
	*index = i;*tmp_dir = dyna->event->EDIR[i];
      }
      break;
    case -1:
      if((tmp0[i] > ZERO) && (ZERO > tmp1[i])) {*index = i;*tmp_dir = dyna->event->EDIR[i];}
      break;
    default:
      std::cerr<<"No EDIR set"<<std::endl;
      std::exit(-1);
    }
  }
}
