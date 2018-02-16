#include"Event.hpp"

Event::Event(double* x0, int* mode){

  *mode = initMode(x0);
  EFNUM = 2;
  EDIR  = new int[EFNUM];
  Eflag = new int[EFNUM];
  EDIR[0] = -1; EDIR[1] = 1;
  Eflag[0] = 1; Eflag[1] = 1;
  EiNUM = -1;
  rev = 1;
}

void Event::eventFunc(double* out, double* xin, double tin){
  out[0] = xin[0];
  out[1] = xin[0];
}

void Event::deventFuncdt(double* out, double* xin, double tin, double* dx){
  out[0] = dx[0];
  out[1] = dx[0];
}

int Event::switchMode(){
  if(EiNUM == 0){
    return (rev==1?2:1);
  }else if(EiNUM == 1){
    return (rev==1?1:2);
  }else{
    fprintf(stderr,"Error: EiNUM\n");
    exit(-1);
  }
}

int Event::initMode(double* xin){
  if(xin[0] > 1e-9) return 1;
  return 2;
}

void Event::reverseEvent(){
  rev *= -1;
  for(unsigned int i = 0; i < EFNUM; i++) EDIR[i] *= -1;
}

