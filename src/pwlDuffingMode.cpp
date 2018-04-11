#include "pwlDuffing.hpp"
// Internal function
unsigned int retNum(unsigned int inID){
  switch(inID){
  case 0: return 1;
  case 1: return 2;
  case 2: return 1;
  case 3: return 2;
  default: 
    fprintf(stderr,"Error: undefined modeID = %d in retNumOfEF\n",inID);
    exit(1);
  }
}

// Class member functions
pwlDuffingMode::pwlDuffingMode(unsigned int inID):ModeProperty(inID,retNum(inID)){
  dyna  = new pwlDuffing(inID);

  for(int i = 0; i < numOfEventFunc; i++) eventFlag[i] = true;
  switch(modeID){
  case 0: eventDir[0] =  1; break;
  case 1: eventDir[0] = -1; eventDir[1] = 1; 
    //eventFlag[1] = false; // For Grazing bifurcation
    break;
  case 2: eventDir[0] = -1; break;
  case 3: eventDir[0] = -1; eventDir[1] = 1; break;
  default: 
    fprintf(stderr,"Error: undefined modeID = %d in pwlDuffingEvent(inID)\n",inID);
    exit(1);
  }
}

pwlDuffingMode::~pwlDuffingMode(){
  delete(dyna);
}

bool pwlDuffingMode::inDomain(const State& state, const Parameter& para){
  double x = state.getX(0);
  double th0, th1, th2, th3;
  th0 = para.getValue(3);
  th1 = para.getValue(4);
  th2 = para.getValue(5);
  th3 = para.getValue(6);


  switch(modeID){
  case 0: 
    if(x - (th1-th0)*3.0/8.0 > ZERO) return false;
    break;
  case 1: 
    if(x - (th1-th0)*3.0/8.0 < ZERO || x - (th1-th2)*3.0/8.0 > ZERO) return false; 
    break;
  case 2:
    if(x - (th3-th2)*3.0/8.0 < ZERO) return false;
    break;
  case 3:
    if(x - (th3-th2)*3.0/8.0 > ZERO || x - (th3-th0)*3.0/8.0 < ZERO) return false;
    break;
  default: 
    fprintf(stderr,"Error: undefined modeID = %d in pwlDuffingEvent::inDomain\n",modeID);
    exit(1);
  }
  return true;
}

void pwlDuffingMode::eventFunction(double* EF, const State& state, const Parameter& para){
  double x = state.getX(0);
  double th0, th1, th2, th3;
  th0 = para.getValue(3);
  th1 = para.getValue(4);
  th2 = para.getValue(5);
  th3 = para.getValue(6);

  switch(modeID){
  case 0: 
    EF[0] = x - (th1-th0)*3.0/8.0;
    break;
  case 1: 
    EF[0] = x - (th1-th0)*3.0/8.0;
    EF[1] = x - (th1-th2)*3.0/8.0;
    break;
  case 2: 
    EF[0] = x - (th3-th2)*3.0/8.0;
    break;
  case 3: 
    EF[0] = x - (th3-th0)*3.0/8.0;
    EF[1] = x - (th3-th2)*3.0/8.0;
    break;
  default: 
    fprintf(stderr,"Error: undefined modeID = %d in pwlDuffingEvent::eventFunction\n",modeID);
    exit(1);
  } 
}

void pwlDuffingMode::dEFdt(double* dEFdt, const State& state, const Parameter& para, double* dxdt){
  switch(modeID){
  case 0: 
    dEFdt[0] = dxdt[0];
    break;
  case 1: 
    dEFdt[0] = dxdt[0];
    dEFdt[1] = dxdt[0];
    break;
  case 2: 
    dEFdt[0] = dxdt[0];
    break;
  case 3: 
    dEFdt[0] = dxdt[0];
    dEFdt[1] = dxdt[0];
    break;
  default: 
    fprintf(stderr,"Error: undefined modeID = %d in pwlDuffingEvent::dEFdt\n",modeID);
    exit(1);
  } 
}

unsigned int pwlDuffingMode::modeDist(int eventIndex){
  switch(modeID){
  case 0: return 1;
  case 1: return (eventIndex==0)?0:2;
  case 2: return 3;
  case 3: return (eventIndex==0)?0:2;
  default: 
    fprintf(stderr,"Error: undefined modeID = %d in pwlDuffingEvent::modeDist\n",modeID);
    exit(1);
  }
}
