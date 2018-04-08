#ifndef _HS_ODE_SOLVER_
#define _HS_ODE_SOLVER_

#include "ODEsolver.hpp"

#define EEPS 1e-10

class StateWithEvent{
public:
  State* state;
  bool eventFlag;
  int  eventIndex;
  int  eventDir;
  
  StateWithEvent(int dim){
    state = new State(dim);
    eventFlag = false;
    eventIndex = -1;
    eventDir = 0;
  }
  ~StateWithEvent(){
    delete state;
  }
  void resetFlags(){
    eventFlag = false;
    eventIndex = -1;
    eventDir = 0;
  }
};

class ModeProperty{
protected:
  unsigned int modeID;
  unsigned int numOfEventFunc;
  int*  eventDir;
  bool* eventFlag;
public:
  Dynamics* dyna;

  ModeProperty(unsigned int inID, unsigned int inNum){
    modeID = inID;
    numOfEventFunc = inNum;
    eventDir  = new int[numOfEventFunc]();
    eventFlag = new bool[numOfEventFunc]();
  }
  virtual ~ModeProperty(){
    delete(eventDir);
    delete(eventFlag);
  }

  virtual void eventFunction(double* out, const State& state, const Parameter& para)=0;
  virtual void dEFdt(double* out, const State& state, const Parameter& para, double* dxdt)=0;
  virtual unsigned int modeDist(int eventIndex)=0;
  virtual bool inDomain(const State& state, const Parameter& para)=0;

  // Accessor
  unsigned int getMode(){return modeID;}
  int  getNumOfEventFunc(){return numOfEventFunc;}
  int  getEventDir(int index){return eventDir[index];}
  bool getEventFlag(int index){return eventFlag[index];}
};

class HSODEsolver:public ODEsolver
{
private:
  bool teventFlag;
    
public:
  HSODEsolver(std::string solverName, double step):ODEsolver(solverName, step){
    teventFlag = false;
  }
  ~HSODEsolver(){};

  bool runHSODEsolver(ModeProperty& mode, const Domain& domain,
		      const State& init, const Parameter& para, double tfinal,
		      StateWithEvent& out, FILE* printDist=NULL, int printDim=-1);
  void eventDetect(ModeProperty& mode, 
		   const State& current, double *ef0, 
		   const State& next,    double *ef1, const Parameter& para,
		   int* index, int* dir);
};

#endif
