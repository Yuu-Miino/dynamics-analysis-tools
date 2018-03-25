#ifndef _HYBRID_SYSTEM_
#define _HYBRID_SYSTEM_

#include "pwlDuffing.hpp"

#define MODE_NUM  4

class HybridSystem{
private:
  ModeProperty* mp[MODE_NUM];
  int mode;
public:
  HybridSystem(int inMode){
    mode = inMode;
    for(int i = 0; i < MODE_NUM; i++) mp[i] = new pwlDuffingMode(i);
  }
  ~HybridSystem(){
    for(int i = 0; i < MODE_NUM; i++) delete mp[i];
  }

  // Poincare map
  void map(const State& inInit, const Parameter& inPara, 
	   const double tfinal, State& dst, FILE* printDist=NULL){
    State init(inInit);
    HSODEsolver odeSolver("RK4",((tfinal-init.getT()) > ZERO ? 1.0:-1.0)*1e-2);
    StateWithEvent swe(init.getDIM());
    
    // If out of domain
    if(!mp[mode]->inDomain(init, inPara)){
      for(int i = 0; i < MODE_NUM; i++){
	if(mp[mode]->inDomain(init, inPara)){ mode = i; break;}
      }
    }

    // Main loop
    while(1){
      odeSolver.runHSODEsolver(*mp[mode], init, inPara, tfinal, swe, printDist);
      init = *swe.state;
      if(swe.eventFlag){
	mode = mp[mode]->modeDist(swe.eventIndex);
	swe.resetFlags();
      }else{
	init.setT(0);
	init.addPhaseDiff(tfinal);
	break;
      }
    }
    dst = init;
  }

  // Accessor
  int getMode() const{return mode;}
};

#endif
