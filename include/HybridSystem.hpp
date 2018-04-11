#ifndef _HYBRID_SYSTEM_
#define _HYBRID_SYSTEM_

#include "pwlDuffing.hpp"
#include <Eigen/Dense>

#include <iostream>

using namespace std;
using namespace Eigen;

class HybridSystem{
private:
  ModeProperty* mp[MODE_NUM];
public:
  HybridSystem(){
    for(int i = 0; i < MODE_NUM; i++) mp[i] = new pwlDuffingMode(i);
  }
  ~HybridSystem(){
    for(int i = 0; i < MODE_NUM; i++) delete mp[i];
  }

  // Additional function
  void getXth(FILE *fp, const Parameter& para){
    double th0, th1, th2, th3;
    th0 = para.getValue(3);
    th1 = para.getValue(4);
    th2 = para.getValue(5);
    th3 = para.getValue(6);
    
    fprintf(fp,"%lf %lf %lf %lf ",
	    (th3-th0)*3.0/8.0,(th1-th0)*3.0/8.0,(th3-th2)*3.0/8.0,(th1-th2)*3.0/8.0);
  };

  void eventFunction(double* EF, const State& state, const Parameter& para, const int mode){
    mp[mode]->eventFunction(EF, state, para);
  }
  void dxdt(double* dxdt, const State& in, const Parameter& para, const int mode){
    mp[mode]->dyna->ode(dxdt, in, para);
  }

  // Poincare map
  bool map(const State& inInit, const Parameter& inPara, int& mode,
	   double tfinal, State& dst, FILE* printDist=NULL){
    // Variables definition
    State init(inInit);
    bool isDivergent = false;
    HSODEsolver odeSolver("RK4",((tfinal-init.getT()) > ZERO ? 1.0:-1.0)*1e-3);
    StateWithEvent swe(init.getDIM());
    Domain domain(2);
    domain.setInterval(0,-5,5);
    domain.setInterval(1,-5,5);
    
    // If out of domain
    if(!mp[mode]->inDomain(init, inPara)){
      for(int i = 0; i < MODE_NUM; i++){
	if(mp[mode]->inDomain(init, inPara)){ mode = i; break;}
      }
    }

    // Main loop
    while(!isDivergent){
      isDivergent = odeSolver.runHSODEsolver(*mp[mode], domain, init, inPara, tfinal, swe, printDist, PRINT_DIM);
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
    return isDivergent;
  }
  void jacobian(MatrixXd& jac, MatrixXd& jacP, int period,
		const State& inInit, const Parameter& inPara, int& mode,
		double tfinal, State& dst){
    // Variables definition    
    bool jacPflag = true;
    MatrixXd tmpJac(JAC_MAT_DIM,JAC_MAT_DIM), tmpJacP(JAC_MAT_DIM,JAC_PARA_DIM);
    
    State init(inInit);
    HSODEsolver odeSolver("RK4",((tfinal-init.getT()) > ZERO ? 1.0:-1.0)*1e-2);
    StateWithEvent swe(init.getDIM());
    Domain domain(2);

    double dxdt[init.getDIM()];
    
    jac = MatrixXd::Identity(JAC_MAT_DIM,JAC_MAT_DIM);

    // If out of domain
    if(!mp[mode]->inDomain(init, inPara)){
      for(int i = 0; i < MODE_NUM; i++){
	if(mp[mode]->inDomain(init, inPara)){ mode = i; break;}
      }
    }

    // Main loop
    for(int periodIter = 0; periodIter < period; periodIter++){
      while(1){
	for(int i = JAC_MAT_DIM; i < init.getDIM(); i++){
	  init.setX(i,0);
	  if(i == JAC_MAT_DIM || i == JAC_MAT_DIM + 4 || i == JAC_MAT_DIM + 8)
	    init.setX(i,1);
	}
	odeSolver.runHSODEsolver(*mp[mode], domain, init, inPara, tfinal, swe);
	init = *swe.state;

	mp[mode]->dyna->ode(dxdt, init, inPara);
	if(swe.eventFlag){
	  // Calc Jacobian matrix
	  for(int row = 0; row < JAC_MAT_DIM; row++)
	    for(int col = 0; col < JAC_MAT_DIM; col++){
	      tmpJac(row,col) = 
		init.getX(JAC_MAT_DIM*(col + 1) + row) - 
		init.getX(JAC_MAT_DIM*(col + 1) + EVENT_STATE_INDEX) * dxdt[row]/dxdt[EVENT_STATE_INDEX];
	    }
	  for(int row = 0; row < JAC_MAT_DIM; row++)
	    for(int col = 0; col < JAC_PARA_DIM; col++)
	      tmpJacP(row,col) = 
		init.getX(JAC_MAT_DIM*(col + JAC_MAT_DIM + 1) + row) - 
		init.getX(JAC_MAT_DIM*(col + JAC_MAT_DIM + 1) + EVENT_STATE_INDEX) * dxdt[row]/dxdt[EVENT_STATE_INDEX];
	  if(jacPflag){ jacP = tmpJacP; jacPflag = false;}
	  else{ jacP = tmpJac*jacP + tmpJacP;}
	  jac = tmpJac * jac;
	  
	  // Mode transition
	  mode = mp[mode]->modeDist(swe.eventIndex);
	  swe.resetFlags();
	}else{
	  // Calc Jacobian matrix
	  for(int row = 0; row < JAC_MAT_DIM; row++)
	    for(int col = 0; col < JAC_MAT_DIM; col++)
	      tmpJac(row,col) = 
		init.getX(JAC_MAT_DIM*(col + 1) + row) - 
		init.getX(JAC_MAT_DIM*(col + 1) + TIME_STATE_INDEX) * dxdt[row]/dxdt[TIME_STATE_INDEX];
	  for(int row = 0; row < JAC_MAT_DIM; row++)
	    for(int col = 0; col < JAC_PARA_DIM; col++)
	      tmpJacP(row,col) = 
		init.getX(JAC_MAT_DIM*(col + JAC_MAT_DIM + 1) + row) - 
		init.getX(JAC_MAT_DIM*(col + JAC_MAT_DIM + 1) + TIME_STATE_INDEX) * dxdt[row]/dxdt[TIME_STATE_INDEX];
	  if(jacPflag){ jacP = tmpJacP; jacPflag = false;}
	  else{ jacP = tmpJac*jacP + tmpJacP;}
	  jac = tmpJac * jac;
	  
	  // Reset
	  init.setT(0);
	  init.setX(TIME_STATE_INDEX,0);
	  break;
	}
      }
      init.setPhaseDiff(0);
      dst = init;
    }
  }

};

#endif
