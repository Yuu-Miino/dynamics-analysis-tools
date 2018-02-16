#ifndef _ODE_
#define _ODE_

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<stdio.h>
#include<math.h>

#include"Dynamics.hpp"
#include"Event.hpp"

#define ZERO 1e-10
#define EEPS 1e-8

#define DIM_MAX 30

class ODE
{
public:
  Dynamics *dyna;
  double h, initStep;

  ODE(Dynamics* in, double step){
    dyna = in;
    h = initStep = step;
  }
  ~ODE(){
  }
  
  void stepRK();
  void stepODE(int printFlag);
  void reverseStep();
};

class ODEwithEvent:public ODE
{
private:
  unsigned int teventFlag;
  unsigned int crossFlag;
public:
  double tevent;
  ODEwithEvent(Dynamics* in, double step):ODE(in, step){
    tevent = 0; teventFlag = 0; crossFlag = 0;
  }
  ~ODEwithEvent(){
  }

  int stepODEwithEvent(int printFlag);
  void eventDetect(double* xold, double told, double* tmp0, double* tmp1, int* index, int* tmp_dir);
};

#endif
