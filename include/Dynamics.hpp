#ifndef _DYNA_
#define _DYNA_

#include "Event.hpp"
#include <math.h>
#include <string>

#define PI 3.14159265359

using namespace std;

class Dynamics{
private:
  std::string name;
public:
  int DIM_state;
  int DIM;
  double *x0, *x, *t;
  double tstart, tfinal, tend;

  // Parameter set
  double k, b, b0, sl, sr;
  double fx1, fx2;
  int state;

  // Event
  Event *event;

  // Constructor
  Dynamics(string namein, char* in);    
  // Destructor
  ~Dynamics(){
    delete x0;
    delete x;
    delete t;
    delete event;
  }   

  // df/dt definition
  double f1(double xin);
  double f2(double xin);

  void func(double* out, double* in, double tin);

  void freshX0();

  // Print profile
  void printProfile(int retFlag);
  void printStd(const char* in=NULL);
};

#endif
