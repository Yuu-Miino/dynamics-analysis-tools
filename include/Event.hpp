#ifndef _EVENT_
#define _EVENT_

#include<iostream>
#include<cstdlib>

class Event{
public:
  unsigned int EFNUM;
  unsigned int EiNUM;
  int *EDIR;
  int *Eflag;
  int EiDIR;
  int rev;
  
  Event(double* x0, int* state);
  ~Event(){
    delete EDIR;
    delete Eflag;
  };
  
  void eventFunc(double* out, double* xin, double tin);
  void deventFuncdt(double* out, double* xin, double tin, double* dx);
  int switchState();
  int initState(double* xin);
  void reverseEvent();
};

#endif
