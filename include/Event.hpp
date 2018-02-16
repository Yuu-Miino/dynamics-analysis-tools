#ifndef _EVENT_
#define _EVENT_

#include<iostream>
#include<cstdlib>

#define EFNUM_MAX 10

class Event{
public:
  unsigned int EFNUM;
  unsigned int EiNUM;
  int *EDIR;
  int *Eflag;
  int EiDIR;
  int rev;
  
  Event(double* x0, int* mode);
  ~Event(){
    delete EDIR;
    delete Eflag;
  };
  
  void eventFunc(double* out, double* xin, double tin);
  void deventFuncdt(double* out, double* xin, double tin, double* dx);
  int switchMode();
  int initMode(double* xin);
  void reverseEvent();
};

#endif
