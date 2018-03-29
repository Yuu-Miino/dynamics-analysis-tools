#ifndef _DYNAMICS_
#define _DYNAMICS_

#include <fstream>
#include <stdio.h>
#include <stdlib.h>

class State{
private:
  int DIM;
  double t, *x;
  double phaseDiff;
public:
  State(int dim){
    DIM = dim;
    x = new double[dim]();
    phaseDiff = 0;
  }
  ~State(){
    delete x;
  }

  // Copy constructor
  State(const State& other){
    DIM = other.DIM;
    x = new double[DIM]();
    setT(other.t);
    for(int i = 0; i < other.DIM; i++) setX(i, other.getX(i));
    setPhaseDiff(other.phaseDiff);
  }
  
  // Operator
  State& operator=(const State& other){
    if(this != &other){
      setT(other.t);
      for(int i = 0; i < other.DIM; i++) setX(i, other.getX(i));
      setPhaseDiff(other.phaseDiff);
    }
    return *this;
  }

  // Accessor
  void setT(double val){t = val;}
  void setPhaseDiff(double val){phaseDiff = val;}
  void setX(int index, double val){x[index] = val;}
  void setX(const double *xin){
    for(int i = 0; i < DIM; i++) x[i] = xin[i];
  }
  void setXfromFile(std::ifstream& ifs, int num){
    for(int i = 0; i < num; i++) ifs>>x[i];
  }
  
  double getT()        const{return t;}
  double getPhaseDiff()const{return phaseDiff;}
  double getX(int index)const{return x[index];}
  void getX(double *xout)const{
    for(int i = 0; i < DIM; i++) xout[i] = x[i];
  }
  int getDIM() const{return DIM;}

  void addX(int index, double val){x[index] += val;}
  void addPhaseDiff(double val){ phaseDiff+= val;}

  // Printer
  void printT(FILE *fp){ fprintf(fp,"%+.12lf ",t+phaseDiff);}
  void printX(FILE *fp, int num){
    for(int i = 0; i < num; i++) fprintf(fp,"%+.12lf ",x[i]);
  }
  
};

class Parameter{
private:
  int DIM;
  double *value;
public:
  Parameter(int dim){
    DIM = dim;
    value = new double[dim]();
  }
  ~Parameter(){
    delete value;
  }
  
  // Copy constructor
  Parameter(const Parameter& other){
    DIM = other.DIM;
    value = new double[DIM]();
    for(int i = 0; i < other.DIM; i++) setValue(i, other.getValue(i));
  }

  // Operator
  Parameter& operator=(const Parameter& other){
    if(this != &other){
      for(int i = 0; i < other.DIM; i++) setValue(i, other.getValue(i));
    }
    return *this;
  }

  // Accessor
  void setValue(int index, double val){value[index] = val;}
  void setValue(const double *vin){
    for(int i = 0; i < DIM; i++) value[i] = vin[i];
  }
  void setValueFromFile(std::ifstream& ifs){
    for(int i = 0; i < DIM; i++) ifs>>value[i];
  }
  double getValue(int index) const{return value[index];}
  void getValue(double *vout) const{
    for(int i = 0; i < DIM; i++) vout[i] = value[i];
  }
  void addValue(int index, double val){ value[index] += val;}

  // Printer
  void printValue(FILE *fp) const{
    for(int i = 0; i < DIM; i++) fprintf(fp,"%+.12lf ",value[i]);
  }
};

class Dynamics{
public:
  Dynamics(){}
  virtual ~Dynamics(){}   

  // df/dt definition
  virtual void ode(double* dxdt, const State& state, const Parameter& para)=0;
};

#endif
