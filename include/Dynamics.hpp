#ifndef _DYNAMICS_
#define _DYNAMICS_

#define ZERO 1e-16

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
    t = other.t;
    for(int i = 0; i < other.DIM; i++) x[i] =  other.x[i];
    phaseDiff = other.phaseDiff;
  }
  
  // Operator
  State& operator=(const State& other){
    if(this != &other){
      t = other.t;
      for(int i = 0; i < other.DIM; i++) x[i] =  other.x[i];
      phaseDiff = other.phaseDiff;
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
  void printT(FILE *fp) const{ fprintf(fp,"%+.12lf ",t+phaseDiff);}
  void printX(FILE *fp, int num) const{
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

class Interval{
private:
  double min, max;
public:
  Interval(){
  }
  Interval(double inMin, double inMax){
    min = inMin;
    max = inMax;
  }
  ~Interval(){
  }

  bool inInterval(double x){
    if(min - x < ZERO && ZERO < max - x) return true;
    return false;
  }

  // Copy constructor
  Interval(const Interval& other){
    min = other.min;
    max = other.max;
  }

  // Operator
  Interval& operator=(const Interval& other){
    if(this != &other){
      min = other.min;
      max = other.max;
    }
    return *this;
  }

  // Accessor
  void setInterval(double inMin, double inMax){min = inMin; max = inMax;};
  void setMin(double inMin){min = inMin;};
  void setMax(double inMax){max = inMax;};
  double getMin()const{return min;};
  double getMax()const{return max;};
};

class Domain{
private:
  Interval* intArray;
  int DIM;
public:
  Domain(){
    DIM = 2;
    intArray = new Interval[DIM]();
    for(int i = 0; i < DIM; i++){
      intArray[i].setMin(-10);
      intArray[i].setMax(+10);
    }
  }
  Domain(int inDim){
    if(inDim == 1){
      fprintf(stderr,"warning: selected dimension is unity. Do you want to declare an interval?\n");
    }
    DIM = inDim;
    intArray = new Interval[DIM]();
    for(int i = 0; i < DIM; i++)
      intArray[i].setInterval(-10,+10);
  }
  ~Domain(){
    delete[] intArray;
  }
  
  bool inDomain (const State& in) const{
    for(int i = 0; i < DIM; i++){
      if(!intArray[i].inInterval(in.getX(i)))
	return false;
    }
    return true;
  }

  // Copy constructor
  Domain(const Domain& other){
    DIM = other.DIM;
    for(int i = 0; i < other.DIM; i++) intArray[i] = other.intArray[i];
  }

  // Operator
  Domain& operator=(const Domain& other){
    if(this != &other){
      for(int i = 0; i < other.DIM; i++) intArray[i] = other.intArray[i];
    }
    return *this;
  }

  // Accessor
  void setInterval(int index, double inMin, double inMax){
    intArray[index].setInterval(inMin,inMax);
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
