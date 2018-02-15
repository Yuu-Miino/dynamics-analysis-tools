#include "Dynamics.hpp"
#include <fstream>
#include <iostream>

Dynamics::Dynamics(string namein, char* in){
  name = namein;    // Name of instance
  ifstream ifs(in); // Name of input file

  // Dimensions
  DIM_state = 2;    // Number of state variables
  DIM       = 2;    // Number of all variables

  x0 = new double[DIM](); // Initial condition
  x  = new double[DIM](); // Current state(s)
  
  // Load parameters and initial conditions from input file
  if(!ifs.fail()){
    ifs>>
      k>>b0>>b>>sl>>sr>>
      x0[0]>>x0[1];
  }else{
    cerr<<"error: file name incorrect ["<<in<<"]"<<endl;
    exit(-1);
  }
  
  // time interval
  t  = new double; // Current time (mod tend)
  tstart = 0.0;    // Start time
  tend   = 2.0*PI; // End time
  tfinal = tstart; // Temporary time to store the total time

  // event detection
  event = new Event(x0, &state); // refer to Event.cpp
}

double Dynamics::f1(double xin){return (xin*sr); }
double Dynamics::f2(double xin){return (xin/sl); }

void Dynamics::func(double* out, double* in, double tin){
  double fx=0, cost = cos(tin);
  switch(state){
  case 1:
    fx = f1(in[0]);
    break;
  case 2:
    fx = f2(in[0]);
    break;
  default:
    break;
  }
  out[0] = in[1];
  out[1] = -k*in[1]- fx + b*cost +b0;
}

void Dynamics::freshX0(){
  for(unsigned int i = 0; i < DIM_state; i++){
    x0[i] = x[i];
  }
}

void Dynamics::printProfile(int retFlag){
  fprintf(stderr,"%10s: %d-dim: [%+.5lf:%+.5lf] [x0, y0] = [%+.5lf %+.5lf], [k b0 b sl sr] = [%+.5lf %+.5lf %+.5lf %+.5lf %+.5lf] ",
	  name.c_str(),DIM,tstart,tend,x0[0],x0[1],k,b0,b,sl,sr);
  if(retFlag == 1) fprintf(stderr,"\n");
}
void Dynamics::printStd(const char* in){
  char name[20] = "dout.dat";
  if(in != NULL){
    strcpy(name,in);
  }
  ofstream ofs(name);
  ofs<<k<<" "<<b0<<" "<<b<<" "<<sl<<" "<<sr<<endl<<x0[0]<<" "<<x0[1]<<endl;
}
