#include "Dynamics.hpp"
#include <fstream>
#include <iostream>

Dynamics::Dynamics(string namein, char* in){
  
  name = namein;
  ifstream ifs(in);
  DIM_state = 2;
  DIM = 2;

  x0 = new double[DIM]();
  x = new double[DIM]();
  t = new double;
  
  if(!ifs.fail()){
    ifs>>k>>b0>>b>>sl>>sr>>x0[0]>>x0[1];
  }else if(in == NULL){
    k = 0.2;  b  = 0.1;  b0 = 0.1;
    x0[0] = 0.3238; x0[1] = -1.292;
    sl = 3; sr = 3;
  }else{
    cerr<<"error: file name incorrect"<<in<<endl;
    exit(-1);
  }
  //s = 3;
  tstart = 0.0; tend = 2.0*PI; tfinal = tstart;

  event = new Event(x0, &state);

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
