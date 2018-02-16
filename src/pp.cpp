#include"Dynamics.hpp"
#include"ODE.hpp"
#include<fstream>

using namespace std;

int main(int argc, char **argv){
  if(argc < 3 || !strcmp("--help",argv[1])){
    fprintf(stderr,"usage: ./pp (counts of maps) (input file name)");
    return -1;
  }  
  FILE *fppoin;
  int mapCount, mapMax, ret;
  Dynamics* dyna;
  ODEwithEvent* ode;

  mapMax = atoi(&argv[1][0]);
  dyna = new Dynamics("pp",argv[2]);
  ode  = new ODEwithEvent(dyna,1e-2);
  fppoin = fopen("pp.poin","w");

  dyna->printProfile(1);

  if(mapMax < 0){
    dyna->tend *= -1;
    dyna->event->reverseEvent();
    ode->reverseStep();
    mapMax *= -1;
  }

  fprintf(stderr,"mapMax = %d\n",mapMax);
  for(mapCount = 0; mapCount <= mapMax; mapCount++){
    fprintf(stderr,"%07.3lf %% (mapCount: %07d )",(double)mapCount*100/(double)mapMax,mapCount);
    dyna->printProfile(0);    
    fprintf(stderr,"\r");
    while(1){
      ret = ode->stepODEwithEvent(1);
      dyna->freshX0();
      
      if(ret){
	dyna->tstart = *(dyna->t);
	dyna->mode = dyna->event->switchMode();
      }else{
	dyna->tstart = 0;
	dyna->tfinal += dyna->tend;
	break;
      }
    }
    
    fprintf(fppoin,"%+015.9lf ",dyna->tfinal+dyna->tend);
    for(int i = 0; i < dyna->DIM; i++) fprintf(fppoin,"%+.10lf ",dyna->x0[i]);
    fprintf(fppoin,"%d \n",dyna->mode);
  }
  fprintf(stderr,"\n");
  
  delete dyna;
  delete ode;
  return 0;
}
