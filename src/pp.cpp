#include"Dynamics.hpp"
#include"ODE.hpp"
#include<fstream>

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cerr<<"error: set the numbers of maps."<<endl;
    return -1;
  }  
  ofstream ofs("poin.dat");
  ofstream ofs2("rMap.dat");
  ofstream ofsP1("rMapP1.dat");
  ofstream ofsP2("rMapP2.dat");
  ofstream ofsP3("rMapP3.dat");
  ofstream ofsP4("rMapP4.dat");
  ofstream ofsP5("rMapP5.dat");
  ofstream ofsP6("rMapP6.dat");

  int mapCount, mapMax;
  mapMax = atoi(&argv[1][0]);

  Dynamics* dyna;
  ODEwithEvent* ode;

  dyna = new Dynamics("pp",argv[2]);
  ode  = new ODEwithEvent(dyna,1e-2);

  int ret, eventCount=0;
  double tmpX;

  dyna->printProfile(1);
  if(mapMax < 0){
    dyna->tend *= -1;
    dyna->event->reverseEvent();
    ode->reverseStep();
    mapMax *= -1;
  }

  cerr<<"mapMax = "<<mapMax<<endl;
  tmpX = dyna->x0[0];
  for(mapCount = 0; mapCount < mapMax; mapCount++){
    fprintf(stderr,"%3.3lf %% (mapCount: %d )",(double)mapCount*100/(double)mapMax,mapCount);
    dyna->printProfile(0);    
    fprintf(stderr,"\r");
    while(1){
      ret = ode->stepODEwithEvent(1);
      dyna->freshX0();
      
      if(ret){
	dyna->tstart = *(dyna->t);
	dyna->state = dyna->event->switchState();
	eventCount++;
      }else{
	dyna->tstart = 0;
	//dyna->x0[2] = dyna->tstart;
	dyna->tfinal += dyna->tend;
	break;
      }
    }
    
    if(mapCount%100 == 0);
    ofs<<dyna->tend<<" "<<dyna->x0[0]<<" "<<dyna->x0[1]<<endl;    
    ofs2<<tmpX<<" "<<dyna->x0[0]<<endl<<dyna->x0[0]<<" "<<dyna->x0[0]<<endl;
    
    switch(eventCount){
    case 0: ofsP1<<tmpX<<" "<<dyna->x0[0]<<endl;
      break;
    case 1: ofsP2<<tmpX<<" "<<dyna->x0[0]<<endl;
      break;
    case 2: ofsP3<<tmpX<<" "<<dyna->x0[0]<<endl;
      break;
    case 3: ofsP4<<tmpX<<" "<<dyna->x0[0]<<endl;
      break;
    case 4: ofsP5<<tmpX<<" "<<dyna->x0[0]<<endl;
      break;
    default:ofsP6<<tmpX<<" "<<dyna->x0[0]<<endl;
      break;
    }
    eventCount = 0;
    tmpX = dyna->x0[0];
  }
  fprintf(stderr,"\n");
  
  if(argc == 4 && argv[3][0] == '1'){
    dyna->printStd(NULL);
  }

  delete dyna;
  delete ode;
  return 0;
}
