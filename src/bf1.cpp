#include "HybridSystem.hpp"
#include "tools.hpp"

#include <iostream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

int main(int argc, char **argv){
  // Help message
  if(argc < 2 || !strcmp("--help",argv[1])){
    fprintf(stderr,"usage:\t./bf1 [-m <counts_of_maps>] [-i <index_of_parameter>] [-e <end_value>] [-r <resolution>] filename\n");
    return 0;
  }

  // Variables definition
  FILE *fpOut;
  int countOfMaps = 100;
  State init(STATE_DIM), dst(STATE_DIM);
  Parameter para(PARA_DIM);
  int mode;

  int paraIndex = 0;
  double startVal, endVal; bool endValFlag = false;
  int resolution = 100;

  // Options
  int opt; opterr = 0;
  while ((opt = getopt(argc, argv, "m:i:e:r:")) != -1) {
    switch(opt){
    case 'm':
      countOfMaps = atoi(optarg);
      break;
    case 'i':
      paraIndex = atoi(optarg);
      break;
    case 'e':
      endVal = atof(optarg);
      endValFlag = true;
      break;
    case 'r':
      resolution = atoi(optarg);
      break;
    default:
      printf("warning: undefined option is selected. %d\n",opt);
      break;
    }
  }

  // Initializations
  string infile = argv[optind];

  
  getFromFileWithMode(init,2,para,mode,infile);
  HybridSystem hs;


  startVal = para.getValue(paraIndex);
  fpOut  = fopen((infile+".bf1").c_str(),"w");

  // Display informations
  fprintf(stderr,"=============== BF1 program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"counts of maps\t: %d\n",countOfMaps);
  fprintf(stderr,"initial values\t: ");  init.printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"initial mode \t: %d\n",mode);
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr); fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"parameter index\t: %d\n",paraIndex);
  fprintf(stderr,"resolution\t: %d\n",resolution);
  fprintf(stderr,"\n");

  if(!endValFlag){
    fprintf(stderr,"Set the end value of parameter continuation: "); std::cin>>endVal;
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"end value of parameter continuation\t: %lf\n",endVal);
  fprintf(stderr,"\n");

  // Print arguments
  printArg(fpOut,argc,argv);

  // Main loop
  for(int resCount = 0; resCount <= resolution; resCount++){
    para.setValue(paraIndex,startVal+(endVal-startVal)/(double)resolution*(double)resCount);
    init.setT(0);
    init.setPhaseDiff(0);

    fprintf(stderr,"%07.3lf %% (resCount: %03d) ", (double)resCount*100/(double)resolution,resCount);
    fprintf(stderr,"mode: %d | %+.12lf | ",mode,para.getValue(paraIndex));
    init.printT(stderr); 
    init.printX(stderr,PRINT_DIM); 
    fprintf(stderr,"\r");
    for(int mapCount = 0; mapCount <= abs(countOfMaps); mapCount++){
      hs.map(init,para,mode,2.0*M_PI,dst,NULL);

      init = dst;

      para.printValue(fpOut);
      init.printT(fpOut);
      init.printX(fpOut,STATE_DIM);
      fprintf(fpOut,"%d \n",mode);
    }
    fprintf(fpOut,"\n");
  }// End of main loop

  fclose(fpOut);
  fprintf(stderr,"\n=============== BF1 program ===============\n");
  return 0;
}
