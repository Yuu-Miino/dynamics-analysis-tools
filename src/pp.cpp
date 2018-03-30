#include "HybridSystem.hpp"
#include "tools.hpp"

#include <stdio.h>
#include <unistd.h>

using namespace std;

int main(int argc, char **argv){
  // Help message
  if(argc < 2 || !strcmp("--help",argv[1])){
    fprintf(stderr,"usage:\t./pp [-m <counts_of_maps>] filename\n");
    return 0;
  }

  // Variables definition
  FILE *fpPoin, *fpOrbit;
  int countOfMaps = 10;
  State init(STATE_DIM), dst(STATE_DIM);
  Parameter para(PARA_DIM);
  int mode;

  // Options
  int opt; opterr = 0;
  while ((opt = getopt(argc, argv, "m:")) != -1) {
    switch(opt){
    case 'm':
      countOfMaps = atoi(optarg);
      break;
    default:
      printf("warning: undefined option is selected. %d\n",opt);
      break;
    }
  }

  // Initializations
  string infile = argv[optind];
  getFromFile(init,2,para,mode,infile);
  HybridSystem hs(mode);

  fpPoin  = fopen((infile+".pp.poin").c_str(),"w");
  fpOrbit = fopen((infile+".pp.orbit").c_str(),"w");

  // Display informations
  fprintf(stderr,"=============== PP program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"counts of maps\t: %d\n",countOfMaps);
  fprintf(stderr,"initial values\t: ");  init.printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"initial mode \t: %d\n",hs.getMode());
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr);  fprintf(stderr,"\n\n");

  // Print arguments
  printArg(fpOrbit,argc,argv);
  printArg(fpPoin,argc,argv);

  // Main loop
  for(int mapCount = 0; mapCount <= abs(countOfMaps); mapCount++){
    fprintf(stderr,"%07.3lf %% (mapCount: %03d) ",
	    (double)mapCount*100/(double)countOfMaps,mapCount);
    fprintf(stderr,"mode: %d | ",hs.getMode()); 
    init.printT(stderr); init.printX(stderr,PRINT_DIM); fprintf(stderr,"\r");

    hs.map(init,para,2.0*M_PI,dst,fpOrbit);

    init = dst;
    init.printT(fpPoin);
    init.printX(fpPoin,PRINT_DIM);
    fprintf(fpPoin,"%d \n",hs.getMode());
  }// End of main loop

  fclose(fpOrbit);
  fclose(fpPoin);
  fprintf(stderr,"\n=============== PP program ===============\n");
  return 0;
}
