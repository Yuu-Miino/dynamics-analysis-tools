#include "pendulum.hpp"
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
  int countOfMaps = 100;
  State init(STATE_DIM), dst(STATE_DIM);
  Parameter para(PARA_DIM);

  pendulum  dyna;
  ODEsolver ode("RK4",1e-2);
  double tfinal = 2.0 * M_PI;
  Domain domain(2);
  domain.setInterval(0,-5,5);
  domain.setInterval(1,-5,5);
  bool divFlag = false;

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
  getFromFile(init,2,para,infile);

  fpPoin  = fopen((infile+".pp.poin").c_str(),"w");
  fpOrbit = fopen((infile+".pp.orbit").c_str(),"w");

  // Display informations
  fprintf(stderr,"=============== PP program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"counts of maps\t: %d\n",countOfMaps);
  fprintf(stderr,"initial values\t: ");  init.printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr);  fprintf(stderr,"\n");
  fprintf(stderr,"\n");

  // Print arguments
  printArg(fpOrbit,argc,argv);
  printArg(fpPoin,argc,argv);

  // Main loop
  for(int mapCount = 0; mapCount <= abs(countOfMaps); mapCount++){
    fprintf(stderr,"%07.3lf %% (mapCount: %03d) ",
	    (double)mapCount*100/(double)countOfMaps,mapCount);
    init.printT(stderr); init.printX(stderr,PRINT_DIM); fprintf(stderr,"\r");

    if(divFlag = ode.runODEsolver(dyna,domain,init,para,tfinal,dst,fpOrbit)){
      fprintf(stderr,"\n\e[41m\e[97mDIVERGENT\e[m\n");
      break;
    }

    init = dst;
    init.setT(0);
    init.addPhaseDiff(tfinal);

    // print
    init.printT(fpPoin);
    init.printX(fpPoin,PRINT_DIM);
    fprintf(fpPoin,"\n");
  }// End of main loop

  fclose(fpOrbit);
  fclose(fpPoin);
  fprintf(stderr,"\n=============== PP program ===============\n");
  return 0;
}
