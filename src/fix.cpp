#include "HybridSystem.hpp"
#include "tools.hpp"

#include <fstream>
#include <iostream>

#include <stdio.h>
#include <unistd.h>
#include <Eigen/Dense>

#define FIX_ITER_MAX 10
#define FIX_MAT_DIM   2

using namespace std;
using namespace Eigen;

int main(int argc, char **argv){
  // Help message
  if(argc < 2 || !strcmp("--help",argv[1])){
    fprintf(stderr,"usage:\t./fix [-p <period>] filename\n");
    return 0;
  }

  // Variables definition
  FILE *fpJac, *fpPt;
  ofstream ofsLog;
  int period = 1;
  State init(STATE_DIM), dst(STATE_DIM);
  Parameter para(PARA_DIM);
  int mode;

  MatrixXd jac, jacP, df(FIX_MAT_DIM,FIX_MAT_DIM);
  VectorXd f(FIX_MAT_DIM), h(FIX_MAT_DIM);

  // Options
  int opt; opterr = 0;
  while ((opt = getopt(argc, argv, "p:")) != -1) {
    switch(opt){
    case 'p':
      period = atoi(optarg);
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
  init.setX(JAC_MAT_DIM,1);
  init.setX(JAC_MAT_DIM+4,1);
  init.setX(JAC_MAT_DIM+8,1);

  fpJac = fopen((infile+".fix.jac").c_str(),"w");
  fpPt  = fopen((infile+".fix.pt").c_str(),"w");
  ofsLog.open(infile+".fix.log");

  // Display informations
  fprintf(stderr,"=============== FIX program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"given period\t: %d\n",period);
  fprintf(stderr,"initial values\t: ");  init.printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"initial mode \t: %d\n",hs.getMode());
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr);  fprintf(stderr,"\n\n");

  // Print arguments
  printArg(fpJac,argc,argv);
  printArg(fpPt,argc,argv);

  // Main loop
  for(int fixIter = 1; fixIter <= FIX_ITER_MAX; fixIter++){
    // Print info
    fprintf(stderr,"fixIter: %03d | ",fixIter);
    fprintf(stderr,"mode: %d | ",hs.getMode()); 
    init.printT(stderr); init.printX(stderr,PRINT_DIM);

    // Get Jacobian matrix
    hs.jacobian(jac,jacP,period,init,para,2.0*M_PI,dst);

    // Newton's method
    df = jac.block<2,2>(0,0) - MatrixXd::Identity(2,2);
    f << dst.getX(0) - init.getX(0), dst.getX(1) - init.getX(1);
    h = df.colPivHouseholderQr().solve(f);
    fprintf(stderr," | %e",f.norm());
    fprintf(stderr,"\r");

    // Log
    ofsLog<<"fixIter: "<<fixIter<<"**************************"<<endl;
    ofsLog<<"init:\n\t"<<init.getX(0)<<"  "<<init.getX(1)<<endl;
    ofsLog<<"dst: \n\t"<<dst.getX(0)<<"  "<<dst.getX(1)<<endl;
    ofsLog<<"jac:\n"<<jac<<endl;
    ofsLog<<"jacP:\n"<<jacP<<endl;
    ofsLog<<"df:\n"<<df<<endl;
    ofsLog<<"f:\n"<<f<<endl;
    ofsLog<<"h:\n"<<h<<endl<<endl;

    // If fixed
    if(f.norm() < 1e-6 && fixIter > 1){
      fprintf(stderr,"\n\nFixed in %2d loops (error: %e) | ",fixIter,f.norm());
      fprintf(stderr,"mode: %d | ",hs.getMode()); 
      init.printT(stderr); init.printX(stderr,PRINT_DIM); fprintf(stderr,"\n");
      
      ComplexEigenSolver<MatrixXd> eig(jac.block<2,2>(0,0));
      cerr<<"Eigenvalues:\n"<<eig.eigenvalues()<<endl;
      
      fprintf(fpJac,"%+.12lf %+.12lf\n%+.12lf %+.12lf",
	      jac(0,0),jac(0,1),jac(1,0),jac(1,1));
      para.printValue(fpPt);
      init.printX(fpPt,PRINT_DIM);
      fprintf(fpPt,"%d",hs.getMode());
      break;
    }
    
    // Renew values 
    for(int i = 0; i < 2; i++) init.addX(i, -h(i));

  }// End of main loop

  fclose(fpJac);
  fclose(fpPt);
  fprintf(stderr,"\n=============== FIX program ===============\n");
  return 0;
}
