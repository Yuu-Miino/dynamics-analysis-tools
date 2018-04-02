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
    fprintf(stderr,"usage:\t./fix [-p <period>] [-C <parameter_index>] filename\n");
    return 0;
  }

  // Variables definition
  FILE *fpJac, *fpPt, *fpCont;
  ofstream ofsLog;
  int period = 1;
  State init(STATE_DIM), dst(STATE_DIM);
  Parameter para(PARA_DIM);
  int mode;

  MatrixXd jac, jacP, df(FIX_MAT_DIM,FIX_MAT_DIM);
  VectorXd f(FIX_MAT_DIM), h(FIX_MAT_DIM);

  bool contFlag = false;
  double endVal = 0, paraStep   = 0;
  int paraIndex = 0, resolution = 0;

  // Options
  int opt; opterr = 0;
  while ((opt = getopt(argc, argv, "p:C:")) != -1) {
    switch(opt){
    case 'p':
      period = atoi(optarg);
      break;
    case 'C': paraIndex = atoi(optarg); contFlag = true; break;
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

  fpJac  = fopen((infile+".fix.jac").c_str(),"w");
  fpPt   = fopen((infile+".fix.pt").c_str(),"w");
  fpCont = fopen((infile+".fix.cont").c_str(),"w");
  ofsLog.open(infile+".fix.log");

  // Display informations
  fprintf(stderr,"=============== FIX program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"given period\t: %d\n",period);
  fprintf(stderr,"initial values\t: ");  init.printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"initial mode \t: %d\n",hs.getMode());
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr);  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  
  if(contFlag){
    fprintf(stderr,"Set the end value  of parameter continuation: "); std::cin>>endVal;
    fprintf(stderr,"Set the resolution of parameter continuation: "); std::cin>>resolution;
    fprintf(stderr,"\n");
    paraStep = (endVal - para.getValue(paraIndex))/(double)resolution;    
  }

  // Print arguments
  printArg(fpJac,argc,argv);
  printArg(fpPt,argc,argv);
  if(contFlag){
    printArg(fpCont,argc,argv);
    fprintf(fpCont,"# k b0 b x0 y0 mode norm(mu1) norm(mu2) arg(mu1) arg(mu2)\n");
  }
  long int posJac = ftell(fpJac);
  long int posPt  = ftell(fpPt);

  // Main loop
  for(int resCount = 0; resCount <= resolution; resCount++){
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

      // If finished
      if(f.norm() < 1e-6 && fixIter > 1){
	ComplexEigenSolver<MatrixXd> eig(jac.block<2,2>(0,0));
	if(norm(eig.eigenvalues()(0)) - 1.0 > ZERO || norm(eig.eigenvalues()(1)) - 1.0 > ZERO) 
	  fprintf(stderr,"\e[41m\e[97m");

	if(!contFlag) cerr<<"\n\n";
	fprintf(stderr,"Finished in %2d loops (error: %e) | ",fixIter,f.norm());
	para.printValue(stderr); cerr<<" | ";
	//fprintf(stderr,"mode: %d | ",hs.getMode());
	//init.printT(stderr); init.printX(stderr,PRINT_DIM); cerr<<" | ";      

	IOFormat oneline(12,DontAlignCols,  "", " , ",  "", "",   "", "");
	cerr<<"MU: "<<eig.eigenvalues().format(oneline);
	
	// Print results
	fseek(fpJac,posJac,SEEK_SET);
	fseek(fpPt,posPt,SEEK_SET);
	fprintf(fpJac,"%+.12lf %+.12lf\n%+.12lf %+.12lf",
		jac(0,0),jac(0,1),jac(1,0),jac(1,1));
	para.printValue(fpPt);
	init.printX(fpPt,PRINT_DIM);
	fprintf(fpPt,"%d",hs.getMode());

	if(contFlag){
	  fprintf(stderr,"\e[m\n");

	  para.printValue(fpCont);
	  init.printX(fpCont,PRINT_DIM);
	  fprintf(fpCont,"%d ",hs.getMode());
	  fprintf(fpCont,"%+.12lf %+.12lf %+.12lf %+.12lf",
		  norm(eig.eigenvalues()(0)),norm(eig.eigenvalues()(1)),arg(eig.eigenvalues()(0)),arg(eig.eigenvalues()(1)));
	  fprintf(fpCont,"\n");
	}
	break;
      }
    
      // Renew values 
      for(int i = 0; i < 2; i++) init.addX(i, -h(i));
      
    }// End of main loop
    para.addValue(paraIndex,paraStep);
  }// End of continuation loop

  fclose(fpJac);
  fclose(fpPt);
  fclose(fpCont);
  fprintf(stderr,"\n=============== FIX program ===============\n");
  return 0;
}
