#include "HybridSystem.hpp"
#include "tools.hpp"

#include <fstream>
#include <iostream>

#include <stdio.h>
#include <unistd.h>
#include <Eigen/Dense>

#define BF_ITER_MAX 10
#define BF_MAT_DIM   3

#define DIFF_EPS 1e-5

using namespace std;
using namespace Eigen;

int main(int argc, char **argv){
  // Help message
  if(argc < 2 || !strcmp("--help",argv[1])){
    fprintf(stderr,"usage:\t./bf2 [-p <period>] [-i <paramter_index>] [-G | -I] [-C <parameter_index> [-s <parameter_step>]] filename\n");
    return 0;
  }

  // Variables definition
  FILE *fpPt, *fpCont;
  ofstream ofsLog;
  int period = 1;
  State init[3]= {STATE_DIM,STATE_DIM,STATE_DIM};
  State dst[3] = {STATE_DIM,STATE_DIM,STATE_DIM};
  Parameter para(PARA_DIM);
  int mode;

  MatrixXd jac[3], jacP[3], df(BF_MAT_DIM,BF_MAT_DIM);
  VectorXd f(BF_MAT_DIM), h(BF_MAT_DIM);

  MatrixXd jacDp(2,2), jacInv(2,2);
  double detJac, dDETdX[2], dDETdP;
  int paraIndex = 1;
  double inputMU = -1;

  bool contFlag = false, finFlag = false, stopFlag = false;
  double paraStep   = 1e-3;
  int contParaIndex= 0;

  // Options
  int opt; opterr = 0;
  while ((opt = getopt(argc, argv, "p:i:GIC:s:")) != -1) {
    switch(opt){
    case 'p':
      period = atoi(optarg);
      break;
    case 'i':
      paraIndex = atoi(optarg);
      break;
    case 'C': contParaIndex = atoi(optarg); contFlag = true;  break;
    case 's': if(contFlag) paraStep = atof(optarg); break;      
    case 'G': inputMU = 1.0; break;
    case 'I': inputMU = -1.0; break;
    default:
      printf("warning: undefined option is selected. %d\n",opt);
      break;
    }
  }

  if(contFlag && (contParaIndex == paraIndex)){
    fprintf(stderr,"error: conflict of paraIndex (%d) and contParaIndex (%d)\n", paraIndex, contParaIndex);
    exit(1);
  }

  // Initializations
  string infile = argv[optind];
  getFromFileWithMode(init[0],2,para,mode,infile);
  HybridSystem hs;
  init[0].setX(JAC_MAT_DIM,1);
  init[0].setX(JAC_MAT_DIM+4,1);
  init[0].setX(JAC_MAT_DIM+8,1);

  fpPt   = fopen((infile+".bf2.pt").c_str(),"w");
  if(contFlag){
    char str[20];
    sprintf(str,".%d-%d(%+1.1e)",paraIndex,contParaIndex,paraStep);
    fpCont = fopen((infile+str+".bf2.cont").c_str(),"w");
  }
  ofsLog.open(infile+".bf2.log");

  // Display informations
  fprintf(stderr,"=============== BF program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"given period\t: %d\n",period);
  fprintf(stderr,"initial values\t: ");  init[0].printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"initial mode \t: %d\n",mode);
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr);  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"parameter index\t: %d\n",paraIndex);
  fprintf(stderr,"inputMU \t: %lf\n",inputMU);
  fprintf(stderr,"\n");

  // Print arguments
  printArg(fpPt,argc,argv);
  if(contFlag){
    printArg(fpCont,argc,argv);
    fprintf(fpCont,"# k b0 b x0 y0 mode norm(mu1) norm(mu2) arg(mu1) arg(mu2)\n");
  }
  long int posPt  = ftell(fpPt);

  // Main loop
  while(!stopFlag){
    for(int bfIter = 1; bfIter <= BF_ITER_MAX; bfIter++){
      // Initialize
      init[1] = init[2] = init[0];
      init[1].addX(0,DIFF_EPS);
      init[2].addX(1,DIFF_EPS);

      // Get Jacobian matrix
      for(int i = 0; i < 3; i++) hs.jacobian(jac[i],jacP[i],period,init[i],para,mode,2.0*M_PI,dst[i]);

      // Derivative of determinant
      for(int i = 1; i < 3; i++){
	jac[i]  = (jac[i] - jac[0])/DIFF_EPS;
	jacP[i] = (jacP[i] - jacP[0])/DIFF_EPS;
      }

      jacDp.block(0,0,2,1) = jacP[1].block(0,paraIndex,2,1);
      jacDp.block(0,1,2,1) = jacP[2].block(0,paraIndex,2,1);
        
      jacInv = (jac[0].block(0,0,2,2) - inputMU * MatrixXd::Identity(2,2)).inverse();
      detJac = (jac[0].block(0,0,2,2) - inputMU * MatrixXd::Identity(2,2)).determinant();

      dDETdX[0] = detJac*((jacInv*jac[1].block(0,0,2,2)).trace());
      dDETdX[1] = detJac*((jacInv*jac[2].block(0,0,2,2)).trace());

      dDETdP = detJac*((jacInv*jacDp).trace());    

      // Newton's method
      df.block(0,0,2,2) = jac[0].block(0,0,2,2) - MatrixXd::Identity(2,2);
      df.block(0,2,2,1) = jacP[0].block(0,paraIndex,2,1);
      df.block(2,0,1,3) << dDETdX[0], dDETdX[1], dDETdP;

      f << dst[0].getX(0) - init[0].getX(0), dst[0].getX(1) - init[0].getX(1), detJac;
      h = df.colPivHouseholderQr().solve(f);

      // Log
      ofsLog<<"bfIter: "<<bfIter<<"**************************"<<endl;
      ofsLog<<"init:\n\t"<<init[0].getX(0)<<"  "<<init[0].getX(1)<<endl;
      ofsLog<<"dst: \n\t"<<dst[0].getX(0)<<"  "<<dst[0].getX(1)<<endl;
      ofsLog<<"detJac: \n\t"<<detJac<<endl;
      ofsLog<<"jac:\n"<<jac[0]<<endl;
      ofsLog<<"jacP:\n"<<jacP[0]<<endl;
      ofsLog<<"df:\n"<<df<<endl;
      ofsLog<<"f:\n"<<f<<endl;
      ofsLog<<"h:\n"<<h<<endl<<endl;

      // If finished
      if(f.norm() < 1e-6 && bfIter > 1){
	finFlag = true;
	if(!contFlag) cerr<<"\n\n";
	fprintf(stderr,"Finished in %2d loops (error: %e) | ",bfIter,f.norm());
	para.printValue(stderr); cerr<<" | ";
      
	ComplexEigenSolver<MatrixXd> eig(jac[0].block<2,2>(0,0));
	IOFormat oneline(12,DontAlignCols,  "", " , ",  "", "",   "", "");
	cerr<<"MU: "<<eig.eigenvalues().format(oneline);
	
	// Print results
	fseek(fpPt,posPt,SEEK_SET);
	para.printValue(fpPt);
	init[0].printX(fpPt,PRINT_DIM);
	fprintf(fpPt,"%d ",mode);
	
	if(contFlag){
	  fprintf(stderr,"\r");
	  
	  para.printValue(fpCont);
	  init[0].printX(fpCont,PRINT_DIM);
	  fprintf(fpCont,"%d ",mode);
	  fprintf(fpCont,"%+.12lf %+.12lf %+.12lf %+.12lf",
		  norm(eig.eigenvalues()(0)),norm(eig.eigenvalues()(1)),arg(eig.eigenvalues()(0)),arg(eig.eigenvalues()(1)));
	  fprintf(fpCont,"\n");
	}
	break;
      }
    
      // Renew values 
      for(int i = 0; i < 2; i++) init[0].addX(i, -h(i));
      para.addValue(paraIndex,-h(2));
    }// End of main loop
    if(!finFlag || !contFlag){
      stopFlag = true;
    }else{
      finFlag = false;
    }
    
    para.addValue(contParaIndex,paraStep);
  }// End of continuation loop

  fclose(fpPt);
  if(contFlag) fclose(fpCont);
  fprintf(stderr,"\n=============== BF program ===============\n");
  return 0;
}

