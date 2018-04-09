#include "HybridSystem.hpp"
#include "tools.hpp"

#include <fstream>
#include <iostream>

#include <stdio.h>
#include <unistd.h>
#include <Eigen/Dense>

#define GZ_ITER_MAX 10
#define GZ_MAT_DIM   4

#define DIFF_EPS 1e-5

using namespace std;
using namespace Eigen;

int main(int argc, char **argv){
  // Help message
  if(argc < 2 || !strcmp("--help",argv[1])){
    fprintf(stderr,"usage:\t./gz2 [-p <period>] [-i <paramter_index>] [-e <event_index>] [-C <parameter_index> [-s <parameter_step>]] filename TAU\n");
    fprintf(stderr,"\t(*TAU is the estimated time when the orbit grazes the boundary.)\n");
    return 0;
  }

  // Variables definition
  FILE *fpPt, *fpCont;
  ofstream ofsLog;
  int period = 1;
  State init(STATE_DIM);
  State dst[2] = {STATE_DIM,STATE_DIM};
  Parameter para(PARA_DIM);
  int initMode, mode[2];

  MatrixXd jac[2], jacP[2], df(GZ_MAT_DIM,GZ_MAT_DIM);
  VectorXd f(GZ_MAT_DIM), h(GZ_MAT_DIM);

  int paraIndex = 1;

  bool isContinuation = false, isFinish = false, isStopped = false;
  double paraStep   = 1e-3;
  int contParaIndex= 0;

  double dxdt[STATE_DIM];
  double TAU;
  double ef[2];
  int eventIndex=0;

  // Options
  int opt; opterr = 0;
  while ((opt = getopt(argc, argv, "p:i:C:s:e:")) != -1) {
    switch(opt){
    case 'p':
      period = atoi(optarg);
      break;
    case 'i':
      paraIndex = atoi(optarg);
      break;
    case 'e':
      eventIndex = atoi(optarg);
      break;
    case 'C': contParaIndex = atoi(optarg); isContinuation = true;  break;
    case 's': if(isContinuation) paraStep = atof(optarg); break;      
    default:
      printf("warning: undefined option is selected. %d\n",opt);
      break;
    }
  }

  if(isContinuation && (contParaIndex == paraIndex)){
    fprintf(stderr,"error: conflict of paraIndex (%d) and contParaIndex (%d)\n", paraIndex, contParaIndex);
    exit(1);
  }

  // Initializations
  string infile = argv[optind];
  TAU = atof(argv[optind+1]);
  getFromFileWithMode(init,2,para,initMode,infile);
  mode[0] = mode[1] = initMode;
  HybridSystem hs;
  init.setX(JAC_MAT_DIM,1);
  init.setX(JAC_MAT_DIM+4,1);
  init.setX(JAC_MAT_DIM+8,1);

  fpPt   = fopen((infile+".gz2.pt").c_str(),"w");
  if(isContinuation){
    char str[20];
    sprintf(str,".%d-%d(%+1.1e)",paraIndex,contParaIndex,paraStep);
    fpCont = fopen((infile+str+".gz2.cont").c_str(),"w");
  }
  ofsLog.open(infile+".gz2.log");

  // Display informations
  fprintf(stderr,"=============== GZ program ===============\n");
  fprintf(stderr,"filename\t: %s\n",argv[optind]);
  fprintf(stderr,"given period\t: %d\n",period);
  fprintf(stderr,"initial values\t: ");  init.printX(stderr,2); fprintf(stderr,"\n");
  fprintf(stderr,"initial mode \t: %d\n",mode[0]);
  fprintf(stderr,"parameters\t: ");  para.printValue(stderr);  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"parameter index\t: %d\n",paraIndex);
  fprintf(stderr,"TAU\t\t: %lf\n",TAU);
  fprintf(stderr,"event index \t: %d\n",eventIndex);
  fprintf(stderr,"\n");

  // Print arguments
  printArg(fpPt,argc,argv);
  if(isContinuation){
    printArg(fpCont,argc,argv);
    fprintf(fpCont,"# k b0 b x0 y0 mode norm(mu1) norm(mu2) arg(mu1) arg(mu2) TAU\n");
  }
  long int posPt  = ftell(fpPt);

  // Main loop
  while(!isStopped){
    for(int gzIter = 1; gzIter <= GZ_ITER_MAX; gzIter++){
      ofsLog<<"gzIter: "<<gzIter<<"**************************"<<endl;
      ofsLog<<"init:\n\t"<<init.getX(0)<<"  "<<init.getX(1)<<" "<<mode[0]<<endl;
      // Get Jacobian matrix
      for(int i = 0; i < 2; i++)
	hs.jacobian(jac[i],jacP[i],period,init,para,mode[i],(i==0)?2.0*M_PI:TAU,dst[i]);
      hs.dxdt(dxdt,dst[1],para,mode[1]);
      hs.eventFunction(ef,dst[1],para,mode[1]);

      // Newton's method
      df.block(0,0,2,2) = jac[0].block(0,0,2,2) - MatrixXd::Identity(2,2);
      df.block(0,2,2,1) = MatrixXd::Zero(2,1);
      df.block(0,3,2,1) = jacP[0].block(0,paraIndex,2,1);
      df.block(2,0,2,2) = jac[1].block(0,0,2,2);
      df.block(2,2,2,1) << dxdt[0], dxdt[1];
      df.block(2,3,2,1) = jacP[1].block(0,paraIndex,2,1);

      f << 
	dst[0].getX(0) - init.getX(0), 
	dst[0].getX(1) - init.getX(1), 
	ef[eventIndex],
	dxdt[0];
      h = df.colPivHouseholderQr().solve(f);

      // Log
      ofsLog<<"dst[0]: \n\t"<<dst[0].getX(0)<<"  "<<dst[0].getX(1)<<" "<<mode[0]<<endl;
      ofsLog<<"dst[1]: \n\t"<<dst[1].getX(0)<<"  "<<dst[1].getX(1)<<" "<<mode[1]<<endl;
      ofsLog<<"jac[0]:\n"<<jac[0]<<endl;
      ofsLog<<"jac[1]:\n"<<jac[1]<<endl;
      ofsLog<<"jacP[0]:\n"<<jacP[0]<<endl;
      ofsLog<<"jacP[1]:\n"<<jacP[1]<<endl;
      ofsLog<<"df:\n"<<df<<endl;
      ofsLog<<"f:\n"<<f<<endl;
      ofsLog<<"h:\n"<<h<<endl<<endl;

      // If finished
      if(f.norm() < 1e-6 && gzIter > 1){
	isFinish = true;
	if(!isContinuation) cerr<<"\n\n";
	fprintf(stderr,"Finished in %2d loops (error: %e) | ",gzIter,f.norm());
	para.printValue(stderr); cerr<<" | ";
	fprintf(stderr,"TAU = %lf | ",TAU);
      
	ComplexEigenSolver<MatrixXd> eig(jac[0].block<2,2>(0,0));
	IOFormat oneline(12,DontAlignCols,  "", " , ",  "", "",   "", "");
	cerr<<"MU: "<<eig.eigenvalues().format(oneline);
	
	// Print results
	fseek(fpPt,posPt,SEEK_SET);
	para.printValue(fpPt);
	init.printX(fpPt,PRINT_DIM);
	fprintf(fpPt,"%d",mode[0]);
	
	if(isContinuation){
	  fprintf(stderr,"\r");
	  
	  para.printValue(fpCont);
	  init.printX(fpCont,PRINT_DIM);
	  fprintf(fpCont,"%d ",mode[0]);
	  fprintf(fpCont,"%+.12lf %+.12lf %+.12lf %+.12lf ",
		  norm(eig.eigenvalues()(0)),norm(eig.eigenvalues()(1)),arg(eig.eigenvalues()(0)),arg(eig.eigenvalues()(1)));
	  fprintf(fpCont,"%+.12lf ",TAU);
	  fprintf(fpCont,"\n");
	}
	break;
      }
    
      // Renew values 
      for(int i = 0; i < 2; i++) init.addX(i, -h(i));
      TAU -= h(2);
      para.addValue(paraIndex,-h(3));
      mode[0] = mode[1] = initMode;

    }// End of main loop
    if(!isFinish || !isContinuation){
      isStopped = true;
    }else{
      isFinish = false;
    }
    
    para.addValue(contParaIndex,paraStep);
  }// End of continuation loop

  fclose(fpPt);
  if(isContinuation) fclose(fpCont);
  fprintf(stderr,"\n=============== GZ program ===============\n");
  return 0;
}

