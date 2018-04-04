#ifndef _TOOLS_
#define _TOOLS_

#include <stdio.h>
#include <fstream>

#include "Dynamics.hpp"

/* 
 * @brief Print the arguments at fp.
 *
 * @param FILE|fp   Out destination.
 * @param  int|argc Number of arguments. 
 * @param  int|argc Number of arguments. 
 * @param char|argv Entity of arguments.
*/

void printArg(FILE *fp, int argc, char** argv){
  fprintf(fp,"# "); 
  for(int i = 0; i < argc; i++) fprintf(fp,"%s ",argv[i]);
  fprintf(fp,"\n");
}

/*
 * @brief Set initial values, parameters and a mode from a file.
 *
 * @param State    |init   Initial value of ODE.
 * @param Parameter|para   Parameters of System equations.
 * @param int      |mode   Current mode of the hybrid system.
 * @param string   |infile The name of data file.
*/
void getFromFileWithMode(State& init, int dim, Parameter& para, int& mode, std::string infile){
  std::ifstream ifs;
  int pos;
  std::string line;

  ifs.open(infile);
  if(!ifs.is_open()){
    fprintf(stderr,"error: cannot open file. [filename = %s]\n",infile.c_str());
    exit(1);
  }
  
  pos = ifs.tellg();
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0] != '#') break;
    pos = ifs.tellg();
  }
  ifs.seekg(pos);

  init.setT(0);
  para.setValueFromFile(ifs);
  init.setXfromFile(ifs,dim);
  ifs>>mode;
}
void getFromFile(State& init, int dim, Parameter& para, std::string infile){
  std::ifstream ifs;
  int pos;
  std::string line;

  ifs.open(infile);
  if(!ifs.is_open()){
    fprintf(stderr,"error: cannot open file. [filename = %s]\n",infile.c_str());
    exit(1);
  }
  
  pos = ifs.tellg();
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0] != '#') break;
    pos = ifs.tellg();
  }
  ifs.seekg(pos);

  init.setT(0);
  para.setValueFromFile(ifs);
  init.setXfromFile(ifs,dim);
}

#endif
