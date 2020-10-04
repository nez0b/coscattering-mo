#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */
extern int access(const char *pathname, int mode);

int nModelParticles=4;
static ModelPrtclsStr ModelPrtcls_[4]=
{
  {"mu-","mu+", 13, "Mmu","Wmu",1,1,0}
, {"tau-","tau+", 15, "Mtau","Wtau",1,1,0}
, {"~chi","~chi~", 1000050, "Mchi","0",1,1,0}
, {"phi","phi", 1000034, "Mphi","Wphi",0,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=13;
int nModelFunc=0;
static char*varNames_[13]={
 "ychis","ychip","yds","ydp","yods","yodp","aEWM1","Mmu","Mtau","Mchi"
,"Mphi","E","Pi"};
char**varNames=varNames_;
static REAL varValues_[13]={
   1.000000E-01,  1.000000E-01,  1.000000E-01,  1.000000E-01,  1.000000E-01,  1.000000E-01,  1.279000E+02,  1.056600E-01,  1.777000E+00,  1.000000E-01
,  1.000000E-01,  2.718282E+00,  3.141593E+00};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)     return 0;
   }
  cErr=1;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   return 0;
}
