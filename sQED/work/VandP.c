#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */
extern int access(const char *pathname, int mode);

int nModelParticles=5;
static ModelPrtclsStr ModelPrtcls_[5]=
{
  {"e-","e+", 11, "Me","Wef",1,1,0}
, {"Ap","Ap", 1000050, "MAp","wAp",2,1,0}
, {"S","S", 1000034, "MS","WS",0,1,0}
, {"~n1","~n1~", 1000035, "Mn1","Wn1",0,1,0}
, {"~n2","~n2~", 1000045, "Mn2","Wn2",0,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=14;
int nModelFunc=0;
static char*varNames_[14]={
 "gD","aEWM1","epsilon","sd","lamPsiS","lamChiPsi","vevH","Me","MAp","MS"
,"Mn1","Mn2","E","Pi"};
char**varNames=varNames_;
static REAL varValues_[14]={
   1.100000E+00,  1.279000E+02,  1.100000E+00,  1.000000E-02,  1.000000E-01,  1.000000E-01,  0.000000E+00,  5.140000E-04,  1.000000E-01,  1.000000E-01
,  1.000000E+00,  1.000000E+00,  2.718282E+00,  3.141593E+00};
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
