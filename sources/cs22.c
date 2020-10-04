#include "micromegas.h"
#include "micromegas_aux.h"
//#include "micromegas_f.h"
#include"../CalcHEP_src/c_source/ntools/include/vegas.h"

double (*sqme22)(int nsub, double GG, REAL *pvect, REAL*cb_coeff, int * err_code)=NULL; 
int  nsub22=0;

/*===========================================================*/
static double Q_ren,Q_fact;
static double PcmOut, totcoef;
static REAL pvect[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static  int PC[5];  
static  int chan=0; 


static double eps=0.0001;

double GGscale=91.187;
/*
double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm, pout;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}
*/          

int  kin22(double PcmIn,REAL * pmass)
{  
   double sqrtS;
   int i;
   for(i=0;i<16;i++) pvect[i]=0;
   sqrtS=sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
//printf(" PcmOut =%E (%E %E %E) \n", PcmOut,sqrtS,pmass[2],pmass[3] );   
   if(PcmOut<sqrtS*1.E-4) return 1;
   totcoef = PcmOut /(32.0*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   return 0;
}

double  dSigma_dCos(double  cos_f)
{
   double  r;
   double sin_f=sqrt(fabs((1-cos_f)*(1+cos_f)));
   int err_code=0;
   
   pvect[11]=PcmOut*cos_f;
   pvect[15]=-pvect[11];
   pvect[10]=PcmOut*sin_f;
   pvect[14]=-pvect[10];
      
   r = (*sqme22)(nsub22,sqrt(4*M_PI*alphaQCD(GGscale)),pvect,NULL,&err_code);
    
   err_code=0;
   return r * totcoef;
}

static double kinematic_22(double PcmIn, double cs, REAL*pmass, REAL*pvect)
{  int i;
   for(i=0;i<16;i++) pvect[i]=0;
   double sqrtS=sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   double PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
//printf(" PcmOut =%E (%E %E %E) \n", PcmOut,sqrtS,pmass[2],pmass[3] );   
   totcoef =  PcmOut /(32*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   pvect[11]=PcmOut*cs;
   pvect[15]=-pvect[11];  
   pvect[10]=sqrt((PcmOut+pvect[11])*(PcmOut-pvect[11]));
   pvect[14]=-pvect[10];
   return totcoef*3.8937966E8;                  
}


double cs22(numout * cc, int nsub, double P, double cos1, double cos2 , int * err) 
{
  int i,k;
  REAL pmass[4];

  passParameters(cc);
  
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);  
  *err=0;
  sqme22=cc->interface->sqme;
  nsub22=nsub; 
  if(kin22(P,pmass)) return 0; else 
  { int err;
    double res= 3.8937966E8*simpson(dSigma_dCos,cos1,cos2,0.3*eps,&err);
    if(err)  printf("error in simpson: cs22.c line 108\n");
    return res;
  }
}



