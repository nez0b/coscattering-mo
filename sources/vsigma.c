#include <sys/utsname.h>
#include "micromegas.h"
#include "micromegas_aux.h"
//#include "micromegas_f.h"

//======================  vSigmaCC =======================

static double T_;
static double M1,M2,sqrtSmin,sqrtSmax;
static CalcHEP_interface * CI;
static int i3=2,i4=3,i5=4,i6=5;
static REAL pmass[6];
static double GG=1.23;

static double sing2(char *s, int nout, double m, double w)
{  int i;
   int nin=2;
   
   double sum0=0,sum1=0;
   if(s[0]<=nin) return 0;
   if(strlen(s)<2) return 0;

   for(i=0;s[i];i++) sum1+=pmass[s[i]-1];
   if(m<=sum1) return 0;
   for(i=nin;i<nin+nout;i++) sum0+=pmass[i];      
   if(m>= sqrtSmax -sum0 +sum1) return 0;   
   return 1/(m*w);     
}

static double s_integrandT_(double  sqrtS )
{  double sv_tot,t,bess, x1,x2,y;
   double ms,md,PcmIn;
   double Rm;
   
   ms = M1 + M2; 
   if(ms>=sqrtS)  return 0;
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;

   if(y-x1-x2>50) return 0;   
      
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   kin22(PcmIn,pmass);
   
   sv_tot=simpson(dSigma_dCos,-1.,1.,1E-3,NULL); 
   bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
   Rm=PcmIn*sqrtS/M1/M2;   
   return  bess*Rm*Rm*sv_tot/T_;    
}   

/*
bessK2(x) = exp(-x)*sqrt(M_PI/2/x)*K2pol(1/x)
bessK1(x) = exp(-x)*sqrt(M_PI/2/x)*K1pol(1/x) 
*/

static double u_integrand_( double u)
{  double z,y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;
   z=1-u*u;
   sqrtS=M1+M2-3*T_*log(z);
   if(sqrtS<=M1+M2 || sqrtS<=pmass[2]+pmass[3]) return 0;
   return s_integrandT_(sqrtS )*6*T_*u/z;
}

double  vcs22(numout * cc,int nsub,int * err)
{
   int i;
   double pcm,r;
   REAL pmass[4], pvect[16];
   
//printf("Mb4a=%E\n", findValW("Mb"));    
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
//printf("Mbb=%E Q=%E\n", findValW("Mb"),  findValW("Q")); 
   if( cc->interface->calcFunc()>0 ) {*err=4;  return 0;}
//printf("Mbb_=%E Q=%E\n", findValW("Mb"),  findValW("Q"));
   *(cc->interface->gtwidth)=0;
   *(cc->interface->twidth)=0;
   *(cc->interface->gswidth)=0;
//printf("Mb4c=%E\n", findValW("Mb"));    
   for(i=0;i<4;i++)  cc->interface->pinf(nsub,1+i,pmass+i,NULL);   
   *err=0;
   if(pmass[0]+pmass[1] <= pmass[2]+pmass[3]) return 0;
   for(i=0;i<16;i++) pvect[i]=0;
//printf("Mb4d=%E\n", findValW("Mb")); 
   pcm= decayPcm(pmass[0]+pmass[1],pmass[2],pmass[3]);
   for(i=0;i<2; i++) pvect[4*i]=pmass[i];
   for(i=2;i<4; i++) pvect[4*i]=sqrt(pmass[i]*pmass[i] +pcm*pcm);
   pvect[8+3]=pcm;
   pvect[12+3]=-pcm;
   r=cc->interface->sqme(nsub,GG,pvect,NULL,err);
   return 3.8937966E8*r*pcm/(16*M_PI*pmass[0]*pmass[1]*(pmass[0]+pmass[1]));  
}

double vSigmaCC(double T,numout* cc, int mode)
{
  int i,err,n,n0,m,w;
  char*s, *pname[6];
  int pdg[6];
  double msum;
  double a=0,factor,dMax=0;
  int spin2,cdim,neutral1,neutral2;
  double oldQ;
  
  double bEps=1.E-4;
  double dI;
  int dmOut;
      
  CI=cc->interface; 
  T_=T;
  


  if(passParameters(cc)) return -1;
  if(Qaddress && CI->nout==2) 
  {  oldQ=*Qaddress;
     for(i=0;i<2;i++) pname[i]=CI->pinf(1,i+1,pmass+i,pdg+i);
     *Qaddress=pmass[0]+pmass[1];
     calcMainFunc();
     if(passParameters(cc)) return -1;
  }
  
  for(i=0;i<2+CI->nout;i++) pname[i]=CI->pinf(1,i+1,pmass+i,pdg+i);  

  M1=pmass[0];
  M2=pmass[1];
  

  if(mode) 
  { if(pname[0][0]!='~' || pname[1][0]!='~') return 0;
    if(T==0 && (M1< Mcdm || M2<Mcdm)) return 0;
    dmOut=0; 
    for(i=2;i<2+CI->nout;i++) if(pname[i][0]=='~') dmOut++;
    if(dmOut==2) return 0; 
  }    
 
  for(i=2,msum=0;i<CI->nout;i++) msum+=pmass[i];
  
  sqrtSmin=M1+M2;
  
  if(msum > sqrtSmin)
  { if(T==0) return 0; else sqrtSmin=msum; }
  sqrtSmax=sqrtSmin-T*log(bEps); 

  n0=0; 
  if(CI->nout>2) for(n=1;(s=CI->den_info(1,n,&m,&w,NULL));n++)
  { double d=sing2(s,CI->nout,CI->va[m],CI->va[w]); 
    if(!isfinite(d)) { printf("non-integrable pole\n"); return 0;}
    if(d>dMax){ dMax=d; n0=n;} 
  }

   if(T==0) a=vcs22(cc,1,&err); else
       {  double eps=1.E-3;
          sqme22=CI->sqme;
          nsub22=1;
          a=simpson(u_integrand_,0.,1.,eps,NULL)*3.8937966E8;
       }  


//   WIDTH_FOR_OMEGA=0;
  
   if(mode)
   {  a*= 1-0.5*dmOut;
      char *p0=pname[0],*p1=pname[1],*c0=NULL,*c1=NULL;
      double s=0,k=M1*M1*M2*M2/sqrt(M1*M2)*K2pol(T/M1)*K2pol(T/M2);            
      for(i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~')
      {  double m=pMass(ModelPrtcls[i].name);
         int dim=ModelPrtcls[i].cdim*(ModelPrtcls[i].spin2+1);
         
              if(strcmp(p0,ModelPrtcls[i].name) ==0){ c0=ModelPrtcls[i].aname;k*=dim;} 
         else if(strcmp(p0,ModelPrtcls[i].aname)==0){ c0=ModelPrtcls[i].name; k*=dim;}
              if(strcmp(p1,ModelPrtcls[i].name) ==0){ c1=ModelPrtcls[i].aname;k*=dim;} 
         else if(strcmp(p1,ModelPrtcls[i].aname)==0){ c1=ModelPrtcls[i].name; k*=dim;}
         if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)!=0) dim*=2;
         if(0.5*(M1+M2)-m >30*T){ k=0;s=1;break;}
         s+=dim*m*m/sqrt(m)*K2pol(T/m)*exp(-(m-0.5*(M1+M2))/T); 
      }
      if(k)      
      {  if(strcmp(p0,p1)) k*=2;
         if(! (  (strcmp(p0,c0)==0 && strcmp(p1,c1)==0)
               ||(strcmp(p0,c1)==0 && strcmp(p1,c0)==0)
              ) ) k*=2;   
      } 
      a*=k/s/s;
   }  

   if(Qaddress && CI->nout==2) { *Qaddress=oldQ; calcMainFunc();}
   
   return a;
}

double vsigmacc_(double *T, int *ccf,int*mode)
{
  numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return vSigmaCC(*T,cc,*mode); 
}
