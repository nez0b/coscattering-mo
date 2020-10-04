/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
//#define MASSES_INFO      
  /* Display information about mass spectrum  */
  
//#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */
   
//#define DECAYS

//#define CROSS_SECTIONS 
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"

//double vsig(double temp) {return vSigmaCC(temp,newProcess("~n1,e-->~n2,e-"),0);}
//double vsigPsiAn(double temp) {return vSigmaCC(temp,newProcess("~n2,~n2~->S,S"),0);}

int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

   double Y[50],dY[50];
   for(int i=0; i<50;i++){ Y[i]=1+i; dY[i]=sqrt(i+1)/2;}
   
  ForceUG=0;  /* to Force Unitary Gauge assign 1 */

  VZdecay=1; VWdecay=1;  

  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
  err=readVar(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
  if(CDM1) 
  { 
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     //printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2) 
  { 
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

  double res1, res2, res3;
  res1 = neqF(0.5,0.01,2);
  res2 = neqF(0.01,0.01,2);
  res3 = neqF(0.002,0.01,2);
  printf("neq(0.1):%3e\tneq(T2):%3e\t neq(T3):%3e\n",res1,res2,res3);

  int a[3]={5,7,9};
  int *ptr = a;
  int i=0;
  /*
  for (i=0;i<3;i++){
      printf("%i\n",*(ptr+i));
  }
  */

  double c[5] = { 5.1, 4.2, 3.3, 2.4, 1.5 }; 

  double *b;
  b = c; 
   //double *d[5];
  for (i = 0; i < 5; i++){ 
     // printf("%d\n", *(*b + i)); 
     //b[i] = &c[i];
     //d[i] = &c[i];
      printf("%f\n", (b[i])); 
      printf("%p\n", (&b+i)); 
  }
  
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif






#ifdef OMEGA
{ int fast=0;
  double Beps=0., cut=0.01;
  double Omega1, Omega2, frac1, frac2;  
  int i,err; 
  //printf("\n==== Calculation of relic density =====\n");   
  double Xf1, Xf2;
     //Omega=darkOmega(&Xf,fast,Beps,&err);
  Omega1=darkOmegaCosc1(&Xf1,vsig);
//  Omega2=darkOmegaCosc2(&Xf2,vsig,vsigPsiAn);
     //printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printf("%.3e %.3e\n",Omega1,Xf1);
  //printf("%.3e %.3e\n",vsig(0.1),vsig(0.01));
 // printf("%.3e %.3e\n",Omega2,Xf2);
     //if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  
}

#endif

#ifdef FREEZEIN
{
  double TR=1E10;
  double omegaFi;  
  toFeebleList(CDM1);
  VWdecay=0; VZdecay=0;
  
  omegaFi=darkOmegaFi(TR,&err);
  printf("omega freeze-in=%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);
}
#endif






  


#ifdef DECAYS
{ char*  pname = pdg2name(25);
  txtList L;
  double width; 
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  } 
  
  pname = pdg2name(24);  
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  } 
}            
#endif

#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;
    
  next=nextOdd(1,&nextM); 
  if(next && nextM<1000)  
  { 
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  }  
}
 
#endif 

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.in  smodels.log  smodels.out  summary.*");  
#endif 



  //killPlots();
  return 0;
}
