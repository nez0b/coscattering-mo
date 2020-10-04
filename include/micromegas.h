#ifndef  __MICROMEGAS__
#define  __MICROMEGAS__

#ifdef __cplusplus
extern "C" {
#endif 

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>

#include"../CalcHEP_src/include/num_out.h"
#include"../CalcHEP_src/c_source/dynamicME/include/dynamic_cs.h"
#include"../CalcHEP_src/c_source/plot/include/plot.h"
#include"../CalcHEP_src/c_source/ntools/include/n_proc.h"
extern char * CDM1, *CDM2, *aCDM1, *aCDM2;

/*
typedef struct numout
{
  void * handle;
  double ** link;
  double *Q,*SC;
  int init;
  CalcHEP_interface * interface; 
} numout;
*/
extern int VVdecay;

extern int sortOddParticles(char * name);

typedef  struct { double par[41]; }  MOcommonSTR;
extern   MOcommonSTR  mocommon_;

#define Mcdm        mocommon_.par[1]
#define ScalarFFPd  mocommon_.par[2]
#define ScalarFFPu  mocommon_.par[3]
#define ScalarFFPs  mocommon_.par[4]
#define pVectorFFPd mocommon_.par[5]
#define pVectorFFPu mocommon_.par[6]
#define pVectorFFPs mocommon_.par[7]
#define SigmaFFPd   mocommon_.par[8]
#define SigmaFFPu   mocommon_.par[9]
#define SigmaFFPs   mocommon_.par[10]

#define ScalarFFNd  mocommon_.par[11]
#define ScalarFFNu  mocommon_.par[12]
#define ScalarFFNs  mocommon_.par[13]
#define pVectorFFNd mocommon_.par[14]
#define pVectorFFNu mocommon_.par[15]
#define pVectorFFNs mocommon_.par[16]
#define SigmaFFNd   mocommon_.par[17]
#define SigmaFFNu   mocommon_.par[18]
#define SigmaFFNs   mocommon_.par[19]

#define Fermi_a     mocommon_.par[20]
#define Fermi_b     mocommon_.par[21]
#define Fermi_c     mocommon_.par[22]

#define Rsun        mocommon_.par[23]
#define rhoDM       mocommon_.par[24]
#define Vearth      mocommon_.par[25]

#define K_dif       mocommon_.par[26]      
#define L_dif       mocommon_.par[27]
#define Delta_dif   mocommon_.par[28]
#define Tau_dif     mocommon_.par[29]
#define Vc_dif      mocommon_.par[30]
#define Rdisk       mocommon_.par[31]

#define deltaY      mocommon_.par[32]
#define dmAsymm     mocommon_.par[33]
#define Vesc        mocommon_.par[34]
#define Vrot        mocommon_.par[35]
#define fracCDM2    mocommon_.par[36]
#define Mcdm1       mocommon_.par[37]
#define Mcdm2       mocommon_.par[38]
#define Tstart      mocommon_.par[39]
#define Tend        mocommon_.par[40]

typedef  struct { int flag[4];}     MOflagsSTR;
extern   MOflagsSTR moflags_;

#define WIMPSIM   moflags_.flag[0] 
#define forRocks  moflags_.flag[1]

typedef  struct { char txt1[20]; char txt2[20];} MoCommonCH; 

extern  MoCommonCH mocommonch_;
         
#define CDM1_ mocommonch_.txt1 
#define CDM2_ mocommonch_.txt2 


extern double sWidth;

extern int ForceUG;
/*============================
     Particles
==============================*/

extern int  pNum(char * name);  /* returns PDG code */
extern double pMass(char * name); /* returns particle mass */
extern char *   pdg2name(int pdg); 
/*======= Subprocesses ===========*/
/*  typedef struct txtListStr
  {  struct txtListStr * next;
     char  *txt;
  } txtListStr;

  typedef txtListStr * txtList;
*/  
extern txtList  makeDecayList(char * pname, int nx);
extern void massFilter(double M, txtList * List);
extern void gammaGluFilter(txtList* List);
extern int process2Lib(char * process,char * lib);
extern void cleanTxtList(txtList L);
extern double findBr(txtList L, char * pattern);
extern void printTxtList(txtList L, FILE *f);
extern void printPartialWidth(double width,txtList l,FILE*f);

/*=============================
     (1)2->2 matrix elements
=============================*/
extern double decay2Info(char * pname, FILE* f);
extern numout *  newProcess(char*Process);
extern double GGscale;
extern double cs22(numout * cc, int nsub,double P,  double cos1, double cos2 , int * err);
extern int  procInfo1(numout*cc, int *nsub, int * nin, int *nout);
extern int procInfo2(numout*cc,int nsub,char**name,REAL*mass);
extern REAL Helicity[2];

extern void PDTList(void);
extern int  setPDT(char*name);
extern char pdfName[50];
void LHAPDFList(void);
extern void  setLHAPDF(int nset, char *name);
extern int restorePDF(char*oldPDF);
extern double hCollider(double Pcm, int pp, int nf, double Qren,double Qfact, char * name1,char *name2,double pTmin,int wrt);
extern double monoJet(void);
double pWidth(char *name, txtList *L);


extern double convStrFun3(double x, double q, int pc1, int pc2, int pp);

/*===================
      Variables 
=====================*/

extern REAL*varAddress(char * name); 
extern int    assignVal(const char * name, double val);
extern int    findVal(char * name,double * val);
extern int    assignValW(const char * name, double val);
extern double findValW(char * name);
extern int    readVar(char *fname);
  
/*===========================================
   Checking of parameters 
  ===========================================*/ 
extern void printVar(FILE *f);
extern void  printMasses(FILE * f, int sort);
extern void  printHiggs(FILE * f);
extern char * nextOdd(int N,double * Mass);

/*=====   Constraints ==== */
extern int Zinvisible(void);
extern int LspNlsp_LEP(double *cs_out);
extern int Zprimelimits(void);

/*=====================================
    Relic density evaluation 
  =====================================*/
  
typedef struct {int err; double weight; char*prtcl[10];} aChannel;
extern aChannel*omegaCh;  
extern aChannel* vSigmaTCh;
extern aChannel*omegaFiCh;

extern double vSigmaCC(double T,numout* cc,int mode);
   
extern int loadHeffGeff(char*fname);
extern double  hEff(double T);
extern double  gEff(double T);
extern double T_s3(double s3);
extern double s3_T(double T);

extern double  h1eff(double x, int eta);
extern double  g1eff(double x, int eta);
extern double Hubble(double T);

extern int toFeebleList(char*pname);
extern int isFeeble(char*name);

extern double hEffLnDiff(double T);
extern double vSigma(double T,double Beps ,int Fast);
extern double darkOmega(double *Xf,int Fast, double Beps,int *err);
extern double darkOmega2(double fast, double Beps0);
extern double darkOmegaExt(double *Xf, double (*f0)(double), double (*f1)(double));
//experimental function
extern double darkOmegaExt1(double *Xf, double (*f0)(double), double (*f1)(double));
extern double darkOmegaCosc1(double *Xf, double (*f0)(double));
extern double darkOmegaCosc2(double * Xf,  double (*f0)(double),double (*f1)(double));
extern double darkOmegaCosc2d(double * Xf,  double (*f0)(double),double (*f1)(double));
extern double neqF( double T, double mf, double g);

extern double vs1120F(double T);       
extern double vs2200F(double T);          
extern double vs1100F(double T);      
extern double vs1210F(double T);       
extern double vs1122F(double T);       
extern double vs2211F(double T);
extern double vs1110F(double T);
extern double vs2220F(double T);
extern double vs1112F(double T);
extern double vs1222F(double T);
extern double vs1220F(double T);
extern double vs2210F(double T);
extern double vs2221F(double T);
extern double vs1211F(double T);

extern double TCoeffF(double T);
extern double dY1F(double T);
extern double dY2F(double T);
extern double Y1F(double T);
extern double Y2F(double T);
extern double YF(double T);

extern double darkOmegaFO(double *Xf,int fast,double Beps);
extern double printChannels(double Xf,double cut,double Beps,int prcnt,FILE *f );   
extern double oneChannel(double Xf,double Beps,char*n1,char*n2,char*n3,char*n4);
extern void improveCrossSection(long n1,long n2,long n3,long n4,double Pcm, 
                                                            double * addr);

extern double Yeq(double T);
extern double Yeq1(double T);
extern double Yeq2(double T);

extern double Beps;
extern int Fast_;
 


#ifdef __cplusplus
}
#endif 

#include"../CalcHEP_src/c_source/SLHAplus/include/SLHAplus.h"

#endif
