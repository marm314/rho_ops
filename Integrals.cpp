#include"Integrals.h"
//Global variables readable everywhere and defining the interval of integration
double theta_inf,theta_sup,phi_inf,phi_sup,r_inf,r_sup,Point_Vr[3],Rotate_grid[3][3];
bool shanr,fishr,inertiar,shanp,fishp,inertiap,nos,R1,R2,RM1,P1,P2,DIPOLE,RHO;
char DIR_POL_HYPER;
MO mos_array[2];
NO nos_array[2];
////Defines for CUBA
void define()
{
#define LAST 4
#define SEED 0
#define NVEC1 1
#define PCORE 10000
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define NNEW 1000
#define FLATNESS 25.
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0
#define KEY 0
#define SPIN NULL
}
////Interval for integration
void define_interval(double Interval[6])
{
 r_inf=Interval[0];
 r_sup=Interval[1];
 theta_inf=Interval[2];
 theta_sup=Interval[3];
 phi_inf=Interval[4];
 phi_sup=Interval[5];
}
//////////////////////////
//Functions description //
//////////////////////////
void integrate_cuba(READ_FCHK_WFN &read_fchk_wfn,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[20], int &fail,double Integrals_interval[6],int &ncores,bool shan,bool fish,bool inertias,
bool r1,bool r2,bool rm1,bool dipole,bool Rho)
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  bool interval=false;
  shanr=shan;
  fishr=fish;
  inertiar=inertias;
  R1=r1;
  R2=r2;
  RM1=rm1;
  DIPOLE=dipole;
  RHO=Rho;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval=false;}
    else{interval=true;}
   }
   else{interval=true;}
  }
  else{interval=true;}
  if(interval)
  {
   define_interval(Integrals_interval);
  }
  else
  {
   theta_inf=ZERO;theta_sup=PI;phi_inf=ZERO;phi_sup=TWO*PI;r_inf=ZERO;
   r_sup=RSUP;
  }
  void * USERDATA=NULL;
  USERDATA=&read_fchk_wfn;
  define();
  cubacores(ncores,PCORE);
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<20;i++)
 {result_integration[i]=integral[i];}
}
void integrate_cubap(READ_FCHK_WFN &read_fchk_wfn,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[20], int &fail,double Integrals_interval[6],int &ncores,bool shan,bool fish,bool inertias,
bool p1,bool p2)
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  bool interval=false;
  shanp=shan;
  fishp=fish;
  inertiap=inertias;
  P1=p1;
  P2=p2;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval=false;}
    else{interval=true;}
   }
   else{interval=true;}
  }
  else{interval=true;}
  if(interval)
  {
   define_interval(Integrals_interval);
  }
  else
  {
   theta_inf=ZERO;theta_sup=PI;phi_inf=ZERO;phi_sup=TWO*PI;r_inf=ZERO;
   r_sup=RSUP;
  }
  void * USERDATA=NULL;
  USERDATA=&read_fchk_wfn;
  define();
  cubacores(ncores,PCORE);
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrandp, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrandp, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrandp, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrandp, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<19;i++)
 {result_integration[i]=integral[i];}
}
void integrate_cuba_sij(READ_FCHK_WFN &read_fchk_wfn,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double **SIJ, int &fail,double Integrals_interval[6],int &ncores,string MOorNO)
{
 int i,j,k,fail_tot=0,verbose=0,nregions,neval;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 bool interval=false;
 nos=false;
 if(MOorNO=="no")
 {nos=true;} //Affects only the fchk files because for wfn only NOs are available
 if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP))
 {
  if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
  {
   if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
   {interval=false;}
   else{interval=true;}
  }
  else{interval=true;}
 }
 else{interval=true;}
 if(interval)
 {
  define_interval(Integrals_interval);
 }
 else
 {
  theta_inf=ZERO;theta_sup=PI;phi_inf=ZERO;phi_sup=TWO*PI;r_inf=ZERO;
  r_sup=RSUP;
 }
 void * USERDATA=NULL;
 USERDATA=&read_fchk_wfn;
 define();
 cubacores(ncores,PCORE);
 if(NCOMP!=read_fchk_wfn.nbasis())
 {
  if(method=="Cuhre")
  {
   Cuhre(NDIM, NCOMP,SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
   MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
  }
  else if(method=="Divone")
  {
   Divonne(NDIM, NCOMP, SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
   MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
   NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
   integral, error, prob);
  }
  else if(method=="Suave")
  {
   Suave(NDIM, NCOMP, SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
   MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
   error, prob);
  }
  else
  {
   Vegas(NDIM, NCOMP, SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
   MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
   error, prob);
  }
  k=0;
  if(read_fchk_wfn.uhf)
  {
   //If we only have alpha NOs in the wfn file(2e- triplet for instance), no_beta_wfn=true. We correct the spin component
   //and do not asume that the first orbital and second orbitals have diff spins (No Sij=zero if i parity is diff  from j).
   if(!read_fchk_wfn.no_beta_wfn)
   {
    for(i=0;i<read_fchk_wfn.nbasis();i++)
    {
     for(j=0;j<=i;j++)
     {
      SIJ[i][j]=integral[k];
      if(j%2!=i%2)
      {SIJ[i][j]=ZERO;}
      k++;
     }
    }
   }
  }
  else
  {
   for(i=0;i<read_fchk_wfn.nbasis();i++)
   {
    for(j=0;j<=i;j++)
    {
     SIJ[i][j]=integral[k];
     k++;
    }
   }
  }
 }
 else
 {
  for(i=0;i<read_fchk_wfn.nbasis();i++)
  {
   read_fchk_wfn.Pair[0]=i;
   if(method=="Cuhre")
   {
    Cuhre(NDIM, NCOMP,SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
    MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
   }
   else if(method=="Divone")
   {
    Divonne(NDIM, NCOMP, SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
    integral, error, prob);
   }
   else if(method=="Suave")
   {
    Suave(NDIM, NCOMP, SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
    error, prob);
   }
   else
   {
    Vegas(NDIM, NCOMP, SIJ_Integrand, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
    MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
    error, prob);
   }
   k=0;
   if(read_fchk_wfn.uhf)
   {
    for(j=i;j<NCOMP;j++)
    {
     SIJ[j][i]=integral[k];
     if(j%2!=i%2)
     {SIJ[j][i]=ZERO;}
     k++;
    }
   }
   else
   {
    for(j=i;j<NCOMP;j++)
    {
     SIJ[j][i]=integral[k];
     k++;
    }
   }
   fail_tot=fail+fail_tot;
  }
  fail=fail_tot;
 }
}
void integrate_tps_fchk(READ_FCHK_WFN &read_fchk_wfn,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[10], int &fail,int &ncores)
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  void * USERDATA=NULL;
  USERDATA=&read_fchk_wfn;
  define();
  cubacores(ncores,PCORE);
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_tps_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_tps_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_tps_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_tps_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<10;i++)
 {
  result_integration[i]=integral[i];
 }
}
void  integrate_dens_sim(N_FCHKS_WFNS two_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double result_integration[7], int &fail,int &ncores,double ROT_MATRIX[3][3])
{
  int i,j,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  void * USERDATA=NULL;
  //Point void pointer to N_FCHKS_WFNS struct direction.
  USERDATA=&two_fchks_wfns;
  define();
  cubacores(ncores,PCORE);
  for(i=0;i<3;i++)
  {
   for(j=0;j<3;j++)
   {
    Rotate_grid[i][j]=ROT_MATRIX[i][j];
   }
  }
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_divergences, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_divergences, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_divergences, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_divergences, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<7;i++)
 {
  result_integration[i]=integral[i];
 }
}

void  integrate_dens_sim2(N_FCHKS_WFNS two_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double result_integration[7], int &fail,int &ncores,double ROT_MATRIX[3][3])
{
  int i,j,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  void * USERDATA=NULL;
  //Point void pointer to N_FCHKS_WFNS struct direction.
  USERDATA=&two_fchks_wfns;
  define();
  cubacores(ncores,PCORE);
  for(i=0;i<3;i++)
  {
   for(j=0;j<3;j++)
   {
    Rotate_grid[i][j]=ROT_MATRIX[i][j];
   }
  }
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_divergences2, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_divergences2, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_divergences2, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_divergences2, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<7;i++)
 {
  result_integration[i]=integral[i];
 }
}

void  integrate_pol_hyperpol(N_FCHKS_WFNS five_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double result_integration[9], int &fail,int &ncores,char dir_pol_hyper)
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  void * USERDATA=NULL;
  //Point void pointer to N_FCHKS_WFNS struct direction.
  USERDATA=&five_fchks_wfns;
  define();
  cubacores(ncores,PCORE);
  if(dir_pol_hyper=='x' || dir_pol_hyper=='X')
  {DIR_POL_HYPER='x';}
  else if(dir_pol_hyper=='y' || dir_pol_hyper=='Y')
  {DIR_POL_HYPER='y';}
  else 
  {DIR_POL_HYPER='z';}
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_pol_hyper, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_pol_hyper, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_pol_hyper, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_pol_hyper, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<9;i++)
 {
  result_integration[i]=integral[i];
 }
}
void integrate_vr_fchk(READ_FCHK_WFN &read_fchk_wfn,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[2], int &fail,int &ncores,double Point_Vr_in[3])
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  void * USERDATA=NULL;
  USERDATA=&read_fchk_wfn;
  define();
  cubacores(ncores,PCORE);
  for(i=0;i<3;i++){Point_Vr[i]=Point_Vr_in[i];}
 //Regard that for cout and to access void * have to cast the
 //pointer to the corresponding type(int, double, string, object...). See e.g.,
 // int i=5;
 // void *a=NULL;
 // a=&i;  //a points to the direction of i
 // cout<<*((int*)a)<<endl; // The first * access to the data and the (int)*
 //casts to int pointer the void pointer.
 //See bellow the assignation of a pointer to object DENSITY in integration function.
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_vr_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_vr_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_vr_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_vr_fchk, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 result_integration[0]=integral[0];
 result_integration[1]=integral[1];
}
///////////////////////////////////////////////////////////////
//Cuba functions                                             //
//Integration of several quantities(density, shannon, etc.)  //
///////////////////////////////////////////////////////////////
int Integrand(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
/* We get 0<=xx[0],xx[1],xx[2]<=1 values and have to tranform them to the correct interval
   defined in spherical coordinates.
   To do so, we first define xx[0]= phi', x[1]=theta' and xx[1]=r'.
   phi in spherical coordinates goes from 0 to 2*PI then it is straitghtforward to use:
                              phi=2*PI*phi',
   whose differential component is d phi=2*PI d phi'clearly this allows to integrate
   from 0 to 1 as from 0 to 2*PI  int _0 to 2 PI  d phi = 2 PI int _0 to 1 phi' d phi'
   then for theta the integral goes from 0 to PI and int _0 to PI sin(theta) d theta=
   2 (remember the sphere) so we need 2= int _0 to 1  d theta'. This theta' is clearly
                           theta' = [1-cos(theta)]/2,
  so that
                       2 d theta'=  sin(theta) d theta
  and finally, r which goes from 0 to infinity int _0 to 1. To do so, we change
                            r= r' /(1-r')
  so that the differential takes the form:
                           dr = dr' /(1-r')^2.                                           */
 #ifndef f
 #define f ff[0]
 #define f1 ff[1]
 #define f2 ff[2]
 #define f3 ff[3]
 #define f4 ff[4]
 #define f5 ff[5]
 #define f6 ff[6]
 #define f7 ff[7]
 #define f8 ff[8]
 #define f9 ff[9]
 #define f10 ff[10]
 #define f11 ff[11]
 #define f12 ff[12]
 #define f13 ff[13]
 #define f14 ff[14]
 #define f15 ff[15]
 #define f16 ff[16]
 #define f17 ff[17]
 #define f18 ff[18]
 #define f19 ff[19]
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 double res,res_S,res_F,factor_jacobian;
 if((theta>=theta_inf && theta<=theta_sup)&&(phi>=phi_inf&&phi<=phi_sup)&&(r>=r_inf&&r<=r_sup))
 {
  double *AUX,*Grad;
  AUX=new double[*ndim];
  Grad=new double[*ndim];
  READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
  read_fchk_wfn_calc=((READ_FCHK_WFN*)userdata);     //Could be avoided *rho_calc by means of the direct use of
  AUX[0]=r*sin(theta)*cos(phi); //((rho*)userdata)[0].rho_eval(AUX,res);
  AUX[1]=r*sin(theta)*sin(phi);
  AUX[2]=r*cos(theta);
  read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
  read_fchk_wfn_calc[0].rho_grad(AUX,Grad);
  read_fchk_wfn_calc=NULL;
  if(res==ZERO)
  {res_F=pow(TEN,-TEN);}
  else
  {res_F=res;}
  res_F=pow(norm3D(Grad),TWO)/res_F;
  delete[] AUX;delete[] Grad;
  AUX=NULL;Grad=NULL;
 }
 else
 {res=ZERO;res_F=ZERO;}
 ///////////////////////////////////////////////////////
 //Example: Sphere radii r=1                          //
 //cout<<"Exact volumen of the sphere"<<4*PI/3<<endl; //
 //  f = 4*PI*pow(xx[0],2);                           //
 ///////////////////////////////////////////////////////
 //Example 2: 3D-Gaussian                            ////////////
 //Exact interal from -inf to inf f(r) d^3r=PI^(3/2)           //
 //cout<<"Gaussian exact: "<<pow(PI,1.5)<<endl;                //
 // f = 4*PI*exp(-r*r)*xx[0]*xx[0]*pow(1/(1-xx[0]),4);         //
 ////////////////////////////////////////////////////////////////
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Density:
  f = factor_jacobian*res;
 // Shannon:
  if(shanr)
  {
   res_S=res;
   if(res_S<pow(TEN,-TWO*SIX))
   {
    res_S=ZERO;
   }
   else
   {
    res_S=-res_S*log(res_S);
   }
   f1 = factor_jacobian*res_S;
  }
  else
  {f1=ZERO;}
 // Fisher:
  if(fishr)
  {
   f2 = factor_jacobian*res_F;
  }
  else
  {f2=ZERO;}
 // T_TF:
  if(res>pow(TEN,-THREE*FIVE))
  {
   f3 = factor_jacobian*pow(res,FIVE/THREE);
  }
  else
  {
   f3=ZERO;
  }
  if(inertiar)
  {
  // ICMx:
   f4 = factor_jacobian*res*r*sin(theta)*cos(phi);
  // ICMy:
   f5 = factor_jacobian*res*r*sin(theta)*sin(phi);
  // ICMz:
   f6 = factor_jacobian*res*r*cos(theta);
  // Ixx:
   f7 = factor_jacobian*res*(pow(r*sin(theta)*sin(phi),TWO)+pow(r*cos(theta),TWO));
  // Ixy:
   f8 = -factor_jacobian*res*(r*sin(theta)*cos(phi))*(r*sin(theta)*sin(phi));
  // Iyy:
   f9 = factor_jacobian*res*(pow(r*sin(theta)*cos(phi),TWO)+pow(r*cos(theta),TWO));
  // Ixz:
   f10 = -factor_jacobian*res*(r*sin(theta)*cos(phi))*(r*cos(theta));
  // Iyz:
   f11 = -factor_jacobian*res*(r*sin(theta)*sin(phi))*(r*cos(theta));
  // Izz:
   f12 = factor_jacobian*res*(pow(r*sin(theta)*cos(phi),TWO)+pow(r*sin(theta)*sin(phi),TWO));
  }
  else
  {
   f4=ZERO;f5=ZERO;f6=ZERO;f7=ZERO;f8=ZERO;f9=ZERO;f10=ZERO;f11=ZERO;f12=ZERO;
  }
 // <r>:
  if(R1)
  {
   f13 = factor_jacobian*r*res;
  }
  else
  {f13=ZERO;}
 // <r2>:
  if(R2)
  {
   f14 = factor_jacobian*r*r*res;
  }
  else
  {f14=ZERO;}
 // mu_dipolar_moment
  if(DIPOLE)
  {
   // mu x
   f15 = factor_jacobian*res*r*sin(theta)*cos(phi);
   // mu y
   f16 = factor_jacobian*res*r*sin(theta)*sin(phi);
   // mu z
   f17 = factor_jacobian*res*r*cos(theta);
  }
  else
  {f15=ZERO;f16=ZERO;f17=ZERO;}
 // rho
  if(RHO)
  {
   f18 = factor_jacobian*res*res;
  }
  else
  {f18=ZERO;}
 // <r-1>:
  if(RM1)
  {
   f19 = factor_jacobian*res/r;
  }
  else
  {f19=ZERO;}
 return 0;
}
int Integrandp(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
 #ifndef fp
 #define fp ff[0]
 #define fp1 ff[1]
 #define fp2 ff[2]
 #define fp3 ff[3]
 #define fp4 ff[4]
 #define fp5 ff[5]
 #define fp6 ff[6]
 #define fp7 ff[7]
 #define fp8 ff[8]
 #define fp9 ff[9]
 #define fp10 ff[10]
 #define fp11 ff[11]
 #define fp12 ff[12]
 #define fp13 ff[13]
 #define fp14 ff[14]
 #define fp15 ff[15]
 #define fp16 ff[16]
 #define fp17 ff[17]
 #define fp18 ff[18]
 #define fp19 ff[19]
 #define p xx[0]/(1-xx[0])
 #define thetap acos(1-2*xx[1])
 #define phip xx[2]*2*PI
 #endif
 double res,res_S,res_F,factor_jacobian;
 if((thetap>=theta_inf && thetap<=theta_sup)&&(phip>=phi_inf&&phip<=phi_sup)&&(p>=r_inf&&p<=r_sup))
 {
  double *AUX,*Grad;
  AUX=new double[*ndim];
  Grad=new double[*ndim];
  READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
  read_fchk_wfn_calc=((READ_FCHK_WFN*)userdata);     //Could be avoided *rho_calc by means of the direct use of
  AUX[0]=p*sin(thetap)*cos(phip); //((rho*)userdata)[0].rho_eval(AUX,res);
  AUX[1]=p*sin(thetap)*sin(phip);
  AUX[2]=p*cos(thetap);
  read_fchk_wfn_calc[0].rho_p_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
  read_fchk_wfn_calc[0].rho_p_grad(AUX,Grad);
  read_fchk_wfn_calc=NULL;
  if(res==ZERO)
  {res_F=pow(TEN,-TEN);}
  else
  {res_F=res;}
  res_F=pow(norm3D(Grad),TWO)/res_F;
  delete[] AUX;delete[] Grad;
  AUX=NULL;Grad=NULL;
 }
 else
 {res=ZERO;res_F=ZERO;}
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Densityp:
  fp = factor_jacobian*res;
 // Shannonp:
  if(shanp)
  {
   res_S=res;
   if(res_S<pow(TEN,-TWO*SIX))
   {
    res_S=ZERO;
   }
   else
   {
    res_S=-res_S*log(res_S);
   }
   fp1 = factor_jacobian*res_S;
  }
  else
  {fp1=ZERO;}
 // Fisherp:
  if(fishp)
  {
   fp2 = factor_jacobian*res_F;
  }
  else
  {fp2=ZERO;}
 // T_TFp:
  fp3 = ZERO;
  if(inertiap)
  {
  // ICMpx:
   fp4 = factor_jacobian*res*p*sin(thetap)*cos(phip);
  // ICMpy:
   fp5 = factor_jacobian*res*p*sin(thetap)*sin(phip);
  // ICMpz:
   fp6 = factor_jacobian*res*p*cos(thetap);
  // Ipxpx:
   fp7 = factor_jacobian*res*(pow(p*sin(thetap)*sin(phip),TWO)+pow(p*cos(thetap),TWO));
  // Ipxpy:
   fp8 = -factor_jacobian*res*(p*sin(thetap)*cos(phip))*(p*sin(thetap)*sin(phip));
  // Ipypy:
   fp9 = factor_jacobian*res*(pow(p*sin(thetap)*cos(phip),TWO)+pow(p*cos(thetap),TWO));
  // Ipxpz:
   fp10 = -factor_jacobian*res*(p*sin(thetap)*cos(phip))*(p*cos(thetap));
  // Ipypz:
   fp11 = -factor_jacobian*res*(p*sin(thetap)*sin(phip))*(p*cos(thetap));
  // Ipzpz:
   fp12 = factor_jacobian*res*(pow(p*sin(thetap)*cos(phip),TWO)+pow(p*sin(thetap)*sin(phip),TWO));
  }
  else
  {
   fp4=ZERO;fp5=ZERO;fp6=ZERO;fp7=ZERO;fp8=ZERO;fp9=ZERO;fp10=ZERO;fp11=ZERO;fp12=ZERO;
  }
 // <p>:
  if(P1)
  {
   fp13 = factor_jacobian*p*res;
  }
  else
  {fp13=ZERO;}
 // <p2>:
  if(P2)
  {
   fp14 = factor_jacobian*p*p*res;
  }
  else
  {fp14=ZERO;}
 //mup is pointelss
  fp15 = ZERO; fp16 = ZERO; fp17 = ZERO;
 // pi^2 is pointless
  fp18 = ZERO;
 // pi^-1 is pointless
  fp19 = ZERO;
 return 0;
}
int SIJ_Integrand(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
 #ifndef fsij
 #define rsij xx[0]/(1-xx[0])
 #define thetasij acos(1-2*xx[1])
 #define phisij xx[2]*2*PI
 #endif
 int i,j,k;
 double factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 double **AOs_grad,*AOs;
 if((thetasij>=theta_inf && thetasij<=theta_sup)&&(phisij>=phi_inf&&phisij<=phi_sup)&&(rsij>=r_inf&&rsij<=r_sup))
 {
  double *AUX;
  AUX=new double[*ndim];
  READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
  read_fchk_wfn_calc=((READ_FCHK_WFN*)userdata);     //Could be avoided *rho_calc by means of the direct use of
  AUX[0]=rsij*sin(thetasij)*cos(phisij); //((rho*)userdata)[0].rho_eval(AUX,res);
  AUX[1]=rsij*sin(thetasij)*sin(phisij);
  AUX[2]=rsij*cos(thetasij);
  if(!read_fchk_wfn_calc[0].wfn)
  {
   AOs=new double[read_fchk_wfn_calc[0].nbasisf];AOs_grad=new double*[3];
   for(i=0;i<3;i++)
   {
    AOs_grad[i]=new double[read_fchk_wfn_calc[0].nbasisf];
   }
   for(i=0;i<read_fchk_wfn_calc[0].nbasisf;i++)
   {
    AOs[i]=ZERO;AOs_grad[0][i]=ZERO;AOs_grad[1][i]=ZERO;AOs_grad[2][i]=ZERO;
   }
   read_fchk_wfn_calc[0].build_AO_AOgrad2(AOs,AOs_grad,AUX);//We work with position [0] as when rho=new rho[1];
   k=0;
   if(nos)
   {
    nos_array[0]=NO(read_fchk_wfn_calc[0],AOs,AOs_grad,read_fchk_wfn_calc[0].Pair[0]);
    for(j=0;j<read_fchk_wfn_calc[0].nbasis();j++){ff[k]=ZERO;}
    for(j=read_fchk_wfn_calc[0].Pair[0];j<*ncomp;j++)
    {
     nos_array[1]=NO(read_fchk_wfn_calc[0],AOs,AOs_grad,j);
     ff[k]=factor_jacobian*nos_array[0].evaluation*nos_array[1].evaluation;
     k++;
    }
   }
   else
   {
    mos_array[0]=MO(read_fchk_wfn_calc[0],AOs,AOs_grad,read_fchk_wfn_calc[0].Pair[0]);
    for(j=0;j<read_fchk_wfn_calc[0].nbasis();j++){ff[k]=ZERO;}
    for(j=read_fchk_wfn_calc[0].Pair[0];j<*ncomp;j++)
    {
     mos_array[1]=MO(read_fchk_wfn_calc[0],AOs,AOs_grad,j);
     ff[k]=factor_jacobian*mos_array[0].evaluation*mos_array[1].evaluation;
     k++;
    }
   }
   for(i=0;i<3;i++)
   {delete[] AOs_grad[i];AOs_grad[i]=NULL;}
   delete[] AOs_grad;AOs_grad=NULL;
   delete[] AOs;AOs=NULL;
  }
  else
  {
   k=0;
   for(i=0;i<read_fchk_wfn_calc[0].nbasis();i++)
   {
    nos_array[0]=NO(read_fchk_wfn_calc[0],AUX,i);
    for(j=0;j<=i;j++)
    {
     nos_array[1]=NO(read_fchk_wfn_calc[0],AUX,j);
     ff[k]=factor_jacobian*nos_array[0].evaluation*nos_array[1].evaluation;
     k++;
    }
   }
  }
  read_fchk_wfn_calc=NULL;
  delete[] AUX;AUX=NULL;
 }
 else
 {
  for(i=0;i<*ncomp;i++)
  {
   ff[i]=ZERO;
  }
 }
 return 0;
}
int Integrand_tps_fchk(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
 #ifndef ftpsfchk
 #define f ff[0]
 #define f1 ff[1]
 #define f2 ff[2]
 #define f3 ff[3]
 #define f4 ff[4]
 #define f5 ff[5]
 #define f6 ff[6]
 #define f7 ff[7]
 #define f8 ff[8]
 #define f9 ff[9]
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 double res,factor_jacobian;
 double *AUX;
 AUX=new double[*ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)userdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r*sin(theta)*cos(phi); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r*sin(theta)*sin(phi);
 AUX[2]=r*cos(theta);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 read_fchk_wfn_calc=NULL;
 delete[] AUX;AUX=NULL;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Density:
  f = factor_jacobian*res;
 // x:
  f1 = factor_jacobian*res*r*sin(theta)*cos(phi);
 // y:
  f2 = factor_jacobian*res*r*sin(theta)*sin(phi);
 // z:
  f3 = factor_jacobian*res*r*cos(theta);
 // xx:
  f4 = factor_jacobian*res*(pow(r*sin(theta)*cos(phi),TWO));
 // xy:
  f5 = factor_jacobian*res*(r*sin(theta)*cos(phi))*(r*sin(theta)*sin(phi));
 // yy:
  f6 = factor_jacobian*res*(pow(r*sin(theta)*sin(phi),TWO));
 // xz:
  f7 = factor_jacobian*res*(r*sin(theta)*cos(phi))*(r*cos(theta));
 // yz:
  f8 = factor_jacobian*res*(r*sin(theta)*sin(phi))*(r*cos(theta));
 // zz:
  f9 = factor_jacobian*res*(pow(r*cos(theta),TWO));
 return 0;
}
int Integrand_vr_fchk(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
 #ifndef fvr
 #define f ff[0]
 #define f1 ff[1]
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 double res,factor_jacobian;
 double *AUX;
 AUX=new double[*ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)userdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r*sin(theta)*cos(phi); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r*sin(theta)*sin(phi);
 AUX[2]=r*cos(theta);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 read_fchk_wfn_calc=NULL;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 AUX[0]=AUX[0]-Point_Vr[0];
 AUX[1]=AUX[1]-Point_Vr[1];
 AUX[2]=AUX[2]-Point_Vr[2];
 // N
  f = factor_jacobian*res;
 // V(r):
  f1 = factor_jacobian*res/(pow(TEN,-SIX*TWO)+norm3D(AUX));
 delete[] AUX;AUX=NULL;
 return 0;
}
int Integrand_divergences(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
/* We get 0<=xx[0],xx[1],xx[2]<=1 values and have to tranform them to the correct interval
   defined in spherical coordinates.
   To do so, we first define xx[0]= phi', x[1]=theta' and xx[1]=r'.
   phi in spherical coordinates goes from 0 to 2*PI then it is straitghtforward to use:
                              phi=2*PI*phi',
   whose differential component is d phi=2*PI d phi'clearly this allows to integrate
   from 0 to 1 as from 0 to 2*PI  int _0 to 2 PI  d phi = 2 PI int _0 to 1 phi' d phi'
   then for theta the integral goes from 0 to PI and int _0 to PI sin(theta) d theta=
   2 (remember the sphere) so we need 2= int _0 to 1  d theta'. This theta' is clearly
                           theta' = [1-cos(theta)]/2,
  so that
                       2 d theta'=  sin(theta) d theta
  and finally, r which goes from 0 to infinity int _0 to 1. To do so, we change
                            r= r' /(1-r')
  so that the differential takes the form:
                           dr = dr' /(1-r')^2.                                           */
 #ifndef f
 #define f ff[0]
 #define f1 ff[1]
 #define f2 ff[2]
 #define f3 ff[3]
 #define f4 ff[4]
 #define f5 ff[5]
 #define f6 ff[6]
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 int i,j;
 double res1,res2,res1_2,res2_2,factor_jacobian,nelec_1,nelec_2;
 double *AUX,*AUX2;
 AUX=new double[*ndim];
 AUX2=new double[*ndim];
 N_FCHKS_WFNS *n_fchk_wfn_calc;                 //Pointer to an object of N_FCHKS_WFNS type.
 n_fchk_wfn_calc=((N_FCHKS_WFNS*)userdata);
 AUX[0]=r*sin(theta)*cos(phi);
 AUX[1]=r*sin(theta)*sin(phi);
 AUX[2]=r*cos(theta);
 //Rotate point for second FCHK/WFN/WFX file. To make the Inertia tensors coincide by means of the grid rotation.
 for(i=0;i<3;i++)
 {
  AUX2[i]=ZERO;
 }
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   AUX2[i]=AUX2[i]+Rotate_grid[i][j]*AUX[j];
  }
 }
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[0])).rho_eval(AUX,res1);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[1])).rho_eval(AUX2,res2);
 nelec_1=(*((n_fchk_wfn_calc[0]).read_fchk_wfn[0])).nelectrons;
 nelec_2=(*((n_fchk_wfn_calc[0]).read_fchk_wfn[1])).nelectrons;
 n_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 delete[] AUX2;
 AUX2=NULL;
 res1_2=res1*res1;
 res2_2=res2*res2;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Density 1:
  f = factor_jacobian*res1;
 // Density 2:
  f1 = factor_jacobian*res2;
 // Density^2 1:
  f2 = factor_jacobian*res1_2;
 // Density^2 2:
  f3 = factor_jacobian*res2_2;
 // Density 1 Density 2:
  f4 = factor_jacobian*res1*res2;
 // (Density_1/N1 - Density_2/N2)^2 :
  f5 = factor_jacobian*pow(res1/(double)nelec_1-res2/(double)nelec_2,TWO);
 //  Density_1 ln [ Density_1 / Density_2 ] :
  if(res2>pow(TEN,-TEN))
  {
   f6 = factor_jacobian*res1*log(res1/res2);
  }
  else
  {
   f6 = ZERO;
  }
 return 0;
}
int Integrand_divergences2(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
/* We get 0<=xx[0],xx[1],xx[2]<=1 values and have to tranform them to the correct interval
   defined in spherical coordinates.
   To do so, we first define xx[0]= phi', x[1]=theta' and xx[1]=r'.
   phi in spherical coordinates goes from 0 to 2*PI then it is straitghtforward to use:
                              phi=2*PI*phi',
   whose differential component is d phi=2*PI d phi'clearly this allows to integrate
   from 0 to 1 as from 0 to 2*PI  int _0 to 2 PI  d phi = 2 PI int _0 to 1 phi' d phi'
   then for theta the integral goes from 0 to PI and int _0 to PI sin(theta) d theta=
   2 (remember the sphere) so we need 2= int _0 to 1  d theta'. This theta' is clearly
                           theta' = [1-cos(theta)]/2,
  so that
                       2 d theta'=  sin(theta) d theta
  and finally, r which goes from 0 to infinity int _0 to 1. To do so, we change
                            r= r' /(1-r')
  so that the differential takes the form:
                           dr = dr' /(1-r')^2.                                           */
 #ifndef f
 #define f ff[0]
 #define f1 ff[1]
 #define f2 ff[2]
 #define f3 ff[3]
 #define f4 ff[4]
 #define f5 ff[5]
 #define f6 ff[6]
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 int i,j;
 double res1,res2,res1_2,res2_2,factor_jacobian,nelec_1,nelec_2;
 double *AUX,*AUX2,*Grad1,*Grad2;
 AUX=new double[*ndim];
 AUX2=new double[*ndim];
 Grad1=new double[*ndim];
 Grad2=new double[*ndim];
 N_FCHKS_WFNS *n_fchk_wfn_calc;                 //Pointer to an object of N_FCHKS_WFNS type.
 n_fchk_wfn_calc=((N_FCHKS_WFNS*)userdata);
 AUX[0]=r*sin(theta)*cos(phi);
 AUX[1]=r*sin(theta)*sin(phi);
 AUX[2]=r*cos(theta);
 //Rotate point for second FCHK/WFN/WFX file. To make the Inertia tensors coincide by means of the grid rotation.
 for(i=0;i<3;i++)
 {
  AUX2[i]=ZERO;
 }
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   AUX2[i]=AUX2[i]+Rotate_grid[i][j]*AUX[j];
  }
 }
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[0])).rho_eval(AUX,res1);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[0])).rho_grad(AUX,Grad1);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[1])).rho_eval(AUX2,res2);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[1])).rho_grad(AUX2,Grad2);
 res1=norm3D(Grad1);
 res2=norm3D(Grad2);
 nelec_1=(*((n_fchk_wfn_calc[0]).read_fchk_wfn[0])).nelectrons;
 nelec_2=(*((n_fchk_wfn_calc[0]).read_fchk_wfn[1])).nelectrons;
 n_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 delete[] AUX2;
 AUX2=NULL;
 delete[] Grad1;
 Grad1=NULL;
 delete[] Grad2;
 Grad2=NULL;
 res1_2=res1*res1;
 res2_2=res2*res2;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Red Grad 1:
  f = factor_jacobian*res1;
 // Red Grad 2:
  f1 = factor_jacobian*res2;
 // Red Grad^2 1:
  f2 = factor_jacobian*res1_2;
 // Red Grad^2 2:
  f3 = factor_jacobian*res2_2;
 // Red Grad 1 Red Grad 2:
  f4 = factor_jacobian*res1*res2;
 // (Red Grad 1*N1^1/3 - Red Grad 2*N2^1/3)^2 :
  f5 = factor_jacobian*pow(res1/((double)nelec_1)-res2/((double)nelec_2),TWO);
 //  Red Grad 1 ln [ Red Grad 1 / Red Grad 2 ] :
  if(res2>pow(TEN,-TEN))
  {
   f6 = factor_jacobian*res1*log(res1/res2);
  }
  else
  {
   f6 = ZERO;
  }
 return 0;
}

int Integrand_pol_hyper(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
/* We get 0<=xx[0],xx[1],xx[2]<=1 values and have to tranform them to the correct interval
   defined in spherical coordinates.
   To do so, we first define xx[0]= phi', x[1]=theta' and xx[1]=r'.
   phi in spherical coordinates goes from 0 to 2*PI then it is straitghtforward to use:
                              phi=2*PI*phi',
   whose differential component is d phi=2*PI d phi'clearly this allows to integrate
   from 0 to 1 as from 0 to 2*PI  int _0 to 2 PI  d phi = 2 PI int _0 to 1 phi' d phi'
   then for theta the integral goes from 0 to PI and int _0 to PI sin(theta) d theta=
   2 (remember the sphere) so we need 2= int _0 to 1  d theta'. This theta' is clearly
                           theta' = [1-cos(theta)]/2,
  so that
                       2 d theta'=  sin(theta) d theta
  and finally, r which goes from 0 to infinity int _0 to 1. To do so, we change
                            r= r' /(1-r')
  so that the differential takes the form:
                           dr = dr' /(1-r')^2.                                           */
 #ifndef f
 #define f ff[0]
 #define f1 ff[1]
 #define f2 ff[2]
 #define f3 ff[3]
 #define f4 ff[4]
 #define f5 ff[5]
 #define f6 ff[6]
 #define f7 ff[7]
 #define f8 ff[8]
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 double res1,res2,res3,res4,res5,alpha,gamma,dir,factor_jacobian;
 double *AUX;
 AUX=new double[*ndim];
 N_FCHKS_WFNS *n_fchk_wfn_calc;                 //Pointer to an object of N_FCHKS_WFNS type.
 n_fchk_wfn_calc=((N_FCHKS_WFNS*)userdata);
 AUX[0]=r*sin(theta)*cos(phi);
 AUX[1]=r*sin(theta)*sin(phi);
 AUX[2]=r*cos(theta);
 if(DIR_POL_HYPER=='x')
 {dir=AUX[0];}
 else if(DIR_POL_HYPER=='y')
 {dir=AUX[1];}
 else
 {dir=AUX[2];}
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[0])).rho_eval(AUX,res1);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[1])).rho_eval(AUX,res2);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[2])).rho_eval(AUX,res3);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[3])).rho_eval(AUX,res4);
 (*((n_fchk_wfn_calc[0]).read_fchk_wfn[4])).rho_eval(AUX,res5);
 n_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Density 1:
  f = factor_jacobian*res1;
 // Density 2:
  f1 = factor_jacobian*res2;
 // Density 3:
  f2 = factor_jacobian*res3;
 // Density 4:
  f3 = factor_jacobian*res4;
 // Density 5:
  f4 = factor_jacobian*res5;
 // alpha:
  alpha=res3-res2; 
  f5 = factor_jacobian*alpha*dir;
 // gamma: 
  gamma=res5-res4-TWO*alpha; 
  f6 = factor_jacobian*gamma*dir;
 // QD1:
  f7 = factor_jacobian*alpha*alpha;
 // QD2: 
  f8 = factor_jacobian*gamma*gamma;
 return 0;
}
///////////////////////
//Cubature           //
///////////////////////
void I_cubature(READ_FCHK_WFN &read_fchk_wfn,string Operation,double presc,double &result_integration,
double &error,double Integrals_interval[6])
{
  double val,quasi_zero,quasi_one;
  void * USERDATA=NULL;
  USERDATA=&read_fchk_wfn;
  quasi_one=1.0e00-1.0e-10;
  quasi_zero=1.0e-10;
  shanr=true;shanp=true;fishp=true;fishr=true;R1=true;R2=true;
  P1=true;P2=true;DIPOLE=true;RHO=true;
  bool interval=false;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval=false;}
    else{interval=true;}
   }
   else{interval=true;}
  }
  else{interval=true;}
  if(interval)
  {
   define_interval(Integrals_interval);
  }
  else
  {
   theta_inf=ZERO;theta_sup=PI;phi_inf=ZERO;phi_sup=TWO*PI;r_inf=ZERO;
   r_sup=RSUP;
  }
  double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
  double xmax[3] = {quasi_one,quasi_one,quasi_one};
  if(Operation=="Density")
  {
   hcubature(1, N_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Densityp")
  {
   hcubature(1, Np_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="R1")
  {
   hcubature(1, r1_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="R2")
  {
   hcubature(1, r2_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="P1")
  {
   hcubature(1, p1_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="P2")
  {
   hcubature(1, p2_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="rho")
  {
   hcubature(1, rho_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Shannon")
  {
   hcubature(1, S_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Shannonp")
  {
   hcubature(1, Sp_Cubature,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Fisher")
  {
   hcubature(1, FI_Cubature,USERDATA,  3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Fisherp")
  {
   hcubature(1, FIp_Cubature,USERDATA,  3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="T_TF")
  {
   hcubature(1, TTF_Cubature,USERDATA,  3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Sij_NO")
  {
   hcubature(1, Sij_NO_Cubature,USERDATA, 3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Sij_MO")
  {
   hcubature(1, Sij_MO_Cubature,USERDATA, 3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else{}
  result_integration=val;
}
//Cubature functions:
int N_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 //Density cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int Np_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_p_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 //Density cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int r1_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!R1){res=ZERO;}
 fval[0] = FOUR*PI*res*r_cubature*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int r2_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!R2){res=ZERO;}
 fval[0] = FOUR*PI*r_cubature*r_cubature*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int p1_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_p_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!P1){res=ZERO;}
 fval[0] = FOUR*PI*r_cubature*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int p2_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_p_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!P2){res=ZERO;}
 fval[0] = FOUR*PI*r_cubature*r_cubature*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int rho_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!RHO){res=ZERO;}
 //Density^2 cubature:
 fval[0] = FOUR*PI*res*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int S_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)&&
 (phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 if(res<pow(TEN,-TWO*SIX)){res=ZERO;}
 if(res>=pow(TEN,-TWO*SIX))
 {
  res=-res*log(res);
 }
 read_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 }
 else{res=ZERO;}
 if(!shanr){res=ZERO;}
 //Shannon cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int Sp_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)&&
 (phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_p_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 if(res<pow(TEN,-TWO*SIX)){res=ZERO;}
 if(res>=pow(TEN,-TWO*SIX))
 {
  res=-res*log(res);
 }
 read_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 }
 else{res=ZERO;}
 if(!shanp){res=ZERO;}
 //ShannonP cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int FI_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)&&
 (phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX,*Grad,eval;
 AUX=new double[ndim];
 Grad=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,eval);//We work with position [0] as when rho=new rho[1];
 read_fchk_wfn_calc[0].rho_grad(AUX,Grad);
 if(eval==ZERO) eval=pow(TEN,-TEN);
 res=pow(norm3D(Grad),TWO)/eval;
 delete[] AUX;
 delete[] Grad;
 AUX=NULL;
 Grad=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!fishr){res=ZERO;}
 //Fisher cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int FIp_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)&&
 (phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX,*Grad,eval;
 AUX=new double[ndim];
 Grad=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_p_eval(AUX,eval);//We work with position [0] as when rho=new rho[1];
 read_fchk_wfn_calc[0].rho_p_grad(AUX,Grad);
 if(eval==ZERO) eval=pow(TEN,-TEN);
 res=pow(norm3D(Grad),TWO)/eval;
 delete[] AUX;
 delete[] Grad;
 AUX=NULL;
 Grad=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!fishp){res=ZERO;}
 //Fisher cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int TTF_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)&&
 (phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX,eval;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,eval);//We work with position [0] as when rho=new rho[1];
 if(eval>pow(TEN,-THREE*FIVE))
 {
  res=pow(eval,FIVE/THREE);
 }
 else
 {
  res=ZERO;
 }
 read_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 }
 else{res=ZERO;}
 //Thomas Fermi cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int Sij_NO_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 nos_array[0]=NO(read_fchk_wfn_calc[0],AUX,read_fchk_wfn_calc[0].Pair[0]);
 nos_array[1]=NO(read_fchk_wfn_calc[0],AUX,read_fchk_wfn_calc[0].Pair[1]);
 res=nos_array[0].evaluation*nos_array[1].evaluation;
 read_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 }
 else{res=ZERO;}
 //Sij cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int Sij_MO_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 mos_array[0]=MO(read_fchk_wfn_calc[0],AUX,read_fchk_wfn_calc[0].Pair[0]);
 mos_array[1]=MO(read_fchk_wfn_calc[0],AUX,read_fchk_wfn_calc[0].Pair[1]);
 res=mos_array[0].evaluation*mos_array[1].evaluation;
 read_fchk_wfn_calc=NULL;
 delete[] AUX;
 AUX=NULL;
 }
 else{res=ZERO;}
 //Sij cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
///////////////////////
//Cubature 2         //
///////////////////////
void I_cubature2(READ_FCHK_WFN &read_fchk_wfn,double presc,double *result_integration,
double error[10],double Integrals_interval[10],bool shan, bool fish, bool r1, bool r2,bool dipolar,bool rho)
{
 int i;
 double val,error_int,quasi_zero,quasi_one;
 READ_FCHK_WFN rho_copies[10]={read_fchk_wfn};
 quasi_one=1.0e00-1.0e-10;
 quasi_zero=1.0e-10;
 R1=r1;R2=r2;
 shanr=shan;
 fishr=fish;
 DIPOLE=dipolar;
 RHO=rho;
 bool interval=false;
 if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP))
 {
  if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
  {
   if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
   {interval=false;}
    else{interval=true;}
  }
  else{interval=true;}
 }
 else{interval=true;}
 if(interval)
 {
  define_interval(Integrals_interval);
 }
 else
 {
  theta_inf=ZERO;theta_sup=PI;phi_inf=ZERO;phi_sup=TWO*PI;r_inf=ZERO;
  r_sup=RSUP;
 }
 double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
 double xmax[3] = {quasi_one,quasi_one,quasi_one};
 omp_set_dynamic(0);     // Explicitly disable dynamic teams
 omp_set_num_threads(10); // Use 10 threads for all consecutive parallel regions
 #pragma omp parallel private(i) firstprivate(xmin,xmax) shared(rho_copies)
 {
  i=omp_get_thread_num();
  if(i==0)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, N_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==1)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, S_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==2)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, FI_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==3)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, TTF_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==4)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, r1_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==5)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, r2_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==6)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, mux_Cubature2,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==7)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, muy_Cubature2,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==8)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, muz_Cubature2,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, rho_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
 }
}
void I_cubature2p(READ_FCHK_WFN &read_fchk_wfn,double presc,double *result_integration,
double error[10],double Integrals_interval[10],bool shan, bool fish, bool p1, bool p2)
{
 int i;
 double val,error_int,quasi_zero,quasi_one;
 READ_FCHK_WFN rho_copies[6]={read_fchk_wfn};
 quasi_one=1.0e00-1.0e-10;
 quasi_zero=1.0e-10;
 P1=p1;P2=p2;
 shanp=shan;
 fishp=fish;
 bool interval=false;
 if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP))
 {
  if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
  {
   if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
   {interval=false;}
   else{interval=true;}
  }
  else{interval=true;}
 }
 else{interval=true;}
 if(interval)
 {
  define_interval(Integrals_interval);
 }
 else
 {
  theta_inf=ZERO;theta_sup=PI;phi_inf=ZERO;phi_sup=TWO*PI;r_inf=ZERO;
  r_sup=RSUP;
 }
 double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
 double xmax[3] = {quasi_one,quasi_one,quasi_one};
 omp_set_dynamic(0);     // Explicitly disable dynamic teams
 omp_set_num_threads(6); // Use 6 threads for all consecutive parallel regions (ony six we don't need 9!)
 #pragma omp parallel private(i) firstprivate(xmin,xmax) shared(rho_copies)
 {
  i=omp_get_thread_num();
  if(i==0)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, Np_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==1)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, Sp_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==2)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, FIp_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else if(i==3)
  {
   #pragma omp critical
   {
   result_integration[i]=ZERO;
   error[i]=ZERO;
   }
  }
  else if(i==4)
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, p1_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
  else
  {
   void * USERDATA=NULL;
   USERDATA=&rho_copies[i];
   hcubature(1, p2_Cubature,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error_int);
   #pragma omp critical
   {
   result_integration[i]=val;
   error[i]=error_int;
   }
  }
 }
}
int mux_Cubature2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!DIPOLE){res=ZERO;}
 fval[0] = FOUR*PI*res*r_cubature*sin(theta_cubature)*cos(phi_cubature)*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int muy_Cubature2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!DIPOLE){res=ZERO;}
 fval[0] = FOUR*PI*res*r_cubature*sin(theta_cubature)*sin(phi_cubature)*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int muz_Cubature2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf && theta_cubature<=theta_sup)
 &&(phi_cubature>=phi_inf&&phi_cubature<=phi_sup)&&(r_cubature>=r_inf&&r_cubature<=r_sup))
 {
 double *AUX;
 AUX=new double[ndim];
 READ_FCHK_WFN *read_fchk_wfn_calc; //Pointer to an object of rho type.
 read_fchk_wfn_calc=((READ_FCHK_WFN*)fdata);     //Could be avoided *rho_calc by means of the direct use of
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((rho*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 read_fchk_wfn_calc[0].rho_eval(AUX,res);//We work with position [0] as when rho=new rho[1];
 delete[] AUX;
 AUX=NULL;
 read_fchk_wfn_calc=NULL;
 }
 else{res=ZERO;}
 if(!DIPOLE){res=ZERO;}
 fval[0] = FOUR*PI*res*r_cubature*cos(theta_cubature)*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
