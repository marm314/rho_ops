#include"Integrals_DMN.h"
#include"legendre_quadrature.h"
#include"sphere_lebedev_rule.h"
//Global variables readable everywhere and defining the interval of integration
int p_ggrid;
double *p_grid,*w_p_grid,*w_ang,*thetau,*phiu;
double Point_intra[3],Point_extra[3];
double theta_inf_dmn,theta_sup_dmn,phi_inf_dmn,phi_sup_dmn,r_inf_dmn,r_sup_dmn,r_intra;
bool interval_dmn,shan_dmn,fishr_dmn,inertiar_dmn,shanp_dmn,fishp_dmn,inertiap_dmn,DIPOLE_dmn;
bool r1_dmn,r2_dmm,p1_dmn,p2_dmn;
////Defines for CUBA
void define_dmn()
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
void define_interval_dmn(double Interval[6])
{
 r_inf_dmn=Interval[0];
 r_sup_dmn=Interval[1];
 theta_inf_dmn=Interval[2];
 theta_sup_dmn=Interval[3];
 phi_inf_dmn=Interval[4];
 phi_sup_dmn=Interval[5];
}
//////////////////////////
//Functions description //
//////////////////////////
void integrate_dmnr(DMN_OPS &DMNr,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[18], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool r1,bool r2,bool dipole)
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  interval_dmn=false;
  shan_dmn=shan;
  fishr_dmn=fish;
  inertiar_dmn=inertias;
  r1_dmn=r1;
  r2_dmm=r2;
  DIPOLE_dmn=dipole;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP_DMN))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval_dmn=false;}
    else{interval_dmn=true;}
   }
   else{interval_dmn=true;}
  }
  else{interval_dmn=true;}
  if(interval_dmn)
  {
   define_interval_dmn(Integrals_interval);
  }
  else
  {
   theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
   r_sup_dmn=RSUP_DMN;
  }
  void * USERDATA=NULL;
  USERDATA=&DMNr;
  define_dmn();
  cubacores(nprocs,PCORE);
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
  Cuhre(NDIM, NCOMP, Integrand_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<18;i++)
 {result_integration[i]=integral[i];}
}
void integrate_dmnp(DMN_P_OPS &DMNp,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[18], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool p1,bool p2)
{
  int i,verbose=0,nregions,neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  interval_dmn=false;
  shanp_dmn=shan;
  fishp_dmn=fish;
  inertiap_dmn=inertias;
  p1_dmn=p1;
  p2_dmn=p2;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP_DMN))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval_dmn=false;}
    else{interval_dmn=true;}
   }
   else{interval_dmn=true;}
  }
  else{interval_dmn=true;}
  if(interval_dmn)
  {
   define_interval_dmn(Integrals_interval);
  }
  else
  {
   theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
   r_sup_dmn=RSUP_DMN;
  }
  void * USERDATA=NULL;
  USERDATA=&DMNp;
  define_dmn();
  cubacores(nprocs,PCORE);
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
  Cuhre(NDIM, NCOMP, Integrandp_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrandp_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrandp_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrandp_dmn, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 for(i=0;i<18;i++)
 {result_integration[i]=integral[i];}
}
void integrate_intra(DMN_OPS &DMNr,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double &result_integration, int &fail,int &nprocs,double POINT_INTRA[3])
{
 int i,verbose=0,nregions,neval;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
 r_sup_dmn=RSUP_DMN;
 for(i=0;i<3;i++){Point_intra[i]=POINT_INTRA[i];}
 void * USERDATA=NULL;
 USERDATA=&DMNr;
 define_dmn();
 cubacores(nprocs,PCORE);
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_intra, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_intra, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_intra, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_intra, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 result_integration=integral[0];
}
void integrate_extra(DMN_OPS &DMNr,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double &result_integration, int &fail,int &nprocs,double POINT_EXTRA[3])
{
 int i,verbose=0,nregions,neval;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
 r_sup_dmn=RSUP_DMN;
 for(i=0;i<3;i++){Point_extra[i]=POINT_EXTRA[i];}
 void * USERDATA=NULL;
 USERDATA=&DMNr;
 define_dmn();
 cubacores(nprocs,PCORE);
 if(method=="Cuhre")
 {
  Cuhre(NDIM, NCOMP, Integrand_extra, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
  MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 }
 else if(method=="Divone")
 {
  Divonne(NDIM, NCOMP, Integrand_extra, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,
  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,BORDER, MAXCHISQ, MINDEVIATION,
  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,STATEFILE,SPIN,&nregions, &neval, &fail,
  integral, error, prob);
 }
 else if(method=="Suave")
 {
  Suave(NDIM, NCOMP, Integrand_extra, USERDATA, NVEC1,EPSREL, EPSABS, verbose | LAST, SEED,
  MINEVAL, MAXEVAL, NNEW, FLATNESS,STATEFILE,SPIN,&nregions, &neval, &fail, integral,
  error, prob);
 }
 else
 {
  Vegas(NDIM, NCOMP, Integrand_extra, USERDATA, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
  MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
  error, prob);
 }
 result_integration=integral[0];
}
///////////////////////////////////////////////////////////////
//Cuba functions                                             //
//Integration of several quantities(density, shannon, etc.)  //
///////////////////////////////////////////////////////////////
int Integrand_dmn(const int *ndim, const double xx[],const int *ncomp, double ff[],
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
 #define r xx[0]/(1-xx[0])
 #define theta acos(1-2*xx[1])
 #define phi xx[2]*2*PI
 #endif
 double res,res_S,res_F,factor_jacobian;
 if((theta>=theta_inf_dmn && theta<=theta_sup_dmn)&&(phi>=phi_inf_dmn&&phi<=phi_sup_dmn)&&(r>=r_inf_dmn&&r<=r_sup_dmn))
 {
  double *AUX,*Grad;
  AUX=new double[*ndim];
  Grad=new double[*ndim];
  DMN_OPS *dmn; //Pointer to an object of DMN_OPS type.
  dmn=((DMN_OPS*)userdata);
  AUX[0]=r*sin(theta)*cos(phi);
  AUX[1]=r*sin(theta)*sin(phi);
  AUX[2]=r*cos(theta);
  dmn[0].grad_rho_r(AUX,Grad,res);//We work with position [0]
  if(res==ZERO)
  {res_F=pow(TEN,-TEN);}
  else
  {res_F=res;}
  res_F=pow(norm3D(Grad),TWO)/res_F;
  delete[] AUX;delete[] Grad;
  Grad=NULL;AUX=NULL;dmn=NULL;
 }
 else
 {res=ZERO;}
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
  // Density:
  f = factor_jacobian*res;
 // Shannon:
  if(shan_dmn)
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
  if(fishr_dmn)
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
   f3 = ZERO;
  }
 if(inertiar_dmn)
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
  if(r1_dmn)
  {
   f13 = factor_jacobian*r*res;
  }
  else
  {f13=ZERO;}
 // <r2>:
  if(r2_dmm)
  {
   f14 = factor_jacobian*r*r*res;
  }
  else
  {f14=ZERO;}
  // mu_dipolar_moment
  if(DIPOLE_dmn)
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
 return 0;
}
int Integrandp_dmn(const int *ndim, const double xx[],const int *ncomp, double ff[],
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
 #define p xx[0]/(1-xx[0])
 #define thetap acos(1-2*xx[1])
 #define phip xx[2]*2*PI
 #endif
 double res,res_S,res_F,factor_jacobian;
 if((thetap>=theta_inf_dmn && thetap<=theta_sup_dmn)&&(phip>=phi_inf_dmn&&phip<=phi_sup_dmn)&&(p>=r_inf_dmn&&p<=r_sup_dmn))
 {
  double *AUX,*Grad;
  AUX=new double[*ndim];
  Grad=new double[*ndim];
  DMN_P_OPS *dmn; //Pointer to an object.
  dmn=((DMN_P_OPS*)userdata);
  AUX[0]=p*sin(thetap)*cos(phip);
  AUX[1]=p*sin(thetap)*sin(phip);
  AUX[2]=p*cos(thetap);
  dmn[0].grad_rho_p(AUX,Grad,res);
  if(res==ZERO)
  {res_F=pow(TEN,-TEN);}
  else
  {res_F=res;}
  res_F=pow(norm3D(Grad),TWO)/res_F;
  delete[] AUX;delete[] Grad;
  AUX=NULL;dmn=NULL;Grad=NULL;
 }
 else
 {res=ZERO;}
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Densityp:
  fp = factor_jacobian*res;
 // Shannonp:
  if(shanp_dmn)
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
  if(fishp_dmn)
  {
   fp2 = factor_jacobian*res_F;
  }
  else
  {fp2=ZERO;}
 // T_TFp:
  fp3 = ZERO;
  if(inertiap_dmn)
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
  if(p1_dmn)
  {
   fp13 = factor_jacobian*p*res;
  }
  else
  {fp13=ZERO;}
 // <p2>:
  if(p2_dmn)
  {
   fp14 = factor_jacobian*p*p*res;
  }
  else
  {fp14=ZERO;}
 //mup is pointelss
  fp15 = ZERO; fp16 = ZERO; fp17 = ZERO;
 return 0;
}
int Integrand_intra(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
 #ifndef fintra
 #define fintra ff[0]
 #define rintra xx[0]/(1-xx[0])
 #define thetaintra acos(1-2*xx[1])
 #define phiintra xx[2]*2*PI
 #endif
 int i;
 double res,factor_jacobian;
 double *AUX,AUX2[3];
 AUX=new double[*ndim];
 DMN_OPS *dmn; //Pointer to an object of DMN_OPS type.
 dmn=((DMN_OPS*)userdata);
 AUX[0]=rintra*sin(thetaintra)*cos(phiintra);
 AUX[1]=rintra*sin(thetaintra)*sin(phiintra);
 AUX[2]=rintra*cos(thetaintra);
 for(i=0;i<3;i++)
 {
  AUX2[i]=AUX[i]+Point_intra[i]/TWO;
  AUX[i]=AUX[i]-Point_intra[i]/TWO;
 }
 res=dmn[0].dm2_coalesc(AUX,AUX2);//We work with position [0]
 delete[] AUX;
 AUX=NULL;dmn=NULL;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Intracule:
 fintra = factor_jacobian*res;
 return 0;
}
int Integrand_extra(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata)
{
 #ifndef fextra
 #define fextra ff[0]
 #define rextra xx[0]/(1-xx[0])
 #define thetaextra acos(1-2*xx[1])
 #define phiextra xx[2]*2*PI
 #endif
 int i;
 double res,factor_jacobian;
 double *AUX,AUX2[3];
 AUX=new double[*ndim];
 DMN_OPS *dmn; //Pointer to an object of DMN_OPS type.
 dmn=((DMN_OPS*)userdata);
 AUX[0]=rextra*sin(thetaextra)*cos(phiextra);
 AUX[1]=rextra*sin(thetaextra)*sin(phiextra);
 AUX[2]=rextra*cos(thetaextra);
 for(i=0;i<3;i++)
 {
  AUX2[i]=Point_extra[i]+AUX[i]/TWO;
  AUX[i]=Point_extra[i]-AUX[i]/TWO;
 }
 res=dmn[0].dm2_coalesc(AUX,AUX2);//We work with position [0]
 delete[] AUX;
 AUX=NULL;dmn=NULL;
 factor_jacobian=FOUR*PI*xx[0]*xx[0]*pow(ONE/(ONE-xx[0]),FOUR);
 // Extracule:
 fextra = factor_jacobian*res;
 return 0;
}
///////////////////////
//Cubature           //
///////////////////////
void I_cubature_dmnr(DMN_OPS &DMNr,string Operation,double presc,double &result_integration,
double &error,double Integrals_interval[6])
{ double val,quasi_zero,quasi_one;
  void * USERDATA=NULL;
  USERDATA=&DMNr;
  quasi_one=1.0e00-1.0e-10;
  quasi_zero=1.0e-10;
  int nelectrons=DMNr.FCHK_for_DMN[0].nelectrons;
  interval_dmn=false;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP_DMN))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval_dmn=false;}
    else{interval_dmn=true;}
   }
   else{interval_dmn=true;}
  }
  else{interval_dmn=true;}
  if(interval_dmn)
  {
   define_interval_dmn(Integrals_interval);
  }
  else
  {
   theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
   r_sup_dmn=RSUP_DMN;
  }
  double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
  double xmax[3] = {quasi_one,quasi_one,quasi_one};
  if(Operation=="Density")
  {
   hcubature(1, N_Cubature_dmnr,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="R1")
  {
   hcubature(1, r1_Cubature_dmnr,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="R2")
  {
   hcubature(1, r2_Cubature_dmnr,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Shannon")
  {
   hcubature(1, S_Cubature_dmnr,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
   val=val/(double)nelectrons+log((double)nelectrons);
  }
  else if(Operation=="Fisher")
  {
   hcubature(1, FI_Cubature_dmnr,USERDATA, 3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
   val=val/(double)nelectrons;
  }
  else if(Operation=="T_TF")
  {
   hcubature(1, TTF_Cubature_dmnr,USERDATA, 3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else{}
  result_integration=val;
}
void I_cubature_dmnp(DMN_P_OPS &DMNp,string Operation,double presc,double &result_integration,
double &error,double Integrals_interval[6])
{ double val,quasi_zero,quasi_one;
  void * USERDATA=NULL;
  USERDATA=&DMNp;
  quasi_one=1.0e00-1.0e-10;
  quasi_zero=1.0e-10;
  int  nelectrons=DMNp.FCHK_for_DMN[0].nelectrons;
  interval_dmn=false;
  if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP_DMN))
  {
   if((Integrals_interval[2]==ZERO && abs(Integrals_interval[3]-PI)<pow(TEN,-EIGHT)))
   {
    if((Integrals_interval[4]==ZERO && abs(Integrals_interval[5]-TWO*PI)<pow(TEN,-EIGHT)))
    {interval_dmn=false;}
    else{interval_dmn=true;}
   }
   else{interval_dmn=true;}
  }
  else{interval_dmn=true;}
  if(interval_dmn)
  {
   define_interval_dmn(Integrals_interval);
  }
  else
  {
   theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
   r_sup_dmn=RSUP_DMN;
  }
  double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
  double xmax[3] = {quasi_one,quasi_one,quasi_one};
  if(Operation=="Densityp")
  {
   hcubature(1, Np_Cubature_dmnp,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="P1")
  {
   hcubature(1, p1_Cubature_dmnp,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="P2")
  {
   hcubature(1, p2_Cubature_dmnp,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
  }
  else if(Operation=="Shannonp")
  {
   hcubature(1, Sp_Cubature_dmnp,USERDATA, 3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
   val=val/(double)nelectrons+log((double)nelectrons);
  }
  else if(Operation=="Fisherp")
  {
   hcubature(1, FIp_Cubature_dmnp,USERDATA, 3, xmin, xmax, 0,1.0e-4,presc, ERROR_INDIVIDUAL,
   &val, &error);
   val=val/(double)nelectrons;
  }
  else{}
  result_integration=val;
}
int N_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_OPS *dmn;
 dmn=((DMN_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation(AUX,AUX);//We work with position [0]
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Density cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int Np_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_P_OPS *dmn;
 dmn=((DMN_P_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation_p(AUX,AUX);//We work with position [0]
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Densityp cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int r1_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_OPS *dmn;
 dmn=((DMN_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation(AUX,AUX);//We work with position [0]
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 fval[0] = FOUR*PI*res*r_cubature*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int r2_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_OPS *dmn;
 dmn=((DMN_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation(AUX,AUX);//We work with position [0]
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 fval[0] = FOUR*PI*r_cubature*r_cubature*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int p1_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_P_OPS *dmn;
 dmn=((DMN_P_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation_p(AUX,AUX);//We work with position [0]
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 fval[0] = FOUR*PI*r_cubature*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int p2_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_P_OPS *dmn;
 dmn=((DMN_P_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation_p(AUX,AUX);//We work with position [0]
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 fval[0] = FOUR*PI*r_cubature*r_cubature*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int S_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)&&
 (phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_OPS *dmn;
 dmn=((DMN_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation(AUX,AUX);//We work with position [0]
 if(res<pow(TEN,-TWO*SIX)){res=ZERO;}
 if(res>=pow(TEN,-TWO*SIX))
 {
  res=-res*log(res);
 }
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Shannon cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int Sp_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)&&
 (phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX;
 AUX=new double[ndim];
 DMN_P_OPS *dmn;
 dmn=((DMN_P_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 res=dmn[0].evaluation_p(AUX,AUX);//We work with position [0] as
 if(res<pow(TEN,-TWO*SIX)){res=ZERO;}
 if(res>=pow(TEN,-TWO*SIX))
 {
  res=-res*log(res);
 }
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Shannonp cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int FI_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)&&
 (phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX,*Grad,eval;
 AUX=new double[ndim];
 Grad=new double[ndim];
 DMN_OPS *dmn;
 dmn=((DMN_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 dmn[0].grad_rho_r(AUX,Grad,eval);
 if(eval==ZERO) eval=pow(TEN,-TEN);
 res=pow(norm3D(Grad),TWO);
 delete[] AUX;
 delete[] Grad;
 Grad=NULL;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Fisher cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int FIp_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)&&
 (phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX,*Grad,eval;
 AUX=new double[ndim];
 Grad=new double[ndim];
 DMN_P_OPS *dmn;
 dmn=((DMN_P_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 dmn[0].grad_rho_p(AUX,Grad,eval);
 if(eval==ZERO) eval=pow(TEN,-TEN);
 res=pow(norm3D(Grad),TWO);
 delete[] Grad;
 delete[] AUX;
 Grad=NULL;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Fisher cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
int TTF_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)&&
 (phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
 double *AUX,eval;
 AUX=new double[ndim];
 DMN_OPS *dmn;
 dmn=((DMN_OPS*)fdata);
 AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature); //((READ_FCHK_WFN*)userdata)[0].rho_eval(AUX,res);
 AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
 AUX[2]=r_cubature*cos(theta_cubature);
 eval=dmn[0].evaluation(AUX,AUX);//We work with position [0] as when Rho=new READ_FCHK_WFN[1];
 if(eval>pow(TEN,-THREE*FIVE))
 {
  res=pow(eval,FIVE/THREE);
 }
 else
 {
  res=ZERO;
 }
 delete[] AUX;
 AUX=NULL;
 dmn=NULL;
 }
 else{res=ZERO;}
 //Thomas Fermi cubature:
 fval[0] = FOUR*PI*res*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 return 0;
}
///////////////////////
//Cubature 2         //
///////////////////////
void I_cubature2_dmnr(DMN_OPS &DMNr,double presc,double *result_integration,
double error[9],double Integrals_interval[9],bool shan, bool fish, bool r1, bool r2,bool dipolar)
{
 int i;
 double val[9],quasi_zero,quasi_one;
 void * USERDATA=NULL;
 USERDATA=&DMNr;
 quasi_one=1.0e00-1.0e-10;
 quasi_zero=1.0e-10;
 r1_dmn=r1;r2_dmm=r2;
 shan_dmn=shan;
 fishr_dmn=fish;
 DIPOLE_dmn=dipolar;
 bool interval=false;
 if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP_DMN))
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
  define_interval_dmn(Integrals_interval);
 }
 else
 {
  theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
  r_sup_dmn=RSUP_DMN;
 }
 double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
 double xmax[3] = {quasi_one,quasi_one,quasi_one};
 hcubature(9, Tot_Cubature_dmnr,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
 val, error);
 for(i=0;i<9;i++)
 {
  result_integration[i]=val[i];
 }
}
void I_cubature2_dmnp(DMN_P_OPS &DMNp,double presc,double *result_integration,
double error[9],double Integrals_interval[9],bool shan, bool fish, bool p1, bool p2)
{
 int i;
 double val[9],quasi_zero,quasi_one;
 void * USERDATA=NULL;
 USERDATA=&DMNp;
 quasi_one=1.0e00-1.0e-10;
 quasi_zero=1.0e-10;
 p1_dmn=p1;p2_dmn=p2;
 shanp_dmn=shan;
 fishp_dmn=fish;
 bool interval=false;
 if((Integrals_interval[0]==ZERO && Integrals_interval[1]==RSUP_DMN))
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
  define_interval_dmn(Integrals_interval);
 }
 else
 {
  theta_inf_dmn=ZERO;theta_sup_dmn=PI;phi_inf_dmn=ZERO;phi_sup_dmn=TWO*PI;r_inf_dmn=ZERO;
  r_sup_dmn=RSUP_DMN;
 }
 double xmin[3] = {quasi_zero,quasi_zero,quasi_zero};
 double xmax[3] = {quasi_one,quasi_one,quasi_one};
 hcubature(9, Tot_Cubature_dmnp,USERDATA,  3, xmin, xmax, 0, 1.0e-4,presc, ERROR_INDIVIDUAL,
 val, error);
 for(i=0;i<9;i++)
 {
  result_integration[i]=val[i];
 }
}
//Cubature functions:
int Tot_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res,res_F,res_S,jacobian;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
  double *AUX,*Grad;
  AUX=new double[ndim];
  Grad=new double[ndim];
  DMN_OPS *dmn; //Pointer to an object of DMN_OPS type.
  dmn=((DMN_OPS*)fdata);
  AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
  AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
  AUX[2]=r_cubature*cos(theta_cubature);
  dmn[0].grad_rho_r(AUX,Grad,res);//We work with position [0]
  dmn=NULL;
  if(res==ZERO)
  {res_F=pow(TEN,-TEN);}
  else
  {res_F=res;}
  res_F=pow(norm3D(Grad),TWO)/res_F;
  delete[] AUX;delete[] Grad;
  AUX=NULL;Grad=NULL;
 }
 else{res=ZERO;res_F=ZERO;}
 jacobian=FOUR*PI*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 //Density cubature:
  fval[0] = jacobian*res;
 // Shannon:
  if(shan_dmn)
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
   fval[1] = jacobian*res_S;
  }
  else
  {fval[1]=ZERO;}
 // Fisher:
  if(fishr_dmn)
  {
   fval[2] = jacobian*res_F;
  }
  else
  {fval[2]=ZERO;}
 // T_TF:
  fval[3] = jacobian*pow(res,FIVE/THREE);
 // <r>:
  if(r1_dmn)
  {
   fval[4] = jacobian*r_cubature*res;
  }
  else
  {fval[4]=ZERO;}
 // <r2>:
  if(r2_dmm)
  {
   fval[5] = jacobian*r_cubature*r_cubature*res;
  }
  else
  {fval[5]=ZERO;}
 // mu:
  if(DIPOLE_dmn)
  {
   fval[6] = jacobian*r_cubature*sin(theta_cubature)*cos(phi_cubature)*res;
   fval[7] = jacobian*r_cubature*sin(theta_cubature)*sin(phi_cubature)*res;
   fval[8] = jacobian*r_cubature*cos(theta_cubature)*res;
  }
  else
  {fval[6]=ZERO;fval[7]=ZERO;fval[8]=ZERO;}
 return 0;
}
int Tot_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
 #ifndef r_cubature
 #define r_cubature x[0]/(1-x[0])
 #define theta_cubature acos(1-2*x[1])
 #define phi_cubature x[2]*2*PI
 #endif
 double res,res_F,res_S,jacobian;
 if((theta_cubature>=theta_inf_dmn && theta_cubature<=theta_sup_dmn)
 &&(phi_cubature>=phi_inf_dmn&&phi_cubature<=phi_sup_dmn)&&(r_cubature>=r_inf_dmn&&r_cubature<=r_sup_dmn))
 {
  double *AUX,*Grad;
  AUX=new double[ndim];
  Grad=new double[ndim];
  DMN_P_OPS *dmn; //Pointer to an object of DMN_OPS type.
  dmn=((DMN_P_OPS*)fdata);
  AUX[0]=r_cubature*sin(theta_cubature)*cos(phi_cubature);
  AUX[1]=r_cubature*sin(theta_cubature)*sin(phi_cubature);
  AUX[2]=r_cubature*cos(theta_cubature);
  dmn[0].grad_rho_p(AUX,Grad,res);//We work with position [0]
  dmn=NULL;
  if(res==ZERO)
  {res_F=pow(TEN,-TEN);}
  else
  {res_F=res;}
  res_F=pow(norm3D(Grad),TWO)/res_F;
  delete[] AUX;delete[] Grad;
  AUX=NULL;Grad=NULL;
 }
 else{res=ZERO;res_F=ZERO;}
 jacobian=FOUR*PI*x[0]*x[0]*pow(ONE/(ONE-x[0]),FOUR);
 //Densityp cubature:
  fval[0] = jacobian*res;
 // Shannonp:
  if(shanp_dmn)
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
   fval[1] = jacobian*res_S;
  }
  else
  {fval[1]=ZERO;}
 // Fisherp:
  if(fishp_dmn)
  {
   fval[2] = jacobian*res_F;
  }
  else
  {fval[2]=ZERO;}
 // T_TFp (no sense):
  fval[3] = ZERO;
 // <p>:
  if(p1_dmn)
  {
   fval[4] = jacobian*r_cubature*res;
  }
  else
  {fval[4]=ZERO;}
 // <p2>:
  if(p2_dmn)
  {
   fval[5] = jacobian*r_cubature*r_cubature*res;
  }
  else
  {fval[5]=ZERO;}
 // mup (no sense):
 fval[6]=ZERO;fval[7]=ZERO;fval[8]=ZERO;
 return 0;
}
