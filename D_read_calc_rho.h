#ifndef _RHO_READ_CALC_H_
#define _RHO_READ_CALC_H_
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<locale>
#include<sstream>
#include<algorithm>
#include"Numbers.h"
#include"Mathematical_Functions.h"
#include"String_ops.h"

using namespace std;
///////////////
//Begin class//
///////////////
class READ_FCHK_WFN
{
 private:
   int counter,counter2,counter3,nu_charge,stype,nprim_shell,contr_coef;
   int spcontr_coef;
   int cartes_coord,nshells,largcontr,npure_d,npure_f,Hi_ang_moment;
   int nalphael,nbetael,nindepbasisf;
   int MO_coef,SCF_rho,MO_beta_coef,spin_SCF_rho,rho_CI,spin_CI_rho;
   double density;
   bool gamess,CI,CAS,relaxed,virtuals;
   bool nu_ch,cart_coor,shtype,nprimsh,satmap,exponents,contract,contractSP,activeSP;
   bool mocoef,scfrho,mo_beta_coef,spinscfrho,rho_CI_bool,spinCIrho,extra1,PS_bool,BETA_MOS,scf_dens_found;
   double *AUX;
   ifstream input_fchk,input_wfn;
   ofstream pseudowfn,pseudowfx;
   string line,name_file;
   /////////////////////////////////////////
   //Declare Functions used by this class //
   /////////////////////////////////////////
   ////////////
   //for fchk//
   ////////////
   //Position space:
   void line_fill(int &number, string line_read);
   void cast_doub(double *filling, string line_read, int &counter, int control);
   void cast_int(int *filling, string line_read, int &counter, int control);
   void build_AO_AOgrad(double *AO,double** AO_grad,double Point[3]);
   void build(int &iprim,int &styp,int &counter,double *C_coef_send,
   double *SPC_coef_send,double *P_exp_send,double Point[3],
   double *AO,double Coord_At[3],bool &activeSP);
   void build_grad(int &iprim,int &styp,int &counter2,double *C_coef_send,
   double *SPC_coef_send,double *P_exp_send,double Point[3],double **AO_grad,
   double Coord_At[3],bool &activeSP);
   double eval(double Coord_At[3],double Point[3],double &expon, int &n, int &l,int &m);
   double eval_grad(int &dir,double Coord_At[3],double Point[3],double &expon, int &n,
   int &l,int &m);
   //Momentum space:
   void build_AOp(double **AOp,double Point[3]);
   void build_AOp_grad(complex<double> **AOps_Grad,double point_p[3]);
   void build_prim_AOp(int iprim,int styp,int &counter,double *C_coef_send,
   double *SPC_coef_send,double *P_exp_send,double Point[3],double **AOp,double Coord_At[3],
   bool &activeSP);
   void eval_p(double Coord_At[3],double point_p[3],double &expon, int &n, int &l,int &m,double &re,double &im);
   void fps(complex<double> &f,double &p,int &quant_num,double &expon);
   ///////////
   //for wfn//
   ///////////
   //Position space:
   void cast_wfn_coord(string in, double **Cartesian_Coor,int &counter);
   void cast_wfn_nuch(string in, double *Nu_charge,int &counter2);
   void cast_wfn_type_asig(string in,int *filling,int &counter);
   void cast_wfn_expon(string in,double *filling,int &counter);
   void cast_wfn_MO(string in,double **filling,int &counter,int &counter2);
   double eval_Primitive_wfn(double Point[3],double pos_nuclei[3],int nlm[3],double exponent);//MO in wfn is NO
   double grad_Primitive_wfn(double Point[3],int dir,double pos_nuclei[3],int nlm[3],
   double exponent);
   //Momentum space:
   void eval_Primitive_p_wfn(double point_p[3],double pos_nuclei[3],int nlm[3],double &expon,double &re,double &im);
   void Center_of_mass();
 public:
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   //Functions with a number 2 are the same as the ones without 2 but either read AOs (AOps) from outside or//
   //create them outside the class                                                                          //
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   READ_FCHK_WFN();
   READ_FCHK_WFN(string,string,bool,bool,bool,int);
   READ_FCHK_WFN(const READ_FCHK_WFN& RHO);
   ~READ_FCHK_WFN();
   void set_BETA_MOS(bool &);
   //Variables
   int natoms,nprimitv,smap,prim_exp,nbasisf,multiplicity,nelectrons,Pair[2];
   int *shell_type, *n_prim_per_shell, *shell_map;
   double *Nu_charge,*Ocupation,**Cartesian_Coor,**Total_rho,**MOcoefA,**MOcoefB;
   double *Prim_exp, *Contr_Coef, *SP_Contr_Coef,Rot_ICM[3][3];
   double  **Spin_rho,**P,**Pbeta,**S,**Sbeta,**Sao;
   bool *SPIN,wfn,open_shell,extra0,rhf,uhf,correlated,error_opens_wfn,overlap,no_beta_wfn,wfx;
   string identity;
   ////////////
   //for both//
   ////////////
   void Quant_fill(int **Quant,int styp);
   //Position space
   void build_AO_AOgrad2(double *AO,double **AO_grad,double Point[3]);
   void rho_eval(double [3],double &);
   void rho_eval_a_b(double [3],double &,double &);
   void rho_grad(double [3],double [3]);
   void rho_grad_a_b(double [3],double [3],double [3]);
   void rho_lapl(double [3],double &);
   void rho_lapl_a_b(double [3],double &,double &);
   void orb_grad(double [3],double **res);
   void build_NO_grad_wfn(double &,double [3],double [3],int &);
   void build_NO_grad_wfn_all(double *,double **,double [3]);
   void build_NO_grad_fchk(double &,double [3],double [3],int &);
   void build_NO_grad_fchk2(double *,double **,double &,double [3],int &);
   void build_MO_grad_fchk(double &,double [3],double [3],int &);
   void build_MO_grad_fchk2(double *, double **,double &,double [3],int &);
   int nbasis(); //Contains total basis (if possible spins separated)
   void muATOMS(double mu[3]);//Get the dipolar moment of the atoms
   void quadrupoleATOMS(double **quadrupole);//Get the dipolar moment of the atoms
   void bring_mo_coefs(double **coefs);//Send MO coefs info
   double Vnuclear(double Point_Vr[3]);//Compute Vn(r)= sum _i ^Natoms Zi/|r-Ri|
   //Momentum space
   void build_AOp2(double **AOp,double point_p[3]);
   void build_AOp_grad2(complex<double> **AOps_Grad,double point_p[3]);
   void rho_p_eval(double [3],double &);
   void rho_p_eval_a_b(double [3],double &,double &);
   void rho_p_grad(double [3], double [3]);
   void rho_p_grad_a_b(double [3], double [3], double [3]);
   void rho_p_lapl(double [3], double &);
   void rho_p_lapl_a_b(double [3], double &,double &);
   void build_NOp_wfn(complex<double> &NOp,double Point[3],int &numMO);
   void build_MOp_fchk(complex<double> &MOp,double Point[3],int &numMO);
   void build_MOp_fchk2(double  **AOp,complex<double> &MOp,int &numMO);
   void grad_MOp_fchk(complex<double> Grad[3],double Point[3],int &numMO);
   void grad_MOp_fchk2(complex<double> **AOps_Grad,complex<double> Grad[3],int &numMO);
};
#endif // _RHO_READ_CALC_H_
