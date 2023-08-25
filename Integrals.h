#ifndef _INTEGRALS_H_
#define _INTEGRALS_H_

#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include<omp.h>
#include"Mathematical_Functions.h"
#include"D_read_calc_rho.h"
#include"NO_class.h"
#include"MO_class.h"
#include"N_fchks_wfns.h"
#include"cuba.h"
#include"cubature.h"
#include"Numbers.h"
#define RSUP 1.0e99

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
void define(); //for CUBA
void define_interval(double Interval[6]);
//Cuba
void integrate_cuba(READ_FCHK_WFN &Rho,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[20], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool r1,bool r2,bool rm1,bool dipole,bool rho);
void integrate_cubap(READ_FCHK_WFN &Rho,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[20], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool p1, bool p2);
void integrate_cuba_sij(READ_FCHK_WFN &Rho,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double **SIJ, int &fail,double Integrals_interval[6],int &nprocs,string MOorNO);
void  integrate_dens_sim(N_FCHKS_WFNS two_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double res_integration[3], int &fail,int &nprocs,double ROT_MATRIX[3][3]);
void  integrate_V_Hartree(N_FCHKS_WFNS two_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double res_integration[7], int &fail,int &nprocs,double displacement[3]);
void  integrate_dens_sim2(N_FCHKS_WFNS two_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double result_integration[7], int &fail,int &nprocs,double ROT_MATRIX[3][3]);
void  integrate_pol_hyperpol(N_FCHKS_WFNS five_fchks_wfns,string method, const int NDIM,
const int NCOMP, const double EPSREL,const double EPSABS, const int MINEVAL, const int MAXEVAL,
double res_integration[9], int &fail,int &nprocs,char dir_pol_hyper);
void integrate_tps_fchk(READ_FCHK_WFN &Rho,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[10], int &fail,int &nprocs);
void integrate_vr_fchk(READ_FCHK_WFN &Rho,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[2], int &fail,int &nprocs,double Point_Vr_in[3]);
//Function redable from cuba.h
int Integrand(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrandp(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int SIJ_Integrand(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrand_tps_fchk(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrand_vr_fchk(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrand_divergences(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrand_divergences2(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrand_pol_hyper(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrand_v_hartree(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
//Cubature
void I_cubature(READ_FCHK_WFN &Rho,string Operation,double presc,double &result_integration,
double &error,double Integrals_interval[6]);
//Function redable from cubature.h
int N_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int Np_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int S_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int Sp_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int FI_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int FIp_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int TTF_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int r1_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int r2_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int p1_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int p2_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int rho_Cubature(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int Sij_NO_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int Sij_MO_Cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
//Cubature 2 (vectorial evaluations of density, shannon, etc.)
void I_cubature2(READ_FCHK_WFN &Rho,double presc,double result_integration[10],
double error[10],double Integrals_interval[9],bool shan,bool fish, bool r1, bool r2,bool dipolar,
bool rho);
void I_cubature2p(READ_FCHK_WFN &Rho,double presc,double result_integration[10],
double error[10],double Integrals_interval[9],bool shan,bool fish, bool p1, bool p2);
int mux_Cubature2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int muy_Cubature2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int muz_Cubature2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
#endif // _INTEGRALS_H_
