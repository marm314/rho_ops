#ifndef _INTEGRALS_DMN_H_
#define _INTEGRALS_DMN_H_

#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<iomanip>
#include<cmath>
#include"Mathematical_Functions.h"
#include"DMN_ops_class.h"
#include"DMN_ops_p_class.h"
#include"cuba.h"
#include"cubature.h"
#include"Numbers.h"
#define RSUP_DMN 1.0e99

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
void define_dmn(); //for CUBA
void define_interval_dmn(double Interval[6]);
//Cuba
void integrate_dmnr(DMN_OPS &DMNr,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[18], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool r1,bool r2,bool dipole);
void integrate_dmnp(DMN_P_OPS &DMNp,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[18], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool p1,bool p2);
void integrate_intrac(DMN_OPS &DMNr,string method,const int NDIM,const int NCOMP,
const double EPSREL,const double EPSABS,const int MINEVAL, const int MAXEVAL,
double result_integration[18], int &fail,double Integrals_interval[6],int &nprocs,
bool shan,bool fish,bool inertias,bool r1,bool r2,bool dipole);
//Function redable from cuba.h
int Integrand_dmn(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Integrandp_dmn(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
//Cubature
void I_cubature_dmnr(DMN_OPS &DMNr,string Operation,double presc,double &result_integration,
double &error,double Integrals_interval[6]);
void I_cubature_dmnp(DMN_P_OPS &DMNp,string Operation,double presc,double &result_integration,
double &error,double Integrals_interval[6]);
//Function redable from cubature.h
int N_Cubature_dmnr(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int Np_Cubature_dmnp(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int S_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int Sp_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int FI_Cubature_dmnr(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int FIp_Cubature_dmnp(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int TTF_Cubature_dmnr(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int r1_Cubature_dmnr(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int r2_Cubature_dmnr(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int p1_Cubature_dmnp(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
int p2_Cubature_dmnp(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval);
//Cubature 2 (vectorial evaluations of density, shannon, etc.)
void I_cubature2_dmnr(DMN_OPS &DMNr,double presc,double result_integration[9],
double error[9],double Integrals_interval[6],bool shan,bool fish, bool r1, bool r2,bool dipolar);
void I_cubature2_dmnp(DMN_P_OPS &DMNp,double presc,double result_integration[9],
double error[9],double Integrals_interval[6],bool shan,bool fish, bool p1, bool p2);
//functions for cubature 2
int Tot_Cubature_dmnr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int Tot_Cubature_dmnp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
#endif // _INTEGRALS_H_

