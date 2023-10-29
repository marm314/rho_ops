#ifndef _INTEGRALS_BECKE_H_
#define _INTEGRALS_BECKE_H_

#include<iostream>
#include<string.h>
#include<fstream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<omp.h>
#include"sphere_lebedev_rule.h"
#include"legendre_quadrature.h"
#include"Numbers.h"
#include"D_read_calc_rho.h"
#ifdef HAVE_LIBXC
#include"xc.h"
#endif

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
void Grid_atomic(READ_FCHK_WFN &Rho,string name,int &natom,int &nradial,int &nang,int &stiff,string partition);
void Integrate_atomic(READ_FCHK_WFN &Rho,double *res_integration);
void Integrate_atomic_paral(vector<READ_FCHK_WFN>Rho,double *res_integration,int &nprocs);
//void Integrate_atomic_sij(READ_FCHK_WFN &Rho,double ***S_Aij,double *Density,int &nbasis);
void Integrate_atomic_sij(vector<READ_FCHK_WFN> Rho,double ***S_Aij,double *Density,int &nbasis,int &nprocs);
void Integrate_atomic_mescal(READ_FCHK_WFN &Rho,double **F_QM,double *V_QM,double &Density,double **Coords_pdb,int &natoms_pdb);
void Clean_quadrature_atomic(string name, int &natoms);
void Grid_avail_atomic(int & Order);
double Xi_XY_table(int &Z1, int &Z2);
void Xi_XY_bcp(READ_FCHK_WFN &Rho, string name, int &iatom, int &jatom, int &natom, double &Xi_rad);
double set_radii(int &Z);
double p_mu(double &mu);
double smooth_stiff(double &mu,int stiff);
#endif // _INTEGRALS_BECKE_H_
