#ifndef _INTEGRALS_BECKE_H_
#define _INTEGRALS_BECKE_H_

#include<iostream>
#include<string.h>
#include<fstream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include"sphere_lebedev_rule.h"
#include"legendre_quadrature.h"
#include"Numbers.h"
#include"D_read_calc_rho.h"

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
void Grid_becke(READ_FCHK_WFN &Rho,string name,int &natom,int &nradial,int &nang,int &stiff);
void Integrate_becke(READ_FCHK_WFN &Rho,double *res_integration);
void clean_quadrature_becke(string name, int &natoms);
void grid_avail_becke(int & Order);
double Xi_XY_val(int &Z1, int &Z2);
double set_radii(int &Z);
double p_mu(double &mu);
double s_mu_stiff(double &mu,int stiff);
#endif // _INTEGRALS_BECKE_H_