#ifndef _INTEGRALS_QUADRATURE_H_
#define _INTEGRALS_QUADRATURE_H_
#include<iostream>
#include<string.h>
#include<fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include"sphere_lebedev_rule.h"
#include"legendre_quadrature.h"
#include"Numbers.h"
#include"D_read_calc_rho.h"
#include"DMN_ops_class.h"
#include"DMN_ops_p_class.h"
#include"MO_class.h"
#include"NO_class.h"
using namespace std;
void integrate_quadrature(void *,string,bool,int ,int ,bool ,double **,double **,int mode=0);
void integrate_quadrature_p(void *,string,bool,int ,int ,bool ,double **,double **);
void build_nos_quadrature(void *,string,double **,int,int);
void build_nos_quadrature_rot(void *data,string name,double **ORBITALS,int order,int order2,double rot);
void build_mos_quadrature(void *,string,double **,int,int);
void build_mos_quadrature_rot(void *data,string name,double **ORBITALS,int order,int order2,double rot);
void integral_calc(string,int &,int &,double **,double **,double [6],double &);
void integral_calc_p(string,int &,int &,double **,double **,double [6],double &);
void calc_sij_mat(double **,double**,double [6],int &,int &,int &);
void calc_mij_mat(double **,double**,double [6],int &,int &,int &);
void integrate_intra_coord(double **,double,double,double,double [3],double[3],int,double*,double*);
void grid_avail(int &);
void clean_quadrature(string,int mode=0);
double theta_rad(double &);
double phi_rad(double &);
void rotate_point(double Point[3],double rot);
#endif //_INTEGRALS_QUADRATURE_H_
