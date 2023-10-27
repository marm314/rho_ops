#ifndef _UTILS_IO_H_
#define _UTILS_IO_H_

#include<iostream>
#include"Mathematical_Functions.h"
#include"D_read_calc_rho.h"
#include"Input_commands.h"
#include"Corr_indicators.h"
#include"MO_class.h"
#include"NO_class.h"
#include"MOp_class.h"
#include"NO_DMN_class.h"
#include"NOp_DMN_class.h"
#include"Numbers.h"
#include"String_ops.h"
#include"gnuplot.h"

using namespace std;

void pre_elf(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &Density_alpha,double &Density_beta,
double &tauW_alpha,double &tauW_beta,double &tau_alpha,double &tau_beta,double &tcurr_alpha,double &tcurr_beta);
void print_int(READ_FCHK_WFN &Read_fchk_wfn,string name_file,double **Sij, int nbasis,double &rho,double &rhoa,double &rhob,
string region);
void mos_to_nos_dmn_sij(double **SIJ,READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands,string name_file,bool &wfn_fchk);
void cube_file(READ_FCHK_WFN &Read_fchk_wfn,string name_file,string op_cube,double cubex,double cubey, double cubez,double stepx,
double stepy,double stepz);
void mos_to_nos_int_fchk_dm1(READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands,string name_file,double **ORBITALS,int &total_grid,
int &nbasis,bool wfn_fchk);

#endif // _UTILS_IO_H_
