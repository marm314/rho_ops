#ifndef _INPUT_COMMANDS_H_
#define _INPUT_COMMANDS_H_

#include<iostream>
#include<algorithm>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include"String_ops.h"
#include"Numbers.h"

using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
class Input
{
 public:
 int i,punctuals_r,punctuals_p,multiplicity,minevals,maxevals,ops,*dmn_order,dmns,nregions,p_rad;
 int *dmn_orderp,dmnsp,ncores,order_grid_r,order_grid_ang,num_plot_ops,extra_lines,orb1,orb2,mo_no_cube;
 bool punctualr,punctualp,scanr,scanp,scanelf,scanindic,cps,indicators,integrals,log,spin_calcs,cas,rotate_grid;
 bool cuba,cubature,quadrature,debug,int_file,dmn,dmnp,dmn_integrals,dmn_thresh,esi_int,symrotdens,symgrad;
 bool dmn_indicators,gnuplot,dmn_plots,nopath,dim3,dim2,Beta_MOs,wfx_print,wfx_print_dmn,store_dmn,print_dm1_fchk;
 bool cubature2,cube,tps,r1_moment,intracule,Vr,scan_localhybs,dens_sim,extracule,int_pol_hyperpol,intra_1rdm,mulliken;
 double **interval_integrals,**coordinates_r,**coordinates_p;
 double init_coord_r[3],init_coord_p[3],init_coord_elf[3],init_coord_lh[3],init_coord_indic[3],interval_integralsCUB[6];
 double step_r,step_p,step_elf,step_lh,step_indic,points_scan_r,points_scan_p,points_scan_elf,points_scan_lh,points_scan_indic,grid_1,grid_2,error_abs,error_rel,dmn_threshold;
 double pointDM1p[3],pointprimeDM1p[3],points_scan_1[2],second_moments_tps[6];
 double pointDM1[3],pointprimeDM1[3],points_scan_2[2],not_used_3d_doub,rotate_angle,mem;
 double cubex,cubey,cubez,step_cube_x,step_cube_y,step_cube_z,Point_intra[3],Point_extra[3],Point_Vr[3],field_F;
 char dir_r,dir_p,dir_elf,dir_lh,dir_indic,scan_1,scan_2,not_used_3d,dir_pol_hyper;
 string name_fchk_wfn,name_log,*name_dmn,*name_dmnp,cuba_cubature,method_cuba,*integral_ops,Sij_region,MOorNO;
 string *dmns_read,name_dm1,name_dm2,path,*extra_lines_plot,*plot_ops,opcube,second_fchk_wfn,third_fchk_wfn;
 string fourth_fchk_wfn,fifth_fchk_wfn;
 Input();
 Input(string);
 ~Input();
};
#endif // _INPUT_COMMANDS_H_



