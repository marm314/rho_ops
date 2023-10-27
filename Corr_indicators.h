#ifndef _CORR_INDICATORS_H_
#define _CORR_INDICATORS_H_

#include<iostream>
#include"Mathematical_Functions.h"
#include"D_read_calc_rho.h"
#include"Input_commands.h"
#include"MO_class.h"
#include"NO_class.h"
#include"MOp_class.h"
#include"NO_DMN_class.h"
#include"NOp_DMN_class.h"
#include"Numbers.h"
#include"String_ops.h"

using namespace std;

double Deviation_idemp(READ_FCHK_WFN &Read_fchk_wfn);
double ID_ni(READ_FCHK_WFN &Read_fchk_wfn);
double IND_ni(READ_FCHK_WFN &Read_fchk_wfn);
double Shannon_ni(READ_FCHK_WFN &Read_fchk_wfn);
double Deviation_idemp_dmn(DMN_OPS &DMN,double &N);
double ID_ni_dmn(DMN_OPS &DMN,double &N);
double IND_ni_dmn(DMN_OPS &DMN,double &N);
void ID_IND_local(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &ID_alpha,double &ID_beta,
double &IND_alpha,double &IND_beta,bool &firstcall);
double Shannon_ni_dmn(DMN_OPS &DMN,double &N);

#endif // _CORR_INDICATORS_H_
