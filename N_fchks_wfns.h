#ifndef _N_FCHKS_WFNS_STRUCT_H_
#define _N_FCHKS_WFNS_STRUCT_H_
#include"D_read_calc_rho.h"
using namespace std;
 struct N_FCHKS_WFNS
{
 //Pointer to two FCHKs or WFNs for more change the 5 to N_files...
 READ_FCHK_WFN *read_fchk_wfn[5];
};
#endif // _N_FCHKS_WFNS_STRUCT_H_
