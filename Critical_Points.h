#ifndef _CRITICAL_POINTS_H_
#define _CRITICAL_POINTS_H_
#include "D_read_calc_rho.h"
#include "Numbers.h"
#include "Mathematical_Functions.h"
using namespace std;
 class CP
 {
  private:
   double direction[3],atom1[3],atom2[3];
   void finesearch(READ_FCHK_WFN &Rho_r,double Point[3],double step);
  public:
   double grad_norm,BCP[3];
   int atomic_num1,atomic_num2;
   CP();
   CP(READ_FCHK_WFN &Rho_r,int AtomA, int AtomB);
   ~CP();
 };
#endif // _CRITICAL_POINTS_H_
