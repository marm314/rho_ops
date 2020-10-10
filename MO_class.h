#ifndef _MO_CLASS_H_
#define _MO_CLASS_H_
#include"D_read_calc_rho.h"
using namespace std;
 class MO
{
 private:
 public:
  //The default constructor is meaningless
  MO();
  MO(READ_FCHK_WFN &,double [3],int);
  MO(READ_FCHK_WFN &,double *,double **,int);
  MO(const MO& MOin);
  ~MO();
  double evaluation,grad[3];
  bool spin;
};
#endif // _MO_CLASS_H_


