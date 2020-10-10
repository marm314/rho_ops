#ifndef _NO_CLASS_H_
#define _NO_CLASS_H_
#include"D_read_calc_rho.h"
using namespace std;
 class NO
{
 private:
 public:
  //The default constructor is meaningless
  NO();
  NO(READ_FCHK_WFN &,double [3],int);
  NO(READ_FCHK_WFN &,double *,double **,int);
  NO(const NO& NOin);
  ~NO();
  double evaluation,grad[3];
  bool spin;
};
#endif // _NO_CLASS_H_

