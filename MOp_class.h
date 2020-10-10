#ifndef _MOp_CLASS_H_
#define _MOp_CLASS_H_
#include"D_read_calc_rho.h"
using namespace std;
 class MOp
{
 private:
 public:
  //The default constructor is meaningless
  MOp();
  MOp(READ_FCHK_WFN &,double [3],int);
  MOp(READ_FCHK_WFN &,double **,complex<double> **,int);
  MOp(const MOp& MOpin);
  ~MOp();
  complex<double> evaluation,grad[3];
  bool spin;
};
#endif // _MOp_CLASS_H_

