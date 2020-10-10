#ifndef _NOp_DMN_CLASS_H_
#define _NOp_DMN_CLASS_H_
#include"DMN_ops_p_class.h"
using namespace std;
 class NOp_DMN
{
 private:
 public:
  //The default constructor is meaningless
  NOp_DMN();
  NOp_DMN(DMN_P_OPS &,bool ,double [3],int);
  ~NOp_DMN();
  complex<double> evaluation;
  double occ;
  bool spin;
};
#endif // _NOp_DMN_CLASS_H_
