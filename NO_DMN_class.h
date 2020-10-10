#ifndef _NO_DMN_CLASS_H_
#define _NO_DMN_CLASS_H_
#include"DMN_ops_class.h"
using namespace std;
 class NO_DMN
{
 private:
 public:
  //The default constructor is meaningless
  NO_DMN();
  NO_DMN(DMN_OPS &,bool ,double [3],int);
  ~NO_DMN();
  double evaluation,occ;
  bool spin;
};
#endif // _NO_DMN_CLASS_H_

