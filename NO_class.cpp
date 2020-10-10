#include"NO_class.h"
NO::NO(){evaluation=ZERO;spin=true;grad[0]=ZERO;grad[1]=ZERO;grad[2]=ZERO;}
NO::NO(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],int numNO)
{
 if(!Read_fchk_wfn.wfn)
 {
  Read_fchk_wfn.build_NO_grad_fchk(evaluation,grad,Point,numNO);
  spin=Read_fchk_wfn.SPIN[numNO];
 }
 else
 {
  Read_fchk_wfn.build_NO_grad_wfn(evaluation,grad,Point,numNO);
  spin=Read_fchk_wfn.SPIN[numNO];
 }
}
NO::NO(READ_FCHK_WFN &Read_fchk_wfn,double *AOs,double **AOs_Grad,int numNO)
{
 if(!Read_fchk_wfn.wfn)
 {
  Read_fchk_wfn.build_NO_grad_fchk2(AOs,AOs_Grad,evaluation,grad,numNO);
  spin=Read_fchk_wfn.SPIN[numNO];
 }
}
NO::~NO()
{}
NO::NO(const NO& NOin)
{
 grad[0]=NOin.grad[0];
 grad[1]=NOin.grad[1];
 grad[2]=NOin.grad[2];
 evaluation=NOin.evaluation;
 spin=NOin.spin;
}
