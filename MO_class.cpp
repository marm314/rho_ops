#include"MO_class.h"
MO::MO(){evaluation=ZERO;spin=true;grad[0]=ZERO;grad[1]=ZERO;grad[2]=ZERO;}
MO::MO(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],int numMO)
{
 if(!Read_fchk_wfn.wfn)
 {
  Read_fchk_wfn.build_MO_grad_fchk(evaluation,grad,Point,numMO);
  spin=Read_fchk_wfn.SPIN[numMO];
 }
 else
 {
  //In the case of wfn we have MOs for single determinat
  //but for multireference calculations we only have NO then
  //we can only use this function:
  Read_fchk_wfn.build_NO_grad_wfn(evaluation,grad,Point,numMO);
  spin=Read_fchk_wfn.SPIN[numMO];
 }
}
MO::MO(READ_FCHK_WFN &Read_fchk_wfn,double *AOs,double **AOs_Grad,int numMO)
{
 if(!Read_fchk_wfn.wfn)
 {
  Read_fchk_wfn.build_MO_grad_fchk2(AOs,AOs_Grad,evaluation,grad,numMO);
  spin=Read_fchk_wfn.SPIN[numMO];
 }
}
MO::~MO()
{}
MO::MO(const MO& MOin)
{
 grad[0]=MOin.grad[0];
 grad[1]=MOin.grad[1];
 grad[2]=MOin.grad[2];
 evaluation=MOin.evaluation;
 spin=MOin.spin;
}
