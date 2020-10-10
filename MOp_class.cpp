#include"MOp_class.h"
MOp::MOp(){evaluation=(ZERO,ZERO);grad[0]=(ZERO,ZERO);grad[1]=(ZERO,ZERO);grad[2]=(ZERO,ZERO);spin=true;}
MOp::MOp(READ_FCHK_WFN &Read_fchk_wfn,double Pointp[3],int numMO)
{
 if(!Read_fchk_wfn.wfn)
 {
  Read_fchk_wfn.build_MOp_fchk(evaluation,Pointp,numMO);
  Read_fchk_wfn.grad_MOp_fchk(grad,Pointp,numMO);
  spin=Read_fchk_wfn.SPIN[numMO];
 }
 else
 {
  //In the case of wfn we have MOs for single determinat
  //but for multireference calculations we only have NO then
  //we can only use this function:
  Read_fchk_wfn.build_NOp_wfn(evaluation,Pointp,numMO);
  spin=Read_fchk_wfn.SPIN[numMO];
 }
}
MOp::MOp(READ_FCHK_WFN &Read_fchk_wfn,double **AOps,complex<double> **AOps_Grad,int numMO)
{
 if(!Read_fchk_wfn.wfn)
 {
  Read_fchk_wfn.build_MOp_fchk2(AOps,evaluation,numMO);
  Read_fchk_wfn.grad_MOp_fchk2(AOps_Grad,grad,numMO);
  spin=Read_fchk_wfn.SPIN[numMO];
 }
}
MOp::~MOp()
{}
MOp::MOp(const MOp& MOpin)
{
 grad[0]=MOpin.grad[0];
 grad[1]=MOpin.grad[1];
 grad[2]=MOpin.grad[2];
 evaluation=MOpin.evaluation;
 spin=MOpin.spin;
}
