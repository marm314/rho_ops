#include"NOp_DMN_class.h"
NOp_DMN::NOp_DMN(){evaluation=ZERO+ZERO*I;spin=true;}
NOp_DMN::NOp_DMN(DMN_P_OPS &dmnp,bool alpha_beta_sep,double Point[3],int numNO)
{
 //If the diagonalization was not done before, diagonalize the DM1
 if(alpha_beta_sep && !dmnp.diag_ab_done)
 {dmnp.diagonalize_ab();}
 if(!alpha_beta_sep && !dmnp.diag_done)
 {dmnp.diagonalize();}
 //If the diagonalization succeed (the fchk and dm1 files were correct) continue
 if(dmnp.diag_done || dmnp.diag_ab_done)
 {
  if(alpha_beta_sep)
  {
   if(numNO%2==0)
   {
    spin=true;
    occ=dmnp.rho_matrixa[numNO/2][numNO/2];
   }
   else
   {
    spin=false;
    occ=dmnp.rho_matrixb[(numNO-1)/2][(numNO-1)/2];
   }
  }
  else
  {
   spin=true;
   occ=dmnp.rho_matrix[numNO][numNO];
  }
  // The total number of NOs for alpha_beta_sep=true is the double of the spins summed case.
  // Max numNO for spins summed is Read_fchk_wfn.nbasisf-1 and for spins separated is
  // 2*Read_fchk_wfn.nbasisf-1 [since we count from 0 in C++]
  dmnp.build_NOp_DMN(evaluation,Point,alpha_beta_sep,numNO);
 }
 else
 {
  cout<<"Could not diagonalize the DM1 since the fchk and/or dm1 file(s) is/are missing"<<endl;
  cout<<"Check the name of the dm1 file and fchk file."<<endl;
  cout<<"The evaluation and occupation of the NOp have been set to ZERO"<<endl;
  //Also check if the set_fchk() function was used in the code.
  evaluation=ZERO+ZERO*I;
  occ=ZERO;
 }
}
NOp_DMN::~NOp_DMN()
{}
