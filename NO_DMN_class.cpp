#include"NO_DMN_class.h"
NO_DMN::NO_DMN(){evaluation=ZERO;spin=true;}
NO_DMN::NO_DMN(DMN_OPS &dmn,bool alpha_beta_sep,double Point[3],int numNO)
{
 //If the diagonalization was not done before, diagonalize the DM1
 if(alpha_beta_sep && !dmn.diag_ab_done)
 {dmn.diagonalize_ab();}
 if(!alpha_beta_sep && !dmn.diag_done)
 {dmn.diagonalize();}
 //If the diagonalization succeed  (the fchk and dm1 files were correct) continue
 if(dmn.diag_done || dmn.diag_ab_done)
 {
  if(alpha_beta_sep)
  {
   if(numNO%2==0)
   {
    spin=true;
    occ=dmn.rho_matrixa[numNO/2][numNO/2];
   }
   else
   {
    spin=false;
    occ=dmn.rho_matrixb[(numNO-1)/2][(numNO-1)/2];
   }
  }
  else
  {
   spin=true;
   occ=dmn.rho_matrix[numNO][numNO];
  }
  // The total number of NOs for alpha_beta_sep=true is the double of the spins summed case.
  // Max numNO for spins summed is Read_fchk_wfn.nbasisf-1 and for spins separated is
  // 2*Read_fchk_wfn.nbasisf-1 [since we count from 0 in C++]
  dmn.build_NO_DMN(evaluation,Point,alpha_beta_sep,numNO);
 }
 else
 {
  cout<<"Could not diagonalize the DM1 since the fchk and/or dm1 file(s) is/are missing"<<endl;
  cout<<"Check the name of the dm1 file and fchk file."<<endl;
  cout<<"The evaluation and occupation of the NO have been set to ZERO"<<endl;
  //Also check if the set_fchk() function was used in the code.
  evaluation=ZERO;
  occ=ZERO;
 }
}
NO_DMN::~NO_DMN()
{}
