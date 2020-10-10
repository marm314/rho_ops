//This class reads .dmn files produced by DMN code of
// Dr. Eduard Matito. & Performs some calculations
// including MOs from fchk files
#ifndef _DMN_OPS_CLASS_H_
#define _DMN_OPS_CLASS_H_
#include<fstream>
#include<omp.h>
#include<vector>
#include"D_read_calc_rho.h"
#include"MO_class.h"
#include"Numbers.h"
using namespace std;
 class DMN_OPS
{
 private:
  int element[4],element_prime[4],RECORD_DELIMITER_LENGTH,nbasisf,nterms;
  double threshold,Dij,Dijkl,Dijklmn,**Array_MOs,**Grad_MOs,**DM1;
  double POINT1[3],POINT_PRIME1[3];
  double Trace,Result,Grad[3];
  bool dmn_file,density,trace,evaluate,grad_density,DMNinCORE;
  bool diagonalize_dm1,diagonalize_dm1_ab,dm2store;
  MO mos_array[4];
  struct indexes_and_elements
  {
   indexes_and_elements(int index[2],int index_prime[2],double dijkl)
   :Dijkl(dijkl)
   {
    indexes[0]=index[0];
    indexes[1]=index[1];
    indexes[2]=index_prime[0];
    indexes[3]=index_prime[1];
   }
   int indexes[4];
   double Dijkl;
  };
  vector <indexes_and_elements> DM2;
  void Read_use_DMN();
 public:
  DMN_OPS();
  DMN_OPS(string, int);
  ~DMN_OPS();
  DMN_OPS(const DMN_OPS &);
  //Set fchk file
  void set_fchk(string name_fchk,string name_log,bool wfn_fchk,bool log,bool cas,int multiplicity);
  //Set threshold
  void set_thershold(double );
  //Set point for intracule
  void set_intra(double [3]);
  //Set store or not the DMN
  void dmnincore(bool,double);
  //Calculate trace
  double calc_tr();
  //Evaluate at r1,r2,...r1',r2'...
  double evaluation(double *,double *);//Density for r1=r1'
  //Density  and gradient at r1=r1'
  void grad_rho_r(double [3],double [3],double &);
  //Diagonalize the DM1 and get NO coefficients (spins summed)
  void diagonalize();
  //Diagonalize the DM1 and get NO coefficients (spins separeted)
  void diagonalize_ab();
  //Build NO(i)
  void build_NO_DMN(double &NO,double point[3],bool &alpha_beta,int &NumNO);
  //Size of the basis
  int nbasis();
  //Fill in AOs Total Density for FCI FCHK
  void tot_dens_fchk(double *tot_density,int &aos);
  //Fill in AOs Spin Density for FCI FCHK
  void spin_dens_fchk(double *tot_density,int &aos);
  //Calculate two particle density at r1=r1' and r2=r2'
  double dm2_coalesc(double Point[3],double Point2[3]);
  READ_FCHK_WFN *FCHK_for_DMN;
  string name_file;
  double **rho_matrix,**NO_coef,**rho_matrixa,**NO_coefa,**rho_matrixb,**NO_coefb,send_threshold;
  bool fchk,dm1,dm2,dm3,dm4,diag_done,diag_ab_done;
};
#endif // _DMN_OPS_CLASS_H_


