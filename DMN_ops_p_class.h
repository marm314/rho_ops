//This class reads .dmn files produced by DMN code of
// Dr. Eduard Matito. & Performs some calculations
// including MOps from fchk files
#ifndef _DMN_OPS_P_CLASS_H_
#define _DMN_OPS_P_CLASS_H_
#include<fstream>
#include<omp.h>
#include"D_read_calc_rho.h"
#include"MOp_class.h"
#include"Numbers.h"
using namespace std;
 class DMN_P_OPS
{
 private:
  int element[4],element_prime[4],RECORD_DELIMITER_LENGTH,nbasisf;
  double threshold,Dij,Dijkl,Dijklmn,**DM1;
  double POINT1[3],POINT_PRIME1[3];
  complex<double> **Array_MOs,**Grad_MOs;
  complex<double> sum,conjugated,conjugated2,Grad[3];
  double Trace;
  bool dmn_file,trace,evaluate,grad_density,DMNinCORE;
  bool diagonalize_dm1,diagonalize_dm1_ab;
  MOp mos_array[4];
  void Read_use_DMN();
 public:
  DMN_P_OPS();
  DMN_P_OPS(string, int);
  ~DMN_P_OPS();
  DMN_P_OPS(const DMN_P_OPS &);
  //Set fchk file
  void set_fchk(string name_fchk,string name_log,bool wfn_fchk,bool log,bool cas,int multiplicity);
  //Set threshold
  void set_thershold(double );
  //Set store or not the DMN
  void dmnincore(bool,double);
  //Calculate trace
  double calc_p_tr();
  //Evaluate at p1,p2,...p1',p2'...
  double evaluation_p(double *,double *);
  //Density  and gradient at p1=p1'
  void grad_rho_p(double [3],double [3],double &);
 //Diagonalize the DM1 and get NO coefficients (spins summed)
  void diagonalize();
  //Diagonalize the DM1 and get NO coefficients (spins separeted)
  void diagonalize_ab();
  //Build NOp(i)
  void build_NOp_DMN(complex<double> &NOp,double point[3],bool &alpha_beta,int &NumNO);
  //Size of the basis
  int nbasis();
  READ_FCHK_WFN *FCHK_for_DMN;
  string name_file;
  double **rho_matrix,**NO_coef,**rho_matrixa,**NO_coefa,**rho_matrixb,**NO_coefb,send_threshold;
  bool fchk,dm1,dm2,dm3,dm4,diag_done,diag_ab_done;
};
#endif // _DMN_OPS_P_CLASS_H_



