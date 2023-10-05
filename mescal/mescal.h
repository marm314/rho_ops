#ifndef _MESCAL_H_
#define _MESCAL_H_
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<algorithm>
#include<cmath>
#define Angs2au 1.8897259886
#define au2Debye 2.541746276  // To transform dipole moments
#define tol4 1e-4
#define tol6 1e-6
#define tol8 1e-8

using namespace std;
class MESCAL
{
 private:
  bool conver_E=false,conver_mu=false,conver_q=false;
  int order[3];
  double Energy_old;
  struct ATOM
  {
   int Z;
   double pos[3],q_perm=0.0e0,q_ind=0.0e0,V_ext=0.0e0,V_q_ind=0.0e0,V_q_perm=0.0e0,V_mu_ind=0.0e0,dipole_ind[3]={0.0e0},pos_wrt_cm[3]={0.0e0};
   double mu_ind[3],F_ext[3]={0.0e0},F_q_ind[3]={0.0e0},F_q_perm[3]={0.0e0},F_mu_ind[3]={0.0e0},alpha[3][3]={0.0e0};
  };
  void Asymbol2Z(int &Z, string symbol);
  double Z2mass(int &Z);
  int Z2val_electrons(int &Z);
  double Z2atomic_pol(int &Z);
  void jacobi(int n, double **m, double **v);
  void read_pdb_file(string name_pbd);
  void read_fragment_file(string name_frag,double **Im,double **Urot,int &ifrag,int &Sum_Val_elect, double &Sum_atomic_pol);
  void Frag_T_inertia(int &ifrag,double Rcm[3],double **Im,double **Urot);
  void set_FV_q_inter_frag(bool &induced);
  void alphaF2mu(int &ifrag, int &iatom, double Field[3]);
  void update_mu_q_ind();
  void set_FV_mu_ind();
  void print_init_sc(string name);
  void print_end_sc(string name);
  void print_iter_info(string name);

 public:
  MESCAL();
  MESCAL(string,string,bool&,bool&);
  ~MESCAL();
  bool perm_q=false,ind_q=false,part_val_e=false;
  int nfragments,maxiter=1000,iter=0;
  double r0=0.0e0,w_mu=0.4,mu_diff_max,q_diff_max,E_diff,threshold_mu,threshold_E,threshold_q,Energy;
  struct FRAGMENT
  {
   string name;
   int natoms;
   vector<ATOM>atoms;
   double **Pi;
  };
  vector<FRAGMENT>fragments;
  void init_output(string name_output);
  void set_FV_ext_punct(double &q_mescal,double Point_mescal[3]);
  void mescal_scs(string name_output);
  void calc_E(string name_output);
  void close_output(string name_output,string sha);

}; 
#endif // _MESCAL_H_
