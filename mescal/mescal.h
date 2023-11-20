#ifndef _MESCAL_H_
#define _MESCAL_H_
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<algorithm>
#include<cmath>
#define Angs2au 1.8897259886  // Transform length units
#define au2eV 27.211399       // Transform energy units
#define au2Debye 2.541746276  // To transform dipole moment units
#define tol4 1e-4
#define tol5 1e-5
#define tol8 1e-8

using namespace std;
class MESCAL
{
 private:
  bool conver_E=false,conver_mu=false,conver_q=false;
  int order[3],nactive=-1;
  double Energy_old;
  string label;
  struct ATOM
  {
   int Z;
   double pos[3],q_perm=0.0e0,q_ind=0.0e0,V_ext=0.0e0,V_q_ind=0.0e0,V_q_perm=0.0e0,V_mu_ind=0.0e0,dipole_ind[3]={0.0e0},pos_wrt_cm[3]={0.0e0};
   double mu_ind[3],F_ext[3]={0.0e0},F_q_ind[3]={0.0e0},F_q_perm[3]={0.0e0},F_mu_ind[3]={0.0e0},alpha[3][3]={0.0e0};
  };
  void Asymbol2Z(int &Z, string symbol);
  int Z2val_electrons(int &Z);
  void Z2label(int Z);
  double Z2mass(int &Z);
  double Z2atomic_pol(int &Z);
  void jacobi(int n, double **m, double **v);
  void init_output(string name_output);
  void close_output(string name_output);
  void read_pdb_file(string name_pbd);
  void read_fragment_file(string name_frag,double **Im,double **Urot,int &ifrag,int &Sum_Val_elect, double &Sum_atomic_pol);
  void Frag_T_inertia(int &ifrag,double Rcm[3],double **Im,double **Urot);
  void Frag_T_inertia_compare(int &ifrag, double **Cartes_coord, int *Zfrag, double Rcm[3],double **Im,double **Urot);
  void set_FV_inter_frag(bool &induced_q,bool &permanent_q);
  void alphaF2mu(int &ifrag, int &iatom, double Field[3]);
  void update_mu_q_ind();
  void print_init_sc(string name);
  void print_end_sc(string name);
  void print_iter_info(string name);
  void print_iactive_info(string name);

 public:
  MESCAL();
  MESCAL(string,string,bool&,bool&);
  ~MESCAL();
  string sha="";
  bool perm_q=false,ind_q=false,part_val_e=false;
  int nfragments,maxiter=1000,iter=0;
  double r0=0.0e0,w_damp=0.4,mu_diff_max,q_diff_max,E_diff,threshold_mu,threshold_E,threshold_q,Energy;
  struct FRAGMENT
  {
   string name;
   int natoms;
   vector<ATOM>atoms;
   double **Pi,Rcm[3],dist_RcmO; // Pi matrix for the atom,atom susceptibility Pi[a_tom][b_atom]
   bool active=true;
  };
  vector<FRAGMENT>fragments;
  int natoms_tot();
  void get_coords(double **Coords);
  void set_FV_ext_qm(double **F_ext,double *V_ext);
  void set_FV_ext_punct(double &q_mescal,double Point_mescal[3]);
  void clean_FV_ext_punct();
  void mescal_scs(string name_output);
  void calc_E(string name_output);
  void deactivate_fragments(double &rad);

}; 
#endif // _MESCAL_H_
