#ifndef _MESCAL_H_
#define _MESCAL_H_
#include<iostream>
#include<iomanip>
#include<omp.h>
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
class Mescal
{
 private:
  bool conver_E=false,conver_mu=false,conver_q=false;
  int order[3],nactive=-1;
  double Energy_old,radius;
  string label,mescal_ofile;
  struct ATOM
  {
   int Z;
   double pos[3],q_perm=0.0e0,q_ind=0.0e0,V_ext=0.0e0,V_q_ind=0.0e0,V_q_perm=0.0e0,V_mu_ind=0.0e0,dipole_ind[3]={0.0e0},pos_wrt_cm[3]={0.0e0};
   double mu_ind[3],F_ext[3]={0.0e0},F_q_ind[3]={0.0e0},F_q_perm[3]={0.0e0},F_mu_ind[3]={0.0e0},alpha[3][3]={0.0e0};
  };
  int Z2val_electrons(int &Z);
  void Z2label(int Z);
  double Z2atomic_pol(int &Z);
  void init_output();
  void close_output();
  void read_pdb_file(string name_pbd);
  void read_fragment_file(string name_frag,double **Im_frag,double **Urot2align,int &ifrag,int &Sum_Val_elect, double &Sum_atomic_pol);
  void Frag_T_inertia(int &ifrag,double Rcm[3],double **Im,double **Urot);
  void Frag_T_inertia_compare(int &ifrag, double **Cartes_coord, int *Zfrag, double Rcm[3],double **Im,double **Urot);
  void set_FV_inter_frag(bool &induced_q,bool &permanent_q);
  void alphaF2mu(int &ifrag, int &iatom, double Field[3]);
  void update_mu_q_ind();
  void print_charges_file();
  void print_init_sc();
  void print_end_sc();
  void print_iter_info();

 public:
  Mescal();
  Mescal(string,string,bool&,bool&,bool&);
  Mescal(const Mescal&Mescal_obj);
  ~Mescal();
  string sha="";
  bool perm_q=false,ind_q=false,part_val_e=false,mute=true,deact_rad=false;
  int nfragments,maxiter=1000,iter=0,ifrac_deact=-1,nthread;
  double r0=0.0e0,w_damp=0.4,mu_diff_max,q_diff_max,E_diff,threshold_mu,threshold_E,threshold_q,Energy;
  struct FRAGMENT
  {
   string name;
   int natoms;
   vector<ATOM>atoms;
   double **Pi,Rcm[3],dist_RcmO; // Pi matrix for the atom,atom susceptibility Pi[a_tom][b_atom]
   bool active=true,i_was_active=false;
  };
  vector<FRAGMENT>fragments;
  int natoms_tot();
  void mescal_get_frag_info();
  void get_coords(double **Coords);
  void set_FV_ext_one(int &ifrag,int &iatom,double F_ext[3],double &V_ext);
  void set_FV_ext_qm(double **F_ext,double *V_ext);
  void set_FV_ext_punct(double &q_mescal,double Point_mescal[3]);
  void get_V_punct(double &V_r,double Point_r[3]);
  void get_ind_q_frag_atom(int &ifrag,int &iatom, double q_charge[2],double coord_q[2][3]);
  void clean();         // It also saves information of fragments that were active before
  void clean_converg(); // Clear only bool variables needed to rerun SCF without modifying the active states
  void mescal_scs();
  void calc_E();
  void deactivate_fragment(int &natoms_in,int *Z,double **coords);
  void deactivate_fragments(double &rad);

};
double Z2mass(int &Z);
void Asymbol2Z(int &Z, string symbol);
void xyz_to_new_xyz(string name_xyz); 
void jacobi_mescal(int n, double **m, double **v);
#endif // _MESCAL_H_
