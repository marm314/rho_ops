#ifndef _MESCAL_H_
#define _MESCAL_H_
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<algorithm>
#include<cmath>
#define Angs2au 1.8897259886
#define tol8 1e-8

using namespace std;
class MESCAL
{
 private:
 int order[3];
 struct ATOM
 {
  int Z;
  double pos[3],charge=0.0e0,charge_ind=0.0e0,V_ext=0.0e0,V_q_ind=0.0e0,V_q_perm=0.0e0,dipole_ind[3]={0.0e0},pos_wrt_cm[3]={0.0e0};
  double F_ext[3]={0.0e0},F_q_ind[3]={0.0e0},F_q_perm[3]={0.0e0},alpha[3][3]={0.0e0};
 };
 void Asymbol2Z(int &Z, string symbol);
 double Z2mass(int &Z);
 int Z2val_electrons(int &Z);
 void jacobi(int n, double **m, double **v);
 void read_pdb_file(string name_pbd);
 void read_fragment_file(string name_frag,double **Im,double **Urot,int &ifrag,int &Sum_Val_elect);
 void Frag_T_inertia(int &ifrag,double Rcm[3],double **Im,double **Urot);

 public:
  MESCAL();
  MESCAL(string,string);
  ~MESCAL();
 int nfragments;
 struct FRAGMENT
 {
  string name;
  int natoms;
  vector<ATOM>atoms;
 };
 vector<FRAGMENT>fragments;
 void init_output(string name_output);
 void set_FV_ext_punct(double &q_mescal,double Point_mescal[3]);
 void set_FV_q_inter_frag(bool &induced);
 void close_output(string name_output);

}; 
#endif // _MESCAL_H_
