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
  double charge,pos[3],dipole[3],pos_wrt_cm[3],alpha[3][3];
 };
 void Asymbol2Z(int &Z, string symbol);
 double Z2mass(int &Z);
 void jacobi(int n, double **m, double **v);
 void read_pdb_file(string name_pbd);
 void read_fragment_file(string name_frag,double **Im,double **Urot,int &ifrag);
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
 void close_output(string name_output);

}; 
#endif // _MESCAL_H_
