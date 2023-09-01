#ifndef _MESCAL_H_
#define _MESCAL_H_
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<algorithm>
#define Angs2au 1.8897259886  

using namespace std;
class MESCAL
{
 private:
 struct ATOM
 {
  int Z;
  double charge,pos[3],dipole[3],alpha[3][3];
 };
 void Asymbol2Z(int &Z, string symbol);
 double Z2mass(int &Z);

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

}; 
#endif // _MESCAL_H_
