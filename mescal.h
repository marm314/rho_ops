#ifndef _MESCAL_H_
#define _MESCAL_H_
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<algorithm>

using namespace std;
class MESCAL
{
 private:
 struct ATOM
 {
  int Z;
  double charge,pos[3],dipole[3];
 };

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
