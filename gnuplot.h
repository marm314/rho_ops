#ifndef _GNUPLOT_H_
#define _GNUPLOT_H_

#include<iostream>
#include<fstream>
#include<string>
#include<stdlib.h>

using namespace std;

class GNUPLOT
{
 private:
  string quotation;
 public:
  GNUPLOT();
  ~GNUPLOT();
  void plot3D(string,string,string,string,string *,int,char,char,double [2],double [2]);
  void plot2D(string,string,string,string,string *,int,char,double [2]);
  string path;
};
#endif // _GNUPLOT_H_
