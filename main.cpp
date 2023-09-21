///////////////////////////////////////////////////////////////////
// Build the density from any fchk, wfn or wfx  files            //
// preform several operations with it                            //
// (e.g., Complexity analysis[Shannon entropy, Fisher Integral,  //
// etc.])                                                        //
// This program belongs to Mauricio Rodriguez Mayorga            //
// Ph.D. student at University of the Basque Country             //
// for support and comments email to: marm3.14@gmail.com         //
///////////////////////////////////////////////////////////////////
#include<iostream>
#include<algorithm>
#include<omp.h>
#include"Mathematical_Functions.h"
#include"D_read_calc_rho.h"
#include"Input_commands.h"
#include"NO_class.h"
#include"MO_class.h"
#include"MOp_class.h"
#include"NO_DMN_class.h"
#include"NOp_DMN_class.h"
#include"Integrals.h"
#include"Integrals_DMN.h"
#include"Integrals_quadrature.h"
#include"Integrals_Becke.h"
#include"gauss_quad.h"
#include"Numbers.h"
#include"String_ops.h"
#include"DMN_ops_class.h"
#include"DMN_ops_p_class.h"
#include"gnuplot.h"
#include"N_fchks_wfns.h"
#include"mescal.h"
#include"gitver.h"
#define Angs2au 1.8897259886
#define au2eV 27.211399
//////////////////////////////
//This class already loads  //
//#include <string>         //
//#include <fstream>        //
//#include <iomanip>        //
//#include <stdio.h>        //
//#include <stdlib.h>       //
//#include <cmath>          //
//#include <locale>         //
//#include <sstream>        //
/////////////////////////////
/////////////////////////////
/////////////////////////////
bool firstcall=true,wfn_fchk=false;
string name_file;
const int RECORD_DELIMITER_LENGTH=4;
/////////////////////////////
using namespace std;
//////////////////////////
//Functions declaration //
//////////////////////////
void punctual(READ_FCHK_WFN &Read_fchk_wfn,double Point[3], double &Density,double &Density_alpha,
double &Density_beta, double Grad[3], double Grad_alpha[3], double Grad_beta[3],
double &Grad_norm,double &Grad_norm_alpha,double &Grad_norm_beta,double &laplacian_r,
double &laplacian_alpha, double &laplacian_beta,double &tauW_alpha,double &tauW_beta,
double &tau_alpha,double &tau_beta,double &k_F_alpha,double &k_F_beta,double &k_s_alpha,
double &k_s_beta,double &s_r_alpha,double &s_r_beta,double &q_red_alpha,double &q_red_beta,
bool wfn_fchk);
void pre_elf(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &Density_alpha,double &Density_beta,
double &tauW_alpha,double &tauW_beta,double &tau_alpha,double &tau_beta,double &tcurr_alpha,double &tcurr_beta);
void mos_to_nos_int_fchk_dm1(READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands,double **ORBITALS,int &total_grid,int &nbasis);
void mos_to_nos_dmn_sij(double **SIJ,READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands);
void print_int(READ_FCHK_WFN &Read_fchk_wfn,string name_file,double **Sij, int nbasis,double &rho,double &rhoa,double &rhob,
string region);
void transform_int(string region,int &nbasis);
void cube_file(READ_FCHK_WFN &Read_fchk_wfn,string name_file,string op_cube,double cubex,double cubey, double cubez,double stepx,
double stepy,double stepz);
double Deviation_idemp(READ_FCHK_WFN &Read_fchk_wfn);
double ID_ni(READ_FCHK_WFN &Read_fchk_wfn);
double IND_ni(READ_FCHK_WFN &Read_fchk_wfn);
double Shannon_ni(READ_FCHK_WFN &Read_fchk_wfn);
double Deviation_idemp_dmn(DMN_OPS &DMN,double &N);
double ID_ni_dmn(DMN_OPS &DMN,double &N);
double IND_ni_dmn(DMN_OPS &DMN,double &N);
void ID_IND_local(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &ID_alpha,double &ID_beta,double &IND_alpha,double &IND_beta);
double Shannon_ni_dmn(DMN_OPS &DMN,double &N);
void calc_time(double DATE[2][4]);
///////////////////////////
// Main Function       ////
// (asks for options   ////
// and writes results) ////
///////////////////////////
int main(int argc, char *argv[])
{
 int i,j,points,num_integrals,nbasis,region_int=1;
 bool fish=false,shan=false,fishp=false,shanp=false,inertiar=false,inertiap=false,pos=false,mom=false,rho=false;
 bool r1=false,r2=false,rm1=false,p1=false,p2=false,not_elf=false,sij=false,r1_moment=false,mos_to_nos_dmn=false,dipolar=false,not_indic=false;
 double counter,Density,Density_alpha,Density_beta,Nelec,tauW_alpha,ELF,ELFa,ELFb,coef_elf;
 double Grad_norm,Grad_norm_alpha,Grad_norm_beta,tauW_beta,tau_alpha,tau_beta,k_F_alpha;
 double k_F_beta,k_s_alpha,k_s_beta,q_red_alpha,q_red_beta,s_r_alpha,s_r_beta,tcurr_alpha,tcurr_beta;
 double laplacian_r,laplacian_alpha,laplacian_beta,laplacian_p,step,I_dyn,I_ndyn,ID_alpha,ID_beta;
 double IND_alpha,IND_beta,Intracule,Extracule,DORI;
 double Point[3]={ZERO,ZERO,ZERO},Grad[3],Grad_alpha[3],Grad_beta[3],RCC[3]={ZERO},mu[3]={ZERO},LOCAL_HYBRIDS_fr[5],DATE[2][4];
 double **Inertia,**Quadrupole,**eigenV,**TPS,shannon=ZERO,shannonp=ZERO,fisher=ZERO,fisherp=ZERO,Tw,Ttf;
 double Integrals_interval[6],Rot_grid_matrix[3][3]={ZERO};
 char direct;
 string method,operation,region_string,line,name_saved,sha;
 ofstream Results;
 ifstream my_chk,date_file,CM_file;
/////////////////////
//Get git version  //
/////////////////////
 gitversion(sha);
/////////////////////
//Begin clock      //
/////////////////////
 system("date +%j' '%H' '%M' '%S>date_RHO.date");
 date_file.open("date_RHO.date");
 date_file>>DATE[0][0]>>DATE[0][1]>>DATE[0][2]>>DATE[0][3];
 date_file.close();
 system("rm date_RHO.date");
/////////////////////
//Ask for commands //
//     file        //
/////////////////////
 if(argc!=2)
 {
  cout<<"Introduce the name of the file with RHO.OPS commands"<<endl;
  do
  {getline(cin,name_file);
  }while(name_file=="");
  name_file.erase(std::remove_if(name_file.begin(),name_file.end(),::isspace),name_file.end());
 }
 else
 {
  string str(argv[1]);
  name_file=str;
 }
 Input Input_commands(name_file);
 name_file=Input_commands.name_fchk_wfn;
 my_chk.open((name_file).c_str());
//Check if the file chk or wfn exist:
 if(my_chk.good()) //Check existence of file
 {
  my_chk.close();//Close it after checked!
  /////////////////////////////////////////////////
  if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
  {
   wfn_fchk=true;//True for wfn
   Results.open((name_file.substr(0,(name_file.length()-4))+"_RHO.out").c_str());
  }
  else
  {Results.open((name_file.substr(0,(name_file.length()-5))+"_RHO.out").c_str());}
  READ_FCHK_WFN Read_fchk_wfn(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,
  Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
  //Check if wfn is correlated and print warning for spin dependant properties
  if((Read_fchk_wfn.wfn)&&(!Read_fchk_wfn.error_opens_wfn))
  {
   Read_fchk_wfn.rho_eval_a_b(Point,Density_alpha,Density_beta);
   if(Read_fchk_wfn.error_opens_wfn)
   {
    if(Input_commands.spin_calcs)
    {
     cout<<Read_fchk_wfn.identity<<" calculation unable to proceed."<<endl;
     Read_fchk_wfn.READ_FCHK_WFN::~READ_FCHK_WFN();
     exit(EXIT_FAILURE);
    }
   }
  }
  // If there are Imag Coefs store them
  if(Input_commands.im_wfn_wfx)
  {
   name_saved=name_file;
   /////////////////////////////////////////////////
   //Store 2nd WFN or WFX
   /////////////////////////////////////////////////
   wfn_fchk=false;
   name_file=Input_commands.second_fchk_wfn;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    wfn_fchk=true;//True for wfn
   }
   READ_FCHK_WFN Read_fchk_wfn_2(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,
   Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
   if(Input_commands.cm)
   {
    if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
    {
     system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
    }
    else
    {
     system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
    }
   }
   name_file=name_saved;
   Read_fchk_wfn.im_wfn_wfx=Input_commands.im_wfn_wfx;
   Read_fchk_wfn.Init_MOim(Read_fchk_wfn_2.MOcoefA);
  }
  //If everything is ok and is an fchk file, we set BETA_MOS variable in the Read_fchk_wfn object.
  //Be aware that set BETA_MOS=true will imply that impaired MOs are built using Beta MOs not Alpha MOs from fchk
  //also notice that this option is pointless when there are no Beta MOs in the fck files (close-shell resctricted).
  //By default BETA_MOS will be false for the Read_fchk_wfn object which means that impaired MOs are also alpha by default.
  if(!Read_fchk_wfn.wfn)
  {Read_fchk_wfn.set_BETA_MOS(Input_commands.Beta_MOs);}
  ////////////////////////////////////
  //Program information printed in  //
  //and Copy Right                  //
  //          FILE_RHO.out          //
  ////////////////////////////////////
  Results<<"#******************************************************************************#";
  Results<<endl;
  Results<<"#******************************************************************************#";
  Results<<endl;
  Results<<"#                           Welcome to RHO.OPS                                 #";
  Results<<endl;
  Results<<"#******************************************************************************#";
  Results<<endl;
  Results<<"#******************************************************************************#";
  Results<<endl;
  Results<<"# Program developed for building the density from fchk, wfn, or wfx files      #";
  Results<<endl;
  Results<<"# then perform several operations with it (e.g., Complexity analysis[Shannon   #";
  Results<<endl;
  Results<<"# entropy, Fisher Integral, etc.]; find bond, ring and cage critical points;   #";
  Results<<endl;
  Results<<"# calculate several quantities for local-hybrid functionals[tw(r),tau(r),etc]; #";
  Results<<endl;
  Results<<"# build and evaluate MOs and NOs as well as calculate overlaps; etc.)          #";
  Results<<endl;
  Results<<"#******************************************************************************#";
  Results<<endl;
  Results<<"# Copyright (C) 2014 M.Sc. Mauricio Rodriguez-Mayorga                          #";
  Results<<endl;
  Results<<"# Ph.D. student at University of the Girona (2018)                             #";
  Results<<endl;
  Results<<"# Postdoctoral reseracher at Vrije Universiteit Amsterdam (2021)               #";
  Results<<endl;
  Results<<"# for support and comments send an email to: marm3.14@gmail.com                #";
  Results<<endl;
  Results<<"#******************************************************************************#";
  Results<<endl;
  Results<<"#********************************************************************************#";
  Results<<endl;
  Results<<"# Remember to cite this program:                                                 #";
  Results<<endl;
  Results<<"# Rodriguez-Mayorga, M.; RHO.OPS: Density Operations Program.                    #";
  Results<<endl;
  Results<<"# Donostia International Physics Center (DIPC)                                   #";
  Results<<endl;
  Results<<"# University of the Basque Country, Donostia, Guipuzkoa, Spain, 2014;            #";
  Results<<endl;
  Results<<"#                                                                                #";
  Results<<endl;
  Results<<"# and please also cite the papers that deal with the most important features of  #";
  Results<<endl;
  Results<<"# the program:                                                                   #";
  Results<<endl;
  Results<<"#· Hahn. T.; Cuba - a library for multidimensional numerical integration.        #";
  Results<<endl;
  Results<<"#  Comput. Phys. Commun. 168, 78-95 (2005) DOI:10.1016/j.cpc.2005.01.010         #";
  Results<<endl;
  Results<<"#· Johnson S. G.; Cubature Package, http://ab-initio.mit.edu/cubature            #";
  Results<<endl;
  Results<<"#· Genz A. C. and Malik A. A.; An adaptive algorithm for numeric integration     #";
  Results<<endl;
  Results<<"#  over an N-dimensional rectangular region.                                     #";
  Results<<endl;
  Results<<"#  J. Comput. Appl. Math. 6 (4), 295-302 (1980).                                 #";
  Results<<endl;
  Results<<"#· Berntsen J., Espelid. T. O., and Genz. A.; An adaptive algorithm for the      #";
  Results<<endl;
  Results<<"#  approximate calculation of multiple integrals.                                #";
  Results<<endl;
  Results<<"#  ACM Trans. Math. Soft. 17 (4), 437-451 (1991).                                #";
  Results<<endl;
  Results<<"#********************************************************************************#";
  Results<<endl;
  Results<<endl;
  ///////////////////////////////////////////////////////////////
  //Print inertia tensor of nuclei CM and coord. reorientation //
  ///////////////////////////////////////////////////////////////
  if(Input_commands.cm)
  { 
   Results<<"Reorientation for file "<<name_file<<endl;
   Results<<endl;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    CM_file.open((name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
    while(getline(CM_file,line))
    {
     Results<<line<<endl;
    }
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
   }
   else
   {
    CM_file.open((name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
    while(getline(CM_file,line))
    {
     Results<<line<<endl;
    }
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
   }
  }
  ////////////////////////////////////////////
  // Deviation from idempotency & other     //
  // correlation indicators                 //
  ////////////////////////////////////////////
  if(Input_commands.indicators)
  {
   double Idem,S;
   I_dyn=ZERO;I_ndyn=ZERO;
   if(Input_commands.dmn_indicators)
   {
    DMN_OPS dmn(Input_commands.name_dm1,1);
    dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
    dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
    if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
    Nelec=dmn.calc_tr();
    dmn.diagonalize_ab();
    I_dyn=ID_ni_dmn(dmn,Nelec);
    I_ndyn=IND_ni_dmn(dmn,Nelec);
    S=Shannon_ni_dmn(dmn,Nelec);
    Idem=Deviation_idemp_dmn(dmn,Nelec);
   }
   else
   {
    I_dyn=ID_ni(Read_fchk_wfn);
    I_ndyn=IND_ni(Read_fchk_wfn);
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#              Deviation from Idempotency, I_T, I_D, I_ND                 #";
   Results<<endl;
   Results<<"#                   and Shannon entropy of ocupations                     #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   if(Input_commands.dmn_indicators)
   {
    Results<<"I_2                : "<<setw(17)<<Idem<<endl;
    Results<<"[I2(ni)=sum_i (ni/Ne)(1-ni)]"<<endl;
   }
   else
   {
    Results<<"I_2                : "<<setw(17)<<Deviation_idemp(Read_fchk_wfn)<<endl;
    Results<<"[I2(ni)=sum_i (ni/Ne)(1-ni)]"<<endl;
   }
   Results<<setprecision(10)<<fixed<<scientific;
   Results<<"I_D                : "<<setw(17)<<I_dyn<<endl;
   Results<<"[ID(ni)= 1/2 sum_i 1/2 sqrt[ni(1-ni)]- ni(1-ni)]"<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Results<<"I_ND               : "<<setw(17)<<I_ndyn<<endl;
   Results<<"[IND(ni)= 1/2 sum_i ni(1-ni)]"<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Results<<"I_T                : "<<setw(17)<<I_dyn+I_ndyn<<endl;
   Results<<"[IT(ni)= 1/4 sum_i sqrt[ni(1-ni)]]"<<endl;
   Results<<"%I_D               : "<<setw(17)<<pow(TEN,TWO)*I_dyn/(I_dyn+I_ndyn+pow(TEN,-TWO*TEN))<<endl;
   Results<<"[%ID(ni)= 100 ID/(ID+IND)]"<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Results<<"%I_ND              : "<<setw(17)<<pow(TEN,TWO)*I_ndyn/(I_dyn+I_ndyn+pow(TEN,-TWO*TEN))<<endl;
   Results<<"[%IND(ni)= 100 IND/(ID+IND)]"<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   if(Input_commands.dmn_indicators)
   {
    Results<<"S(ni)              : "<<setw(17)<<S<<endl;
   }
   else
   {
    Results<<"S(ni)              : "<<setw(17)<<Shannon_ni(Read_fchk_wfn)<<endl;
   }
   Results<<"[S(ni)=-sum_i (ni/Ne)ln(ni/Ne)]"<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ///////////////////////////////////////////////////
  // Prepare .basis file, copute intracule 1rdm    //
  // and overlap matrix in primitives              //
  ///////////////////////////////////////////////////
  if(Input_commands.intra_1rdm_sij)
  {
   if(name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')
   {
    Results<<"#*************************************************************************#";
    Results<<endl;
    Results<<"#            Preparing basis and dm1 files in primitives,                 #";
    Results<<endl;
    Results<<"#   computing intracule-like integrals, and primitives overlap matrix     #";
    Results<<endl;
    Results<<"#*************************************************************************#";
    Results<<endl;
    Results<<"#*************************************************************************#";
    Results<<endl;
    string name_basis,name_dm1,name_in;
    bool overlap=false,last=false;
    int mode=1,nucleous,nr,nang,elements[2];
    int **Quant,Lmax=0,Nroot_Lmax_plus_1,nx_exp[2],ny_exp[2],nz_exp[2];
    double nx,ny,nz,Exp,Dij,alpha,a,b;
    double Atom_coord[3]={ZERO},Atom_coord2[3]={ZERO};
    double **Coef,**CoefT,**Temp,**dm1_prim,**dm1,**Smunu,**Intra_1,*r_gauss,*w_gauss;
    Quant=new int*[35];
    for(i=0;i<35;i++)
    {
     Quant[i]=new int[3];
    }
    Read_fchk_wfn.Quant_fill(Quant,0);
    name_in=name_file.substr(0,name_file.length()-3); 
    name_basis=name_in+"basis";
    name_dm1=name_in+"dm1";
    ofstream basis_out(name_basis.c_str());
    // Basis file and primitives quadrature preparation
    basis_out<<"# Z  \t\t X_A \t\t Y_A \t\t Z_A \t\t Exp \t\t Quant(nx,ny,nz)\n";
    CoefT=new double*[Read_fchk_wfn.nprimitv];
    dm1_prim=new double*[Read_fchk_wfn.nprimitv];
    for(i=0;i<Read_fchk_wfn.nprimitv;i++)
    {
     CoefT[i]=new double[Read_fchk_wfn.nbasis()];
     dm1_prim[i]=new double[Read_fchk_wfn.nprimitv];
     nucleous=(int)Read_fchk_wfn.Nu_charge[Read_fchk_wfn.shell_map[i]-1];
     nx=(double)Quant[Read_fchk_wfn.shell_type[i]-1][0];
     ny=(double)Quant[Read_fchk_wfn.shell_type[i]-1][1];
     nz=(double)Quant[Read_fchk_wfn.shell_type[i]-1][2];
     if((Read_fchk_wfn.shell_type[i]-1)>Lmax){Lmax=Read_fchk_wfn.shell_type[i]-1;}
     Exp=Read_fchk_wfn.Prim_exp[i];
     Atom_coord[0]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[i]-1][0];
     Atom_coord[1]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[i]-1][1];
     Atom_coord[2]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[i]-1][2];
     basis_out<<setw(3)<<nucleous<<"\t"<<setprecision(8)<<fixed<<scientific<<Atom_coord[0]<<"\t"<<Atom_coord[1]<<"\t"<<Atom_coord[2];
     basis_out<<"\t"<<Exp<<"\t\t"<<(int)nx<<"\t"<<(int)ny<<"\t"<<(int)nz<<endl;
    }    
    basis_out.close();
    if(Lmax==3){Lmax=1;}
    else if(Lmax==9){Lmax=2;}
    else if(Lmax==19){Lmax=3;}
    else if(Lmax==34){Lmax=4;}
    else {Lmax=0;}
    Nroot_Lmax_plus_1=Lmax+1;
    Results<<"Lmax                  :"<<setw(10)<<Lmax<<endl;
    Results<<"Quadrature rule order :"<<setw(10)<<Nroot_Lmax_plus_1<<endl;
    Results<<"[Quadrature for primitive integrals. Order = Lmax + 1 ]"<<endl;
    alpha=ZERO;
    a=ZERO;
    b=ONE;
    gauss_hermite_rule(name_in,alpha,a,b,Nroot_Lmax_plus_1);
    r_gauss=new double[Nroot_Lmax_plus_1];
    w_gauss=new double[Nroot_Lmax_plus_1];
    //Read quadrature info
    ifstream read_quad;
    // Read weights
    read_quad.open((name_in+"_w.txt").c_str());
    for(i=0;i<Nroot_Lmax_plus_1;i++)
    {
     read_quad>>w_gauss[i];
    }
    read_quad.close();
    // Read roots
    read_quad.open((name_in+"_x.txt").c_str());
    for(i=0;i<Nroot_Lmax_plus_1;i++)
    {
     read_quad>>r_gauss[i];
    }
    read_quad.close();
    // Remove quadrature files
    system(("rm "+name_in+"_r.txt").c_str());
    system(("rm "+name_in+"_w.txt").c_str());
    system(("rm "+name_in+"_x.txt").c_str());
    // DM1 file and scan 'intracule like' coordinate
    dm1=new double*[Read_fchk_wfn.nbasis()];
    Coef=new double*[Read_fchk_wfn.nbasis()];
    Temp=new double*[Read_fchk_wfn.nbasis()];
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     dm1[i]=new double[Read_fchk_wfn.nbasis()];
     for(j=0;j<Read_fchk_wfn.nbasis();j++)
     {
      dm1[i][j]=ZERO;
      if(i==j)
      {
       dm1[i][i]=Read_fchk_wfn.Ocupation[i];
      }
     }
     Coef[i]=new double[Read_fchk_wfn.nprimitv];
     Temp[i]=new double[Read_fchk_wfn.nprimitv];
     for(j=0;j<Read_fchk_wfn.nprimitv;j++)
     {
      Coef[i][j]=Read_fchk_wfn.MOcoefA[i][j];
      CoefT[j][i]=Coef[i][j];
     }
    }
    matmul_full(Read_fchk_wfn.nbasis(),Read_fchk_wfn.nbasis(),Read_fchk_wfn.nprimitv,dm1,Coef,Temp);
    matmul_full(Read_fchk_wfn.nprimitv,Read_fchk_wfn.nbasis(),Read_fchk_wfn.nprimitv,CoefT,Temp,dm1_prim);
    nr=Input_commands.order_grid_r;
    nang=Input_commands.order_grid_ang;
    if(nr>0 && nang>0)
    {
     grid_avail(nang);
    }
    else
    {
     if(nr==-1 && nang==-1){overlap=true;}
     nr=1;nang=1;mode=3;
    }
    Intra_1=new double*[nr];
    for(i=0;i<nr;i++)
    {
     Intra_1[i]=new double[nang];
     for(j=0;j<nang;j++)
     {
      Intra_1[i][j]=ZERO;
     }
    }
    int k=0;
    void *data=NULL; // This is not used to we can passed it as NULL
    double **nth1=NULL,**nth2=NULL;
    integrate_quadrature(data,name_in,false,nr,nang,false,nth1,nth2,mode);
    ofstream dm1_out(name_dm1.c_str(),ios::out | ios::binary);
    ofstream overlap_out((name_in+"overlap").c_str());
    overlap_out<<setprecision(10)<<fixed<<scientific;
    for(i=0;i<Read_fchk_wfn.nprimitv;i++)
    {
     Atom_coord[0]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[i]-1][0];
     Atom_coord[1]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[i]-1][1];
     Atom_coord[2]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[i]-1][2];
     nx_exp[0]=Quant[Read_fchk_wfn.shell_type[i]-1][0];
     ny_exp[0]=Quant[Read_fchk_wfn.shell_type[i]-1][1];
     nz_exp[0]=Quant[Read_fchk_wfn.shell_type[i]-1][2];
     for(j=0;j<Read_fchk_wfn.nprimitv;j++)
     {
      if(!overlap)
      {
       if(abs(dm1_prim[i][j])>pow(TEN,-TEN)||(i==j && i==(Read_fchk_wfn.nprimitv-1))) 
       {
        if(i==j && i==(Read_fchk_wfn.nprimitv-1))
        {
         last=true;
        }
        Atom_coord2[0]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[j]-1][0];
        Atom_coord2[1]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[j]-1][1];
        Atom_coord2[2]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[j]-1][2];
        nx_exp[1]=Quant[Read_fchk_wfn.shell_type[j]-1][0];
        ny_exp[1]=Quant[Read_fchk_wfn.shell_type[j]-1][1];
        nz_exp[1]=Quant[Read_fchk_wfn.shell_type[j]-1][2];
        Dij=dm1_prim[i][j];
        integrate_intra_coord(Intra_1,Dij,Read_fchk_wfn.Prim_exp[i],Read_fchk_wfn.Prim_exp[j],Atom_coord,Atom_coord2,nx_exp,ny_exp,nz_exp,
        Nroot_Lmax_plus_1,r_gauss,w_gauss,nr,nang,last,overlap);
        elements[0]=i+1;
        elements[1]=j+1;
        dm1_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        dm1_out.write((char*) &elements[0], sizeof(elements[0]));
        dm1_out.write((char*) &elements[1], sizeof(elements[1]));
        dm1_out.write((char*) &Dij, sizeof(Dij));
        dm1_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       }
      }
      else
      {
       if(abs(dm1_prim[i][j])>pow(TEN,-TEN))
       {
        Dij=dm1_prim[i][j];
        elements[0]=i+1;
        elements[1]=j+1;
        dm1_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
        dm1_out.write((char*) &elements[0], sizeof(elements[0]));
        dm1_out.write((char*) &elements[1], sizeof(elements[1]));
        dm1_out.write((char*) &Dij, sizeof(Dij));
        dm1_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
       } 
       Atom_coord2[0]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[j]-1][0];
       Atom_coord2[1]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[j]-1][1];
       Atom_coord2[2]=Read_fchk_wfn.Cartesian_Coor[Read_fchk_wfn.shell_map[j]-1][2];
       nx_exp[1]=Quant[Read_fchk_wfn.shell_type[j]-1][0];
       ny_exp[1]=Quant[Read_fchk_wfn.shell_type[j]-1][1];
       nz_exp[1]=Quant[Read_fchk_wfn.shell_type[j]-1][2];
       Dij=ONE;
       if(i==j && i==(Read_fchk_wfn.nprimitv-1))
       {
        last=true;
       }
       integrate_intra_coord(Intra_1,Dij,Read_fchk_wfn.Prim_exp[i],Read_fchk_wfn.Prim_exp[j],Atom_coord,Atom_coord2,nx_exp,ny_exp,nz_exp,
       Nroot_Lmax_plus_1,r_gauss,w_gauss,nr,nang,last,overlap);
       if(abs(Intra_1[0][0])<pow(TEN,-TEN)){Intra_1[0][0]=ZERO;}
       overlap_out<<setw(20)<<Intra_1[0][0]; // Stored in there the S_mu,nu value
       k++;
       if(k==5){overlap_out<<endl;k=0;}
      }
     }
    }
    Dij=ZERO;
    elements[0]=0;
    elements[1]=0;
    dm1_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
    dm1_out.write((char*) &elements[0], sizeof(elements[0]));
    dm1_out.write((char*) &elements[1], sizeof(elements[1]));
    dm1_out.write((char*) &Dij, sizeof(Dij));
    dm1_out.seekp(RECORD_DELIMITER_LENGTH, ios::cur);
    dm1_out.close();
    overlap_out.close();
    if(!overlap)
    {
     system(("/bin/rm -rf "+name_in+"overlap").c_str()); // Delete it, if it is empty.
     Results<<"# In the following: I(u) is spherically-averaged [ I(u) / 4*PI ]"<<endl;
     Results<<"#      u                  I(u)                I(u)u**2"<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     for(i=0;i<nr;i++)
     {
      Results<<setw(20)<<Intra_1[i][0]<<setw(20)<<Intra_1[i][1]<<setw(20)<<Intra_1[i][2]<<endl;
     }
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      delete[] dm1[i];dm1[i]=NULL;
      delete[] Coef[i];Coef[i]=NULL;
      delete[] Temp[i];Temp[i]=NULL;
     }
    }
    else
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      delete[] Temp[i];Temp[i]=NULL;
      delete[] dm1[i];dm1[i]=NULL;
     }
     delete[] Temp;Temp=NULL;
     Smunu=new double*[Read_fchk_wfn.nprimitv];
     Temp=new double*[Read_fchk_wfn.nprimitv];
     for(i=0;i<Read_fchk_wfn.nprimitv;i++)
     {
      Smunu[i]=new double[Read_fchk_wfn.nprimitv];
      Temp[i]=new double[Read_fchk_wfn.nprimitv];
      for(j=0;j<Read_fchk_wfn.nprimitv;j++)
      {
       Temp[i][j]=ZERO;
       Smunu[i][j]=ZERO;
      }
     }
     ifstream read_Smunu((name_in+"overlap").c_str());
     for(i=0;i<Read_fchk_wfn.nprimitv;i++)
     {
      for(j=0;j<Read_fchk_wfn.nprimitv;j++)
      {
       read_Smunu>>Smunu[i][j];
      }
     }
     read_Smunu.close();
     matmul(Read_fchk_wfn.nprimitv,Smunu,dm1_prim,Temp);
     // Use this double precision variable to store Tr [1D S] = N electrons
     Intra_1[0][0]=ZERO;
     for(i=0;i<Read_fchk_wfn.nprimitv;i++)
     {
      Intra_1[0][0]=Intra_1[0][0]+Temp[i][i];
      for(j=0;j<Read_fchk_wfn.nprimitv;j++){Temp[i][j]=ZERO;}
     }
     Results<<endl;
     Results<<" N electrons:"<<setprecision(6)<<fixed<<setw(12)<<Intra_1[0][0]<<endl;
     Results<<endl;
     // S^MO = C S^Prim C^T
     matmul_full(Read_fchk_wfn.nbasis(),Read_fchk_wfn.nprimitv,Read_fchk_wfn.nprimitv,Coef,Smunu,Temp);
     matmul_full(Read_fchk_wfn.nprimitv,Read_fchk_wfn.nprimitv,Read_fchk_wfn.nbasis(),Temp,CoefT,Smunu);
     ofstream overlap_out2((name_in+"overlap_mo").c_str());
     overlap_out2<<setprecision(10)<<fixed<<scientific;
     // Now delete Coef and Temp
     k=0;
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      delete[] Coef[i];Coef[i]=NULL;
      for(j=0;j<Read_fchk_wfn.nbasis();j++)
      {
       if(abs(Smunu[i][j])<pow(TEN,-EIGHT)){Smunu[i][j]=ZERO;}
       //if(i!=j && abs(Smunu[i][j])>pow(TEN,-SEVEN)){Results<<" Warning orbitals "<<setw(5)<<i+1<<setw(5)<<j+1<<" are not orthonormal."<<endl;}
       //if(i==j && abs(Smunu[i][i]-ONE)>pow(TEN,-SEVEN)){Results<<" Warning orbital  "<<setw(5)<<i+1<<setw(5)<<i+1<<" is not orthonormal."<<endl;}
       overlap_out2<<setw(20)<<Smunu[i][j]; // Stored in there the S_mu,nu value
       k++;
       if(k==5){overlap_out2<<endl;k=0;}
      }
     }
     overlap_out2.close();
     for(i=0;i<Read_fchk_wfn.nprimitv;i++)
     {
      delete[] Temp[i];Temp[i]=NULL;
     }
    } 
    // Deallocate arrays
    clean_quadrature(name_in,mode);
    for(i=0;i<nr;i++)
    {
     delete[] Intra_1[i];Intra_1[i]=NULL;
    }
    for(i=0;i<Read_fchk_wfn.nprimitv;i++)
    {
     delete[] CoefT[i];CoefT[i]=NULL;
     delete[] dm1_prim[i];dm1_prim[i]=NULL;
    }
    for(i=0;i<35;i++)
    {
     delete[] Quant[i];Quant[i]=NULL;
    }
    delete[] Intra_1;Intra_1=NULL;
    delete[] Quant;Quant=NULL;
    delete[] dm1;dm1=NULL;
    delete[] Coef;Coef=NULL;
    delete[] CoefT;CoefT=NULL;
    delete[] Temp;Temp=NULL;
    delete[] dm1_prim;dm1_prim=NULL;
    delete[] r_gauss;r_gauss=NULL;
    delete[] w_gauss;w_gauss=NULL;
    Results<<endl;
    if(!overlap){Results<<" See the files "<<name_basis<<" "<<name_dm1<<endl;}
    else{Results<<" See the files "<<name_in+"overlap"<<" "<<name_in+"overlap_mo"<<endl;}
    Results<<endl;
    Results<<"#*************************************************************************#";
    Results<<endl;
   }
  }	  
  ///////////////////////////////////////////////
  // Punctual evals and scans for WFN and FCHK //
  ///////////////////////////////////////////////
  if(Input_commands.punctualr)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#                       Punctual Evaluation (r)                           #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   for(i=0;i<Input_commands.punctuals_r;i++)
   {
    Point[0]=Input_commands.coordinates_r[i][0];
    Point[1]=Input_commands.coordinates_r[i][1];
    Point[2]=Input_commands.coordinates_r[i][2];
    Results<<"Coordinates of the point: \t";
    Results<<setprecision(10)<<fixed<<scientific;
    Results<<setw(17)<<Point[0]<<"\t"<<setw(17)<<Point[1]<<"\t"<<setw(17)<<Point[2]<<endl;
    punctual(Read_fchk_wfn,Point,Density,Density_alpha,Density_beta,Grad,Grad_alpha,Grad_beta,Grad_norm,
    Grad_norm_alpha,Grad_norm_beta,laplacian_r,laplacian_alpha,laplacian_beta,tauW_alpha,
    tauW_beta,tau_alpha,tau_beta,k_F_alpha,k_F_beta,k_s_alpha,k_s_beta,s_r_alpha,s_r_beta,
    q_red_alpha,q_red_beta,wfn_fchk);
    Results<<"Rho(r)             : "<<setw(17)<<Density<<endl;
    Results<<"Rho(r)  alpha      : "<<setw(17)<<Density_alpha<<endl;
    Results<<"Rho(r)  beta       : "<<setw(17)<<Density_beta<<endl;
    Results<<"|Grad[rho(r)]|     : "<<setw(17)<<Grad_norm<<endl;
    Results<<"Gradx[rho(r)]      : "<<setw(17)<<Grad[0]<<endl;
    Results<<"Grady[rho(r)]      : "<<setw(17)<<Grad[1]<<endl;
    Results<<"Gradz[rho(r)]      : "<<setw(17)<<Grad[2]<<endl;
    Results<<"|Grad[rho(r)]|alpha: "<<setw(17)<<Grad_norm_alpha<<endl;
    Results<<"Gradx[rho(r)] alpha: "<<setw(17)<<Grad_alpha[0]<<endl;
    Results<<"Grady[rho(r)] alpha: "<<setw(17)<<Grad_alpha[1]<<endl;
    Results<<"Gradz[rho(r)] alpha: "<<setw(17)<<Grad_alpha[2]<<endl;
    Results<<"|Grad[rho(r)]| beta: "<<setw(17)<<Grad_norm_beta<<endl;
    Results<<"Gradx[rho(r)]  beta: "<<setw(17)<<Grad_beta[0]<<endl;
    Results<<"Grady[rho(r)]  beta: "<<setw(17)<<Grad_beta[1]<<endl;
    Results<<"Gradz[rho(r)]  beta: "<<setw(17)<<Grad_beta[2]<<endl;
    Results<<"Laplacian[rho(r)]  : "<<setw(17)<<laplacian_r<<endl;
    Results<<"Laplacian[rho(r)] a: "<<setw(17)<<laplacian_alpha<<endl;
    Results<<"Laplacian[rho(r)] b: "<<setw(17)<<laplacian_beta<<endl;
    Results<<"tauW[rho(r)]       : "<<setw(17)<<tauW_alpha+tauW_beta<<endl;
    Results<<"tauW[rho(r)] alpha : "<<setw(17)<<tauW_alpha<<endl;
    Results<<"tauW[rho(r)] beta  : "<<setw(17)<<tauW_beta<<endl;
    Results<<"tauW is the Weizsacker kinetic energy"<<endl;
    Results<<"tau_total          : "<<setw(17)<<tau_alpha+tau_beta<<endl;
    Results<<"tau_alpha          : "<<setw(17)<<tau_alpha<<endl;
    Results<<"tau_beta           : "<<setw(17)<<tau_beta<<endl;
    Results<<"tau alpha or beta = 1/2 sum n_i | Grad phi_alpha or beta |^2 "<<endl;
    ID_IND_local(Read_fchk_wfn,Point,ID_alpha,ID_beta,IND_alpha,IND_beta);
    Results<<"ID_total           : "<<setw(17)<<ID_alpha+ID_beta<<endl;
    Results<<"ID_alpha           : "<<setw(17)<<ID_alpha<<endl;
    Results<<"ID_beta            : "<<setw(17)<<ID_beta<<endl;
    Results<<"ID alpha or beta = Corr. Indicator by Matito et al."<<endl;
    Results<<"IND_total          : "<<setw(17)<<IND_alpha+IND_beta<<endl;
    Results<<"IND_alpha          : "<<setw(17)<<IND_alpha<<endl;
    Results<<"IND_beta           : "<<setw(17)<<IND_beta<<endl;
    Results<<"IND alpha or beta = Corr. Indicator by Matito et al."<<endl;
    Results<<"k_F[rho(r)]  alpha : "<<setw(17)<<k_F_alpha<<endl;
    Results<<"k_F[rho(r)]  beta  : "<<setw(17)<<k_F_beta<<endl;
    Results<<"k_F[rho(r)] = (3 PI^2)^1/3 (rho(r))^1/3"<<endl;
    Results<<"k_s[rho(r)]  alpha : "<<setw(17)<<k_s_alpha<<endl;
    Results<<"k_s[rho(r)]  beta  : "<<setw(17)<<k_s_beta<<endl;
    Results<<"(k_s[rho(r)] = sqrt( (4*k_F[rho(r)])/(PI*a_o) )"<<endl;
    Results<<"a_o is the Bohr radius in a.u. is equal to 1)"<<endl;
    Results<<"s[rho(r)]  alpha   : "<<setw(17)<<s_r_alpha<<endl;
    Results<<"s[rho(r)]  beta    : "<<setw(17)<<s_r_beta<<endl;
    Results<<"s[rho(r)]   = |Grad [rho(r)]|/ (2*k_F*rho(r))"<<endl;
    Results<<"q_red[rho(r)] alpha: "<<setw(17)<<q_red_alpha<<endl;
    Results<<"q_red[rho(r)] beta : "<<setw(17)<<q_red_beta<<endl;
    Results<<"q_red[rho(r)] = Laplacian[rho(r) alpha or beta] / (k_F[rho(r)] alpha or beta k_F[rho(r)] alpha or beta rho(r) alpha or beta)"<<endl;
    Results<<endl;
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  if(Input_commands.punctualp)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#                       Punctual Evaluation (p)                           #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   for(i=0;i<Input_commands.punctuals_p;i++)
   {
    Point[0]=Input_commands.coordinates_p[i][0];
    Point[1]=Input_commands.coordinates_p[i][1];
    Point[2]=Input_commands.coordinates_p[i][2];
    Results<<"Coordinates of the point: \t";
    Results<<setprecision(10)<<fixed<<scientific;
    Results<<setw(17)<<Point[0]<<"\t"<<setw(17)<<Point[1]<<"\t"<<setw(17)<<Point[2]<<endl;
    Read_fchk_wfn.rho_p_eval(Point,Density);
    Read_fchk_wfn.rho_p_grad(Point,Grad);
    Read_fchk_wfn.rho_p_lapl(Point,laplacian_p);
    Read_fchk_wfn.rho_p_eval_a_b(Point,Density_alpha,Density_beta);
    Read_fchk_wfn.rho_p_grad_a_b(Point,Grad_alpha,Grad_beta);
    Read_fchk_wfn.rho_p_lapl_a_b(Point,laplacian_alpha,laplacian_beta);
    Results<<"Rho(p)             : "<<setw(17)<<Density<<endl;
    Results<<"Rho(p)  alpha      : "<<setw(17)<<Density_alpha<<endl;
    Results<<"Rho(p)  beta       : "<<setw(17)<<Density_beta<<endl;
    Results<<"|Grad[rho(p)]|     : "<<setw(17)<<norm3D(Grad)<<endl;
    Results<<"Gradx[rho(p)]      : "<<setw(17)<<Grad[0]<<endl;
    Results<<"Grady[rho(p)]      : "<<setw(17)<<Grad[1]<<endl;
    Results<<"Gradz[rho(p)]      : "<<setw(17)<<Grad[2]<<endl;
    Results<<"|Grad[rho(p)]|alpha: "<<setw(17)<<norm3D(Grad_alpha)<<endl;
    Results<<"Gradx[rho(p)] alpha: "<<setw(17)<<Grad_alpha[0]<<endl;
    Results<<"Grady[rho(p)] alpha: "<<setw(17)<<Grad_alpha[1]<<endl;
    Results<<"Gradz[rho(p)] alpha: "<<setw(17)<<Grad_alpha[2]<<endl;
    Results<<"|Grad[rho(p)]| beta: "<<setw(17)<<norm3D(Grad_beta)<<endl;
    Results<<"Gradx[rho(p)]  beta: "<<setw(17)<<Grad_beta[0]<<endl;
    Results<<"Grady[rho(p)]  beta: "<<setw(17)<<Grad_beta[1]<<endl;
    Results<<"Gradz[rho(p)]  beta: "<<setw(17)<<Grad_beta[2]<<endl;
    Results<<"Laplacian[rho(p)]  : "<<setw(17)<<laplacian_p<<endl;
    Results<<"Laplacian[rho(p)] a: "<<setw(17)<<laplacian_alpha<<endl;
    Results<<"Laplacian[rho(p)] b: "<<setw(17)<<laplacian_beta<<endl;
    Results<<endl;
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  if(Input_commands.scanr)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#Coord. of the point: \t\t\t\t\t  Density \t\t Grad \t\t Laplacian #";
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Point[0]=Input_commands.init_coord_r[0];
   Point[1]=Input_commands.init_coord_r[1];
   Point[2]=Input_commands.init_coord_r[2];
   Read_fchk_wfn.rho_eval(Point,Density);
   Read_fchk_wfn.rho_grad(Point,Grad);
   Read_fchk_wfn.rho_lapl(Point,laplacian_r);
   Grad_norm=norm3D(Grad);
   step=Input_commands.step_r;
   points=Input_commands.points_scan_r;
   direct=Input_commands.dir_r;
   for(i=0;i<points;i++)
   {
    Results<<setw(17)<<Point[0]<<" "<<setw(17)<<Point[1]<<" "<<setw(17)<<Point[2]<<" "<<setw(17)<<Density;
    Results<<" "<<setw(17)<<Grad_norm<<" "<<setw(17)<<laplacian_r<<endl;
    if(direct=='x' || direct=='X')
    {Point[0]=Point[0]+step;}
    else if(direct=='y' || direct=='Y')
    {Point[1]=Point[1]+step;}
    else
    {Point[2]=Point[2]+step;}
    Read_fchk_wfn.rho_eval(Point,Density);
    Read_fchk_wfn.rho_grad(Point,Grad);
    Read_fchk_wfn.rho_lapl(Point,laplacian_r);
    Grad_norm=norm3D(Grad);
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  if(Input_commands.scanp)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#Coord. of the moment: \t\t\t\t\t  Density(p) \t\t Grad(p)    Laplacian(p) #";
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Point[0]=Input_commands.init_coord_p[0];
   Point[1]=Input_commands.init_coord_p[1];
   Point[2]=Input_commands.init_coord_p[2];
   Read_fchk_wfn.rho_p_eval(Point,Density);
   Read_fchk_wfn.rho_p_grad(Point,Grad);
   Read_fchk_wfn.rho_p_lapl(Point,laplacian_p);
   step=Input_commands.step_p;
   points=Input_commands.points_scan_p;
   direct=Input_commands.dir_p;
   for(i=0;i<points;i++)
   {
    Grad_norm=norm3D(Grad);
    Results<<setw(17)<<Point[0]<<" "<<setw(17)<<Point[1]<<" "<<setw(17)<<Point[2]<<" "<<setw(17)<<Density;
    Results<<" "<<setw(17)<<Grad_norm<<" "<<setw(17)<<laplacian_p<<endl;
    if(direct=='x' || direct=='X')
    {Point[0]=Point[0]+step;}
    else if(direct=='y' || direct=='Y')
    {Point[1]=Point[1]+step;}
    else
    {Point[2]=Point[2]+step;}
    Read_fchk_wfn.rho_p_eval(Point,Density);
    Read_fchk_wfn.rho_p_grad(Point,Grad);
    Read_fchk_wfn.rho_p_lapl(Point,laplacian_p);
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  if(Input_commands.scanelf)
  {
   coef_elf=(THREE/FIVE)*pow(SIX*PI*PI,TWO/THREE);
   if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
   {
    not_elf=true;
    cout<<"Unable to compute the ELF from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
    cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
    cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
   }
   if(!not_elf)
   {
    Results<<"#*************************************************************************#";
    Results<<endl;
    Results<<"#Coord. of the point: \t\t\t\t  Density \t\t Grad \t\t ELF \t\t ELFa \t ELFb"<<endl;
    Results<<setprecision(8)<<fixed<<scientific;
    Point[0]=Input_commands.init_coord_elf[0];
    Point[1]=Input_commands.init_coord_elf[1];
    Point[2]=Input_commands.init_coord_elf[2];
    pre_elf(Read_fchk_wfn,Point,Density_alpha,Density_beta,tauW_alpha,tauW_beta,tau_alpha,tau_beta,tcurr_alpha,tcurr_beta);
    Density=Density_alpha+Density_beta;
    step=Input_commands.step_elf;
    points=Input_commands.points_scan_elf;
    direct=Input_commands.dir_elf;
    for(i=0;i<points;i++)
    {
     if(pow(Density,FIVE/THREE)>pow(TEN,-TEN))
     {
      ELF=TWO*(tau_alpha+tau_beta-tauW_alpha-tauW_beta-tcurr_alpha-tcurr_beta)/(coef_elf*pow(Density,FIVE/THREE));
      ELF=ONE/(ONE+ELF*ELF);
     }
     else
     {ELF=ZERO;}
     if(pow(Density_alpha,FIVE/THREE)>pow(TEN,-TEN))
     {
      ELFa=TWO*(tau_alpha-tauW_alpha-tcurr_alpha)/(coef_elf*pow(Density_alpha,FIVE/THREE));
      ELFa=ONE/(ONE+ELFa*ELFa);
     }
     else
     {ELFa=ZERO;}
     if(pow(Density_beta,FIVE/THREE)>pow(TEN,-TEN))
     {
      ELFb=TWO*(tau_beta-tauW_beta-tcurr_beta)/(coef_elf*pow(Density_beta,FIVE/THREE));
      ELFb=ONE/(ONE+ELFb*ELFb);
     }
     else
     {ELFb=ZERO;}
     Results<<setw(15)<<Point[0]<<" "<<setw(15)<<Point[1]<<" "<<setw(15)<<Point[2]<<" "<<setw(15)<<Density;
     Results<<" "<<setw(15)<<Grad_norm<<" "<<setw(15)<<ELF<<" "<<setw(15)<<ELFa<<" "<<setw(15)<<ELFb<<endl;
     if(direct=='x' || direct=='X')
     {Point[0]=Point[0]+step;}
     else if(direct=='y' || direct=='Y')
     {Point[1]=Point[1]+step;}
     else
     {Point[2]=Point[2]+step;}
     pre_elf(Read_fchk_wfn,Point,Density_alpha,Density_beta,tauW_alpha,tauW_beta,tau_alpha,tau_beta,tcurr_alpha,tcurr_beta);
     Density=Density_alpha+Density_beta;
    }
    Results<<"#[JCP,92,5397(1990)]"<<endl;
    Results<<"#[Notice that a represents alpha and b beta]"<<endl;
    Results<<"#*************************************************************************#";
    Results<<endl;
   }
  }
  if(Input_commands.scanindic)
  {
   if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
   {
    not_indic=true;
    cout<<"Unable to compute the INDICATORS from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
    cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
    cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
   }
   if(!not_indic)
   {
    Results<<"#*************************************************************************#";
    Results<<endl;
    Results<<"#Coord. of the point: \t\t  Density \t\t IDa \t\t IDb \t\t INDa \t INDb"<<endl;
    Results<<setprecision(8)<<fixed<<scientific;
    Point[0]=Input_commands.init_coord_indic[0];
    Point[1]=Input_commands.init_coord_indic[1];
    Point[2]=Input_commands.init_coord_indic[2];
    Read_fchk_wfn.rho_eval(Point,Density);
    ID_IND_local(Read_fchk_wfn,Point,ID_alpha,ID_beta,IND_alpha,IND_beta);
    step=Input_commands.step_indic;
    points=Input_commands.points_scan_indic;
    direct=Input_commands.dir_indic;
    for(i=0;i<points;i++)
    {
     Results<<setw(15)<<Point[0]<<" "<<setw(15)<<Point[1]<<" "<<setw(15)<<Point[2]<<" "<<setw(15)<<Density;
     Results<<" "<<setw(15)<<ID_alpha<<" "<<setw(15)<<ID_beta<<" "<<setw(15)<<IND_alpha<<" "<<setw(15)<<IND_beta<<endl;
     if(direct=='x' || direct=='X')
     {Point[0]=Point[0]+step;}
     else if(direct=='y' || direct=='Y')
     {Point[1]=Point[1]+step;}
     else
     {Point[2]=Point[2]+step;}
     Read_fchk_wfn.rho_eval(Point,Density);
     ID_IND_local(Read_fchk_wfn,Point,ID_alpha,ID_beta,IND_alpha,IND_beta);
    }
    Results<<"#[Notice that a represents alpha and b beta]"<<endl;
    Results<<"#*************************************************************************#";
    Results<<endl;
   }
  }
  if(Input_commands.scan_localhybs)
  {
   double tw_div_t,t2,q_red_total,s_red_total;
   if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
   {
    cout<<"Unable to compute tau(r) from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
    cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
    cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"# tauW(r) is the Weizsacker kinetic energy"<<endl;
   Results<<"# tau(r) = 1/2 sum n_i | Grad phi|^2 "<<endl;
   Results<<"# s[rho(r)]   = |Grad [rho(r)]|/ (2*k_F*rho(r))"<<endl;
   Results<<"# q_red[rho(r)] = Laplacian[rho(r)] / (4 k_F[rho(r)] k_F[rho(r)] rho(r))"<<endl;
   Results<<"# t^2[rho(r)] = (pi/3)^1/3 (a_o |Grad rho(r)|^2)/(16 rho^{7/3}(r))"<<endl;
   Results<<"# theta(x) = 4(1- 2 qred/ (s_red)^2 + (q_red)^2/(s_red)^4) -> DORI JCP, 142, 074112 (2015)"<<endl;
   Results<<"# Coord. of the point: \t t(w)/t(r) \t 1/[1+0.5 t^2(r)] \t 0.676/[1+8.908 theta(x)] \t [s(r)/(0.73+s(r))]^2 \t erf[0.20 s(r)]"<<endl;
   Results<<setprecision(8)<<fixed<<scientific;
   Point[0]=Input_commands.init_coord_lh[0];
   Point[1]=Input_commands.init_coord_lh[1];
   Point[2]=Input_commands.init_coord_lh[2];
   punctual(Read_fchk_wfn,Point,Density,Density_alpha,Density_beta,Grad,Grad_alpha,Grad_beta,Grad_norm,
   Grad_norm_alpha,Grad_norm_beta,laplacian_r,laplacian_alpha,laplacian_beta,tauW_alpha,
   tauW_beta,tau_alpha,tau_beta,k_F_alpha,k_F_beta,k_s_alpha,k_s_beta,s_r_alpha,s_r_beta,
   q_red_alpha,q_red_beta,wfn_fchk);
   step=Input_commands.step_lh;
   points=Input_commands.points_scan_lh;
   direct=Input_commands.dir_lh;
   for(i=0;i<points;i++)
   {
    tw_div_t=(tauW_alpha+tauW_beta)/(tau_alpha+tau_beta+pow(TEN,-TWO*TEN));
    t2=pow(PI/THREE,ONE_THIRD)*a_o*pow(Grad_norm,TWO)/(TWO*EIGHT*pow(Density,SEVEN*ONE_THIRD)+pow(TEN,-TWO*TEN));
    q_red_total=q_red_alpha+q_red_beta;
    s_red_total=s_r_alpha+s_r_beta;
    DORI=FOUR*(ONE-TWO*q_red_total/(s_red_total*s_red_total+pow(TEN,-TWO*TEN))+q_red_total*q_red_total/(pow(s_red_total,FOUR)+pow(TEN,-TWO*TEN)));
    // See  JCP, 142, 074112 (2015)
    //      JCP, 140, 18A510 (2014)
    //      JCP, 131, 154112 (2009)
    //      JPCA, 113, 11898 (2009)
    LOCAL_HYBRIDS_fr[0]=tw_div_t;
    LOCAL_HYBRIDS_fr[1]=ONE/(ONE+HALF*t2+pow(TEN,-TWO*TEN));
    LOCAL_HYBRIDS_fr[2]=0.676/(ONE+8.908*DORI);
    LOCAL_HYBRIDS_fr[3]=s_red_total*s_red_total/(pow(0.73+s_red_total,TWO)+pow(TEN,-TWO*TEN));
    LOCAL_HYBRIDS_fr[4]=erf(0.20*s_red_total);
    Results<<setw(15)<<Point[0]<<" "<<setw(15)<<Point[1]<<" "<<setw(15)<<Point[2]<<" ";
    for(j=0;j<5;j++)
    {
     Results<<setw(15)<<LOCAL_HYBRIDS_fr[j]<<" ";
    }
    Results<<endl;
    if(direct=='x' || direct=='X')
    {Point[0]=Point[0]+step;}
    else if(direct=='y' || direct=='Y')
    {Point[1]=Point[1]+step;}
    else
    {Point[2]=Point[2]+step;}
    punctual(Read_fchk_wfn,Point,Density,Density_alpha,Density_beta,Grad,Grad_alpha,Grad_beta,Grad_norm,
    Grad_norm_alpha,Grad_norm_beta,laplacian_r,laplacian_alpha,laplacian_beta,tauW_alpha,
    tauW_beta,tau_alpha,tau_beta,k_F_alpha,k_F_beta,k_s_alpha,k_s_beta,s_r_alpha,s_r_beta,
    q_red_alpha,q_red_beta,wfn_fchk);
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ////////////////////////////////////
  //  DMN punctual evaluations      //
  ////////////////////////////////////
  if(Input_commands.dmn && !wfn_fchk)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#          Read the DMN output file and evaluate n-Order RDM              #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   for(i=0;i<Input_commands.dmns;i++)
   {
    DMN_OPS dmn(Input_commands.name_dmn[i],Input_commands.dmn_order[i]);
    dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
    if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
    if(dmn.dm1)
    {
     Results<<"DM1"<<endl;
     Results<<"****"<<endl;
     Results<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     Results<<"Tr 1               : "<<setw(17)<<dmn.calc_tr()<<endl;
     Results<<"RHO(r1,r1')        : ";
     Results<<setw(17)<<dmn.evaluation(Input_commands.pointDM1,Input_commands.pointprimeDM1)<<endl;
     Results<<"r1                 : "<<endl;
     Results<<" ";
     if(Input_commands.pointDM1[0]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointDM1[0];
     Results<<" ";
     if(Input_commands.pointDM1[1]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointDM1[1];
     Results<<" ";
     if(Input_commands.pointDM1[2]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointDM1[2];
     Results<<endl;
     Results<<"r1'                : "<<endl;
     Results<<" ";
     if(Input_commands.pointprimeDM1[0]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointprimeDM1[0];
     Results<<" ";
     if(Input_commands.pointprimeDM1[1]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointprimeDM1[1];
     Results<<" ";
     if(Input_commands.pointprimeDM1[2]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointprimeDM1[2];
     Results<<endl;
     if(Input_commands.pointDM1[0]==Input_commands.pointprimeDM1[0] && Input_commands.pointDM1[1]==Input_commands.pointprimeDM1[1])
     {
      if(Input_commands.pointDM1[2]==Input_commands.pointprimeDM1[2])
      {
       dmn.grad_rho_r(Input_commands.pointDM1,Grad,Density);
       Results<<"Rho(r)             : "<<setw(17)<<Density<<endl;
       Results<<"|Grad RHO(r)|      : "<<setw(17)<<norm3D(Grad)<<endl;
       Results<<"Gradx RHO(r)       : "<<setw(17)<<Grad[0]<<endl;
       Results<<"Grady RHO(r)       : "<<setw(17)<<Grad[1]<<endl;
       Results<<"Gradz RHO(r)       : "<<setw(17)<<Grad[2]<<endl;
      }
     }
     Results<<endl;
     Results<<"NO occupations (spins summed):"<<endl;
     Results<<endl;
     dmn.diagonalize();
     counter=ZERO;
     for(j=0;j<Read_fchk_wfn.nbasisf;j++)
     {
      counter=counter+ONE;
      Results<<" ";
      if(abs(dmn.rho_matrix[j][j])>=pow(TEN,-SIX))
      {
       if(dmn.rho_matrix[j][j]>=ZERO) Results<<" ";
       Results<<dmn.rho_matrix[j][j];
      }
      else{Results<<" "<<ZERO;}
      if((counter/FIVE)-(int)(counter/FIVE) ==ZERO) Results<<endl;
     }
     Results<<endl;
     Results<<endl;
     dmn.diagonalize_ab();
     Results<<"NO occupations (alpha):"<<endl;
     Results<<endl;
     counter=ZERO;
     for(j=0;j<Read_fchk_wfn.nbasisf;j++)
     {
      counter=counter+ONE;
      Results<<" ";
      if(abs(dmn.rho_matrixa[j][j])>=pow(TEN,-SIX))
      {
       if(dmn.rho_matrixa[j][j]>=ZERO) Results<<" ";
       Results<<dmn.rho_matrixa[j][j];
      }
      else{Results<<" "<<ZERO;}
      if((counter/FIVE)-(int)(counter/FIVE) ==ZERO) Results<<endl;
     }
     Results<<endl;
     Results<<endl;
     Results<<"NO occupations (beta):"<<endl;
     Results<<endl;
     counter=0;
     for(j=0;j<Read_fchk_wfn.nbasisf;j++)
     {
      counter=counter+ONE;
      Results<<" ";
      if(abs(dmn.rho_matrixb[j][j])>=pow(TEN,-SIX))
      {
       if(dmn.rho_matrixb[j][j]>=ZERO) Results<<" ";
       Results<<dmn.rho_matrixb[j][j];
      }
      else{Results<<" "<<ZERO;}
      if((counter/FIVE)-(int)(counter/FIVE) ==ZERO) Results<<endl;
     }
     Results<<endl;
    }
    else if(dmn.dm2)
    {
     Results<<"DM2"<<endl;
     Results<<"****"<<endl;
     Results<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     Results<<"Tr 2               : "<<setw(17)<<dmn.calc_tr()<<endl;
     Results<<endl;
    }
    else if(dmn.dm3)
    {
     Results<<"DM3"<<endl;
     Results<<"****"<<endl;
     Results<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     Results<<"Tr 3                : "<<setw(17)<<dmn.calc_tr()<<endl;
     Results<<endl;
    }
    else
    {}
   }
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ////////////////////////////////////
  //  DMNp punctual evaluations     //
  ////////////////////////////////////
  if(Input_commands.dmnp && !wfn_fchk)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#          Read the DMN output file and evaluate n-Order RDMp             #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   for(i=0;i<Input_commands.dmnsp;i++)
   {
    DMN_P_OPS dmnp(Input_commands.name_dmnp[i],Input_commands.dmn_orderp[i]);
    dmnp.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
    if(Input_commands.dmn_thresh){dmnp.set_thershold(Input_commands.dmn_threshold);}
    if(dmnp.dm1)
    {
     Results<<"DM1p"<<endl;
     Results<<"****"<<endl;
     Results<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     Results<<"Tr 1               : "<<setw(17)<<dmnp.calc_p_tr()<<endl;
     Results<<"RHO(p1,p1')        : ";
     Results<<setw(17)<<dmnp.evaluation_p(Input_commands.pointDM1p,Input_commands.pointprimeDM1p)<<endl;
     Results<<"p1                 : "<<endl;
     Results<<" ";
     if(Input_commands.pointDM1p[0]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointDM1p[0];
     Results<<" ";
     if(Input_commands.pointDM1p[1]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointDM1p[1];
     Results<<" ";
     if(Input_commands.pointDM1p[2]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointDM1p[2];
     Results<<endl;
     Results<<"p1'                : "<<endl;
     Results<<" ";
     if(Input_commands.pointprimeDM1p[0]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointprimeDM1p[0];
     Results<<" ";
     if(Input_commands.pointprimeDM1p[1]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointprimeDM1p[1];
     Results<<" ";
     if(Input_commands.pointprimeDM1p[2]>=ZERO){Results<<" ";}
     Results<<Input_commands.pointprimeDM1p[2];
     Results<<endl;
     if(Input_commands.pointDM1p[0]==Input_commands.pointprimeDM1p[0] && Input_commands.pointDM1p[1]==Input_commands.pointprimeDM1p[1])
     {
      if(Input_commands.pointDM1p[2]==Input_commands.pointprimeDM1p[2])
      {
       dmnp.grad_rho_p(Input_commands.pointDM1p,Grad,Density);
       Results<<"Rho(p)             : "<<setw(17)<<Density<<endl;
       Results<<"|Grad RHO(p)|      : "<<setw(17)<<norm3D(Grad)<<endl;
       Results<<"Gradpx RHO(p)      : "<<setw(17)<<Grad[0]<<endl;
       Results<<"Gradpy RHO(p)      : "<<setw(17)<<Grad[1]<<endl;
       Results<<"Gradpz RHO(p)      : "<<setw(17)<<Grad[2]<<endl;
      }
     }
     Results<<endl;
     Results<<"NOp occupations (spins summed):"<<endl;
     Results<<endl;
     dmnp.diagonalize();
     counter=ZERO;
     for(j=0;j<Read_fchk_wfn.nbasisf;j++)
     {
      counter=counter+ONE;
      Results<<" ";
      if(abs(dmnp.rho_matrix[j][j])>=pow(TEN,-SIX))
      {
       if(dmnp.rho_matrix[j][j]>=ZERO) Results<<" ";
       Results<<dmnp.rho_matrix[j][j];
      }
      else{Results<<" "<<ZERO;}
      if((counter/FIVE)-(int)(counter/FIVE) ==ZERO) Results<<endl;
     }
     Results<<endl;
     Results<<endl;
     dmnp.diagonalize_ab();
     Results<<"NOp occupations (alpha):"<<endl;
     Results<<endl;
     counter=ZERO;
     for(j=0;j<Read_fchk_wfn.nbasisf;j++)
     {
      counter=counter+ONE;
      Results<<" ";
      if(abs(dmnp.rho_matrixa[j][j])>=pow(TEN,-SIX))
      {
       if(dmnp.rho_matrixa[j][j]>=ZERO) Results<<" ";
       Results<<dmnp.rho_matrixa[j][j];
      }
      else{Results<<" "<<ZERO;}
      if((counter/FIVE)-(int)(counter/FIVE) ==ZERO) Results<<endl;
     }
     Results<<endl;
     Results<<endl;
     Results<<"NOp occupations (beta):"<<endl;
     Results<<endl;
     counter=0;
     for(j=0;j<Read_fchk_wfn.nbasisf;j++)
     {
      counter=counter+ONE;
      Results<<" ";
      if(abs(dmnp.rho_matrixb[j][j])>=pow(TEN,-SIX))
      {
       if(dmnp.rho_matrixb[j][j]>=ZERO) Results<<" ";
       Results<<dmnp.rho_matrixb[j][j];
      }
      else{Results<<" "<<ZERO;}
      if((counter/FIVE)-(int)(counter/FIVE) ==ZERO) Results<<endl;
     }
     Results<<endl;
    }
    else if(dmnp.dm2)
    {
     Results<<"DM2p"<<endl;
     Results<<"****"<<endl;
     Results<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     Results<<"Tr 2               : "<<setw(17)<<dmnp.calc_p_tr()<<endl;
     Results<<endl;
    }
    else if(dmnp.dm3)
    {
     Results<<"DM3p"<<endl;
     Results<<"****"<<endl;
     Results<<endl;
     Results<<setprecision(10)<<fixed<<scientific;
     Results<<"Tr 3              : "<<setw(17)<<dmnp.calc_p_tr()<<endl;
     Results<<endl;
    }
    else
    {}
   }
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  /////////////////////////
  // Integration options //
  /////////////////////////
  if(Input_commands.integrals)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#                        Integral Options                                 #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   if(Input_commands.cuba) //Cuba
   {
    int counter_int,Tot_Ovrlps,MIN_EVALS,MAX_EVALS,fail;
    double **SIJ,rhoa=ZERO,rhob=ZERO,error_abs, error_rel,res_integration[20];
    double Im[3][3],Quad_rho[3][3];
    Density=ZERO;
    fail=10000;
    if(Input_commands.debug)
    {
     cout<<"The estimated error of Cuba is equal to the max(Eabs,Erel|Integration|)"<<endl;
     cout<<"Absolute error in the accuracy (Eabs) of the integration (e.g. 1e-4)"<<endl;
     cout<<"Relative error in the accuracy (Erel) of the integration (e.g. 1e-3)"<<endl;
    }
    error_abs=Input_commands.error_abs;
    error_rel=Input_commands.error_rel;
    MIN_EVALS=Input_commands.minevals;
    MAX_EVALS=Input_commands.maxevals;
    method=Input_commands.method_cuba;
    if(method=="vegas"){method="Vegas";}
    if(method=="suave"){method="Suave";}
    if(method=="divonne"){method="Divonne";}
    if(method!="Divonne" && method!="Vegas" && method!="Suave"){method="Cuhre";}
    for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
    {
     operation=Input_commands.integral_ops[num_integrals];
     if(operation=="density"){pos=true;}
     if(operation=="densityp"){mom=true;}
     if(operation=="inertia"){pos=true;inertiar=true;}
     if(operation=="quadrupole"){pos=true;inertiar=true;}
     if(operation=="inertiap"){mom=true;inertiap=true;}
     if(operation=="shannon"){shan=true;pos=true;}
     if(operation=="fisher"){fish=true;pos=true;}
     if(operation=="shannonp"){shanp=true;mom=true;}
     if(operation=="fisherp"){fishp=true;mom=true;}
     if(operation=="ttf"){pos=true;}
     if(operation=="r1"){pos=true;r1=true;}
     if(operation=="r2"){pos=true;r2=true;}
     if(operation=="rm1"){pos=true;rm1=true;}
     if(operation=="p1"){mom=true;p1=true;}
     if(operation=="p2"){mom=true;p2=true;}
     if(operation=="tw"){pos=true;fish=true;}
     if(operation=="dipolar"){pos=true;dipolar=true;}
     if(operation=="rho"){pos=true;rho=true;}
     if(operation=="sij"){sij=true;}
    }
    if((pos || mom) && sij)
    {
     cout<<"Evaluation of density rho(r) or pi(p) properties is not compatible with Sij ints creation."<<endl;
     cout<<"resubmit the calculation separating these evaluations"<<endl;
     mom=false;pos=false;sij=false;fish=false;shan=false;shanp=false;fishp=false;
    }
    for(i=0;i<6;i++)
    {
     Integrals_interval[i]=Input_commands.interval_integralsCUB[i];
     //Functions borrowed from Quadrature
     //Pass from -180 to 180 and 0 to 180 to rad and 0 to 360 and 0 to 180
     if(i==2 || i==3)
     {Integrals_interval[i]=theta_rad(Integrals_interval[i]);}
     if(i==4 || i==5)
     {Integrals_interval[i]=phi_rad(Integrals_interval[i]);}
    }
    if(pos)
    {
     if(Input_commands.dmn_integrals && !wfn_fchk)
     {
      if(rho)
      {
       cout<<"Warning, The < rho > integral is not implemented for FCHK + DM1 files!"<<endl;
       cout<<"Convert the DM1 file + FCHK file into a DM1_FCHK file"<<endl;
      }
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      integrate_dmnr(dmn,method,3,18,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,
      fail,Integrals_interval,Input_commands.nprocs,shan,fish,inertiar,r1,r2,dipolar);
     }
     else
     {
      integrate_cuba(Read_fchk_wfn,method,3,20,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,
      fail,Integrals_interval,Input_commands.nprocs,shan,fish,inertiar,r1,r2,rm1,dipolar,rho);
     }
     Results<<"The result of the integration  N  = "<<setw(17)<<res_integration[0]<<endl;
     Nelec=(double)Read_fchk_wfn.nelectrons;
     if(shan)
     {
      shannon=res_integration[1]/Nelec+log(Nelec);
      Results<<"The result of the integration S_r = "<<setw(17)<<shannon<<endl;
      Results<<"The Shannon entropy power     J_r = ";
      Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon);
      Results<<endl;
      Results<<"[J_r=(1/(2*pi*e)) * e^(2/3 *S_r)]"<<endl;
     }
     if(fish)
     {
      fisher=res_integration[2]/Nelec;
      Results<<"The result of the integration I_r = "<<setw(17)<<fisher<<endl;
      Tw=fisher/EIGHT;
      Results<<"The result of the integration T_w = "<<setw(17)<<Tw<<endl;
      Results<<"The result of the integration T_W = "<<setw(17)<<Nelec*Tw<<endl;
      Results<<"(T_w is per electron and T_W is the total)"<<endl;
     }
     Ttf=res_integration[3]*pow(THREE*PI*PI,TWO/THREE)*THREE/TEN;
     Results<<"The result of the integration T_TF= "<<setw(17)<<Ttf<<endl;
     if(inertiar)
     {
      if(abs(res_integration[0]-(int)res_integration[0])>=HALF)
      {Density=ceil(res_integration[0]);}
      else
      {Density=floor(res_integration[0]);}
      if(Density==ZERO){Density=ONE;}
      for(i=0;i<3;i++)
      {
       RCC[i]=res_integration[i+4]/Density;
       if(abs(RCC[i])<pow(TEN,-SIX)){RCC[i]=ZERO;}
      }
      Im[0][0]=pow(RCC[1],TWO)+pow(RCC[2],TWO);
      Im[0][1]=-RCC[0]*RCC[1];
      Im[0][2]=-RCC[0]*RCC[2];
      Im[1][0]=Im[0][1];
      Im[1][1]=pow(RCC[0],TWO)+pow(RCC[2],TWO);
      Im[1][2]=-RCC[1]*RCC[2];
      Im[2][0]=Im[0][2];
      Im[2][1]=Im[1][2];
      Im[2][2]=pow(RCC[0],TWO)+pow(RCC[1],TWO);
      Results<<endl;
      Results<<"Coordinates of the CC: \t";
      Results<<setprecision(10)<<fixed<<scientific;
      Results<<setw(17)<<RCC[0]<<" "<<setw(17)<<RCC[1]<<" "<<setw(17)<<RCC[2]<<endl;
      Results<<endl;
      Inertia=new double*[3];
      eigenV=new double*[3];
      Quadrupole=new double*[3];
      for(i=0;i<3;i++)
      {
       Inertia[i]=new double[3];
       eigenV[i]=new double[3];
       Quadrupole[i]=new double[3];
      }
      Read_fchk_wfn.quadrupoleATOMS(Quadrupole);
      counter_int=0;
      for(i=0;i<3;i++)
      {
       for(j=0;j<=i;j++)
       {
        Inertia[i][j]=res_integration[counter_int+7]/Density-Im[i][j];
        Inertia[j][i]=Inertia[i][j];
        Quad_rho[i][j]=res_integration[counter_int+7];
        if(i!=j)
        {
         Quad_rho[i][j]=-Quad_rho[i][j];
        }
        Quad_rho[j][i]=Quad_rho[i][j];
        counter_int++;
       }
      }
      Results<<"The inertia(r) tensor ( "<<setw(17)<<Inertia[0][0]<<" ";//see below this function
      Results<<setw(17)<<Inertia[0][1]<<" "<<setw(17)<<Inertia[0][2]<<" )"<<endl;
      Results<<"                      ( "<<setw(17)<<Inertia[1][0]<<" ";
      Results<<setw(17)<<Inertia[1][1]<<" "<<setw(17)<<Inertia[1][2]<<" )"<<endl;
      Results<<"                      ( "<<setw(17)<<Inertia[2][0]<<" ";
      Results<<setw(17)<<Inertia[2][1]<<" "<<setw(17)<<Inertia[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Diagonalized tensor:"<<endl;
      Results<<endl;
      i=3;
      jacobi(i,Inertia,eigenV); //Diagonalization
      for(i=0;i<3;i++)
      {
       for(j=0;j<3;j++)
       {
        if(abs(Inertia[i][j])<pow(TEN,-EIGHT) && i!=j){Inertia[i][j]=ZERO;}
        if(abs(eigenV[i][j])<pow(TEN,-SEVEN)){eigenV[i][j]=ZERO;}
       }
      }
      Results<<" ( "<<setw(17)<<Inertia[0][0]<<" "<<setw(17)<<Inertia[0][1]<<" ";
      Results<<setw(17)<<Inertia[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Inertia[1][0]<<" "<<setw(17)<<Inertia[1][1]<<" ";
      Results<<setw(17)<<Inertia[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Inertia[2][0]<<" "<<setw(17)<<Inertia[2][1]<<" ";
      Results<<setw(17)<<Inertia[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<endl;
      Results<<"Eigenvectors (columns):"<<endl;
      Results<<endl;
      Results<<" ( "<<setw(17)<<eigenV[0][0]<<" "<<setw(17)<<eigenV[0][1]<<" ";
      Results<<setw(17)<<eigenV[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<eigenV[1][0]<<" "<<setw(17)<<eigenV[1][1]<<" ";
      Results<<setw(17)<<eigenV[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<eigenV[2][0]<<" "<<setw(17)<<eigenV[2][1]<<" ";
      Results<<setw(17)<<eigenV[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Trace of I_diag [tr(I_diag)]      = "<<setw(17)<<Inertia[0][0]+Inertia[1][1]
      +Inertia[2][2]<<endl;
      Results<<endl;
      Quadrupole[0][0]=Quadrupole[0][0]-HALF*(Quad_rho[1][1]+Quad_rho[2][2]-Quad_rho[0][0]);
      Quadrupole[0][1]=Quadrupole[0][1]-Quad_rho[0][1];
      Quadrupole[0][2]=Quadrupole[0][2]-Quad_rho[0][2];
      Quadrupole[1][0]=Quadrupole[1][0]-Quad_rho[1][0];
      Quadrupole[1][1]=Quadrupole[1][1]-HALF*(Quad_rho[0][0]+Quad_rho[2][2]-Quad_rho[1][1]);
      Quadrupole[1][2]=Quadrupole[1][2]-Quad_rho[1][2];
      Quadrupole[2][0]=Quadrupole[2][0]-Quad_rho[2][0];
      Quadrupole[2][1]=Quadrupole[2][1]-Quad_rho[2][1];
      Quadrupole[2][2]=Quadrupole[2][2]-HALF*(Quad_rho[0][0]+Quad_rho[1][1]-Quad_rho[2][2]);
      for(i=0;i<3;i++)
      {
       for(j=0;j<3;j++)
       {
        Quadrupole[i][j]=AUtoDAng*Quadrupole[i][j];
        if(abs(Quadrupole[i][j])<pow(TEN,-SIX))
        {Quadrupole[i][j]=ZERO;}
       }
      }
      Results<<"Primitive Quadrupole Moment (Debye-Angstrom):"<<endl;
      Results<<endl;
      Results<<" ( "<<setw(17)<<Quadrupole[0][0]<<" "<<setw(17)<<Quadrupole[0][1]<<" "<<setw(17)<<Quadrupole[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Quadrupole[1][0]<<" "<<setw(17)<<Quadrupole[1][1]<<" "<<setw(17)<<Quadrupole[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Quadrupole[2][0]<<" "<<setw(17)<<Quadrupole[2][1]<<" "<<setw(17)<<Quadrupole[2][2]<<" )"<<endl;
      Results<<endl;
      Quad_rho[0][0]=HALF*(TWO*Quadrupole[0][0]-Quadrupole[1][1]-Quadrupole[2][2]);
      Quad_rho[1][1]=HALF*(TWO*Quadrupole[1][1]-Quadrupole[0][0]-Quadrupole[2][2]);
      Quad_rho[2][2]=HALF*(TWO*Quadrupole[2][2]-Quadrupole[1][1]-Quadrupole[0][0]);
      Quad_rho[0][1]=THREE*HALF*Quadrupole[0][1];Quad_rho[1][0]=THREE*HALF*Quadrupole[1][0];
      Quad_rho[0][2]=THREE*HALF*Quadrupole[0][2];Quad_rho[2][0]=THREE*HALF*Quadrupole[2][0];
      Quad_rho[1][2]=THREE*HALF*Quadrupole[1][2];Quad_rho[2][1]=THREE*HALF*Quadrupole[2][1];
      Results<<"Traceless Quadrupole Moment (Debye-Angstrom) [Gamess, Molpro]:"<<endl;
      Results<<endl;
      Results<<" ( "<<setw(17)<<Quad_rho[0][0]<<" "<<setw(17)<<Quad_rho[0][1]<<" "<<setw(17)<<Quad_rho[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Quad_rho[1][0]<<" "<<setw(17)<<Quad_rho[1][1]<<" "<<setw(17)<<Quad_rho[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Quad_rho[2][0]<<" "<<setw(17)<<Quad_rho[2][1]<<" "<<setw(17)<<Quad_rho[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Trace of Quadrupole               = "<<Quad_rho[0][0]+Quad_rho[1][1]+Quad_rho[2][2]<<endl;
      Results<<endl;
      Quad_rho[0][0]=(Quadrupole[0][0]-(Quadrupole[0][0]+Quadrupole[1][1]+Quadrupole[2][2])/THREE);
      Quad_rho[1][1]=(Quadrupole[1][1]-(Quadrupole[0][0]+Quadrupole[1][1]+Quadrupole[2][2])/THREE);
      Quad_rho[2][2]=(Quadrupole[2][2]-(Quadrupole[0][0]+Quadrupole[1][1]+Quadrupole[2][2])/THREE);
      Quad_rho[0][1]=Quadrupole[0][1];Quad_rho[1][0]=Quadrupole[1][0];
      Quad_rho[0][2]=Quadrupole[0][2];Quad_rho[2][0]=Quadrupole[2][0];
      Quad_rho[1][2]=Quadrupole[1][2];Quad_rho[2][1]=Quadrupole[2][1];
      Results<<"Traceless Quadrupole Moment (Debye-Angstrom) [Gaussian]:"<<endl;
      Results<<endl;
      Results<<" ( "<<setw(17)<<Quad_rho[0][0]<<" "<<setw(17)<<Quad_rho[0][1]<<" "<<setw(17)<<Quad_rho[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Quad_rho[1][0]<<" "<<setw(17)<<Quad_rho[1][1]<<" "<<setw(17)<<Quad_rho[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Quad_rho[2][0]<<" "<<setw(17)<<Quad_rho[2][1]<<" "<<setw(17)<<Quad_rho[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Trace of Quadrupole               = "<<Quad_rho[0][0]+Quad_rho[1][1]+Quad_rho[2][2]<<endl;
      Results<<endl;
      for(i=0;i<3;i++)
      {
       delete[] Inertia[i];Inertia[i]=NULL;
       delete[] eigenV[i];eigenV[i]=NULL;
       delete[] Quadrupole[i];Quadrupole[i]=NULL;
      }
      delete[] Inertia;Inertia=NULL;
      delete[] eigenV;eigenV=NULL;
      delete[] Quadrupole;Quadrupole=NULL;
     }
     if(r1){Results<<"The result of the integration <r> = "<<setw(17)<<res_integration[13]<<endl;}
     if(r2){Results<<"The result of the integration <r2>= "<<setw(17)<<res_integration[14]<<endl;}
     if(rm1){Results<<"The result of the integ.     <r-1>= "<<setw(17)<<res_integration[19]<<endl;}
     if(dipolar)
     {
      Read_fchk_wfn.muATOMS(mu);
      Results<<"Atomic contribution to mux        = "<<setw(17)<<mu[0]<<endl;
      Results<<"Atomic contribution to muy        = "<<setw(17)<<mu[1]<<endl;
      Results<<"Atomic contribution to muz        = "<<setw(17)<<mu[2]<<endl;
      mu[0]=-res_integration[15]+mu[0];mu[1]=-res_integration[16]+mu[1];mu[2]=-res_integration[17]+mu[2];
      Results<<"The result of the integration mux = "<<setw(17)<<mu[0]<<endl;
      Results<<"The result of the integration muy = "<<setw(17)<<mu[1]<<endl;
      Results<<"The result of the integration muz = "<<setw(17)<<mu[2]<<endl;
      Results<<"The norm of the dipolar mom. |mu| = "<<setw(17)<<norm3D(mu)<<endl;
     }
     if(rho && !Input_commands.dmn_integrals)
     {
      Results<<"The result of < rho >             = "<<setw(17)<<res_integration[18]<<endl;
     }
     Results<<endl;
     Results<<"\t \t Obtained with "<<method<<", fail report "<<fail<<endl;
     Results<<endl;
    }
    if(mom)
    {
     if(Input_commands.dmn_integrals && !wfn_fchk)
     {
      DMN_P_OPS dmnp(Input_commands.name_dm1,1);
      dmnp.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmnp.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmnp.set_thershold(Input_commands.dmn_threshold);}
      integrate_dmnp(dmnp,method,3,18,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,
      fail,Integrals_interval,Input_commands.nprocs,shanp,fishp,inertiap,p1,p2);
     }
     else
     {
      integrate_cubap(Read_fchk_wfn,method,3,20,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,
      fail,Integrals_interval,Input_commands.nprocs,shanp,fishp,inertiap,p1,p2);
     }
     Results<<"The result of the integration  Np = "<<setw(17)<<res_integration[0]<<endl;
     Nelec=(double)Read_fchk_wfn.nelectrons;
     if(shanp)
     {
      shannonp=res_integration[1]/Nelec+log(Nelec);
      Results<<"The result of the integration S_p = "<<setw(17)<<shannonp<<endl;
      Results<<"The Shannon entropy power     J_p = ";
      Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp);
      Results<<endl;
      Results<<"[J_p=(1/(2*pi*e)) * e^(2/3 *S_p)]"<<endl;
     }
     if(fishp)
     {
      fisherp=res_integration[2]/Nelec;
      Results<<"The result of the integration I_p = "<<setw(17)<<fisherp<<endl;
     }
     if(inertiap)
     {
      if(abs(res_integration[0]-(int)res_integration[0])>=HALF)
      {Density=ceil(res_integration[0]);}
      else
      {Density=floor(res_integration[0]);}
      if(Density==ZERO){Density=ONE;}
      for(i=0;i<3;i++){RCC[i]=res_integration[i+4]/Density;}
      Im[0][0]=pow(RCC[1],TWO)+pow(RCC[2],TWO);
      Im[0][1]=-RCC[0]*RCC[1];
      Im[0][2]=-RCC[0]*RCC[2];
      Im[1][0]=Im[0][1];
      Im[1][1]=pow(RCC[0],TWO)+pow(RCC[2],TWO);
      Im[1][2]=-RCC[1]*RCC[2];
      Im[2][0]=Im[0][2];
      Im[2][1]=Im[1][2];
      Im[2][2]=pow(RCC[0],TWO)+pow(RCC[1],TWO);
      for(i=0;i<3;i++){if(RCC[i]<pow(TEN,-SIX)){RCC[i]=ZERO;}}
      Results<<endl;
      Results<<"Coordinates of the CCp: ";
      Results<<setprecision(10)<<fixed<<scientific;
      Results<<setw(17)<<RCC[0]<<" "<<setw(17)<<RCC[1]<<" "<<setw(17)<<RCC[2]<<endl;
      Results<<endl;
      Inertia=new double*[3];
      for(i=0;i<3;i++){Inertia[i]=new double[3];}
      eigenV=new double*[3];
      for(i=0;i<3;i++){eigenV[i]=new double[3];}
      counter_int=0;
      for(i=0;i<3;i++)
      {
       for(j=0;j<=i;j++)
       {
        Inertia[i][j]=res_integration[counter_int+7]/Density-Im[i][j];
        Inertia[j][i]=Inertia[i][j];
        counter_int++;
       }
      }
      Results<<"The inertia(p) tensor ( "<<setw(17)<<Inertia[0][0]<<" ";//see below this function
      Results<<setw(17)<<Inertia[0][1]<<" "<<setw(17)<<Inertia[0][2]<<" )"<<endl;
      Results<<"                      ( "<<setw(17)<<Inertia[1][0]<<" ";
      Results<<setw(17)<<Inertia[1][1]<<" "<<setw(17)<<Inertia[1][2]<<" )"<<endl;
      Results<<"                      ( "<<setw(17)<<Inertia[2][0]<<" ";
      Results<<setw(17)<<Inertia[2][1]<<" "<<setw(17)<<Inertia[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Diagonalized tensor:"<<endl;
      Results<<endl;
      i=3;
      jacobi(i,Inertia,eigenV); //Diagonalization
      for(i=0;i<3;i++)
      {
       for(j=0;j<3;j++)
       {
        if(abs(Inertia[i][j])<pow(TEN,-EIGHT) && i!=j){Inertia[i][j]=ZERO;}
        if(abs(eigenV[i][j])<pow(TEN,-SEVEN)){eigenV[i][j]=ZERO;}
       }
      }
      Results<<" ( "<<setw(17)<<Inertia[0][0]<<" "<<setw(17)<<Inertia[0][1]<<" ";
      Results<<setw(17)<<Inertia[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Inertia[1][0]<<" "<<setw(17)<<Inertia[1][1]<<" ";
      Results<<setw(17)<<Inertia[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Inertia[2][0]<<" "<<setw(17)<<Inertia[2][1]<<" ";
      Results<<setw(17)<<Inertia[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<endl;
      Results<<"Eigenvectors (columns):"<<endl;
      Results<<endl;
      Results<<" ( "<<setw(17)<<eigenV[0][0]<<" "<<setw(17)<<eigenV[0][1]<<" ";
      Results<<setw(17)<<eigenV[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<eigenV[1][0]<<" "<<setw(17)<<eigenV[1][1]<<" ";
      Results<<setw(17)<<eigenV[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<eigenV[2][0]<<" "<<setw(17)<<eigenV[2][1]<<" ";
      Results<<setw(17)<<eigenV[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Trace of I_diag [tr(I_diag)]      = "<<setw(17)<<Inertia[0][0]+Inertia[1][1]
      +Inertia[2][2]<<endl;
      Results<<endl;
      for(i=0;i<3;i++){delete[] Inertia[i];} delete[] Inertia;
      for(i=0;i<3;i++){delete[] eigenV[i];} delete[] eigenV;
     }
     if(p1){Results<<"The result of the integration <p> = "<<setw(17)<<res_integration[13]<<endl;}
     if(p2)
     {
      Results<<"The result of the integration <p2>= "<<setw(17)<<res_integration[14]<<endl;
      Results<<"The kinetic energy [from PI(p)] T = "<<setw(17)<<res_integration[14]/TWO<<endl;
     }
     Results<<endl;
     Results<<"\t \t Obtained with "<<method<<", fail report "<<fail<<endl;
     Results<<endl;
    }
    if(fish==true && shan==true)
    {
     Results<<"The Fisher-Shannon product    P_r = ";
     Results<<setw(17)<<ONE/THREE*fisher*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon)<<endl;
     Results<<"[P_r = (1/3)*I_r*J_r]"<<endl;
    }
    if(fishp==true && shanp==true)
    {
     Results<<"The Fisher-Shannon product    P_p = ";
     Results<<setw(17)<<ONE/THREE*fisherp*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp)<<endl;
     Results<<"[P_p = (1/3)*I_p*J_p]"<<endl;
    }
    if(shanp==true && shan==true)
    {
     Results<<"The Shannon-Shannonp  (S_r + S_p) = ";
     Results<<setw(17)<<shannon+shannonp<<endl;
     Results<<"[S_Total = S_r + S_p]"<<endl;
     Results<<"[S_r + S_p >= 3(1+ln pi) ~ 6.434189658]"<<endl;
     Results<<"The Shannon-Shannonp  (J_r x J_p) = ";
     Results<<setw(17)<<(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp))*(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon))<<endl;
     Results<<"[J_r*J_p >= 1/4]"<<endl;
    }
    if(fish==true && fishp==true)
    {
     Results<<"The Fisher-Fisherp    (I_r x I_p) = ";
     Results<<setw(17)<<fisher*fisherp<<endl;
     Results<<"[I_r*I_p >= 36]"<<endl;
    }
    if(sij)
    {
     if(Read_fchk_wfn.wfn && Input_commands.MOorNO!="no")
     {
      cout<<"Warning! wfn files only contain MOs for DFT and HF calcs! (NOs are used in wfn files)"<<endl;
      Input_commands.MOorNO="no";
     }
     if(!Read_fchk_wfn.wfn && Input_commands.MOorNO=="no")
     {
      if(!Input_commands.dmn_integrals)
      {
       if(!Read_fchk_wfn.overlap)
       {
        cout<<"The PS matrix is not available. The Sij matrix will be computed in MOs"<<endl;
        Input_commands.MOorNO="mo";
       }
      }
      else
      {
       mos_to_nos_dmn=true;
       Input_commands.MOorNO="mo";
      }
     }
     Results<<endl;
     Results<<"List of int files generated:"<<endl;
     Results<<endl;
     nbasis=Read_fchk_wfn.nbasis();
     SIJ=new double*[nbasis];
     for(i=0;i<nbasis;i++)
     {SIJ[i]=new double[i+1];}
     if(Read_fchk_wfn.wfn)
     {
      // N(N+1)/2 =SUM_i ^N i
      Tot_Ovrlps=(nbasis*nbasis+nbasis)/2;
      if(Input_commands.debug && Tot_Ovrlps>10000)
      {cout<<"Posibly the number of NOs pairs "<<Tot_Ovrlps<<" is greater than 10000 and Cuba vectorial evals is not allowed";}
     }
     else
     {Tot_Ovrlps=Read_fchk_wfn.nbasis();}
     for(i=0;i<Input_commands.nregions;i++)
     {
      ostringstream convert;
      convert<<region_int;
      region_string=convert.str();
      region_int++;
      for(j=0;j<6;j++)
      {
       Integrals_interval[j]=Input_commands.interval_integrals[i][j];
       //Functions borrowed from Quadrature
       //Pass from -180 to 180 and 0 to 180 to rad and 0 to 360 and 0 to 180
       if(j==2 || j==3)
       {Integrals_interval[j]=theta_rad(Integrals_interval[j]);}
       if(j==4 || j==5)
       {Integrals_interval[j]=phi_rad(Integrals_interval[j]);}
      }
      integrate_cuba_sij(Read_fchk_wfn,method,3,Tot_Ovrlps,error_rel,error_abs,MIN_EVALS,MAX_EVALS,SIJ,
      fail,Integrals_interval,Input_commands.nprocs,Input_commands.MOorNO);
      if(mos_to_nos_dmn)
      {
       mos_to_nos_dmn_sij(SIJ,Read_fchk_wfn,Input_commands);
      }
      print_int(Read_fchk_wfn,name_file,SIJ,Read_fchk_wfn.nbasis(),Density,rhoa,rhob,region_string);
      if(Read_fchk_wfn.wfn)
      {
       Results<<name_file.substr(0,(name_file.length()-4))+"_X"+region_string+".int"<<endl;
      }
      else
      {
       Results<<name_file.substr(0,(name_file.length()-5))+"_X"+region_string+".int"<<endl;
      }
     }
     for(i=0;i<nbasis;i++)
     {delete[] SIJ[i];SIJ[i]=NULL;}
     delete[] SIJ;
     SIJ=NULL;
     Results<<endl;
     Results<<"\t \t Obtained with "<<method<<", fail report "<<fail<<endl;
     Results<<endl;
    }
    Results<<endl;
    Results<<"(fail 0 means finish without problem, -1 error, >1 accuracy";
    Results<<" not achived and 10000 not calculated)"<<endl;
   }
   else if(Input_commands.cubature) //Cubature
   {
    double result_integration;
    Results<<setprecision(10)<<fixed<<scientific;
    double error,error_abs;
    if(Input_commands.debug)
    {
     cout<<"Absolute error in the accuracy of the integration (e.g. 1e-4)"<<endl;
    }
    error_abs=Input_commands.error_abs;
    method="Cubature";
    for(i=0;i<6;i++)
    {
     Integrals_interval[i]=Input_commands.interval_integralsCUB[i];
     //Functions borrowed from Quadrature
     //Pass from -180 to 180 and 0 to 180 to rad and 0 to 360 and 0 to 180
     if(i==2 || i==3)
     {Integrals_interval[i]=theta_rad(Integrals_interval[i]);}
     if(i==4 || i==5)
     {Integrals_interval[i]=phi_rad(Integrals_interval[i]);}
    }
    for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
    {
     error=pow(TEN,FOUR);
     operation=Input_commands.integral_ops[num_integrals];
     if(operation=="density")
     {
      operation="Density";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnr(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of the integration  N  = "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="rho")
     {
      operation="rho";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       cout<<"Warning, The < rho > integral is not implemented for FCHK + DM1 files!"<<endl;
       cout<<"Convert the DM1 file + FCHK file into a DM1_FCHK file"<<endl;
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of < rho >             = "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="densityp")
     {
      operation="Densityp";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of the integration  Np = "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="r1")
     {
      operation="R1";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of the integration <r> = "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="r2")
     {
      operation="R2";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of the integration <r2>= "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="p1")
     {
      operation="P1";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of the integration <p> = "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="p2")
     {
      operation="P2";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Results<<"The result of the integration <p2>= "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
      Results<<"The kinetic energy [from PI(p)] T = "<<setw(17)<<result_integration/TWO<<endl;
     }
     else if(operation=="shannon")
     {
      operation="Shannon";
      shan=true;
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnr(dmn,operation,error_abs,result_integration,error,Integrals_interval);
       shannon=result_integration;
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
       Nelec=(double)Read_fchk_wfn.nelectrons;
       shannon=result_integration/Nelec+log(Nelec);
      }
      Results<<"The result of the integration S_r = "<<setw(17)<<shannon;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
      Results<<"The Shannon entropy power     J_r = ";
      Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon);
      Results<<endl;
      Results<<"[J_r=(1/(2*pi*e)) * e^(2/3 *S_r)]"<<endl;
     }
     else if(operation=="shannonp")
     {
      operation="Shannonp";
      shanp=true;
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
       shannonp=result_integration;
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
       Nelec=(double)Read_fchk_wfn.nelectrons;
       shannonp=result_integration/Nelec+log(Nelec);
      }
      Results<<"The result of the integration S_p = "<<setw(17)<<shannonp;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
      Results<<"The Shannon entropy power     J_p = ";
      Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp);
      Results<<endl;
      Results<<"[J_p=(1/(2*pi*e)) * e^(2/3 *S_p)]"<<endl;
     }
     else if(operation=="fisher")
     {
      operation="Fisher";
      fish=true;
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnr(dmn,operation,error_abs,result_integration,error,Integrals_interval);
       fisher=result_integration;
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
       Nelec=(double)Read_fchk_wfn.nelectrons;
       fisher=result_integration/Nelec;
      }
      Results<<"The result of the integration I_r = "<<setw(17)<<fisher;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error/Nelec<<endl;
     }
     else if(operation=="fisherp")
     {
      operation="Fisherp";
      fishp=true;
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_P_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnp(dmn,operation,error_abs,result_integration,error,Integrals_interval);
       fisherp=result_integration;
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
       Nelec=(double)Read_fchk_wfn.nelectrons;
       fisherp=result_integration/Nelec;
      }
      Results<<"The result of the integration I_p = "<<setw(17)<<fisherp;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error/Nelec<<endl;
     }
     else if(operation=="tw")
     {
      operation="Fisher";//Since Tw[rho]=1/8  Fisher[rho]
      cout<<"(Fisher integral will be calculated since Tw[rho]=1/8 Fisher[rho])"<<endl;
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      I_cubature_dmnr(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      Tw=result_integration/EIGHT;
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
       Nelec=(double)Read_fchk_wfn.nelectrons;
       Tw=result_integration/(Nelec*EIGHT);
      }
      Results<<"The result of the integration T_w = "<<setw(17)<<Tw;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error/Nelec<<endl;
      Results<<"The result of the integration T_W = "<<setw(17)<<Nelec*Tw;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
      Results<<"(T_w is per electron and T_W is the total)"<<endl;
     }
     else if(operation=="ttf")
     {
      operation="T_TF";
      if(Input_commands.dmn_integrals && !wfn_fchk)
      {
       DMN_OPS dmn(Input_commands.name_dm1,1);
       dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
       dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
       if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
       I_cubature_dmnr(dmn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      else
      {
       I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      }
      Ttf=result_integration*pow(THREE*PI*PI,TWO/THREE)*THREE/TEN;
      Results<<"The result of the integration T_TF= "<<setw(17)<<Ttf;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
     }
     else if(operation=="sij")
     {
      if(Input_commands.MOorNO=="no"){operation="Sij_NO";}
      else{operation="Sij_MO";}
      if((!Read_fchk_wfn.overlap && !Read_fchk_wfn.wfn) && operation=="Sij_NO")
      {
       cout<<"The PS matrix is not available. The Sij matrix will be computed in MOs "<<endl;
       operation="Sij_MO";
      }
      Read_fchk_wfn.Pair[0]=Input_commands.orb1-1;
      Read_fchk_wfn.Pair[1]=Input_commands.orb2-1;
      I_cubature(Read_fchk_wfn,operation,error_abs,result_integration,error,Integrals_interval);
      Results<<"The result of the integration Sij = "<<setw(17)<<result_integration;
      Results<<"\t obtained with "<<method<<", +/- "<<setw(17)<<error<<endl;
      if(operation=="Sij_NO")
      {
       Results<<"for NOs "<<Input_commands.orb1<<" "<<Input_commands.orb2<<endl;
      }
      else
      {
       Results<<"for MOs "<<Input_commands.orb1<<" "<<Input_commands.orb2<<endl;
      }
      Results<<endl;
      Results<<"\t obtained with "<<method<<endl;
      Results<<endl;
     }
     else{}
    }
    if(fish==true && shan==true)
    {
     Results<<"The Fisher-Shannon product    P_r = ";
     Results<<setw(17)<<ONE/THREE*fisher*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon)<<endl;
     Results<<"[P_r = (1/3)*I_r*J_r]"<<endl;
    }
    if(fishp==true && shanp==true)
    {
     Results<<"The Fisher-Shannon product    P_p = ";
     Results<<setw(17)<<ONE/THREE*fisherp*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp)<<endl;
     Results<<"[P_p = (1/3)*I_p*J_p]"<<endl;
    }
    if(shanp==true && shan==true)
    {
     Results<<"The Shannon-Shannonp  (S_r + S_p) = ";
     Results<<setw(17)<<shannon+shannonp<<endl;
     Results<<"[S_Total = S_r + S_p]"<<endl;
     Results<<"[S_r + S_p >= 3(1+ln pi) ~ 6.434189658]"<<endl;
     Results<<"The Shannon-Shannonp  (J_r x J_p) = ";
     Results<<setw(17)<<(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp))*(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon))<<endl;
     Results<<"[J_r*J_p >= 1/4]"<<endl;
    }
    if(fish==true && fishp==true)
    {
     Results<<"The Fisher-Fisherp    (I_r x I_p) = ";
     Results<<setw(17)<<fisher*fisherp<<endl;
     Results<<"[I_r*I_p >= 36]"<<endl;
    }
    Results<<"+/- 10000 means that the quantity was not calculated"<<endl;
   }
   else if(Input_commands.cubature2) //Cubature vectorial for DMN and parallel for WFN & FCHK.
   {
    double res_integration[10];
    Results<<setprecision(10)<<fixed<<scientific;
    double error[10],error_abs;
    if(Input_commands.debug)
    {
     cout<<"Absolute error in the accuracy of the integration (e.g. 1e-4)"<<endl;
    }
    error_abs=Input_commands.error_abs;
    method="Cubature 2";
    for(i=0;i<9;i++)
    {
     Integrals_interval[i]=Input_commands.interval_integralsCUB[i];
     //Functions borrowed from Quadrature
     //Pass from -180 to 180 and 0 to 180 to rad and 0 to 360 and 0 to 180
     if(i==2 || i==3)
     {Integrals_interval[i]=theta_rad(Integrals_interval[i]);}
     if(i==4 || i==5)
     {Integrals_interval[i]=phi_rad(Integrals_interval[i]);}
    }
    for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
    {
     operation=Input_commands.integral_ops[num_integrals];
     if(operation=="density"){pos=true;}
     if(operation=="densityp"){mom=true;}
     if(operation=="shannon"){shan=true;pos=true;}
     if(operation=="fisher"){fish=true;pos=true;}
     if(operation=="shannonp"){shanp=true;mom=true;}
     if(operation=="fisherp"){fishp=true;mom=true;}
     if(operation=="ttf"){pos=true;}
     if(operation=="r1"){pos=true;r1=true;}
     if(operation=="r2"){pos=true;r2=true;}
     if(operation=="p1"){mom=true;p1=true;}
     if(operation=="p2"){mom=true;p2=true;}
     if(operation=="tw"){pos=true;fish=true;}
     if(operation=="dipolar"){pos=true;dipolar=true;}
     if(operation=="rho"){pos=true;rho=true;}
     if(operation=="sij"){cout<<"Use Cuba, Quadrature or Cubature for Sij calculations"<<endl;}
    }
    for(i=0;i<10;i++){error[i]=pow(TEN,FOUR);}
    if(pos)
    {
     if(Input_commands.dmn_integrals && !wfn_fchk)
     {
      if(rho)
      {
       cout<<"Warning, The < rho > integral is not implemented for FCHK + DM1 files!"<<endl;
       cout<<"Convert the DM1 file + FCHK file into a DM1_FCHK file"<<endl;
      }
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      I_cubature2_dmnr(dmn,error_abs,res_integration,error,Integrals_interval,shan,fish,r1,r2,dipolar);
     }
     else
     {
      I_cubature2(Read_fchk_wfn,error_abs,res_integration,error,Integrals_interval,shan,fish,r1,r2,dipolar,rho);
     }
     Results<<"The result of the integration  N  = "<<setw(17)<<res_integration[0]<<endl;
     Nelec=(double)Read_fchk_wfn.nelectrons;
     if(shan)
     {
      shannon=res_integration[1]/Nelec+log(Nelec);
      Results<<"The result of the integration S_r = "<<setw(17)<<shannon<<endl;
      Results<<"The Shannon entropy power     J_r = ";
      Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon);
      Results<<endl;
      Results<<"[J_r=(1/(2*pi*e)) * e^(2/3 *S_r)]"<<endl;
     }
     if(fish)
     {
      fisher=res_integration[2]/Nelec;
      Results<<"The result of the integration I_r = "<<setw(17)<<fisher<<endl;
      Tw=fisher/EIGHT;
      Results<<"The result of the integration T_w = "<<setw(17)<<Tw<<endl;
      Results<<"The result of the integration T_W = "<<setw(17)<<Nelec*Tw<<endl;
      Results<<"(T_w is per electron and T_W is the total)"<<endl;
     }
     Ttf=res_integration[3]*pow(THREE*PI*PI,TWO/THREE)*THREE/TEN;
     Results<<"The result of the integration T_TF= "<<setw(17)<<Ttf<<endl;
     if(r1){Results<<"The result of the integration <r> = "<<setw(17)<<res_integration[4]<<endl;}
     if(r2){Results<<"The result of the integration <r2>= "<<setw(17)<<res_integration[5]<<endl;}
     if(dipolar)
     {
      Read_fchk_wfn.muATOMS(mu);
      mu[0]=-res_integration[6]+mu[0];mu[1]=-res_integration[7]+mu[1];mu[2]=-res_integration[8]+mu[2];
      Results<<"The norm of the dipolar mom. |mu| = "<<setw(17)<<norm3D(mu)<<endl;
      Results<<"The result of the integration mux = "<<setw(17)<<mu[0]<<endl;
      Results<<"The result of the integration muy = "<<setw(17)<<mu[1]<<endl;
      Results<<"The result of the integration muz = "<<setw(17)<<mu[2]<<endl;
     }
     if(rho)
     {
      Results<<"The result of < rho >             = "<<setw(17)<<res_integration[9]<<endl;
     }
     Results<<endl;
     Results<<"\t \t Obtained with "<<method<<endl;
     Results<<endl;
    }
    if(mom)
    {
     if(Input_commands.dmn_integrals && !wfn_fchk)
     {
      DMN_P_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      I_cubature2_dmnp(dmn,error_abs,res_integration,error,Integrals_interval,shanp,fishp,p1,p2);
     }
     else
     {
      I_cubature2p(Read_fchk_wfn,error_abs,res_integration,error,Integrals_interval,shanp,fishp,p1,p2);
     }
     Results<<"The result of the integration  Np = "<<setw(17)<<res_integration[0]<<endl;
     Nelec=(double)Read_fchk_wfn.nelectrons;
     if(shanp)
     {
      shannonp=res_integration[1]/Nelec+log(Nelec);
      Results<<"The result of the integration S_p = "<<setw(17)<<shannonp<<endl;
      Results<<"The Shannon entropy power     J_p = ";
      Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp);
      Results<<endl;
      Results<<"[J_p=(1/(2*pi*e)) * e^(2/3 *S_p)]"<<endl;
     }
     if(fishp)
     {
      fisherp=res_integration[2]/Nelec;
      Results<<"The result of the integration I_p = "<<setw(17)<<fisherp<<endl;
     }
     if(p1){Results<<"The result of the integration <p> = "<<setw(17)<<res_integration[4]<<endl;}
     if(p2)
     {
      Results<<"The result of the integration <p2>= "<<setw(17)<<res_integration[5]<<endl;
      Results<<"The kinetic energy [from PI(p)] T = "<<setw(17)<<res_integration[5]/TWO<<endl;
     }
     Results<<endl;
     Results<<"\t \t Obtained with "<<method<<endl;
     Results<<endl;
    }
    if(fish==true && shan==true)
    {
     Results<<"The Fisher-Shannon product    P_r = ";
     Results<<setw(17)<<ONE/THREE*fisher*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon)<<endl;
     Results<<"[P_r = (1/3)*I_r*J_r]"<<endl;
    }
    if(fishp==true && shanp==true)
    {
     Results<<"The Fisher-Shannon product    P_p = ";
     Results<<setw(17)<<ONE/THREE*fisherp*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp)<<endl;
     Results<<"[P_p = (1/3)*I_p*J_p]"<<endl;
    }
    if(shanp==true && shan==true)
    {
     Results<<"The Shannon-Shannonp  (S_r + S_p) = ";
     Results<<setw(17)<<shannon+shannonp<<endl;
     Results<<"[S_Total = S_r + S_p]"<<endl;
     Results<<"[S_r + S_p >= 3(1+ln pi) ~ 6.434189658]"<<endl;
     Results<<"The Shannon-Shannonp  (J_r x J_p) = ";
     Results<<setw(17)<<(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp))*(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon))<<endl;
     Results<<"[J_r*J_p >= 1/4]"<<endl;
    }
    if(fish==true && fishp==true)
    {
     Results<<"The Fisher-Fisherp    (I_r x I_p) = ";
     Results<<setw(17)<<fisher*fisherp<<endl;
     Results<<"[I_r*I_p >= 36]"<<endl;
    }
    Results<<endl;
   }
   else if(Input_commands.quadrature) 
   {
    int grid_theta_phi;
    double **evals_density,**evals_grad_mod,res;
    double Interval[6];
    method="Quadrature";
    void *data=NULL,*data2=NULL;
    grid_theta_phi=Input_commands.order_grid_ang;
    grid_avail(grid_theta_phi);
    if(!Input_commands.dmn_integrals)
    {
     data=&Read_fchk_wfn;
    }
    for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
    {
     if(Input_commands.integral_ops[num_integrals]=="density"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="densityp"){mom=true;}
     if(Input_commands.integral_ops[num_integrals]=="r1"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="p1"){mom=true;}
     if(Input_commands.integral_ops[num_integrals]=="r2"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="p2"){mom=true;}
     if(Input_commands.integral_ops[num_integrals]=="ttf"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="shannon"){pos=true;shan=true;}
     if(Input_commands.integral_ops[num_integrals]=="shannonp"){mom=true;shanp=true;}
     if(Input_commands.integral_ops[num_integrals]=="fisher" || Input_commands.integral_ops[num_integrals]=="tw" )
     {fish=true;pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="fisherp")
     {fishp=true;mom=true;}
     if(Input_commands.integral_ops[num_integrals]=="inertia"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="quadrupole"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="inertiap"){mom=true;}
     if(Input_commands.integral_ops[num_integrals]=="dipolar"){pos=true;}
     if(Input_commands.integral_ops[num_integrals]=="rho"){pos=true;rho=true;}
     if(Input_commands.integral_ops[num_integrals]=="sij"){sij=true;}
     if(Input_commands.integral_ops[num_integrals]=="r1_moment"){r1_moment=true;}
    }
    if(pos)
    {
     evals_density=new double*[Input_commands.order_grid_r];
     for(i=0;i<Input_commands.order_grid_r;i++)
     {
      evals_density[i]=new double[grid_theta_phi];
     }
     if(fish)
     {
      evals_grad_mod=new double*[Input_commands.order_grid_r];
      for(i=0;i<Input_commands.order_grid_r;i++)
      {
       evals_grad_mod[i]=new double[grid_theta_phi];
      }
     }
     //Build all required evaluations of the Density and gradients
     if(!Input_commands.dmn_integrals)
     {
      integrate_quadrature(data,name_file,Input_commands.dmn_integrals,Input_commands.order_grid_r,Input_commands.order_grid_ang,
      fish,evals_density,evals_grad_mod);
     }
     else
     {
      if(rho)
      {
       cout<<"Warning, The < rho > integral is not implemented for FCHK + DM1 files!"<<endl;
       cout<<"Convert the DM1 file + FCHK file into a DM1_FCHK file"<<endl;
      }
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      data=&dmn;
      integrate_quadrature(data,name_file,Input_commands.dmn_integrals,Input_commands.order_grid_r,Input_commands.order_grid_ang,
      fish,evals_density,evals_grad_mod);
     }
     for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
     {
      for(i=0;i<6;i++){Interval[i]=Input_commands.interval_integrals[num_integrals][i];}
      integral_calc(Input_commands.integral_ops[num_integrals],Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      if(Input_commands.integral_ops[num_integrals]=="density")
      {
       Density=res;
       Results<<"The result of the integration  N  = "<<setw(17)<<Density;
       Results<<"\t obtained with "<<method<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="rho")
      {
       Results<<"The result of < rho >             = "<<setw(17)<<res;
       Results<<"\t obtained with "<<method<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="shannon")
      {
       shannon=res/(double)Read_fchk_wfn.nelectrons+log((double)Read_fchk_wfn.nelectrons);
       Results<<"The result of the integration S_r = "<<setw(17)<<shannon;
       Results<<"\t obtained with "<<method<<endl;
       shan=true;
       Results<<"The Shannon entropy power     J_r = ";
       Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon);
       Results<<endl;
       Results<<"[J_r=(1/(2*pi*e)) * e^(2/3 *S_r)]"<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="fisher")
      {
       fisher=res/(double) Read_fchk_wfn.nelectrons;
       Results<<"The result of the integration I_r = "<<setw(17)<<fisher;
       Results<<"\t obtained with "<<method<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="tw")
      {
       Tw=res/EIGHT;
       fisher=res/(double) Read_fchk_wfn.nelectrons;
       Results<<"The result of the integration T_w = "<<setw(17)<<Tw/(double) Read_fchk_wfn.nelectrons;
       Results<<"\t obtained with "<<method<<endl;
       Results<<"The result of the integration T_W = "<<setw(17)<<Tw;
       Results<<"\t obtained with "<<method<<endl;
       Results<<"(T_w is per electron and T_W is the total)"<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="ttf")
      {
       Ttf=res*pow(THREE*PI*PI,TWO/THREE)*THREE/TEN;
       Results<<"The result of the integration T_TF= "<<setw(17)<<Ttf;
       Results<<"\t obtained with "<<method<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="r1")
      {
       Results<<"The result of the integration <r> = "<<setw(17)<<res;
       Results<<"\t obtained with "<<method<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="r2")
      {
       Results<<"The result of the integration <r2>= "<<setw(17)<<res;
       Results<<"\t obtained with "<<method<<endl;
      }
      if((Input_commands.integral_ops[num_integrals]=="inertia" || Input_commands.integral_ops[num_integrals]=="quadrupole"))
      {
       int j;
       double RCC[3];
       double Im[3][3],Quad_rho[3][3];
       Quadrupole=new double*[3];
       Inertia=new double*[3];
       eigenV=new double*[3];
       for(i=0;i<3;i++)
       {
        Quadrupole[i]=new double[3];
        Inertia[i]=new double[3];
        eigenV[i]=new double[3];
       }
       Read_fchk_wfn.quadrupoleATOMS(Quadrupole);
       string operation="density";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       if(abs(res-(int)res)>=HALF)
       {Density=ceil(res);}
       else
       {Density=floor(res);}
       if(Density==ZERO){Density=ONE;}
       operation="Icmx";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       RCC[0]=res/Density;
       operation="Icmy";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       RCC[1]=res/Density;
       operation="Icmz";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       RCC[2]=res/Density;
       Im[0][0]=pow(RCC[1],TWO)+pow(RCC[2],TWO);
       Im[0][1]=-RCC[0]*RCC[1];
       Im[0][2]=-RCC[0]*RCC[2];
       Im[1][0]=Im[0][1];
       Im[1][1]=pow(RCC[0],TWO)+pow(RCC[2],TWO);
       Im[1][2]=-RCC[1]*RCC[2];
       Im[2][0]=Im[0][2];
       Im[2][1]=Im[1][2];
       Im[2][2]=pow(RCC[0],TWO)+pow(RCC[1],TWO);
       operation="Ixx";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       Quad_rho[0][0]=res;
       Inertia[0][0]=res/Density-Im[0][0];
       operation="Ixy";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       Quad_rho[1][0]=-res;
       Quad_rho[0][1]=-res;
       Inertia[1][0]=res/Density-Im[1][0];
       Inertia[0][1]=Inertia[1][0];
       operation="Iyy";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       Quad_rho[1][1]=res;
       Inertia[1][1]=res/Density-Im[1][1];
       operation="Ixz";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       Quad_rho[2][0]=-res;
       Quad_rho[0][2]=-res;
       Inertia[2][0]=res/Density-Im[2][0];
       Inertia[0][2]=Inertia[2][0];
       operation="Iyz";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       Quad_rho[2][1]=-res;
       Quad_rho[1][2]=-res;
       Inertia[2][1]=res/Density-Im[2][1];
       Inertia[1][2]=Inertia[2][1];
       operation="Izz";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       Quad_rho[2][2]=res;
       Inertia[2][2]=res/Density-Im[2][2];
       Results<<endl;
       Results<<"Coordinates of the CC : ";
       Results<<setprecision(10)<<fixed<<scientific;
       for(i=0;i<3;i++){if(RCC[i]<pow(TEN,-SIX)){RCC[i]=ZERO;}}
       Results<<setw(17)<<RCC[0]<<" "<<setw(17)<<RCC[1]<<" "<<setw(17)<<RCC[2]<<endl;
       Results<<endl;
       Results<<"The inertia(r) tensor ( "<<setw(17)<<Inertia[0][0]<<" ";//see below this function
       Results<<setw(17)<<Inertia[0][1]<<" "<<setw(17)<<Inertia[0][2]<<" )"<<endl;
       Results<<"                      ( "<<Inertia[1][0]<<" ";
       Results<<setw(17)<<Inertia[1][1]<<"\t"<<setw(17)<<Inertia[1][2]<<" )"<<endl;
       Results<<"                      ( "<<setw(17)<<Inertia[2][0]<<" ";
       Results<<setw(17)<<Inertia[2][1]<<" "<<setw(17)<<Inertia[2][2]<<" )"<<endl;
       Results<<endl;
       Results<<"Diagonalized tensor:"<<endl;
       Results<<endl;
       i=3;
       jacobi(i,Inertia,eigenV); //Diagonalization
       for(i=0;i<3;i++)
       {
        for(j=0;j<3;j++)
        {
         if(Inertia[i][j]<pow(TEN,-EIGHT) && i!=j){Inertia[i][j]=ZERO;}
         if(eigenV[i][j]<pow(TEN,-SEVEN)){eigenV[i][j]=ZERO;}
        }
       }
       Results<<" ( "<<setw(17)<<Inertia[0][0]<<" "<<setw(17)<<Inertia[0][1]<<" ";
       Results<<setw(17)<<Inertia[0][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Inertia[1][0]<<" "<<setw(17)<<Inertia[1][1]<<" ";
       Results<<setw(17)<<Inertia[1][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Inertia[2][0]<<" "<<setw(17)<<Inertia[2][1]<<" ";
       Results<<setw(17)<<Inertia[2][2]<<" )"<<endl;
       Results<<endl;
       Results<<endl;
       Results<<"Eigenvectors (columns):"<<endl;
       Results<<endl;
       Results<<" ( "<<setw(17)<<eigenV[0][0]<<" "<<setw(17)<<eigenV[0][1]<<" ";
       Results<<setw(17)<<eigenV[0][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<eigenV[1][0]<<" "<<setw(17)<<eigenV[1][1]<<" ";
       Results<<setw(17)<<eigenV[1][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<eigenV[2][0]<<" "<<setw(17)<<eigenV[2][1]<<" ";
       Results<<setw(17)<<eigenV[2][2]<<" )"<<endl;
       Results<<endl;
       Results<<"Trace of Ir_diag [tr(Ir_diag)]    = "<<setw(17)<<Inertia[0][0]+Inertia[1][1]
       +Inertia[2][2]<<endl;
       Results<<endl;
       Quadrupole[0][0]=Quadrupole[0][0]-HALF*(Quad_rho[1][1]+Quad_rho[2][2]-Quad_rho[0][0]);
       Quadrupole[0][1]=Quadrupole[0][1]-Quad_rho[0][1];
       Quadrupole[0][2]=Quadrupole[0][2]-Quad_rho[0][2];
       Quadrupole[1][0]=Quadrupole[1][0]-Quad_rho[1][0];
       Quadrupole[1][1]=Quadrupole[1][1]-HALF*(Quad_rho[0][0]+Quad_rho[2][2]-Quad_rho[1][1]);
       Quadrupole[1][2]=Quadrupole[1][2]-Quad_rho[1][2];
       Quadrupole[2][0]=Quadrupole[2][0]-Quad_rho[2][0];
       Quadrupole[2][1]=Quadrupole[2][1]-Quad_rho[2][1];
       Quadrupole[2][2]=Quadrupole[2][2]-HALF*(Quad_rho[0][0]+Quad_rho[1][1]-Quad_rho[2][2]);
       for(i=0;i<3;i++)
       {
        for(j=0;j<3;j++)
        {
         Quadrupole[i][j]=AUtoDAng*Quadrupole[i][j];
         if(abs(Quadrupole[i][j])<pow(TEN,-SIX))
         {Quadrupole[i][j]=ZERO;}
        }
       }
       Results<<"Primitive Quadrupole Moment (Debye-Angstrom):"<<endl;
       Results<<endl;
       Results<<" ( "<<setw(17)<<Quadrupole[0][0]<<" "<<setw(17)<<Quadrupole[0][1]<<" "<<setw(17)<<Quadrupole[0][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Quadrupole[1][0]<<" "<<setw(17)<<Quadrupole[1][1]<<" "<<setw(17)<<Quadrupole[1][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Quadrupole[2][0]<<" "<<setw(17)<<Quadrupole[2][1]<<" "<<setw(17)<<Quadrupole[2][2]<<" )"<<endl;
       Results<<endl;
       Quad_rho[0][0]=HALF*(TWO*Quadrupole[0][0]-Quadrupole[1][1]-Quadrupole[2][2]);
       Quad_rho[1][1]=HALF*(TWO*Quadrupole[1][1]-Quadrupole[0][0]-Quadrupole[2][2]);
       Quad_rho[2][2]=HALF*(TWO*Quadrupole[2][2]-Quadrupole[1][1]-Quadrupole[0][0]);
       Quad_rho[0][1]=THREE*HALF*Quadrupole[0][1];Quad_rho[1][0]=THREE*HALF*Quadrupole[1][0];
       Quad_rho[0][2]=THREE*HALF*Quadrupole[0][2];Quad_rho[2][0]=THREE*HALF*Quadrupole[2][0];
       Quad_rho[1][2]=THREE*HALF*Quadrupole[1][2];Quad_rho[2][1]=THREE*HALF*Quadrupole[2][1];
       Results<<"Traceless Quadrupole Moment (Debye-Angstrom) [Gamess, Molpro]:"<<endl;
       Results<<endl;
       Results<<" ( "<<setw(17)<<Quad_rho[0][0]<<" "<<setw(17)<<Quad_rho[0][1]<<" "<<setw(17)<<Quad_rho[0][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Quad_rho[1][0]<<" "<<setw(17)<<Quad_rho[1][1]<<" "<<setw(17)<<Quad_rho[1][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Quad_rho[2][0]<<" "<<setw(17)<<Quad_rho[2][1]<<" "<<setw(17)<<Quad_rho[2][2]<<" )"<<endl;
       Results<<endl;
       Results<<"Trace of Quadrupole               = "<<Quad_rho[0][0]+Quad_rho[1][1]+Quad_rho[2][2]<<endl;
       Results<<endl;
       Quad_rho[0][0]=(Quadrupole[0][0]-(Quadrupole[0][0]+Quadrupole[1][1]+Quadrupole[2][2])/THREE);
       Quad_rho[1][1]=(Quadrupole[1][1]-(Quadrupole[0][0]+Quadrupole[1][1]+Quadrupole[2][2])/THREE);
       Quad_rho[2][2]=(Quadrupole[2][2]-(Quadrupole[0][0]+Quadrupole[1][1]+Quadrupole[2][2])/THREE);
       Quad_rho[0][1]=Quadrupole[0][1];Quad_rho[1][0]=Quadrupole[1][0];
       Quad_rho[0][2]=Quadrupole[0][2];Quad_rho[2][0]=Quadrupole[2][0];
       Quad_rho[1][2]=Quadrupole[1][2];Quad_rho[2][1]=Quadrupole[2][1];
       Results<<"Traceless Quadrupole Moment (Debye-Angstrom) [Gaussian]:"<<endl;
       Results<<endl;
       Results<<" ( "<<setw(17)<<Quad_rho[0][0]<<" "<<setw(17)<<Quad_rho[0][1]<<" "<<setw(17)<<Quad_rho[0][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Quad_rho[1][0]<<" "<<setw(17)<<Quad_rho[1][1]<<" "<<setw(17)<<Quad_rho[1][2]<<" )"<<endl;
       Results<<" ( "<<setw(17)<<Quad_rho[2][0]<<" "<<setw(17)<<Quad_rho[2][1]<<" "<<setw(17)<<Quad_rho[2][2]<<" )"<<endl;
       Results<<endl;
       Results<<"Trace of Quadrupole               = "<<Quad_rho[0][0]+Quad_rho[1][1]+Quad_rho[2][2]<<endl;
       Results<<"\t obtained with "<<method<<endl;
       Results<<endl;
       for(i=0;i<3;i++)
       {
        delete[] eigenV[i];eigenV[i]=NULL;
        delete[] Inertia[i];Inertia[i]=NULL;
        delete[] Quadrupole[i];Quadrupole[i]=NULL;
       }
       delete[] eigenV;eigenV=NULL;
       delete[] Inertia;Inertia=NULL;
       delete[] Quadrupole;Quadrupole=NULL;
      }
      if(Input_commands.integral_ops[num_integrals]=="dipolar")
      {
       double mu[3];
       Read_fchk_wfn.muATOMS(mu);
       //Icmx, Icmy and Icm already evaluate the needed integrals for mu!
       operation="Icmx";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       mu[0]=mu[0]-res;
       operation="Icmy";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       mu[1]=mu[1]-res;
       operation="Icmz";
       integral_calc(operation,Input_commands.order_grid_r,grid_theta_phi,
       evals_density,evals_grad_mod,Interval,res);
       mu[2]=mu[2]-res;
       Results<<"The norm of the dipolar mom. |mu| = "<<setw(17)<<norm3D(mu)<<endl;
       Results<<"The result of the integration mux = "<<setw(17)<<mu[0]<<endl;
       Results<<"The result of the integration muy = "<<setw(17)<<mu[1]<<endl;
       Results<<"The result of the integration muz = "<<setw(17)<<mu[2]<<endl;
      }
     }
     if(fish==true && shan==true)
     {
      Results<<"The Fisher-Shannon product    P_r = ";
      Results<<setw(17)<<ONE/THREE*fisher*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon)<<endl;
      Results<<"[P_r = (1/3)*I_r*J_r]"<<endl;
     }
     for(i=0;i<Input_commands.order_grid_r;i++)
     {
      delete[] evals_density[i];
      evals_density[i]= NULL;
      if(fish)
      {
       delete[] evals_grad_mod[i];
       evals_grad_mod[i]=NULL;
      }
     }
     delete[] evals_density;
     evals_density=NULL;
     if(fish){delete[] evals_grad_mod;evals_grad_mod=NULL;}
     clean_quadrature(name_file);
    }
    if(mom)
    {
     evals_density=new double*[Input_commands.order_grid_r];
     for(i=0;i<Input_commands.order_grid_r;i++)
     {
      evals_density[i]=new double[grid_theta_phi];
     }
     if(fishp)
     {
      evals_grad_mod=new double*[Input_commands.order_grid_r];
      for(i=0;i<Input_commands.order_grid_r;i++)
      {
       evals_grad_mod[i]=new double[grid_theta_phi];
      }
     }
     //Build all required evaluations of the Density and gradients
     if(!Input_commands.dmn_integrals)
     {
      integrate_quadrature_p(data,name_file,Input_commands.dmn_integrals,Input_commands.order_grid_r,Input_commands.order_grid_ang,
      fishp,evals_density,evals_grad_mod);
     }
     else
     {
      DMN_P_OPS dmn1(Input_commands.name_dm1,1);
      dmn1.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn1.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn1.set_thershold(Input_commands.dmn_threshold);}
      data2=&dmn1;
      integrate_quadrature_p(data2,name_file,Input_commands.dmn_integrals,Input_commands.order_grid_r,Input_commands.order_grid_ang,
      fishp,evals_density,evals_grad_mod);
     }
     for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
     {
      for(i=0;i<6;i++){Interval[i]=Input_commands.interval_integrals[num_integrals][i];}
      integral_calc_p(Input_commands.integral_ops[num_integrals],Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      if(Input_commands.integral_ops[num_integrals]=="densityp")
      {
       Density=res;
       Results<<"The result of the integration  Np = "<<setw(17)<<Density;
       Results<<"\t obtained with "<<method<<endl;
      }
      if(Input_commands.integral_ops[num_integrals]=="shannonp")
      {
       shannonp=res/(double)Read_fchk_wfn.nelectrons+log((double)Read_fchk_wfn.nelectrons);
       Results<<"The result of the integration S_p = "<<setw(17)<<shannonp;
       Results<<"\t obtained with "<<method<<endl;
       Results<<"The Shannon entropy power     J_p = ";
       Results<<setw(17)<<ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp);
       Results<<endl;
       Results<<"[J_p=(1/(2*pi*e)) * e^(2/3 *S_p)]"<<endl;
       shanp=true;
     }
     if(Input_commands.integral_ops[num_integrals]=="fisherp")
     {
      fisherp=res/(double) Read_fchk_wfn.nelectrons;
      Results<<"The result of the integration I_p = "<<setw(17)<<fisherp;
      Results<<"\t obtained with "<<method<<endl;
     }
     if(Input_commands.integral_ops[num_integrals]=="p1")
     {
      Results<<"The result of the integration <p> = "<<setw(17)<<res;
      Results<<"\t obtained with "<<method<<endl;
     }
     if(Input_commands.integral_ops[num_integrals]=="p2")
     {
      Results<<"The result of the integration <p2>= "<<setw(17)<<res;
      Results<<"\t obtained with "<<method<<endl;
      Results<<"The kinetic energy [from PI(p)] T = "<<setw(17)<<res/TWO<<endl;
     }
     if(Input_commands.integral_ops[num_integrals]=="inertiap")
     {
      int j;
      double RCC[3];
      double Im[3][3];
      Inertia=new double*[3];
      for(i=0;i<3;i++){Inertia[i]=new double[3];}
      eigenV=new double*[3];
      for(i=0;i<3;i++)
      {eigenV[i]=new double[3];}
      string operation="densityp";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      if(abs(res-(int)res)>=HALF)
      {Density=ceil(res);}
      else
      {Density=floor(res);}
      if(Density==ZERO){Density=ONE;}
      operation="Icmx";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      RCC[0]=res/Density;
      operation="Icmy";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      RCC[1]=res/Density;
      operation="Icmz";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      RCC[2]=res/Density;
      Im[0][0]=pow(RCC[1],TWO)+pow(RCC[2],TWO);
      Im[0][1]=-RCC[0]*RCC[1];
      Im[0][2]=-RCC[0]*RCC[2];
      Im[1][0]=Im[0][1];
      Im[1][1]=pow(RCC[0],TWO)+pow(RCC[2],TWO);
      Im[1][2]=-RCC[1]*RCC[2];
      Im[2][0]=Im[0][2];
      Im[2][1]=Im[1][2];
      Im[2][2]=pow(RCC[0],TWO)+pow(RCC[1],TWO);
      operation="Ixx";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      Inertia[0][0]=res/Density-Im[0][0];
      operation="Ixy";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      Inertia[1][0]=res/Density-Im[1][0];
      Inertia[0][1]=Inertia[1][0];
      operation="Iyy";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      Inertia[1][1]=res/Density-Im[1][1];
      operation="Ixz";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      Inertia[2][0]=res/Density-Im[2][0];
      Inertia[0][2]=Inertia[2][0];
      operation="Iyz";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      Inertia[2][1]=res/Density-Im[2][1];
      Inertia[1][2]=Inertia[2][1];
      operation="Izz";
      integral_calc_p(operation,Input_commands.order_grid_r,grid_theta_phi,
      evals_density,evals_grad_mod,Interval,res);
      Inertia[2][2]=res/Density-Im[2][2];
      Results<<endl;
      Results<<"Coordinates of the CCp: ";
      Results<<setprecision(10)<<fixed<<scientific;
      for(i=0;i<3;i++){if(RCC[i]<pow(TEN,-SIX)){RCC[i]=ZERO;}}
      Results<<setw(17)<<RCC[0]<<" "<<setw(17)<<RCC[1]<<" "<<setw(17)<<RCC[2]<<endl;
      Results<<endl;
      Results<<"The inertia(p) tensor ( "<<setw(17)<<Inertia[0][0]<<" ";//see below this function
      Results<<setw(17)<<Inertia[0][1]<<" "<<setw(17)<<Inertia[0][2]<<" )"<<endl;
      Results<<"                      ( "<<setw(17)<<Inertia[1][0]<<" ";
      Results<<setw(17)<<Inertia[1][1]<<" "<<setw(17)<<Inertia[1][2]<<" )"<<endl;
      Results<<"                      ( "<<setw(17)<<Inertia[2][0]<<" ";
      Results<<setw(17)<<Inertia[2][1]<<" "<<setw(17)<<Inertia[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Diagonalized tensor:"<<endl;
      Results<<endl;
      i=3;
      jacobi(i,Inertia,eigenV); //Diagonalization
      for(i=0;i<3;i++)
      {
       for(j=0;j<3;j++)
       {
        if(Inertia[i][j]<pow(TEN,-EIGHT) && i!=j){Inertia[i][j]=ZERO;}
        if(eigenV[i][j]<pow(TEN,-SEVEN)){eigenV[i][j]=ZERO;}
       }
      }
      Results<<" ( "<<setw(17)<<Inertia[0][0]<<" "<<setw(17)<<Inertia[0][1]<<" ";
      Results<<setw(17)<<Inertia[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Inertia[1][0]<<" "<<setw(17)<<Inertia[1][1]<<" ";
      Results<<setw(17)<<Inertia[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<Inertia[2][0]<<" "<<setw(17)<<Inertia[2][1]<<" ";
      Results<<setw(17)<<Inertia[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<endl;
      Results<<"Eigenvectors (columns):"<<endl;
      Results<<endl;
      Results<<" ( "<<setw(17)<<eigenV[0][0]<<" "<<setw(17)<<eigenV[0][1]<<" ";
      Results<<setw(17)<<eigenV[0][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<eigenV[1][0]<<" "<<setw(17)<<eigenV[1][1]<<" ";
      Results<<setw(17)<<eigenV[1][2]<<" )"<<endl;
      Results<<" ( "<<setw(17)<<eigenV[2][0]<<" "<<setw(17)<<eigenV[2][1]<<" ";
      Results<<setw(17)<<eigenV[2][2]<<" )"<<endl;
      Results<<endl;
      Results<<"Trace of Ip_diag [tr(Ip_diag)]    = "<<setw(17)<<Inertia[0][0]+Inertia[1][1]
      +Inertia[2][2];
      Results<<"\t obtained with "<<method<<endl;
      Results<<endl;
      for(i=0;i<3;i++)
      {delete[] eigenV[i];eigenV[i]=NULL;}
      delete[] eigenV;eigenV=NULL;
      for(i=0;i<3;i++)
      {delete[] Inertia[i];Inertia[i]=NULL;}
      delete[] Inertia;Inertia=NULL;
     }
    }
    if(fishp==true && shanp==true)
    {
     Results<<"The Fisher-Shannon product    P_p = ";
     Results<<setw(17)<<ONE/THREE*fisherp*ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp)<<endl;
     Results<<"[P_p = (1/3)*I_p*J_p]"<<endl;
    }
    for(i=0;i<Input_commands.order_grid_r;i++)
    {
     delete[] evals_density[i];
     evals_density[i]= NULL;
     if(fishp)
     {
      delete[] evals_grad_mod[i];
     }
    }
    delete[] evals_density;
    evals_density=NULL;
    if(fishp){delete[] evals_grad_mod;evals_grad_mod=NULL;}
    clean_quadrature(name_file);
    }
    if(shanp==true && shan==true)
    {
     Results<<"The Shannon-Shannonp  (S_r + S_p) = ";
     Results<<setw(17)<<shannon+shannonp<<endl;
     Results<<"[S_Total = S_r + S_p]"<<endl;
     Results<<"[S_r + S_p >= 3(1+ln pi) ~ 6.434189658]"<<endl;
     Results<<"The Shannon-Shannonp  (J_r x J_p) = ";
     Results<<setw(17)<<(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannonp))*(ONE/(TWO*PI*eul)*pow(eul,(TWO/THREE)*shannon))<<endl;
     Results<<"[J_r*J_p >= 1/4]"<<endl;
    }
    if(fish==true && fishp==true)
    {
     Results<<"The Fisher-Fisherp    (I_r x I_p) = ";
     Results<<setw(17)<<fisher*fisherp<<endl;
     Results<<"[I_r*I_p >= 36]"<<endl;
    }
    if(sij || r1_moment)
    {
     Results<<endl;
     Results<<"List of int files generated:"<<endl;
     Results<<endl;
     int j,total_grid=grid_theta_phi*Input_commands.order_grid_r;
     double **ORBITALS,**Sij,rhoa=ZERO,rhob=ZERO;//Total rows=Number of total evaluations & Columns=Number of total MOs
     ORBITALS=new double*[total_grid];           //for r1_moment, Sij will contain moments instead of overlaps of orbitals!
     Sij=new double*[Read_fchk_wfn.nbasis()];
     data=&Read_fchk_wfn;
     for(i=0;i<total_grid;i++)
     {ORBITALS[i]=new double[Read_fchk_wfn.nbasis()];}
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {Sij[i]=new double[i+1];}
     nbasis=Read_fchk_wfn.nbasis();
     if(Read_fchk_wfn.wfn){Input_commands.MOorNO="no";}
     if(Input_commands.MOorNO=="no" && (!Read_fchk_wfn.overlap && !Read_fchk_wfn.wfn))
     {
      if(!Input_commands.dmn_integrals)
      {
       //Only FCHK files calcs when no LOG was given and NOs were asked but NO DMN integrals were required will
       //enter here!!! But be aware that MOorNO will be set to MO since we will the use diagonalization of DM1
       //and the MOs four fuilding the NOs
       cout<<"The PS matrix is not available. The Sij matrix will be computed in MOs"<<endl;
      }
      else
      {mos_to_nos_dmn=true;}
      Input_commands.MOorNO="mo";
     }
     if(!Input_commands.rotate_grid)
     {
      if(Input_commands.MOorNO=="no")
      {
       build_nos_quadrature(data,name_file,ORBITALS,Input_commands.order_grid_r,Input_commands.order_grid_ang);
      }
      else
      {
       build_mos_quadrature(data,name_file,ORBITALS,Input_commands.order_grid_r,Input_commands.order_grid_ang);
      }
     }
     else
     {
      if(Input_commands.MOorNO=="no")
      {
       build_nos_quadrature_rot(data,name_file,ORBITALS,Input_commands.order_grid_r,Input_commands.order_grid_ang,Input_commands.rotate_angle);
      }
      else
      {
       build_mos_quadrature_rot(data,name_file,ORBITALS,Input_commands.order_grid_r,Input_commands.order_grid_ang,Input_commands.rotate_angle);
      }
     }
     if(Input_commands.dmn_integrals && mos_to_nos_dmn)
     {
      mos_to_nos_int_fchk_dm1(Read_fchk_wfn,Input_commands,ORBITALS,total_grid,nbasis);
     }
     for(num_integrals=0;num_integrals<Input_commands.ops;num_integrals++)
     {
      if(Input_commands.integral_ops[num_integrals]=="sij")
      {
       ostringstream convert;
       convert<<region_int;
       region_string=convert.str();
       for(i=0;i<6;i++){Interval[i]=Input_commands.interval_integrals[num_integrals][i];}
       calc_sij_mat(Sij,ORBITALS,Interval,Input_commands.order_grid_r,grid_theta_phi,nbasis);
       if(Read_fchk_wfn.uhf)
       {
        //If we only have alpha NOs in the wfn file(2e- triplet for instance), no_beta_wfn=true. We correct the spin component
        //and do not asume that the first orbital and second orbitals have diff spins (No Sij=zero if i parity is diff from j).
        if(!Read_fchk_wfn.no_beta_wfn)
        {
         for(i=0;i<Read_fchk_wfn.nbasis();i++)
         {
          for(j=0;j<=i;j++)
          {
           if(i%2!=j%2)
           {Sij[i][j]=ZERO;}
          }
         }
        }
       }
       print_int(Read_fchk_wfn,name_file,Sij,nbasis,Density,rhoa,rhob,region_string);
       if(Read_fchk_wfn.wfn)
       {
        Results<<name_file.substr(0,(name_file.length()-4))+"_X"+region_string+".int"<<endl;
       }
       else
       {
        Results<<name_file.substr(0,(name_file.length()-5))+"_X"+region_string+".int"<<endl;
       }
       region_int++;
      }
      if(Input_commands.integral_ops[num_integrals]=="r1_moment")
      {
       ostringstream convert;
       convert<<region_int;
       region_string=convert.str();
       for(i=0;i<6;i++){Interval[i]=Input_commands.interval_integrals[num_integrals][i];}
       calc_mij_mat(Sij,ORBITALS,Interval,Input_commands.order_grid_r,grid_theta_phi,nbasis);
       if(Read_fchk_wfn.uhf)
       {
        //If we only have alpha NOs in the wfn file(2e- triplet for instance), no_beta_wfn=true. We correct the spin component
        //and do not asume that the first orbital and second orbitals have diff spins (No Sij=zero if i parity is diff from j).
        if(!Read_fchk_wfn.no_beta_wfn)
        {
         for(i=0;i<Read_fchk_wfn.nbasis();i++)
         {
          for(j=0;j<=i;j++)
          {
           if(i%2!=j%2)
           {Sij[i][j]=ZERO;}
          }
         }
        }
       }
       print_int(Read_fchk_wfn,name_file,Sij,nbasis,Density,rhoa,rhob,region_string);
       if(Read_fchk_wfn.wfn)
       {
        Results<<name_file.substr(0,(name_file.length()-4))+"_X"+region_string+".int"<<endl;
       }
       else
       {
        Results<<name_file.substr(0,(name_file.length()-5))+"_X"+region_string+".int"<<endl;
       }
       region_int++;
      }
     }
     Results<<endl;
     clean_quadrature(name_file);
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {delete[] Sij[i]; Sij[i]=NULL;}
     for(i=0;i<total_grid;i++)
     {delete[] ORBITALS[i];ORBITALS[i]=NULL;}
     delete[] ORBITALS; ORBITALS=NULL;
     delete[] Sij; Sij=NULL;
    }
   }
   else // Becke
   {
    int grid_theta_phi,nprops=7; // Note: set nprops to numbers of properties !
    double *res_integration,mu[3]={ZERO};
    double T_TF,fact_TF=pow(THREE*PI*PI,TWO/THREE)*THREE/TEN;
#ifdef HAVE_LIBXC
    nprops++; // Add E_xc using LIBXC as an extra integrated property
#endif
    res_integration=new double[nprops*Read_fchk_wfn.natoms+nprops]; // N, mu_x, mu_y, mu_z, r, r^2 (per atom) + molecular N, mu_x, mu_y and mu_z, r, r^2
    if(Input_commands.nprocs>omp_get_max_threads()){Input_commands.nprocs=omp_get_max_threads();}  // nprocs < max OMP threads avail
    if(Input_commands.nprocs>Read_fchk_wfn.natoms){Input_commands.nprocs=Read_fchk_wfn.natoms;}    // nprocs < natoms
    vector<READ_FCHK_WFN>Read_fchk_wfn_th;
    if(wfn_fchk) // WFN/WFX
    {
#ifdef HAVE_LIBXC
    cout<<"Warning! LibXC OMP integrations are only available in serial mode."<<endl; 
    Input_commands.nprocs=1;
#else 
     for(i=0;i<Input_commands.nprocs;i++)
     {
      Read_fchk_wfn_th.push_back(Read_fchk_wfn); // Generate as many copies as nprocs
     }
#endif
    }
    else         // FCHK
    {
     if(Input_commands.nprocs>1)
     {
      cout<<"Parallelization for integrating with Becke/TFVC is currently available only for WFN/WFX files"<<endl;
      Input_commands.nprocs=1;
     }
    }
    method="Becke/TFVC quadrature";
    grid_theta_phi=Input_commands.order_grid_ang;
    grid_avail_becke(grid_theta_phi);
    Grid_becke(Read_fchk_wfn,name_file,Read_fchk_wfn.natoms,Input_commands.order_grid_r,grid_theta_phi,Input_commands.stiff);
    if(wfn_fchk)  // WFN/WFX
    {
#ifdef HAVE_LIBXC
     Integrate_becke(Read_fchk_wfn,res_integration);
#else 
     Integrate_becke_paral(Read_fchk_wfn_th,res_integration,Input_commands.nprocs);
#endif
    }
    else          // FCHK
    {
     Integrate_becke(Read_fchk_wfn,res_integration);
    }
    Results<<endl;
    for(i=0;i<Read_fchk_wfn.natoms;i++)
    {
     Results<<" Atom "<<setw(4)<<i+1<<endl;
     for(j=0;j<nprops;j++)
     {
      if(abs(res_integration[i*nprops+j])<pow(TEN,-EIGHT)){res_integration[i*nprops+j]=ZERO;}
     }
     Results<<" N electrons   = "<<setw(17)<<res_integration[i*nprops]<<endl;
     Results<<" mux           = "<<setw(17)<<res_integration[i*nprops+1]<<endl;;
     Results<<" muy           = "<<setw(17)<<res_integration[i*nprops+2]<<endl;
     Results<<" muz           = "<<setw(17)<<res_integration[i*nprops+3]<<endl;
     Results<<" <r>           = "<<setw(17)<<res_integration[i*nprops+4]<<endl;
     Results<<" <r^2>         = "<<setw(17)<<res_integration[i*nprops+5]<<endl;
     Results<<" T_TF          = "<<setw(17)<<fact_TF*res_integration[i*nprops+6]<<endl;
#ifdef HAVE_LIBXC
     Results<<" E_xc          = "<<setw(17)<<res_integration[i*nprops+7]<<endl;
#endif
    Results<<endl;
    }
    Results<<endl;
    Density=res_integration[nprops*Read_fchk_wfn.natoms];
    Read_fchk_wfn.muATOMS(mu);
    for(i=0;i<3;i++)
    {
     if(abs(mu[i])<pow(TEN,-EIGHT)){mu[i]=ZERO;}
    }
    Results<<"Atomic contribution to mux        = "<<setw(17)<<mu[0]<<endl;
    Results<<"Atomic contribution to muy        = "<<setw(17)<<mu[1]<<endl;
    Results<<"Atomic contribution to muz        = "<<setw(17)<<mu[2]<<endl;
    for(i=0;i<3;i++)
    {
     mu[i]=-res_integration[nprops*Read_fchk_wfn.natoms+i+1]+mu[i];
     if(abs(mu[i])<pow(TEN,-EIGHT)){mu[i]=ZERO;}
    }
    Results<<"The result of the integration mux = "<<setw(17)<<mu[0]<<endl;
    Results<<"The result of the integration muy = "<<setw(17)<<mu[1]<<endl;
    Results<<"The result of the integration muz = "<<setw(17)<<mu[2]<<endl;
    Results<<"The norm of the dipolar mom. |mu| = "<<setw(17)<<norm3D(mu)<<endl;
    Results<<"The result of the integration <r> = "<<setw(17)<<res_integration[nprops*Read_fchk_wfn.natoms+4]<<endl;
    Results<<"The result of the integration <r2>= "<<setw(17)<<res_integration[nprops*Read_fchk_wfn.natoms+5]<<endl;
    T_TF=fact_TF*res_integration[nprops*Read_fchk_wfn.natoms+6];
    Results<<"The result of the integration T_TF= "<<setw(17)<<T_TF<<endl;
#ifdef HAVE_LIBXC
    Results<<"The result of the integration E_xc= "<<setw(17)<<res_integration[nprops*Read_fchk_wfn.natoms+7]<<endl; // MRM: add one to 7 if there is a new
#endif                                                                                                               // prop. after T_TF
    Results<<"The result of the integration  N  = "<<setw(17)<<Density;
    Results<<"\t obtained with "<<method<<endl;
    Results<<endl;
    Results<<" Using a radial grid of  "<<setw(5)<<Input_commands.order_grid_r<<" points"<<endl;
    Results<<" Using an agular grid of "<<setw(5)<<grid_theta_phi<<" points"<<endl;
    Results<<" Using OMP running on    "<<setw(5)<<Input_commands.nprocs<<" threads"<<endl;
    Results<<endl;
    clean_quadrature_becke(name_file,Read_fchk_wfn.natoms);
    delete[] res_integration; res_integration=NULL;
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  //////////////////////////////////////////////////////
  // Divergences as dissimilarities between densities //
  // Quadratic distance, Quantum similarity and       //
  // Kullback-Leibler distance                        //
  // &                                                //
  // V_Hartree                                        //
  //////////////////////////////////////////////////////
  if(Input_commands.dens_sim || Input_commands.v_hartree)
  {
   if(Input_commands.dens_sim)
   {
    Results<<"#*************************************************************************#";
    Results<<endl;
    Results<<"#      Use two densities to quantify the dissimilarities:                 #";
    Results<<endl;
    Results<<"#      Quadratic distance (QD), Quantum similarity index (QSI) and        #";
    Results<<endl;
    Results<<"#      Kullback-Leibler distance (KL).                                    #";
    Results<<endl;
    Results<<"#      See Chapter 9 in Information-theoretic measures of atomic and      #";
    Results<<endl;
    Results<<"#      molecular systems. PhD Thesis of Sheila Lopez Rosa. UGR 2010       #";
    Results<<endl;
    Results<<"#*************************************************************************#";
   }
   else
   {
    Results<<"#*************************************************************************#";
    Results<<endl;
    Results<<"#      Use two densities to compute V_Hartree                             #";
    Results<<endl;
    Results<<"#*************************************************************************#";
   }
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   bool wfn_fchk_in=false;
   int fail=1000,k;
   double res_integration[19]={ZERO},error_abs,error_rel,MIN_EVALS,MAX_EVALS;
   string name_file_saved;
   error_abs=Input_commands.error_abs;
   error_rel=Input_commands.error_rel;
   MIN_EVALS=Input_commands.minevals;
   MAX_EVALS=Input_commands.maxevals;
   method=Input_commands.method_cuba;
   if(method=="vegas"){method="Vegas";}
   if(method=="suave"){method="Suave";}
   if(method=="divonne"){method="Divonne";}
   if(method!="Divonne" && method!="Vegas" && method!="Suave"){method="Cuhre";}
   /////////////////////////////////////////////////
   //Store 2nd FCHK or WFN
   /////////////////////////////////////////////////
   Results<<"Storing second FCHK/WFN/WFX file..."<<endl;
   Results<<endl;
   name_file_saved=name_file;
   name_file=Input_commands.second_fchk_wfn;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    wfn_fchk_in=true;//True for wfn
   }
   READ_FCHK_WFN Read_fchk_wfn_2(name_file,Input_commands.name_log,wfn_fchk_in,Input_commands.log,Input_commands.cas,
   Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
   //Print inertia tensor of nuclei CM and coord. reorientation
   if(Input_commands.cm) // Include $cm KEYWORD
   {
    if(Input_commands.dens_sim)
    {
     Results<<"Reorientation for file "<<name_file<<endl;
     Results<<endl;
    }
    if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
    {
     CM_file.open((name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
     while(getline(CM_file,line))
     {
      Results<<line<<endl;
     }
     CM_file.close();
     system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
    }
    else
    {
     CM_file.open((name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
     while(getline(CM_file,line))
     {
      Results<<line<<endl;
     }
     CM_file.close();
     system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
    }
    Results<<"Note that if the LOG file was included or the MULTIPLICITY."<<endl;
    Results<<"The second FCHK/WFN/WFX will use the same information as for the first FCHK/WFN/WFX."<<endl;
    Results<<endl;
   }
   /////////////////////////////////////////////////
   //Point to two FCHKs or WFNs
   /////////////////////////////////////////////////
   N_FCHKS_WFNS two_fchks_wfns;
   two_fchks_wfns.read_fchk_wfn[0]=&Read_fchk_wfn;
   two_fchks_wfns.read_fchk_wfn[1]=&Read_fchk_wfn_2;
   //How densities are defined:
   /*
    We use the pointer [0] or [1] which belongs to the struct two_fchks_wfns
    Point[0]=0;Point[1]=0;Point[2]=0;
    (*(two_fchks_wfns.read_fchk_wfn[0])).rho_eval(Point,Density);
    cout<<Density<<endl;
    (*(two_fchks_wfns.read_fchk_wfn[1])).rho_eval(Point,Density);
    cout<<Density<<endl;
   */
   //Obtain the rotation matrix between the two systems:
   for(i=0;i<3;i++)
   {
    for(j=0;j<3;j++)
    {
     for(k=0;k<3;k++)
     {
      Rot_grid_matrix[i][j]=Rot_grid_matrix[i][j]+Read_fchk_wfn_2.Rot_ICM[k][i]*Read_fchk_wfn.Rot_ICM[k][j];
     }
    }
   }
   if(Input_commands.symrotdens)
   {
    int counter_int=0;
    double aux,aux2;
    double **Im,**Rot_dens_a,**Rot_dens_b,**Idens_a,**Idens_b,**Eigenv_a,**Eigenv_b;
    Rot_dens_a=new double*[3];
    Rot_dens_b=new double*[3];
    Idens_a=new double*[3];
    Idens_b=new double*[3];
    Eigenv_a=new double*[3];
    Eigenv_b=new double*[3];
    Im=new double*[3];
    for(i=0;i<3;i++)
    {
     Rot_dens_a[i]=new double[3];
     Rot_dens_b[i]=new double[3];
     Idens_a[i]=new double[3];
     Idens_b[i]=new double[3];
     Eigenv_a[i]=new double[3];
     Eigenv_b[i]=new double[3];
     Im[i]=new double[3];
     for(j=0;j<3;j++)
     {
      Rot_grid_matrix[i][j]=ZERO;
      Rot_dens_a[i][j]=ZERO;
      Rot_dens_b[i][j]=ZERO;
      Idens_a[i][j]=ZERO;
      Idens_b[i][j]=ZERO;
      Eigenv_a[i][j]=ZERO;
      Eigenv_b[i][j]=ZERO;
      Im[i][j]=ZERO;
     }
    }
    Integrals_interval[0]=ZERO;
    Integrals_interval[1]=1e99;
    //Functions borrowed from Quadrature
    //Pass from -180 to 180 and 0 to 180 to rad and 0 to 360 and 0 to 180
    aux=ZERO;aux2=TEN*NINE*TWO;
    Integrals_interval[2]=theta_rad(aux);
    Integrals_interval[3]=theta_rad(aux2);
    Integrals_interval[5]=phi_rad(aux2);
    aux2=-TEN*NINE*TWO;
    Integrals_interval[4]=phi_rad(aux2);
    Results<<"Computing Inertia Tensor for "<<name_file_saved<<" density"<<endl;
    Results<<endl;
    if(!Input_commands.symrot_no)
    {
     integrate_cuba(Read_fchk_wfn,method,3,20,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,
     fail,Integrals_interval,Input_commands.nprocs,false,false,true,false,false,false,false,false);
    }
    Results<<"The result of the integration N_1 = "<<setw(17)<<res_integration[0]<<endl;
    if(abs(res_integration[0]-(int)res_integration[0])>=HALF)
    {Density=ceil(res_integration[0]);}
    else
    {Density=floor(res_integration[0]);}
    if(Density==ZERO){Density=ONE;}
    for(i=0;i<3;i++){RCC[i]=res_integration[i+4]/Density;}
    Im[0][0]=pow(RCC[1],TWO)+pow(RCC[2],TWO);
    Im[0][1]=-RCC[0]*RCC[1];
    Im[0][2]=-RCC[0]*RCC[2];
    Im[1][0]=Im[0][1];
    Im[1][1]=pow(RCC[0],TWO)+pow(RCC[2],TWO);
    Im[1][2]=-RCC[1]*RCC[2];
    Im[2][0]=Im[0][2];
    Im[2][1]=Im[1][2];
    Im[2][2]=pow(RCC[0],TWO)+pow(RCC[1],TWO);
    for(i=0;i<3;i++){if(RCC[i]<pow(TEN,-SIX)){RCC[i]=ZERO;}}
    Results<<endl;
    Results<<"Coordinates of the CC: \t";
    Results<<setprecision(10)<<fixed<<scientific;
    Results<<setw(17)<<RCC[0]<<" "<<setw(17)<<RCC[1]<<" "<<setw(17)<<RCC[2]<<endl;
    Results<<endl;
    for(i=0;i<Read_fchk_wfn.natoms;i++)
    {
     for(j=0;j<3;j++)
     {
      Read_fchk_wfn.Cartesian_Coor[i][j]=Read_fchk_wfn.Cartesian_Coor[i][j]-RCC[j];
     }
    }
    counter_int=0;
    for(i=0;i<3;i++)
    {
     for(j=0;j<=i;j++)
     {
      Idens_a[i][j]=res_integration[counter_int+7]/Density-Im[i][j];
      Idens_a[j][i]=Idens_a[i][j];
      counter_int++;
     }
    }
    Results<<"The inertia(r) tensor ( "<<setw(17)<<Idens_a[0][0]<<" ";//see below this function
    Results<<setw(17)<<Idens_a[0][1]<<" "<<setw(17)<<Idens_a[0][2]<<" )"<<endl;
    Results<<"                      ( "<<setw(17)<<Idens_a[1][0]<<" ";
    Results<<setw(17)<<Idens_a[1][1]<<" "<<setw(17)<<Idens_a[1][2]<<" )"<<endl;
    Results<<"                      ( "<<setw(17)<<Idens_a[2][0]<<" ";
    Results<<setw(17)<<Idens_a[2][1]<<" "<<setw(17)<<Idens_a[2][2]<<" )"<<endl;
    Results<<endl;
    Results<<"Diagonalized and ordered tensor:"<<endl;
    Results<<endl;
    i=3;
    jacobi(i,Idens_a,Eigenv_a); //Diagonalization
    for(i=0;i<3;i++)
    {
     for(j=0;j<3;j++)
     {
      if(abs(Idens_a[i][j])<pow(TEN,-EIGHT) && i!=j){Idens_a[i][j]=ZERO;}
      if(abs(Eigenv_a[i][j])<pow(TEN,-SEVEN)){Eigenv_a[i][j]=ZERO;}
     }
    }
    for(i=0;i<3;i++)
    {
     for(j=i;j<3;j++)
     {
      if(abs(Idens_a[i][i])<abs(Idens_a[j][j]))
      {
       aux=Idens_a[i][i];
       Idens_a[i][i]=Idens_a[j][j];
       Idens_a[j][j]=aux;
       for(k=0;k<3;k++)
       {
        aux=Eigenv_a[k][i];
        Eigenv_a[k][i]=Eigenv_a[k][j];
        Eigenv_a[k][j]=aux;
       }
      }
     }
    }
    aux=ZERO;
    Results<<" ( "<<setw(17)<<Idens_a[0][0]<<" "<<setw(17)<<Idens_a[0][1]<<" ";
    Results<<setw(17)<<Idens_a[0][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Idens_a[1][0]<<" "<<setw(17)<<Idens_a[1][1]<<" ";
    Results<<setw(17)<<Idens_a[1][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Idens_a[2][0]<<" "<<setw(17)<<Idens_a[2][1]<<" ";
    Results<<setw(17)<<Idens_a[2][2]<<" )"<<endl;
    Results<<endl;
    Results<<endl;
    Results<<"Eigenvectors (columns):"<<endl;
    Results<<endl;
    Results<<" ( "<<setw(17)<<Eigenv_a[0][0]<<" "<<setw(17)<<Eigenv_a[0][1]<<" ";
    Results<<setw(17)<<Eigenv_a[0][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Eigenv_a[1][0]<<" "<<setw(17)<<Eigenv_a[1][1]<<" ";
    Results<<setw(17)<<Eigenv_a[1][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Eigenv_a[2][0]<<" "<<setw(17)<<Eigenv_a[2][1]<<" ";
    Results<<setw(17)<<Eigenv_a[2][2]<<" )"<<endl;
    Results<<endl;
    Results<<"Trace of I_diag [tr(I_diag)]      = "<<setw(17)<<Idens_a[0][0]+Idens_a[1][1]
    +Idens_a[2][2]<<endl;
    Results<<endl;
    for(i=0;i<19;i++){res_integration[i]=ZERO;}
    Results<<"Computing Inertia Tensor for "<<name_file<<" density"<<endl;
    Results<<endl;
    if(!Input_commands.symrot_no)
    {
     integrate_cuba(Read_fchk_wfn_2,method,3,20,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,
     fail,Integrals_interval,Input_commands.nprocs,false,false,true,false,false,false,false,false);
    }
    Results<<"The result of the integration N_2 = "<<setw(17)<<res_integration[0]<<endl;
    if(abs(res_integration[0]-(int)res_integration[0])>=HALF)
    {Density=ceil(res_integration[0]);}
    else
    {Density=floor(res_integration[0]);}
    if(Density==ZERO){Density=ONE;}
    for(i=0;i<3;i++){RCC[i]=res_integration[i+4]/Density;}
    Im[0][0]=pow(RCC[1],TWO)+pow(RCC[2],TWO);
    Im[0][1]=-RCC[0]*RCC[1];
    Im[0][2]=-RCC[0]*RCC[2];
    Im[1][0]=Im[0][1];
    Im[1][1]=pow(RCC[0],TWO)+pow(RCC[2],TWO);
    Im[1][2]=-RCC[1]*RCC[2];
    Im[2][0]=Im[0][2];
    Im[2][1]=Im[1][2];
    Im[2][2]=pow(RCC[0],TWO)+pow(RCC[1],TWO);
    for(i=0;i<3;i++){if(RCC[i]<pow(TEN,-SIX)){RCC[i]=ZERO;}}
    Results<<endl;
    Results<<"Coordinates of the CC: \t";
    Results<<setprecision(10)<<fixed<<scientific;
    Results<<setw(17)<<RCC[0]<<" "<<setw(17)<<RCC[1]<<" "<<setw(17)<<RCC[2]<<endl;
    Results<<endl;
    for(i=0;i<Read_fchk_wfn_2.natoms;i++)
    {
     for(j=0;j<3;j++)
     {
      Read_fchk_wfn_2.Cartesian_Coor[i][j]=Read_fchk_wfn_2.Cartesian_Coor[i][j]-RCC[j];
     }
    }
    counter_int=0;
    for(i=0;i<3;i++)
    {
     for(j=0;j<=i;j++)
     {
      Idens_b[i][j]=res_integration[counter_int+7]/Density-Im[i][j];
      Idens_b[j][i]=Idens_b[i][j];
      counter_int++;
     }
    }
    Results<<"The inertia(r) tensor ( "<<setw(17)<<Idens_b[0][0]<<" ";//see below this function
    Results<<setw(17)<<Idens_b[0][1]<<" "<<setw(17)<<Idens_b[0][2]<<" )"<<endl;
    Results<<"                      ( "<<setw(17)<<Idens_b[1][0]<<" ";
    Results<<setw(17)<<Idens_b[1][1]<<" "<<setw(17)<<Idens_b[1][2]<<" )"<<endl;
    Results<<"                      ( "<<setw(17)<<Idens_b[2][0]<<" ";
    Results<<setw(17)<<Idens_b[2][1]<<" "<<setw(17)<<Idens_b[2][2]<<" )"<<endl;
    Results<<endl;
    Results<<"Diagonalized and ordered tensor:"<<endl;
    Results<<endl;
    i=3;
    jacobi(i,Idens_b,Eigenv_b); //Diagonalization
    for(i=0;i<3;i++)
    {
     for(j=0;j<3;j++)
     {
      if(abs(Idens_b[i][j])<pow(TEN,-EIGHT) && i!=j){Idens_b[i][j]=ZERO;}
      if(abs(Eigenv_b[i][j])<pow(TEN,-SEVEN)){Eigenv_b[i][j]=ZERO;}
     }
    }
    for(i=0;i<3;i++)
    {
     for(j=i;j<3;j++)
     {
      if(abs(Idens_b[i][i])<abs(Idens_b[j][j]))
      {
       aux=Idens_b[i][i];
       Idens_b[i][i]=Idens_b[j][j];
       Idens_b[j][j]=aux;
       for(k=0;k<3;k++)
       {
        aux=Eigenv_b[k][i];
        Eigenv_b[k][i]=Eigenv_b[k][j];
        Eigenv_b[k][j]=aux;
       }
      }
     }
    }
    aux=ZERO;
    Results<<" ( "<<setw(17)<<Idens_b[0][0]<<" "<<setw(17)<<Idens_b[0][1]<<" ";
    Results<<setw(17)<<Idens_b[0][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Idens_b[1][0]<<" "<<setw(17)<<Idens_b[1][1]<<" ";
    Results<<setw(17)<<Idens_b[1][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Idens_b[2][0]<<" "<<setw(17)<<Idens_b[2][1]<<" ";
    Results<<setw(17)<<Idens_b[2][2]<<" )"<<endl;
    Results<<endl;
    Results<<endl;
    Results<<"Eigenvectors (columns):"<<endl;
    Results<<endl;
    Results<<" ( "<<setw(17)<<Eigenv_b[0][0]<<" "<<setw(17)<<Eigenv_b[0][1]<<" ";
    Results<<setw(17)<<Eigenv_b[0][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Eigenv_b[1][0]<<" "<<setw(17)<<Eigenv_b[1][1]<<" ";
    Results<<setw(17)<<Eigenv_b[1][2]<<" )"<<endl;
    Results<<" ( "<<setw(17)<<Eigenv_b[2][0]<<" "<<setw(17)<<Eigenv_b[2][1]<<" ";
    Results<<setw(17)<<Eigenv_b[2][2]<<" )"<<endl;
    Results<<endl;
    Results<<"Trace of I_diag [tr(I_diag)]      = "<<setw(17)<<Idens_b[0][0]+Idens_b[1][1]
    +Idens_b[2][2]<<endl;
    Results<<endl;
    for(i=0;i<3;i++)
    {
     for(j=0;j<3;j++)
     {
      for(k=0;k<3;k++)
      {
       Rot_grid_matrix[i][j]=Rot_grid_matrix[i][j]+Eigenv_b[i][k]*Eigenv_a[j][k];
      }
     }
    }
   }
   name_file=name_file_saved;
   if(Input_commands.dens_sim)
   {
    Results<<"The final rotation matrix for "<<Input_commands.second_fchk_wfn<<" file is"<<endl;
    if(Input_commands.symrot_no)
    {
     for(i=0;i<3;i++)
     {
      for(j=0;j<3;j++)
      {
       Rot_grid_matrix[i][j]=ZERO;
       if(i==j) Rot_grid_matrix[i][j]=ONE;
      }
     }
    }
    Results<<endl;
    for(i=0;i<3;i++)
    {
     for(j=0;j<3;j++)
     {
      if(abs(Rot_grid_matrix[i][j])<pow(TEN,-TEN)){Rot_grid_matrix[i][j]=ZERO;}
      Results<<setw(17)<<Rot_grid_matrix[i][j];
     }
     Results<<endl;
    }
    Results<<endl;
    integrate_dens_sim(two_fchks_wfns,method,3,7,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,fail,Input_commands.nprocs,Rot_grid_matrix);
   }
   else
   {
    for(i=0;i<3;i++){Input_commands.init_coord_r[i]=Input_commands.init_coord_r[i]*Angs2au;}
    integrate_V_Hartree(two_fchks_wfns,method,6,2,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,fail,Input_commands.nprocs,
    Input_commands.init_coord_r);
   }
   Density_alpha=res_integration[0];
   Density_beta=res_integration[1];
   if(Input_commands.dens_sim)
   {
    Results<<"The result of the integration N_1 = "<<setw(17)<<res_integration[0]<<endl;
    Results<<"The result of the integration N_2 = "<<setw(17)<<res_integration[1]<<endl;
    Results<<"The result of < rho_1 >           = "<<setw(17)<<res_integration[2]<<endl;
    Results<<"The result of < rho_2 >           = "<<setw(17)<<res_integration[3]<<endl;
    Results<<"The result of int rho_1 rho_2 dr  = "<<setw(17)<<res_integration[4]<<endl;
    Results<<"The result of the integration QD  = "<<setw(17)<<pow(res_integration[5],HALF)<<endl;
    Results<<"QD[rho_1,rho_2] = (int [rho_1-rho_2]^2 dr)^1/2"<<endl;
    Results<<"The result of the integration QSI = "<<setw(17)<<res_integration[4]/sqrt(pow(TEN,-TEN)+abs(res_integration[2]*res_integration[3]))<<endl;
    Results<<"QSI[rho_1,rho_2] = int rho_1 rho_2 dr/sqrt[< rho_1 > < rho_2 >]"<<endl;
    Results<<"The result of the integration KL  = "<<setw(17)<<(res_integration[6]/(double)Read_fchk_wfn.nelectrons)
                                                   +log((round(res_integration[1])+pow(TEN,-TEN))/(round(res_integration[0])+pow(TEN,-TEN)))<<endl;
    Results<<"KL[rho_1,rho_2] = int rho_1 ln [rho_1/rho_2]"<<endl;
    Results<<endl;
   }
   else
   {
    Results<<"The result of N_1 x N_2           = "<<setw(17)<<res_integration[0]<<endl;
    Results<<"The result of V_Hartree (au)      = "<<setw(17)<<HALF*res_integration[1]<<endl;
    Results<<"The result of V_Hartree (eV)      = "<<setw(17)<<HALF*res_integration[1]*au2eV<<endl;
    Results<<"V_Hartree = 1/2 int rho(r) rho(r') / | r - r' | dr dr' "<<endl;
    Results<<endl;
   }
   Results<<endl;
   Results<<" obtained with "<<Input_commands.method_cuba<<" fail report: "<<fail<<endl;
   Results<<endl;
   Results<<endl;
   if(Input_commands.symgrad)
   {
    integrate_dens_sim2(two_fchks_wfns,method,3,7,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,fail,Input_commands.nprocs,Rot_grid_matrix);
    Results<<"The result of the integration G_1 = "<<setw(17)<<res_integration[0]<<endl;
    Results<<"The result of the integration G_2 = "<<setw(17)<<res_integration[1]<<endl;
    Results<<"The result of int [g(r)_1]^2      = "<<setw(17)<<res_integration[2]<<endl;
    Results<<"The result of int [g(r)_2]^2      = "<<setw(17)<<res_integration[3]<<endl;
    Results<<"The result of int g(r)_1 g(r)_2 dr= "<<setw(17)<<res_integration[4]<<endl;
    Results<<"The result of the integration QDg = "<<setw(17)<<pow(res_integration[5],HALF)<<endl;
    Results<<"QDg[g(r)_1,g(r)_2] = (int [g(r)_1-g(r)_1]^2 dr)^1/2"<<endl;
    Results<<"The result of the integration QSIg= "<<setw(17)<<res_integration[4]/sqrt(pow(TEN,-TEN)+abs(res_integration[2]*res_integration[3]))<<endl;
    Results<<"QSIg[g(r)_2,g(r)_2] = int g(r)_1 g(r)_2 dr/sqrt[int [g(r)_1]^2 dr int [g(r)_1]^2 dr]"<<endl;
    Results<<"The result of the integration KLg = "<<setw(17)<<(res_integration[6]/((double)Read_fchk_wfn.nelectrons))
                                                   +log((round(Density_beta)+pow(TEN,-TEN))/(round(Density_alpha)+pow(TEN,-TEN)))<<endl;
    Results<<"KLs[g(r)_1,g(r)_2] = int g(r)_1 ln [g(r)_1/g(r)_2]"<<endl;
    Results<<"where g(r) = |Grad rho(r)|"<<endl;
    Density_alpha=ZERO;
    Density_beta=ZERO;
    Results<<endl;
    Results<<endl;
    Results<<" obtained with "<<Input_commands.method_cuba<<" fail report: "<<fail<<endl;
    Results<<endl;
   }
   fail=1000;
   Results<<"(fail 0 means finish without problem, -1 error, >1 accuracy";
   Results<<" not achived and 10000 not calculated)"<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  //////////////////////////////////////////////////////
  // Polarizabilities and Hyperpolarizabilities from  //
  // densities. See Nakano PRA, 55, 1503, 1997        //
  //////////////////////////////////////////////////////
  if(Input_commands.int_pol_hyperpol)
  {
   cout<<"Computing integrated polarizabilities"<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#      Use field densities to compute polarizabilities and                #";
   Results<<endl;
   Results<<"#      hyperpolarizabilities using Eqs. 22 and 23 in 18 and 19 of         #";
   Results<<endl;
   Results<<"#      Nakano et al. PRA, 55, 1503 (1997).                                #";
   Results<<endl;
   Results<<"#      We also include new indicators simmilar to QD for Eqs. 22 and 23   #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   int fail=1000;
   double res_integration[9]={ZERO},field,error_abs,error_rel,MIN_EVALS,MAX_EVALS;
   Nelec=Read_fchk_wfn.nelectrons;
   string name_file_saved;
   field=Input_commands.field_F;
   error_abs=Input_commands.error_abs;
   error_rel=Input_commands.error_rel;
   MIN_EVALS=Input_commands.minevals;
   MAX_EVALS=Input_commands.maxevals;
   method=Input_commands.method_cuba;
   if(method=="vegas"){method="Vegas";}
   if(method=="suave"){method="Suave";}
   if(method=="divonne"){method="Divonne";}
   if(method!="Divonne" && method!="Vegas" && method!="Suave"){method="Cuhre";}
   Results<<"Storing the rest of FCHK/WFN/WFX files..."<<endl;
   Results<<endl;
   name_file_saved=name_file;
   /////////////////////////////////////////////////
   //Store 2nd FCHK or WFN or WFX
   /////////////////////////////////////////////////
   wfn_fchk=false;
   name_file=Input_commands.second_fchk_wfn;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    wfn_fchk=true;//True for wfn
   }
   READ_FCHK_WFN Read_fchk_wfn_2(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,
   Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
   }
   else
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
   }
   /////////////////////////////////////////////////
   //Store 3rd FCHK or WFN or WFX
   /////////////////////////////////////////////////
   wfn_fchk=false;
   name_file=Input_commands.third_fchk_wfn;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    wfn_fchk=true;//True for wfn
   }
   READ_FCHK_WFN Read_fchk_wfn_3(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,
   Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
   }
   else
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
   }
   /////////////////////////////////////////////////
   //Store 4th FCHK or WFN or WFX
   /////////////////////////////////////////////////
   wfn_fchk=false;
   name_file=Input_commands.fourth_fchk_wfn;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    wfn_fchk=true;//True for wfn
   }
   READ_FCHK_WFN Read_fchk_wfn_4(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,
   Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
   }
   else
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
   }
   /////////////////////////////////////////////////
   //Store 5th FCHK or WFN or WFX
   /////////////////////////////////////////////////
   wfn_fchk=false;
   name_file=Input_commands.fifth_fchk_wfn;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    wfn_fchk=true;//True for wfn
   }
   READ_FCHK_WFN Read_fchk_wfn_5(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,
   Input_commands.cm,Input_commands.multiplicity);//Construct with parametric construction.
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
   }
   else
   {
    system(("/bin/rm  "+name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
   }
   Results<<"Note that if the LOG file was included or the MULTIPLICITY,"<<endl;
   Results<<"the 2nd, 3rd, 4th, and 5th FCHK/WFN/WFX will use the same information as for the first FCHK/WFN/WFX."<<endl;
   Results<<endl;
   name_file=name_file_saved;
   /////////////////////////////////////////////////
   //Point to five FCHKs or WFNs
   /////////////////////////////////////////////////
   N_FCHKS_WFNS five_fchks_wfns;
   five_fchks_wfns.read_fchk_wfn[0]=&Read_fchk_wfn;
   five_fchks_wfns.read_fchk_wfn[1]=&Read_fchk_wfn_2;
   five_fchks_wfns.read_fchk_wfn[2]=&Read_fchk_wfn_3;
   five_fchks_wfns.read_fchk_wfn[3]=&Read_fchk_wfn_4;
   five_fchks_wfns.read_fchk_wfn[4]=&Read_fchk_wfn_5;
   //How densities are defined:
   /*
    We use the pointer [0] or [1] which belongs to the struct two_fchks_wfns
    Point[0]=0;Point[1]=0;Point[2]=0;
    (*(two_fchks_wfns.read_fchk_wfn[0])).rho_eval(Point,Density);
    cout<<Density<<endl;
    (*(two_fchks_wfns.read_fchk_wfn[1])).rho_eval(Point,Density);
    cout<<Density<<endl;
   */
   Results<<endl;
   integrate_pol_hyperpol(five_fchks_wfns,method,3,9,error_rel,error_abs,MIN_EVALS,MAX_EVALS,res_integration,fail,Input_commands.nprocs,Input_commands.dir_pol_hyper);
   Results<<"The result of the integration N_1 = "<<setw(17)<<res_integration[0]<<endl;
   Results<<"The result of the integration N_2 = "<<setw(17)<<res_integration[1]<<endl;
   Results<<"The result of the integration N_3 = "<<setw(17)<<res_integration[2]<<endl;
   Results<<"The result of the integration N_4 = "<<setw(17)<<res_integration[3]<<endl;
   Results<<"The result of the integration N_5 = "<<setw(17)<<res_integration[4]<<endl;
   Results<<"The result of alpha_ii            = "<<setw(17)<<-res_integration[5]/(TWO*field)<<endl;
   Results<<"The result of gamma_iiii          = "<<setw(17)<<-res_integration[6]/(TWO*pow(field,THREE))<<endl;
   Results<<"[i = coordinate used x, y or z]"<<endl;
   Results<<"The result of the integration QD1 = "<<setw(17)<<sqrt(res_integration[7])/(TWO*field*Nelec)<<endl; 
   Results<<"The result of the integration QD2 = "<<setw(17)<<sqrt(res_integration[8])/(TWO*pow(field,THREE)*Nelec)<<endl;
   Results<<"QD1[rho(r,F)] = 1/2FN (int [rho(r,F)-rho(r,-F)]^2 dr)^1/2"<<endl;    
   Results<<"QD2[rho(r,F)] = 1/2F^3 N (int [rho(r,2F)-rho(r,-2F)-2(rho(r,F)-rho(r,-F))]^2 dr)^1/2"<<endl;    
   Results<<"Note: The 1/2! and 1/3! terms are not included in alpha_ii and gamma_iiii."<<endl;    
   Results<<endl;
   Results<<endl;
   Results<<" obtained with "<<Input_commands.method_cuba<<" fail report: "<<fail<<endl;
   Results<<endl;
   fail=1000;
   Results<<"(fail 0 means finish without problem, -1 error, >1 accuracy";
   Results<<" not achived and 10000 not calculated)"<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  /////////////////////////////////////////////
  // Total position spread tensor evaluation //
  /////////////////////////////////////////////
  if(Input_commands.tps)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#      Use the density to evaluate the Total Position Spread (TPS)        #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   int fail=1000;
   double res_integration_fchk[10]={0};
   TPS=new double*[3];
   for(i=0;i<3;i++){TPS[i]=new double[3];}
   Results<<setprecision(10)<<fixed<<scientific;
   integrate_tps_fchk(Read_fchk_wfn,Input_commands.method_cuba,3,10,Input_commands.error_rel,Input_commands.error_abs,
   Input_commands.minevals,Input_commands.maxevals,res_integration_fchk,fail,Input_commands.nprocs);
   Results<<" N  =  Int Rho(r) dr       = "<<setw(17)<<res_integration_fchk[0]<<endl;
   Results<<"   <x>               <y>               <z>"<<endl;
   for(i=1;i<4;i++)
   {
    if(abs(res_integration_fchk[i])<pow(TEN,-EIGHT)){res_integration_fchk[i]=ZERO;}
   }
   Results<<setw(17)<<res_integration_fchk[1]<<" "<<setw(17)<<res_integration_fchk[2];
   Results<<" "<<setw(17)<<res_integration_fchk[3]<<endl;
   Results<<endl;
   Results<<"TPS approximated: "<<endl;
   Results<<endl;
   //xx
   TPS[0][0]=(res_integration_fchk[4]-res_integration_fchk[1]*res_integration_fchk[1]);
   //yx
   TPS[0][1]=(res_integration_fchk[5]-res_integration_fchk[1]*res_integration_fchk[2]);
   //xz
   TPS[0][2]=(res_integration_fchk[7]-res_integration_fchk[1]*res_integration_fchk[3]);
   //xy
   TPS[1][0]=(res_integration_fchk[5]-res_integration_fchk[1]*res_integration_fchk[2]);
   //yy
   TPS[1][1]=(res_integration_fchk[6]-res_integration_fchk[2]*res_integration_fchk[2]);
   //yz
   TPS[1][2]=(res_integration_fchk[8]-res_integration_fchk[2]*res_integration_fchk[3]);
   //zx
   TPS[2][0]=(res_integration_fchk[7]-res_integration_fchk[1]*res_integration_fchk[3]);
   //zy
   TPS[2][1]=(res_integration_fchk[8]-res_integration_fchk[2]*res_integration_fchk[3]);
   //zz
   TPS[2][2]=(res_integration_fchk[9]-res_integration_fchk[3]*res_integration_fchk[3]);
   for(i=0;i<3;i++)
   {
    for(j=0;j<3;j++)
    {
     if(abs(TPS[i][j])<pow(TEN,-EIGHT)){TPS[i][j]=ZERO;}
    }
   }
   //Print appproximated TPS
   Results<<" ( "<<setw(17)<<TPS[0][0]<<" "<<setw(17)<<TPS[0][1]<<" ";
   Results<<setw(17)<<TPS[0][2]<<" )"<<endl;
   Results<<" ( "<<setw(17)<<TPS[1][0]<<" "<<setw(17)<<TPS[1][1]<<" ";
   Results<<setw(17)<<TPS[1][2]<<" )"<<endl;
   Results<<" ( "<<setw(17)<<TPS[2][0]<<" "<<setw(17)<<TPS[2][1]<<" ";
   Results<<setw(17)<<TPS[2][2]<<" )"<<endl;
   Results<<endl;
   Results<<"TPS exact       : "<<endl;
   Results<<endl;
   //xx
   TPS[0][0]=TPS[0][0]+Input_commands.second_moments_tps[0];
   //yx
   TPS[0][1]=TPS[0][1]+Input_commands.second_moments_tps[1];
   //xz
   TPS[0][2]=TPS[0][2]+Input_commands.second_moments_tps[2];
   //xy
   TPS[1][0]=TPS[1][0]+Input_commands.second_moments_tps[1];
   //yy
   TPS[1][1]=TPS[1][1]+Input_commands.second_moments_tps[3];
   //yz
   TPS[1][2]=TPS[1][2]+Input_commands.second_moments_tps[4];
   //zx
   TPS[2][0]=TPS[2][0]+Input_commands.second_moments_tps[2];
   //zy
   TPS[2][1]=TPS[2][1]+Input_commands.second_moments_tps[4];
   //zz
   TPS[2][2]=TPS[2][2]+Input_commands.second_moments_tps[5];
   //Print exact TPS
   Results<<" ( "<<setw(17)<<TPS[0][0]<<" "<<setw(17)<<TPS[0][1]<<" ";
   Results<<setw(17)<<TPS[0][2]<<" )"<<endl;
   Results<<" ( "<<setw(17)<<TPS[1][0]<<" "<<setw(17)<<TPS[1][1]<<" ";
   Results<<setw(17)<<TPS[1][2]<<" )"<<endl;
   Results<<" ( "<<setw(17)<<TPS[2][0]<<" "<<setw(17)<<TPS[2][1]<<" ";
   Results<<setw(17)<<TPS[2][2]<<" )"<<endl;
   Results<<endl;
   Results<<endl;
   Results<<" obtained with "<<Input_commands.method_cuba<<" fail report: "<<fail<<endl;
   Results<<endl;
   fail=1000;
   Results<<"(fail 0 means finish without problem, -1 error, >1 accuracy";
   Results<<" not achived and 10000 not calculated)"<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   for(i=0;i<3;i++)
   {
    delete[] TPS[i];TPS[i]=NULL;
   }
   delete[] TPS;TPS=NULL;
  }
  ////////////////////////////////////
  //  Intracule evaluations         //
  ////////////////////////////////////
  if(Input_commands.intracule && !wfn_fchk)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#          Intracule density with normalization N(N-1)/2                  #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   int MIN_EVALS,MAX_EVALS,fail=10000;
   double error_abs,error_rel;
   string method;
   DMN_OPS dmn(Input_commands.name_dm2,2);
   dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
   dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
   if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
   error_abs=Input_commands.error_abs;
   error_rel=Input_commands.error_rel;
   MIN_EVALS=Input_commands.minevals;
   MAX_EVALS=Input_commands.maxevals;
   method=Input_commands.method_cuba;
   if(method=="vegas"){method="Vegas";}
   if(method=="suave"){method="Suave";}
   if(method=="divonne"){method="Divonne";}
   if(method!="Divonne" && method!="Vegas" && method!="Suave"){method="Cuhre";}
   Results<<"Intracule density"<<endl;
   Results<<"*****************"<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Results<<"Tr DM2             : "<<setw(17)<<HALF*dmn.calc_tr()<<endl;
   Results<<"r12                : ";
   if(Input_commands.Point_intra[0]>=ZERO){Results<<" ";}
   Results<<Input_commands.Point_intra[0];
   Results<<" ";
   if(Input_commands.Point_intra[1]>=ZERO){Results<<" ";}
   Results<<Input_commands.Point_intra[1];
   Results<<" ";
   if(Input_commands.Point_intra[2]>=ZERO){Results<<" ";}
   Results<<Input_commands.Point_intra[2]<<endl;
   Intracule=ZERO;
   integrate_intra(dmn,method,3,1,error_rel,error_abs,MIN_EVALS,MAX_EVALS,Intracule,fail,
   Input_commands.nprocs,Input_commands.Point_intra);
   Results<<"I(x12,y12,z12)     : "<<setw(17)<<HALF*Intracule<<endl;
   Results<<endl;
   Results<<"\t \t Obtained with "<<method<<", fail report "<<fail<<endl;
   Results<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ////////////////////////////////////
  //  Extracule evaluations         //
  ////////////////////////////////////
  if(Input_commands.extracule && !wfn_fchk)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#          Extracule density with normalization N(N-1)/2                  #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   int MIN_EVALS,MAX_EVALS,fail=10000;
   double error_abs,error_rel;
   string method;
   DMN_OPS dmn(Input_commands.name_dm2,2);
   dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
   dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
   if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
   error_abs=Input_commands.error_abs;
   error_rel=Input_commands.error_rel;
   MIN_EVALS=Input_commands.minevals;
   MAX_EVALS=Input_commands.maxevals;
   method=Input_commands.method_cuba;
   if(method=="vegas"){method="Vegas";}
   if(method=="suave"){method="Suave";}
   if(method=="divonne"){method="Divonne";}
   if(method!="Divonne" && method!="Vegas" && method!="Suave"){method="Cuhre";}
   Results<<"Extracule density"<<endl;
   Results<<"*****************"<<endl;
   Results<<endl;
   Results<<setprecision(10)<<fixed<<scientific;
   Results<<"Tr DM2             : "<<setw(17)<<HALF*dmn.calc_tr()<<endl;
   Results<<"R12                : ";
   if(Input_commands.Point_extra[0]>=ZERO){Results<<" ";}
   Results<<Input_commands.Point_extra[0];
   Results<<" ";
   if(Input_commands.Point_extra[1]>=ZERO){Results<<" ";}
   Results<<Input_commands.Point_extra[1];
   Results<<" ";
   if(Input_commands.Point_extra[2]>=ZERO){Results<<" ";}
   Results<<Input_commands.Point_extra[2]<<endl;
   Extracule=ZERO;
   integrate_extra(dmn,method,3,1,error_rel,error_abs,MIN_EVALS,MAX_EVALS,Extracule,fail,
   Input_commands.nprocs,Input_commands.Point_extra);
   Results<<"E(X12,Y12,Z12)     : "<<setw(17)<<HALF*Extracule<<endl;
   Results<<endl;
   Results<<"\t \t Obtained with "<<method<<", fail report "<<fail<<endl;
   Results<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ///////////////////////////////////
  //  Potential V(r) integrated    //
  ///////////////////////////////////
  if(Input_commands.Vr)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#           Use the density to evaluate the Potential V(r)                #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   int fail=1000;
   double res_integration_fchk[2]={ZERO},Vn,Ve;
   Results<<setprecision(10)<<fixed<<scientific;
   integrate_vr_fchk(Read_fchk_wfn,Input_commands.method_cuba,3,2,Input_commands.error_rel,Input_commands.error_abs,
   Input_commands.minevals,Input_commands.maxevals,res_integration_fchk,fail,Input_commands.nprocs,Input_commands.Point_Vr);
   Vn=Read_fchk_wfn.Vnuclear(Input_commands.Point_Vr);
   Ve=res_integration_fchk[1];
   Results<<" N  =  Int Rho(r) dr       = "<<setw(17)<<res_integration_fchk[0]<<endl;
   Results<<" Ve(r)                     = "<<setw(17)<<Ve<<endl;
   Results<<"[Ve(r) = Int Rho(r')/|r-r'| dr']"<<endl;
   Results<<" Vn(r)                     = "<<setw(17)<<Vn<<endl;
   Results<<"[Vn(r) = sum_i ^Natoms Zi/|r-Ri|]"<<endl;
   Results<<" Vn(r) - Ve(r)             = "<<setw(17)<<Vn-Ve<<endl;
   Results<<"    x                 y                 z "<<endl;
   Results<<setw(17)<<Input_commands.Point_Vr[0]<<" "<<setw(17)<<Input_commands.Point_Vr[1];
   Results<<" "<<setw(17)<<Input_commands.Point_Vr[2]<<endl;
   Results<<endl;
   Results<<endl;
   Results<<" obtained with "<<Input_commands.method_cuba<<" fail report: "<<fail<<endl;
   Results<<endl;
   fail=1000;
   Results<<"(fail 0 means finish without problem, -1 error, >1 accuracy";
   Results<<" not achived and 10000 not calculated)"<<endl;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ////////////////////////
  // Plot using gnuplot //
  ////////////////////////
  if(Input_commands.gnuplot)
  {
   bool exchange,dens,densp,lhybds_plot;
   double tw_div_t,t2,q_red_total,s_red_total,var1,var2;
   ofstream files_data,files_data2,files_data3,files_data4,files_data5,files_data6;
   string data,name_gnuplot_commands,name_gnuplot_commands2,name_eps;
   fish=false;fishp=false;shan=false;shanp=false;
   dens=false;densp=false;lhybds_plot=false;exchange=false;
   if(Read_fchk_wfn.wfn)
   {
    name_gnuplot_commands=name_file.substr(0,(name_file.length()-4));
   }
   else
   {
    name_gnuplot_commands=name_file.substr(0,(name_file.length()-5));
   }
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#                        Plotting Options                                 #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   Results<<"See files:"<<endl;
   Results<<endl;
   for(i=0;i<Input_commands.num_plot_ops;i++)
   {
    if(Input_commands.plot_ops[i]=="density"){dens=true;}
    if(Input_commands.plot_ops[i]=="shannon"){dens=true;shan=true;}
    if(Input_commands.plot_ops[i]=="fisher"){dens=true;fish=true;}
    if(Input_commands.plot_ops[i]=="local_hybrids")
    {
     lhybds_plot=true;
     dens=false;
     shan=false;
     fish=false;
    }
    if(Input_commands.plot_ops[i]=="densityp"){densp=true;}
    if(Input_commands.plot_ops[i]=="shannonp"){densp=true;shanp=true;}
    if(Input_commands.plot_ops[i]=="fisherp"){densp=true;fishp=true;}
   }
   if(Input_commands.dim2)
   {
    GNUPLOT plot;
    if(!Input_commands.nopath)
    {plot.path=Input_commands.path;}
    if(dens)
    {
     files_data.open((name_gnuplot_commands+"_density.dat").c_str());
     if(shan)
     {files_data2.open((name_gnuplot_commands+"_shannon.dat").c_str());}
     if(fish)
     {files_data3.open((name_gnuplot_commands+"_fisher.dat").c_str());}
     if(Input_commands.dmn_plots)
     {
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       if(Input_commands.scan_1=='x' || Input_commands.scan_1=='X')
       {
        Point[1]=ZERO;
        Point[2]=ZERO;
        Point[0]=var1;
       }
       if(Input_commands.scan_1=='y' || Input_commands.scan_1=='Y')
       {
        Point[0]=ZERO;
        Point[2]=ZERO;
        Point[1]=var1;
       }
       if(Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')
       {
        Point[0]=ZERO;
        Point[1]=ZERO;
        Point[2]=var1;
       }
       Density=dmn.evaluation(Point,Point);
       files_data<<var1<<" "<<Density<<endl;
       if(shan)
       {
        if(Density>=pow(TEN,-TWO*TEN))
        {files_data2<<var1<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
        else{files_data2<<var1<<" "<<Density<<endl;}
       }
       if(fish)
       {
        if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
        dmn.grad_rho_r(Point,Grad,Density);
        files_data3<<var1<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
       }
      }
     }
     else
     {
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       if(Input_commands.scan_1=='x' || Input_commands.scan_1=='X')
       {
        Point[1]=ZERO;
        Point[2]=ZERO;
        Point[0]=var1;
       }
       if(Input_commands.scan_1=='y' || Input_commands.scan_1=='Y')
       {
        Point[0]=ZERO;
        Point[2]=ZERO;
        Point[1]=var1;
       }
       if(Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')
       {
        Point[0]=ZERO;
        Point[1]=ZERO;
        Point[2]=var1;
       }
       Read_fchk_wfn.rho_eval(Point,Density);
       files_data<<var1<<" "<<Density<<endl;
       if(shan)
       {
        if(Density>=pow(TEN,-TWO*TEN))
        {files_data2<<var1<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
        else{files_data2<<var1<<" "<<Density<<endl;}
       }
       if(fish)
       {
        if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
        Read_fchk_wfn.rho_grad(Point,Grad);
        files_data3<<var1<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
       }
      }
     }
     files_data.close();
     if(shan)
     {files_data2.close();}
     if(fish)
     {files_data3.close();}
    }
    if(densp)
    {
     files_data.open((name_gnuplot_commands+"_densityp.dat").c_str());
     if(shanp)
     {files_data2.open((name_gnuplot_commands+"_shannonp.dat").c_str());}
     if(fishp)
     {files_data3.open((name_gnuplot_commands+"_fisherp.dat").c_str());}
     if(Input_commands.dmn_plots)
     {
      DMN_P_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       if(Input_commands.scan_1=='x' || Input_commands.scan_1=='X')
       {
        Point[1]=ZERO;
        Point[2]=ZERO;
        Point[0]=var1;
       }
       if(Input_commands.scan_1=='y' || Input_commands.scan_1=='Y')
       {
        Point[0]=ZERO;
        Point[2]=ZERO;
        Point[1]=var1;
       }
       if(Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')
       {
        Point[0]=ZERO;
        Point[1]=ZERO;
        Point[2]=var1;
       }
       Density=dmn.evaluation_p(Point,Point);
       files_data<<var1<<" "<<Density<<endl;
       if(shanp)
       {
        if(Density>=pow(TEN,-TWO*TEN))
        {files_data2<<var1<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
        else{files_data2<<var1<<" "<<Density<<endl;}
       }
       if(fishp)
       {
        if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
        dmn.grad_rho_p(Point,Grad,Density);
        files_data3<<var1<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
       }
      }
     }
     else
     {
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       if(Input_commands.scan_1=='x' || Input_commands.scan_1=='X')
       {
        Point[1]=ZERO;
        Point[2]=ZERO;
        Point[0]=var1;
       }
       if(Input_commands.scan_1=='y' || Input_commands.scan_1=='Y')
       {
        Point[0]=ZERO;
        Point[2]=ZERO;
        Point[1]=var1;
       }
       if(Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')
       {
        Point[0]=ZERO;
        Point[1]=ZERO;
        Point[2]=var1;
       }
       Read_fchk_wfn.rho_p_eval(Point,Density);
       files_data<<var1<<" "<<Density<<endl;
       if(shanp)
       {
        if(Density>=pow(TEN,-TWO*TEN))
        {files_data2<<var1<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
        else{files_data2<<var1<<" "<<Density<<endl;}
       }
       if(fishp)
       {
        if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
        Read_fchk_wfn.rho_p_grad(Point,Grad);
        files_data3<<var1<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
       }
      }
     }
    }
    files_data.close();
    if(shanp)
    {files_data2.close();}
    if(fishp)
    {files_data3.close();}
    if(lhybds_plot)
    {
     files_data.open((name_gnuplot_commands+"_density.dat").c_str());
     files_data2.open((name_gnuplot_commands+"_twt.dat").c_str());
     files_data3.open((name_gnuplot_commands+"_t2.dat").c_str());
     files_data4.open((name_gnuplot_commands+"_dori.dat").c_str());
     files_data5.open((name_gnuplot_commands+"_grad_red.dat").c_str());
     files_data6.open((name_gnuplot_commands+"_erf.dat").c_str());
     for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
     {
      if(Input_commands.scan_1=='x' || Input_commands.scan_1=='X')
      {
       Point[1]=ZERO;
       Point[2]=ZERO;
       Point[0]=var1;
      }
      if(Input_commands.scan_1=='y' || Input_commands.scan_1=='Y')
      {
       Point[0]=ZERO;
       Point[2]=ZERO;
       Point[1]=var1;
      }
      if(Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')
      {
       Point[0]=ZERO;
       Point[1]=ZERO;
       Point[2]=var1;
      }
      punctual(Read_fchk_wfn,Point,Density,Density_alpha,Density_beta,Grad,Grad_alpha,Grad_beta,Grad_norm,
      Grad_norm_alpha,Grad_norm_beta,laplacian_r,laplacian_alpha,laplacian_beta,tauW_alpha,
      tauW_beta,tau_alpha,tau_beta,k_F_alpha,k_F_beta,k_s_alpha,k_s_beta,s_r_alpha,s_r_beta,
      q_red_alpha,q_red_beta,wfn_fchk);
      tw_div_t=(tauW_alpha+tauW_beta)/(tau_alpha+tau_beta+pow(TEN,-TWO*TEN));
      t2=pow(PI/THREE,ONE_THIRD)*a_o*pow(Grad_norm,TWO)/(TWO*EIGHT*pow(Density,SEVEN*ONE_THIRD)+pow(TEN,-TWO*TEN));
      q_red_total=q_red_alpha+q_red_beta;
      s_red_total=s_r_alpha+s_r_beta;
      DORI=FOUR*(ONE-TWO*q_red_total/(s_red_total*s_red_total+pow(TEN,-TWO*TEN))+q_red_total*q_red_total/(pow(s_red_total,FOUR)+pow(TEN,-TWO*TEN)));
      // See  JCP, 142, 074112 (2015)
      //      JCP, 140, 18A510 (2014)
      //      JCP, 131, 154112 (2009)
      //      JPCA, 113, 11898 (2009)
      LOCAL_HYBRIDS_fr[0]=tw_div_t;
      LOCAL_HYBRIDS_fr[1]=ONE/(ONE+HALF*t2+pow(TEN,-TWO*TEN));
      LOCAL_HYBRIDS_fr[2]=0.676/(ONE+8.908*DORI);
      LOCAL_HYBRIDS_fr[3]=s_red_total*s_red_total/(pow(0.73+s_red_total,TWO)+pow(TEN,-TWO*TEN));
      LOCAL_HYBRIDS_fr[4]=erf(0.20*s_red_total);
      files_data<<var1<<" "<<Density<<endl;
      files_data2<<var1<<" "<<LOCAL_HYBRIDS_fr[0]<<endl;
      files_data3<<var1<<" "<<LOCAL_HYBRIDS_fr[1]<<endl;
      files_data4<<var1<<" "<<LOCAL_HYBRIDS_fr[2]<<endl;
      files_data5<<var1<<" "<<LOCAL_HYBRIDS_fr[3]<<endl;
      files_data6<<var1<<" "<<LOCAL_HYBRIDS_fr[4]<<endl;
     }
     files_data.close();
     files_data2.close();
     files_data3.close();
     files_data4.close();
     files_data5.close();
     files_data6.close();
    }
    //Send for plotting (except for Local hybrids)
    for(i=0;i<Input_commands.num_plot_ops;i++)
    {
     if(Input_commands.plot_ops[i]=="density" || Input_commands.plot_ops[i]=="local_hybrids")
     {
      if(Input_commands.plot_ops[i]=="local_hybrids")
      {
       Input_commands.plot_ops[i]="density";
       cout<<"Notice that the .eps and the .plot2 files were only generated for the density"<<endl;
      }
      operation="{/Symbol r}(r)";
      data=name_gnuplot_commands+"_density.dat";
     }
     if(Input_commands.plot_ops[i]=="shannon" && shan)
     {
      operation="S[{/Symbol r}(r)]";
      data=name_gnuplot_commands+"_shannon.dat";
     }
     if(Input_commands.plot_ops[i]=="fisher"  && fish)
     {
      operation="I[{/Symbol r}(r)]";
      data=name_gnuplot_commands+"_fisher.dat";
     }
     if(Input_commands.plot_ops[i]=="densityp")
     {
      operation="{/Symbol p}(p)";
      data=name_gnuplot_commands+"_densityp.dat";
     }
     if(Input_commands.plot_ops[i]=="shannonp" && shanp)
     {
      operation="S[{/Symbol p}(p)]";
      data=name_gnuplot_commands+"_shannonp.dat";
     }
     if(Input_commands.plot_ops[i]=="fisherp" && fishp)
     {
      operation="I[{/Symbol p}(p)]";
      data=name_gnuplot_commands+"_fisherp.dat";
     }
     name_eps=name_gnuplot_commands+"_"+Input_commands.plot_ops[i];
     name_gnuplot_commands2=name_gnuplot_commands+"_"+Input_commands.plot_ops[i];
     plot.plot2D(data,name_gnuplot_commands2,name_eps,operation,Input_commands.extra_lines_plot,Input_commands.extra_lines,
     Input_commands.scan_1,Input_commands.points_scan_1);
     Results<<"    "+name_eps+".eps"<<endl;
    }
   }
   if(Input_commands.dim3)
   {
    GNUPLOT plot;
    if(!Input_commands.nopath)
    {plot.path=Input_commands.path;}
    char not_used,aux;
    if(((Input_commands.scan_1=='x' || Input_commands.scan_1=='X') && (Input_commands.scan_2=='y' || Input_commands.scan_2=='Y')) ||
      ((Input_commands.scan_2=='x' || Input_commands.scan_2=='X') && (Input_commands.scan_1=='y' || Input_commands.scan_1=='Y')))
       {
        not_used='z';
        Point[2]=Input_commands.not_used_3d_doub;
        if(!(Input_commands.scan_1=='x' || Input_commands.scan_1=='X'))
        {
         exchange=true;
         aux=Input_commands.scan_1;
         Input_commands.scan_1=Input_commands.scan_2;
         Input_commands.scan_2=aux;
        }
       }
    if(((Input_commands.scan_1=='x' || Input_commands.scan_1=='X') && (Input_commands.scan_2=='z' || Input_commands.scan_2=='Z')) ||
       ((Input_commands.scan_2=='x' || Input_commands.scan_2=='X') && (Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')))
       {
        not_used='y';
        Point[1]=Input_commands.not_used_3d_doub;
        if(!(Input_commands.scan_1=='x' || Input_commands.scan_1=='X'))
        {
         exchange=true;
         aux=Input_commands.scan_1;
         Input_commands.scan_1=Input_commands.scan_2;
         Input_commands.scan_2=aux;
        }
       }
     if((((Input_commands.scan_1=='y' || Input_commands.scan_1=='Y') && (Input_commands.scan_2=='z' || Input_commands.scan_2=='Z'))) ||
        ((Input_commands.scan_2=='y' || Input_commands.scan_2=='Y') && (Input_commands.scan_1=='z' || Input_commands.scan_1=='Z')))
       {
        not_used='x';
        Point[0]=Input_commands.not_used_3d_doub;
        if(!(Input_commands.scan_1=='y' || Input_commands.scan_1=='Y'))
        {
         exchange=true;
         aux=Input_commands.scan_1;
         Input_commands.scan_1=Input_commands.scan_2;
         Input_commands.scan_2=aux;
        }
       }
    //Leave one space when var1 changes value for scan in table
    if(dens)
    {
     files_data.open((name_gnuplot_commands+"_density.dat").c_str());
     if(shan)
     {files_data2.open((name_gnuplot_commands+"_shannon.dat").c_str());}
     if(fish)
     {files_data3.open((name_gnuplot_commands+"_fisher.dat").c_str());}
     if(Input_commands.dmn_plots)
     {
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       for(var2=Input_commands.points_scan_1[0];var2<=Input_commands.points_scan_1[1];var2=var2+Input_commands.grid_1)
       {
        if(not_used=='z')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[1]=var1;
         }
         else
         {
          Point[0]=var1;
          Point[1]=var2;
         }
         Density=dmn.evaluation(Point,Point);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shan)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
          }
         if(fish)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          dmn.grad_rho_r(Point,Grad,Density);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='y')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[2]=var1;
         }
         else
         {
          Point[0]=var1;
          Point[2]=var2;
         }
         Density=dmn.evaluation(Point,Point);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shan)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fish)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          dmn.grad_rho_r(Point,Grad,Density);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='x')
        {
         if(exchange)
         {
          Point[1]=var2;
          Point[2]=var1;
         }
         else
         {
          Point[1]=var1;
          Point[2]=var2;
         }
         Density=dmn.evaluation(Point,Point);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shan)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fish)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          dmn.grad_rho_r(Point,Grad,Density);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
       }
       files_data<<endl;
       if(shan){files_data2<<endl;}
       if(fish){files_data3<<endl;}
      }
     }
     else
     {
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       for(var2=Input_commands.points_scan_1[0];var2<=Input_commands.points_scan_1[1];var2=var2+Input_commands.grid_1)
       {
        if(not_used=='z')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[1]=var1;
          aux=Input_commands.scan_1;
          Input_commands.scan_1=Input_commands.scan_2;
          Input_commands.scan_2=aux;
         }
         else
         {
          Point[0]=var1;
          Point[1]=var2;
         }
         Read_fchk_wfn.rho_eval(Point,Density);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shan)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
          }
         if(fish)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          Read_fchk_wfn.rho_grad(Point,Grad);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='y')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[2]=var1;
          aux=Input_commands.scan_1;
          Input_commands.scan_1=Input_commands.scan_2;
          Input_commands.scan_2=aux;
         }
         else
         {
          Point[0]=var1;
          Point[2]=var2;
         }
         Read_fchk_wfn.rho_eval(Point,Density);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shan)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fish)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          Read_fchk_wfn.rho_grad(Point,Grad);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='x')
        {
         if(exchange)
         {
          Point[1]=var2;
          Point[2]=var1;
          aux=Input_commands.scan_1;
          Input_commands.scan_1=Input_commands.scan_2;
          Input_commands.scan_2=aux;
         }
         else
         {
          Point[1]=var1;
          Point[2]=var2;
         }
         Read_fchk_wfn.rho_eval(Point,Density);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shan)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fish)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          Read_fchk_wfn.rho_grad(Point,Grad);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
       }
       files_data<<endl;
       if(shan){files_data2<<endl;}
       if(fish){files_data3<<endl;}
      }
     }
     files_data.close();
     if(shan)
     {files_data2.close();}
     if(fish)
     {files_data3.close();}
    }
    if(densp)
    {
     files_data.open((name_gnuplot_commands+"_densityp.dat").c_str());
     if(shanp)
     {files_data2.open((name_gnuplot_commands+"_shannonp.dat").c_str());}
     if(fishp)
     {files_data3.open((name_gnuplot_commands+"_fisherp.dat").c_str());}
     if(Input_commands.dmn_plots)
     {
      DMN_P_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       for(var2=Input_commands.points_scan_1[0];var2<=Input_commands.points_scan_1[1];var2=var2+Input_commands.grid_1)
       {
        if(not_used=='z')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[1]=var1;
         }
         else
         {
          Point[0]=var1;
          Point[1]=var2;
         }
         Density=dmn.evaluation_p(Point,Point);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shanp)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fishp)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          dmn.grad_rho_p(Point,Grad,Density);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='y')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[2]=var1;
         }
         else
         {
          Point[0]=var1;
          Point[2]=var2;
         }
         Density=dmn.evaluation_p(Point,Point);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shanp)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fishp)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          dmn.grad_rho_p(Point,Grad,Density);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='x')
        {
         if(exchange)
         {
          Point[1]=var2;
          Point[2]=var1;
         }
         else
         {
          Point[1]=var1;
          Point[2]=var2;
         }
         Density=dmn.evaluation_p(Point,Point);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shanp)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fishp)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          dmn.grad_rho_p(Point,Grad,Density);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
       }
       files_data<<endl;
       if(shanp){files_data2<<endl;}
       if(fishp){files_data3<<endl;}
      }
     }
     else
     {
      for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
      {
       for(var2=Input_commands.points_scan_1[0];var2<=Input_commands.points_scan_1[1];var2=var2+Input_commands.grid_1)
       {
        if(not_used=='z')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[1]=var1;
         }
         else
         {
          Point[0]=var1;
          Point[1]=var2;
         }
         Read_fchk_wfn.rho_p_eval(Point,Density);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shanp)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fishp)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          Read_fchk_wfn.rho_p_grad(Point,Grad);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='y')
        {
         if(exchange)
         {
          Point[0]=var2;
          Point[2]=var1;
         }
         else
         {
          Point[0]=var1;
          Point[2]=var2;
         }
         Read_fchk_wfn.rho_p_eval(Point,Density);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shanp)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fishp)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          Read_fchk_wfn.rho_p_grad(Point,Grad);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
        if(not_used=='x')
        {
         if(exchange)
         {
          Point[1]=var2;
          Point[2]=var1;
         }
         else
         {
          Point[1]=var1;
          Point[2]=var2;
         }
         Read_fchk_wfn.rho_p_eval(Point,Density);
         files_data<<var1<<" "<<var2<<" "<<Density<<endl;
         if(shanp)
         {
          if(Density>=pow(TEN,-TWO*TEN))
          {files_data2<<var1<<" "<<var2<<" "<<-(Density/(double)Read_fchk_wfn.nelectrons)*log((Density/(double)Read_fchk_wfn.nelectrons))<<endl;}
          else{files_data2<<var1<<" "<<var2<<" "<<Density<<endl;}
         }
         if(fishp)
         {
          if(Density==ZERO){Density=pow(TEN,-THREE*TEN);}
          Read_fchk_wfn.rho_p_grad(Point,Grad);
          files_data3<<var1<<" "<<var2<<" "<<pow(norm3D(Grad),TWO)/(Density*(double)Read_fchk_wfn.nelectrons)<<endl;
         }
        }
       }
       files_data<<endl;
       if(shanp){files_data2<<endl;}
       if(fishp){files_data3<<endl;}
      }
     }
     files_data.close();
     if(shanp)
     {files_data2.close();}
     if(fishp)
     {files_data3.close();}
    }
    if(lhybds_plot)
    {
     files_data.open((name_gnuplot_commands+"_density.dat").c_str());
     files_data2.open((name_gnuplot_commands+"_twt.dat").c_str());
     files_data3.open((name_gnuplot_commands+"_t2.dat").c_str());
     files_data4.open((name_gnuplot_commands+"_dori.dat").c_str());
     files_data5.open((name_gnuplot_commands+"_grad_red.dat").c_str());
     files_data6.open((name_gnuplot_commands+"_erf.dat").c_str());
     for(var1=Input_commands.points_scan_1[0];var1<=Input_commands.points_scan_1[1];var1=var1+Input_commands.grid_1)
     {
      for(var2=Input_commands.points_scan_1[0];var2<=Input_commands.points_scan_1[1];var2=var2+Input_commands.grid_1)
      {
       if(not_used=='z')
       {
        if(exchange)
        {
         Point[0]=var2;
         Point[1]=var1;
        }
        else
        {
         Point[0]=var1;
         Point[1]=var2;
        }
       }
       if(not_used=='y')
       {
        if(exchange)
        {
         Point[0]=var2;
         Point[2]=var1;
        }
        else
        {
         Point[0]=var1;
         Point[2]=var2;
        }
       }
       if(not_used=='x')
       {
        if(exchange)
        {
         Point[1]=var2;
         Point[2]=var1;
        }
        else
        {
         Point[1]=var1;
         Point[2]=var2;
        }
       }
       punctual(Read_fchk_wfn,Point,Density,Density_alpha,Density_beta,Grad,Grad_alpha,Grad_beta,Grad_norm,
       Grad_norm_alpha,Grad_norm_beta,laplacian_r,laplacian_alpha,laplacian_beta,tauW_alpha,
       tauW_beta,tau_alpha,tau_beta,k_F_alpha,k_F_beta,k_s_alpha,k_s_beta,s_r_alpha,s_r_beta,
       q_red_alpha,q_red_beta,wfn_fchk);
       tw_div_t=(tauW_alpha+tauW_beta)/(tau_alpha+tau_beta+pow(TEN,-TWO*TEN));
       t2=pow(PI/THREE,ONE_THIRD)*a_o*pow(Grad_norm,TWO)/(TWO*EIGHT*pow(Density,SEVEN*ONE_THIRD)+pow(TEN,-TWO*TEN));
       q_red_total=q_red_alpha+q_red_beta;
       s_red_total=s_r_alpha+s_r_beta;
       DORI=FOUR*(ONE-TWO*q_red_total/(s_red_total*s_red_total+pow(TEN,-TWO*TEN))+q_red_total*q_red_total/(pow(s_red_total,FOUR)+pow(TEN,-TWO*TEN)));
       // See  JCP, 142, 074112 (2015)
       //      JCP, 140, 18A510 (2014)
       //      JCP, 131, 154112 (2009)
       //      JPCA, 113, 11898 (2009)
       LOCAL_HYBRIDS_fr[0]=tw_div_t;
       LOCAL_HYBRIDS_fr[1]=ONE/(ONE+HALF*t2+pow(TEN,-TWO*TEN));
       LOCAL_HYBRIDS_fr[2]=0.676/(ONE+8.908*DORI);
       LOCAL_HYBRIDS_fr[3]=s_red_total*s_red_total/(pow(0.73+s_red_total,TWO)+pow(TEN,-TWO*TEN));
       LOCAL_HYBRIDS_fr[4]=erf(0.20*s_red_total);
       files_data<<var1<<" "<<var2<<" "<<Density<<endl;
       files_data2<<var1<<" "<<var2<<" "<<LOCAL_HYBRIDS_fr[0]<<endl;
       files_data3<<var1<<" "<<var2<<" "<<LOCAL_HYBRIDS_fr[1]<<endl;
       files_data4<<var1<<" "<<var2<<" "<<LOCAL_HYBRIDS_fr[2]<<endl;
       files_data5<<var1<<" "<<var2<<" "<<LOCAL_HYBRIDS_fr[3]<<endl;
       files_data6<<var1<<" "<<var2<<" "<<LOCAL_HYBRIDS_fr[4]<<endl;
      }
      files_data<<endl;
      files_data2<<endl;
      files_data3<<endl;
      files_data4<<endl;
      files_data5<<endl;
      files_data6<<endl;
     }
     files_data.close();
     files_data2.close();
     files_data3.close();
     files_data4.close();
     files_data5.close();
     files_data6.close();
    }
    //Send for plotting
    for(i=0;i<Input_commands.num_plot_ops;i++)
    {
     if(Input_commands.plot_ops[i]=="density" || Input_commands.plot_ops[i]=="local_hybrids")
     {
      if(Input_commands.plot_ops[i]=="local_hybrids")
      {
       Input_commands.plot_ops[i]="density";
       cout<<"Notice that the .eps and the .plot3 files were only generated for the density"<<endl;
      }
      operation="{/Symbol r}(r)";
      data=name_gnuplot_commands+"_density.dat";
     }
     if(Input_commands.plot_ops[i]=="shannon")
     {
      operation="S[{/Symbol r}(r)]";
      data=name_gnuplot_commands+"_shannon.dat";
     }
     if(Input_commands.plot_ops[i]=="fisher")
     {
      operation="I[{/Symbol r}(r)]";
      data=name_gnuplot_commands+"_fisher.dat";
     }
     if(Input_commands.plot_ops[i]=="densityp")
     {
      operation="{/Symbol p}(p)";
      data=name_gnuplot_commands+"_densityp.dat";
     }
     if(Input_commands.plot_ops[i]=="shannonp")
     {
      operation="S[{/Symbol p}(p)]";
      data=name_gnuplot_commands+"_shannonp.dat";
     }
     if(Input_commands.plot_ops[i]=="fisherp")
     {
      operation="I[{/Symbol p}(p)]";
      data=name_gnuplot_commands+"_fisherp.dat";
     }
     name_eps=name_gnuplot_commands+"_"+Input_commands.plot_ops[i];
     name_gnuplot_commands2=name_gnuplot_commands+"_"+Input_commands.plot_ops[i];
     plot.plot3D(data,name_gnuplot_commands2,name_eps,operation,Input_commands.extra_lines_plot,Input_commands.extra_lines,
     Input_commands.scan_1,Input_commands.scan_2,Input_commands.points_scan_1,Input_commands.points_scan_2);
     Results<<"    "+name_eps+".eps"<<endl;
     Results<<"    "+name_eps+"_map.eps"<<endl;
    }
   }
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  /////////////////////////////////////
  // Perform a MESCAL calculation    //
  /////////////////////////////////////
  if(Input_commands.mescal)
  {
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#    MESCAL response to the QM electronic density (F_elec) + F_nuclei     #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   string mescal_file;
   if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
   {
    mescal_file=name_file.substr(0,(name_file.length()-4))+"_MESCAL.out";
   }
   else
   {
    mescal_file=name_file.substr(0,(name_file.length()-5))+"_MESCAL.out";
   }
   MESCAL mescal(mescal_file,Input_commands.mescal_pdb);
   // Init F_ext (here we could have an 'else if' to send info QM -> MM integrating the density)
   if(Input_commands.mescal_punctual)
   {
    double Point_mescal[3];
    for(i=0;i<Input_commands.npoints_mescal;i++)
    {
     for(j=0;j<3;j++){Point_mescal[j]=Input_commands.Point_mescal[i][j];}
     mescal.set_F_ext_punct(Input_commands.q_mescal,Point_mescal);
    }
   }
   else
   {
    Results<<"Comment: No external punctual charge read. Thus, no F_ext employed in MESCAL"<<endl;
   }

   // call mescal SCF for mu.

   mescal.close_output(mescal_file);
   Results<<endl;
  }
  /////////////////////////////////////
  // Transform int files for ESI-3c  //
  /////////////////////////////////////
  if(Input_commands.esi_int)
  {
   nbasis=Read_fchk_wfn.nbasis();
   transform_int(Input_commands.Sij_region,nbasis);
  }
  ///////////////////////////
  // Mulliken analysis     //
  ///////////////////////////
  if(Input_commands.mulliken && !wfn_fchk)
  {
   int k;
   double **P_matrix,**NO_coef,**Mulliken_matrix,*Gross_occ,*Nu_charges;
   bool more=true;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<"#                    Mulliken Population Analysis                         #";
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
   Results<<endl;
   nbasis=Read_fchk_wfn.nbasisf;
   P_matrix=new double*[nbasis];
   NO_coef=new double*[nbasis];
   Mulliken_matrix=new double*[nbasis];
   Gross_occ=new double[nbasis];
   Nu_charges=new double[Read_fchk_wfn.natoms];
   for(i=0;i<Read_fchk_wfn.natoms;i++)
   {
    Nu_charges[i]=Read_fchk_wfn.Nu_charge[i];
   }
   for(i=0;i<nbasis;i++)
   {
    P_matrix[i]=new double[nbasis];
    NO_coef[i]=new double[nbasis];
    Mulliken_matrix[i]=new double[nbasis];
   }
   for(i=0;i<nbasis;i++)
   {
    for(j=0;j<nbasis;j++)
    {
     Mulliken_matrix[i][j]=ZERO;
     NO_coef[i][j]=Read_fchk_wfn.S[j][i];
     if(j<=i)
     {
      P_matrix[i][j]=Read_fchk_wfn.P[i][j];
      P_matrix[j][i]=Read_fchk_wfn.P[i][j];
     }
    }  
   }
   for(i=0;i<nbasis;i++)
   {
    for(j=0;j<nbasis;j++)
    {
     for(k=0;k<nbasis;k++)
     {
      if(i==j)
      {
       Mulliken_matrix[i][i]=Mulliken_matrix[i][i]+P_matrix[k][k]*pow(NO_coef[k][i],TWO);
      }
      else
      {
       Mulliken_matrix[i][j]=Mulliken_matrix[i][j]+P_matrix[k][k]*NO_coef[k][i]*NO_coef[k][j]*Read_fchk_wfn.Sao[i][j];
      }
     }
    }
   }
   Results<<endl;
   Results<<" Full Mulliken population analysis:"<<endl;
   Results<<endl;
   Results<<setprecision(14)<<fixed<<scientific;
   k=0;
   do
   {
    for(i=0+k;i<nbasis;i++)
    {
     for(j=0+k;j<=i;j++)
     {
      if(j<k+5)
      {
       if(Mulliken_matrix[i][j]>=ZERO){Results<<"  ";}
       else{Results<<" ";}
       Results<<Mulliken_matrix[i][j];
       more=false;
      }
      else
      {
       j=i+1;
       more=true;
      }
     }
     Results<<endl;
    }
    Results<<endl;
    k=k+5;
   }while(more);
   Density=ZERO;
   for(i=0;i<nbasis;i++)
   {
    Gross_occ[i]=ZERO;
    for(j=0;j<nbasis;j++)
    {
     Gross_occ[i]=Gross_occ[i]+Mulliken_matrix[j][i];
    }
    if(abs(Gross_occ[i])<pow(TEN,-TEN)){Gross_occ[i]=ZERO;}
    Density=Density+Gross_occ[i];
   }
   Results<<endl;
   Results<<" Gross orbital populations:"<<endl;
   Results<<endl;
   ifstream log_file(Input_commands.name_log.c_str());
   more=true;
   while(getline(log_file,line) && more)
   {
    if(line=="     Gross orbital populations:")
    {
     more=false;
    }
   }
   j=0;
   for(i=0;i<nbasis;i++)
   {
    if(!more)
    {
     getline(log_file,line);
     line=line.substr(0,23);
     if(line[9]!=' ' && i!=0){j++;}
    }
    else
    {
     line="";
    }
    Nu_charges[j]=Nu_charges[j]-Gross_occ[i];
    Results<<line<<setw(21)<<Gross_occ[i]<<endl;
   }
   if(!more)
   {
    Results<<endl;
    while(getline(log_file,line))
    {
     if(line==" Mulliken atomic charges:" || line==" Mulliken charges:")
     {
      Results<<line<<endl;
      Results<<endl;
      getline(log_file,line);
      for(i=0;i<Read_fchk_wfn.natoms;i++)
      {
       getline(log_file,line);
       Results<<line.substr(0,11)<<setw(21)<<Nu_charges[i]<<endl;
      }
     }
    }
   }
   log_file.close();
   Results<<endl;
   Results<<endl;
   Results<<" N electrons:"<<setprecision(6)<<fixed<<setw(12)<<Density<<endl;
   Results<<endl;
   for(i=0;i<nbasis;i++)
   {
    delete[] P_matrix[i];P_matrix[i]=NULL;
    delete[] NO_coef[i];NO_coef[i]=NULL;
    delete[] Mulliken_matrix[i];Mulliken_matrix[i]=NULL;
   }
   delete[] P_matrix;P_matrix=NULL;
   delete[] NO_coef;NO_coef=NULL;
   delete[] Mulliken_matrix;Mulliken_matrix=NULL;
   delete[] Gross_occ;Gross_occ=NULL;
   delete[] Nu_charges;Nu_charges=NULL;
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  ///////////////////////////////
  // Print a DM1 fchk file     //
  ///////////////////////////////
  if((Input_commands.print_dm1_fchk && !wfn_fchk) && Input_commands.store_dmn)
  {
   int i,aos,counteraos=0;
   double *tot_density;
   string aux;
   ofstream dm1fchk;
   ifstream read_hf_fchk,formatted_data;
   read_hf_fchk.open((name_file).c_str());
   dm1fchk.open((name_file.substr(0,(name_file.length()-5))+"_dm1.fchk").c_str());
   Results<<endl;
   Results<<"See file "<<name_file.substr(0,(name_file.length()-5))+"_dm1.fchk"<<endl;
   Results<<"Note: 'Total SCF Density' will contain the FCI density obtained using the DM1 matrix"<<endl;
   if(Read_fchk_wfn.multiplicity!=1)
   {
    Results<<"Note2: 'Spin SCF Density' will contain the FCI spin density obtained using the DM1 matrix"<<endl;
   }
   while(getline(read_hf_fchk,aux))
   {
    if(aux.length()>17)
    {
     if(aux.substr(0,17)=="Total SCF Density")
     {
      dm1fchk<<aux<<endl;
      aux=aux.substr(56,aux.length());
      stringstream ss(aux); //string to int in ss
      ss>>aos;
      tot_density=new double [aos];
      for(i=0;i<aos;i++)
      {
       read_hf_fchk>>tot_density[i];
       tot_density[i]=ZERO;
      }
      cout<<"See file "<<name_file.substr(0,(name_file.length()-5))+"_dm1.fchk"<<endl;
      DMN_OPS dmn(Input_commands.name_dm1,1);
      dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
      dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
      if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
      dmn.tot_dens_fchk(tot_density,aos);
      FILE *pFile;
      pFile=fopen((name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str(),"w");
      for(i=0;i<aos;i++)
      {
       if(tot_density[i]>=ZERO)
       {
        fprintf(pFile,"  %14.8E",tot_density[i]);
       }
       else
       {
        fprintf(pFile," %14.8E",tot_density[i]);
       }
       counteraos++;
       if(counteraos==5)
       {
        fprintf(pFile,"\n");
        counteraos=0;
       }
      }
      fclose(pFile);
      formatted_data.open((name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str());
      while(getline(formatted_data,aux))
      {
       dm1fchk<<aux<<endl;
      }
      formatted_data.close();
      system(("rm "+name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str());
      if(Read_fchk_wfn.multiplicity!=1)
      {
       getline(read_hf_fchk,aux);
       dm1fchk<<aux<<endl;
       while(aux.substr(0,16)!="Spin SCF Density")
       {
        getline(read_hf_fchk,aux);
        dm1fchk<<aux<<endl;
       }
       for(i=0;i<aos;i++)
       {
        read_hf_fchk>>tot_density[i];
        tot_density[i]=ZERO;
       }
       dmn.spin_dens_fchk(tot_density,aos);
       FILE *pFile;
       pFile=fopen((name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str(),"w");
       counteraos=0;
       for(i=0;i<aos;i++)
       {
        if(tot_density[i]>=ZERO)
        {
         fprintf(pFile,"  %14.8E",tot_density[i]);
        }
        else
        {
         fprintf(pFile," %14.8E",tot_density[i]);
        }
        counteraos++;
        if(counteraos==5)
        {
         fprintf(pFile,"\n");
         counteraos=0;
        }
       }
       fclose(pFile);
       formatted_data.open((name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str());
       while(getline(formatted_data,aux))
       {
        dm1fchk<<aux<<endl;
       }
       formatted_data.close();
       system(("rm "+name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str());
      }
      delete[] tot_density;
      tot_density=NULL;
     }
     else
     {dm1fchk<<aux<<endl;}
    }
    else
    {dm1fchk<<aux<<endl;}
   }
   dm1fchk.close();
   read_hf_fchk.close();
   system(("sed s/'Beta'/'B3TA'/g "+name_file.substr(0,(name_file.length()-5))+"_dm1.fchk >"+name_file.substr(0,(name_file.length()-5))+"newdens.rho").c_str());
   system(("mv "+name_file.substr(0,(name_file.length()-5))+"newdens.rho "+name_file.substr(0,(name_file.length()-5))+"_dm1.fchk").c_str());
   if(Input_commands.debug){cout<<"Be aware that Beta was changed to B3TA in the fchk file!"<<endl;}
   if(!Input_commands.cas)
   {
    formatted_data.open((name_file.substr(0,(name_file.length()-5))+"_dm1.fchk").c_str());
    dm1fchk.open((name_file.substr(0,(name_file.length()-5))+"_dm1.fchk2").c_str());
    i=0;
    aux="0";
    while(getline(formatted_data,aux))
    {
     if(i==1)
     {dm1fchk<<"SP        CASSCF                                                      gen"<<endl;}
     else
     {
      if(aux.length()>27)
      {
       if(aux.substr(0,27)=="Number of contracted shells")
       {
        dm1fchk<<"Number of CAS Electrons                    I"<<setw(17)<<Read_fchk_wfn.nelectrons<<endl;
        dm1fchk<<"Number of CAS Orbitals                     I"<<setw(17)<<Read_fchk_wfn.nbasisf<<endl;
        dm1fchk<<aux<<endl;
       }
       else
       {dm1fchk<<aux<<endl;}
      }
      else
      {dm1fchk<<aux<<endl;}
     }
     i++;
    }
    formatted_data.close();
    system(("mv "+name_file.substr(0,(name_file.length()-5))+"_dm1.fchk2 "+name_file.substr(0,(name_file.length()-5))+"_dm1.fchk").c_str());
   }
  }
  ////////////////////////////////////////
  // Print a pseudo wfx file for ESI-3c //
  ////////////////////////////////////////
  if(Input_commands.wfx_print)
  {
   ofstream wfx_file_out;
   ifstream wfx_file_rm;
   wfx_file_rm.open(("pseudo_"+name_file.substr(0,(name_file.length()-5))+".wfx").c_str());
   //Delete it if it was created using only the fchk and log files
   if(wfx_file_rm.good() && Input_commands.wfx_print_dmn)
   {system(("rm pseudo_"+name_file.substr(0,(name_file.length()-5))+".wfx").c_str());}
   wfx_file_rm.close();
   if(!wfn_fchk) //fchk
   {
    if(Input_commands.wfx_print_dmn) //Open only for dmn
    {
     wfx_file_out.open(("pseudo_"+name_file.substr(0,(name_file.length()-5))+".wfx").c_str());
    }
   }
   else
   {
    wfx_file_out.open(("pseudo_"+name_file.substr(0,(name_file.length()-4))+".wfx").c_str());
   }
   if(Input_commands.wfx_print_dmn)
   {
    wfx_file_out<<"<Number of Occupied Molecular Orbitals>"<<endl;
    //Always alpha and beta separated for dmn
    DMN_OPS dmn(Input_commands.name_dm1,1);
    dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
    dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
    if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
    dmn.diagonalize_ab();
    wfx_file_out<<2*Read_fchk_wfn.nbasisf<<endl;
    wfx_file_out<<"</Number of Occupied Molecular Orbitals>"<<endl;
    wfx_file_out<<"<Molecular Orbital Occupation Numbers>"<<endl;
    wfx_file_out<<setprecision(10)<<fixed<<scientific;
    for(i=0;i<Read_fchk_wfn.nbasisf;i++)
    {
     if(dmn.diag_ab_done)
     {
      wfx_file_out<<setw(17)<<dmn.rho_matrixa[i][i]<<endl;
      wfx_file_out<<setw(17)<<dmn.rho_matrixb[i][i]<<endl;
     }
     else
     {
      wfx_file_out<<ZERO<<endl;
     }
    }
    wfx_file_out<<"</Molecular Orbital Occupation Numbers>                                       "<<endl;
    wfx_file_out.close();
   }
   if(!wfn_fchk && !Input_commands.wfx_print_dmn) //fchk file and no dmn
   {
    if(!Input_commands.log)
    {
     cout<<"Warning! Include the $LOG keyword and the .log file name (see the manual)"<<endl;
    }
   }
   if(wfn_fchk) //wfn file
   {
    wfx_file_out<<"<Number of Occupied Molecular Orbitals>"<<endl;
    //As they are in the wfn file (unless virtual NOs where built or NOs where permuted)
    wfx_file_out<<Read_fchk_wfn.nbasis()<<endl;
    wfx_file_out<<"</Number of Occupied Molecular Orbitals>"<<endl;
    wfx_file_out<<"<Molecular Orbital Occupation Numbers>"<<endl;
    wfx_file_out<<setprecision(10)<<fixed<<scientific;
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     wfx_file_out<<setw(17)<<Read_fchk_wfn.Ocupation[i]<<endl;
    }
    wfx_file_out<<"</Molecular Orbital Occupation Numbers>"<<endl;
    wfx_file_out.close();
   }
  }
  ///////////////////////
  // Print a cube file //
  ///////////////////////
  if(Input_commands.cube)
  {
   Results<<endl;
   if(Read_fchk_wfn.wfn)
   {
    Results<<"See file "<<name_file.substr(0,(name_file.length()-4))+".cube"<<endl;
   }
   else
   {
    Results<<"See file "<<name_file.substr(0,(name_file.length()-5))+".cube"<<endl;
   }
   if(Input_commands.opcube=="mo" || Input_commands.opcube=="no"){Read_fchk_wfn.Pair[0]=Input_commands.mo_no_cube-1;}
   if(Input_commands.opcube=="no" && !Input_commands.log && !Read_fchk_wfn.wfn)
   {
    cout<<"Unable to compute NOs using an FCHK file without the LOG file"<<endl;
   }
   else
   {
    cube_file(Read_fchk_wfn,name_file,Input_commands.opcube,Input_commands.cubex,Input_commands.cubey,Input_commands.cubez,
    Input_commands.step_cube_x,Input_commands.step_cube_y,Input_commands.step_cube_z);
   }
   Results<<endl;
   Results<<"#*************************************************************************#";
   Results<<endl;
  }
  /////////////////////
  // End clock       //
  /////////////////////
  system("date +%j' '%H' '%M' '%S>date_RHO.date");
  date_file.open("date_RHO.date");
  date_file>>DATE[1][0]>>DATE[1][1]>>DATE[1][2]>>DATE[1][3];
  date_file.close();
  system("rm date_RHO.date");
  calc_time(DATE);
  Results<<"#*************************************************************************#";
  Results<<endl;
  Results<<"#          .--.                  Try not!                                 #"<<endl;
  Results<<"#::\\`--._,'.::.`._.--'/::     Do or do not.                               #"<<endl;
  Results<<"#::::.  ` __::__ '  .::::    There is no try.                             #"<<endl;
  Results<<"#::::::-:.`'..`'.:-::::::                                                 #"<<endl;
  Results<<"#::::::::\\ `--' /::::::::              -Yoda                              #"<<endl;
  Results<<"#*************************************************************************#";
  Results<<endl;
  Results<<endl;
  Results<<"Git sha: "<<sha<<endl;
  Results<<endl;
  Results<<setprecision(0)<<fixed;
  Results<<"# It took "<<DATE[0][0]<<" days "<<DATE[0][1]<<" hours "<<DATE[0][2]<<" min. ";
  Results<<DATE[0][3]<<" secs. "<<endl;
  Results<<"#*************************************************************************#";
  Results<<endl;
  Results<<"#                        End of job!                                      #";
  Results<<endl;
  Results<<"#*************************************************************************#";
  Results<<endl;
  Results.close();
  if(!wfn_fchk)
  {
   cout<<"See "<<name_file.substr(0,(name_file.length()-5))<<"_RHO.out for results.";
   cout<<endl;
  }
  else
  {
   cout<<"See "<<name_file.substr(0,(name_file.length()-4))<<"_RHO.out for results.";
   cout<<endl;
  }
 }
 else{cout<<"File "<<name_file<<" not found!"<<endl;}
////////////////////////////////////////////////////////////////////////
 return 0;
}
//////////////////////////
//Functions description //
//////////////////////////
void punctual(READ_FCHK_WFN &Read_fchk_wfn,double Point[3], double &Density,double &Density_alpha,
double &Density_beta, double Grad[3], double Grad_alpha[3], double Grad_beta[3],
double &Grad_norm,double &Grad_norm_alpha,double &Grad_norm_beta,double &laplacian_r,
double &laplacian_alpha, double &laplacian_beta,double &tauW_alpha,double &tauW_beta,
double &tau_alpha,double &tau_beta,double &k_F_alpha,double &k_F_beta,double &k_s_alpha,
double &k_s_beta,double &s_r_alpha,double &s_r_beta,double &q_red_alpha,double &q_red_beta,
bool wfn_fchk)
{
 int nbasis,i,j;
 double **NO_orb_grad,AUX[3]={ZERO},shift_min=pow(TEN,-TWO*TEN);
 nbasis=Read_fchk_wfn.nbasis();
 NO_orb_grad=new double*[4];
 for(i=0;i<4;i++)
 {NO_orb_grad[i]=new double[nbasis];}
 for(i=0;i<4;i++)
 {
  for(j=0;j<nbasis;j++)
  {NO_orb_grad[i][j]=ZERO;}
 }
///////////////////////////////////////////
 if(Read_fchk_wfn.overlap || Read_fchk_wfn.wfn){Read_fchk_wfn.orb_grad(Point,NO_orb_grad);}//Get NOs and NO gradients
///////////////////////////////////////////
 Read_fchk_wfn.rho_eval(Point,Density);
 Read_fchk_wfn.rho_eval_a_b(Point,Density_alpha,Density_beta);
 Read_fchk_wfn.rho_grad(Point,Grad);
 Read_fchk_wfn.rho_lapl(Point,laplacian_r);
 Grad_norm=norm3D(Grad);//Norm of the gradient
 Read_fchk_wfn.rho_grad_a_b(Point,Grad_alpha,Grad_beta);
 Grad_norm_alpha=norm3D(Grad_alpha);
 Grad_norm_beta=norm3D(Grad_beta);
 Read_fchk_wfn.rho_lapl_a_b(Point,laplacian_alpha,laplacian_beta);
 tauW_alpha=pow(Grad_norm_alpha,TWO)/(EIGHT*Density_alpha+shift_min);
 tauW_beta=pow(Grad_norm_beta,TWO)/(EIGHT*Density_beta+shift_min);
 k_F_alpha=pow(THREE*PI*PI*Density_alpha+shift_min,(ONE_THIRD));
 k_F_beta=pow(THREE*PI*PI*Density_beta+shift_min,(ONE_THIRD));
 k_s_alpha=pow(FOUR*k_F_alpha/(PI*a_o),(HALF));
 k_s_beta=pow(FOUR*k_F_beta/(PI*a_o),(HALF));
 s_r_alpha=Grad_norm_alpha/(TWO*k_F_alpha*Density_alpha+shift_min);
 s_r_beta=Grad_norm_beta/(TWO*k_F_beta*Density_beta+shift_min);
 q_red_alpha=laplacian_alpha/(FOUR*k_F_alpha*k_F_alpha*Density_alpha+shift_min);
 q_red_beta=laplacian_beta/(FOUR*k_F_beta*k_F_beta*Density_beta+shift_min);
 if((Read_fchk_wfn.correlated && Read_fchk_wfn.open_shell) && (Read_fchk_wfn.wfn && firstcall))
 {cout<<"Warning open-shell correlated wfn file!"<<endl;}
 tau_alpha=ZERO;
 tau_beta=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.correlated)
  {
   if(Read_fchk_wfn.open_shell)
   {
    for(i=0;i<nbasis;i++)
    {
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     if(i%2==0)
     {
      tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
     }
     else
     {
      tau_beta=tau_beta+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
     }
    }
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO)/TWO;
    }
    tau_beta=tau_alpha;
   }
  }
  else
  {
   if(Read_fchk_wfn.open_shell)
   {
    if(firstcall)
    {
     cout<<"Warning! Unable to compute tau_alpha and tau_beta with open-shell"<<endl;
     cout<<"correlated wfn file!"<<endl;
     firstcall=false;
    }
    tau_alpha=ZERO;
    tau_beta=ZERO;
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO)/TWO;
    }
    tau_beta=tau_alpha;
   }
  }
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<nbasis;i++)
    {
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     if(i%2==0)
     {
      tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
     }
     else
     {
      tau_beta=tau_beta+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
     }
    }
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
    }
    tau_alpha=tau_alpha/TWO; //Correct the ocupation!
    tau_beta=tau_alpha;
   }
  }
  else if(firstcall)
  {
   cout<<"Warning! Unable to compute tau_alpha and tau_beta without the"<<endl;
   cout<<"S^1/2 P S^1/2 matrix (tau_alpha and tau_beta will be 0.0e0)."<<endl;
   firstcall=false;
   tau_alpha=ZERO;
   tau_beta=ZERO;
  }
  else
  {
   firstcall=false;
   tau_alpha=ZERO;
   tau_beta=ZERO;
  }
 }
 tau_alpha=HALF*tau_alpha;
 tau_beta=HALF*tau_beta;
//delete dynamic arrays
 for(i=0;i<4;i++)
 {delete[] NO_orb_grad[i];NO_orb_grad[i]=NULL;}
 delete[] NO_orb_grad;NO_orb_grad=NULL;
}
//Compute intermediate quantities used by ELF
void pre_elf(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &Density_alpha,double &Density_beta,
double &tauW_alpha,double &tauW_beta,double &tau_alpha,double &tau_beta,double &tcurr_alpha,double &tcurr_beta)
{
 if(Read_fchk_wfn.wfn)
 {
  int nbasis,i,j,k;
  double Grad_norm_alpha,shift_min=pow(TEN,-TWO*TEN);
  nbasis=Read_fchk_wfn.nbasis();
  Density_alpha=ZERO;Density_beta=ZERO;tauW_alpha=ZERO;tauW_beta=ZERO,tau_alpha=ZERO;tau_beta=ZERO;tcurr_alpha=ZERO;tcurr_beta=ZERO;
  if(Read_fchk_wfn.im_wfn_wfx)
  {
   complex<double>rhoa(ZERO,ZERO);
   complex<double>ztmp0(ZERO,ZERO);
   complex<double>ztmpI(ZERO,ONE);
   complex<double> **NO_orb_grad,AUX[3],gradA[3],grad_currA[3];
   NO_orb_grad=new complex<double>*[4];
   for(i=0;i<3;i++)
   {
    gradA[i]=ztmp0;grad_currA[i]=ztmp0;
   }
   for(i=0;i<4;i++)
   {
    NO_orb_grad[i]=new complex<double>[nbasis];
    for(j=0;j<nbasis;j++)
    {NO_orb_grad[i][j]=ztmp0;}
   }
   Read_fchk_wfn.orb_gradCC(Point,NO_orb_grad);
   if(Read_fchk_wfn.open_shell)
   {
    cout<<"Warning! pre_elf not prepared for open shell"<<endl;
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     rhoa=rhoa+(HALF*Read_fchk_wfn.Ocupation[i])*conj(NO_orb_grad[0][i])*NO_orb_grad[0][i];
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     tau_alpha=tau_alpha+(HALF*Read_fchk_wfn.Ocupation[i])*pow(norm3DCC(AUX),TWO);
     for(j=0;j<3;j++)
     {
      gradA[j]=gradA[j]+(HALF*Read_fchk_wfn.Ocupation[i])
              *(conj(NO_orb_grad[0][i])*NO_orb_grad[j+1][i]+conj(NO_orb_grad[j+1][i])*NO_orb_grad[0][i]);
      grad_currA[j]=grad_currA[j]-ztmpI*HALF*(HALF*Read_fchk_wfn.Ocupation[i])
              *(conj(NO_orb_grad[0][i])*NO_orb_grad[j+1][i]-conj(NO_orb_grad[j+1][i])*NO_orb_grad[0][i]);
     }
    }
    Grad_norm_alpha=norm3DCC(gradA);
    Density_alpha=real(rhoa);Density_beta=Density_alpha;
    tau_alpha=HALF*tau_alpha;tau_beta=tau_alpha;
    tauW_alpha=pow(Grad_norm_alpha,TWO)/(EIGHT*Density_alpha+shift_min);tauW_beta=tauW_alpha;
    tcurr_alpha=HALF*pow(norm3DCC(grad_currA),TWO)/(Density_alpha+shift_min);tcurr_beta=tcurr_alpha;
   }
   //delete dynamic arrays
   for(i=0;i<4;i++)
   {delete[] NO_orb_grad[i];NO_orb_grad[i]=NULL;}
   delete[] NO_orb_grad;NO_orb_grad=NULL;
  }
  else
  {
   double Grad_alpha[3]={ZERO},Grad_beta[3]={ZERO};
   double Grad_norm_beta;
   double **NO_orb_grad,AUX[3]={ZERO};
   NO_orb_grad=new double*[4];
   for(i=0;i<4;i++)
   {
    NO_orb_grad[i]=new double[nbasis];
    for(j=0;j<nbasis;j++)
    {NO_orb_grad[i][j]=ZERO;}
   }
   if(Read_fchk_wfn.overlap || Read_fchk_wfn.wfn){Read_fchk_wfn.orb_grad(Point,NO_orb_grad);}//Get NOs and NO gradients
   if(!Read_fchk_wfn.correlated)
   {
    if(Read_fchk_wfn.open_shell)
    {
     for(i=0;i<nbasis;i++)
     {
      for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
      if(i%2==0)
      {
       for(k=0;k<3;k++)
       {
        Grad_alpha[k]=Grad_alpha[k]+TWO*NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       }
       Density_alpha=Density_alpha+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i];
       tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
      }
      else
      {
       for(k=0;k<3;k++)
       {
        Grad_beta[k]=Grad_beta[k]+TWO*NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       }
       Density_beta=Density_beta+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i];
       tau_beta=tau_beta+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
      }
     }
    }
    else
    {
     for(i=0;i<nbasis;i++)
     {
      for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
      for(k=0;k<3;k++)
      {
       Grad_alpha[k]=Grad_alpha[k]+NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       Grad_beta[k]=Grad_alpha[k];
      }
      Density_alpha=Density_alpha+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
      tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO)/TWO;
     }
     tau_beta=tau_alpha;
     Density_beta=Density_alpha;
    }
   }
   else
   {
    if(Read_fchk_wfn.open_shell)
    {
     Density_alpha=ZERO;
     Density_beta=ZERO;
     tau_alpha=ZERO;
     tau_beta=ZERO;
    }
    else
    {
     for(i=0;i<nbasis;i++)
     {
      for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
      for(k=0;k<3;k++)
      {
       Grad_alpha[k]=Grad_alpha[k]+NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       Grad_beta[k]=Grad_alpha[k];
      }
      Density_alpha=Density_alpha+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
      tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO)/TWO;
     }
     tau_beta=tau_alpha;
     Density_beta=Density_alpha;
    }
   }
   Grad_norm_alpha=norm3D(Grad_alpha);
   Grad_norm_beta=norm3D(Grad_beta);
   tau_alpha=HALF*tau_alpha;
   tau_beta=HALF*tau_beta;
   tauW_alpha=pow(Grad_norm_alpha,TWO)/(EIGHT*Density_alpha+shift_min);
   tauW_beta=pow(Grad_norm_beta,TWO)/(EIGHT*Density_beta+shift_min);
   //delete dynamic arrays
   for(i=0;i<4;i++)
   {delete[] NO_orb_grad[i];NO_orb_grad[i]=NULL;}
   delete[] NO_orb_grad;NO_orb_grad=NULL;
  }
 }
 else
 {
  cout<<"pre_elf function not implemented for FCHK files"<<endl;
 }
}
//Change MOs to NOs for quadrature using the DMN information
void mos_to_nos_int_fchk_dm1(READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands,double **ORBITALS,int &total_grid,int &nbasis)
{
 //Change MOs to NOs in ORBITALS using the DM1 file.
 int i,j,k,counter;
 double *NO_at_r;
 NO_at_r=new double[nbasis];
 DMN_OPS dmn(Input_commands.name_dm1,1);
 dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
 dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
 if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
 dmn.diagonalize_ab();
 for(i=0;i<total_grid;i++)
 {
  counter=0;
  if(Read_fchk_wfn.uhf)
  {
   for(j=0;j<nbasis;j=j+2)
   {
    NO_at_r[j]=ZERO;
    NO_at_r[j+1]=ZERO;
    for(k=0;k<Read_fchk_wfn.nbasisf;k++)
    {
     NO_at_r[j]=NO_at_r[j]+dmn.NO_coefa[k][counter]*ORBITALS[i][2*k];
     NO_at_r[j+1]=NO_at_r[j+1]+dmn.NO_coefb[k][counter]*ORBITALS[i][2*k+1];
    }
    counter++;
   }
  }
  else
  {
   for(j=0;j<nbasis;j++)
   {
    NO_at_r[j]=ZERO;
    for(k=0;k<Read_fchk_wfn.nbasisf;k++)
    {
     NO_at_r[j]=NO_at_r[j]+dmn.NO_coefa[k][j]*ORBITALS[i][k];
    }
   }
  }
  for(j=0;j<nbasis;j++)
  {
   ORBITALS[i][j]=NO_at_r[j];
  }
 }
 delete[] NO_at_r; NO_at_r=NULL;
}
//Change SIJ in MOs to SIJ in NOs using DMN information
void mos_to_nos_dmn_sij(double **SIJ,READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands)
{
 int i,j,nbasis;
 double **aux,**aux1,**NO_coef;
 nbasis=Read_fchk_wfn.nbasis();
 aux=new double*[nbasis];
 aux1=new double*[nbasis];
 NO_coef=new double*[nbasis];
 for(i=0;i<nbasis;i++)
 {
  aux[i]=new double[nbasis];
  aux1[i]=new double[nbasis];
  NO_coef[i]=new double[nbasis];
 }
 for(i=0;i<nbasis;i++)
 {
  for(j=0;j<=i;j++)
  {
   aux[i][j]=ZERO;
   aux[j][i]=ZERO;
   NO_coef[i][j]=ZERO;
   NO_coef[j][i]=ZERO;
   aux1[i][j]=SIJ[i][j];
   aux1[j][i]=aux1[i][j];
  }
 }
 DMN_OPS dmn(Input_commands.name_dm1,1);
 dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
 dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
 if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
 if(Read_fchk_wfn.uhf)
 {
  for(i=0;i<nbasis;i++)
  {
   for(j=0;j<nbasis/2;j++)
   {
    if(i%2==0)
    {
     NO_coef[2*j][i]=dmn.NO_coefa[j][i/2];
    }
    else
    {
     NO_coef[2*j+1][i]=dmn.NO_coefb[j][(i-1)/2];
    }
   }
  }
  dmn.diagonalize_ab();
  matmul(nbasis,aux1,NO_coef,aux);
  mat_transpose(nbasis,NO_coef,aux1);
  mat_equal(nbasis,aux1,NO_coef);
  matmul(nbasis,NO_coef,aux,aux1);
 }
 else
 {
  dmn.diagonalize();
  matmul(nbasis,aux1,dmn.NO_coef,aux);
  mat_transpose(nbasis,dmn.NO_coef,NO_coef);
  matmul(nbasis,NO_coef,aux,aux1);
 }
 for(i=0;i<nbasis;i++)
 {
  for(j=0;j<=i;j++)
  {
   SIJ[i][j]=aux1[i][j];
  }
 }
 for(i=0;i<nbasis;i++)
 {
  delete[] NO_coef[i];NO_coef[i]=NULL;
  delete[] aux[i];aux[i]=NULL;
  delete[] aux1[i];aux1[i]=NULL;
 }
 delete[] NO_coef;NO_coef=NULL;
 delete[] aux;aux=NULL;
 delete[] aux1;aux1=NULL;
}
//Print INT file
void print_int(READ_FCHK_WFN &Read_fchk_wfn,string name_file,double **Sij,int nbasis,double &rho,double &rhoa,
double &rhob,string region)
{
 int i,j,k;
 bool more=false;
 ofstream file_int;
 if(Read_fchk_wfn.wfn)
 {
  file_int.open((name_file.substr(0,(name_file.length()-4))+"_X"+region+".int").c_str());
 }
 else
 {
  file_int.open((name_file.substr(0,(name_file.length()-5))+"_X"+region+".int").c_str());
 }
 file_int<<setprecision(14)<<fixed<<scientific;
 file_int<<" MOLECULAR SCF ENERGY (AU)  =        0.00000000000"<<endl;
 file_int<<endl;
 file_int<<" INTEGRATION IS OVER ATOM   X    "<<region<<endl;
 file_int<<" RESULTS OF THE INTEGRATION"<<endl;
 file_int<<"              N   "<<rho<<"  NET CHARGE"<<endl;
 file_int<<"              G   "<<endl;
 file_int<<"              K   0.00000000000000E+00        E(ATOM)   0.00000000000000E+00"<<endl;
 file_int<<"              L   0.00000000000000E+00"<<endl;
 file_int<<endl;
 file_int<<"          The Atomic Overlap Matrix"<<endl;
 file_int<<endl;
 if(!Read_fchk_wfn.wfn)
 {
  if(Read_fchk_wfn.uhf)
  {
   file_int<<"Unrestricted ";
   if(Read_fchk_wfn.open_shell)
   {file_int<<"Open-Shell Wavefunction"<<endl;}
   else
   {file_int<<"Closed-Shell Wavefunction"<<endl;}
  }
  else
  {
   file_int<<"Restricted Closed-Shell Wavefunction"<<endl;
  }
 }
 else
 {
  if(!Read_fchk_wfn.correlated)
  {
   if(Read_fchk_wfn.uhf)
   {
    file_int<<"Unrestricted ";
    if(Read_fchk_wfn.open_shell && Read_fchk_wfn.multiplicity!=1)
    {
     file_int<<"Open-Shell Wavefunction";
    }
    else
    {
     file_int<<"Closed-Shell Wavefunction";
    }
    file_int<<endl;
   }
   else
   {
    file_int<<"Restricted Closed-Shell Wavefunction"<<endl;
   }
  file_int<<endl;
  }
  else
  {
   file_int<<"Restricted Closed-Shell Wavefunction"<<endl; //Asummed! .wfn information
   file_int<<endl;                                         //is missing.
  }
 }
 file_int<<endl;
 k=0;
 do
 {
  for(i=0+k;i<nbasis;i++)
  {
   for(j=0+k;j<=i;j++)
   {
    if(j<k+5)
    {
     if(Sij[i][j]>=ZERO){file_int<<"  ";}
     else{file_int<<" ";}
     file_int<<Sij[i][j];
     more=false;
    }
    else
    {
     j=i+1;
     more=true;
    }
   }
  file_int<<endl;
  }
 file_int<<endl;
 k=k+5;
 }while(more);
 file_int<<" ALPHA ELECTRONS (NA)   "<<rhoa<<endl;
 file_int<<" BETA ELECTRONS (NB)   "<<rhob<<endl;
 file_int<<endl;
 file_int<<" NORMAL TERMINATION OF PROAIMV"<<endl;
 file_int.close();
}
//Transform the INT file to be readable by ESI-3c
void transform_int(string region,int &nbasis)
{
 int i,j,k,counter,total;
 bool alpha=false,copying=false,more=true;
 double *sij,**SIJ;
 string name,str;
 if(wfn_fchk)
 {name=name_file.substr(0,(name_file.length()-4))+"_X"+region+".int";}
 else
 {name=name_file.substr(0,(name_file.length()-5))+"_X"+region+".int";}
 system(("cp "+name+" tmp87.txt").c_str());
 ifstream input;
 input.open("tmp87.txt");
 ofstream output;
 //Create tmp88.txt with overlap matrix information
 output.open("tmp88.txt");
 total=0;
 for(i=0;i<=nbasis;i++)
 {total=i+total;}
 sij=new double[total];
 SIJ=new double*[nbasis];
 for(i=0;i<nbasis;i++)
 {SIJ[i]=new double[nbasis];}
 while(getline(input,str) && !alpha)
 {
  if(" ALPHA"==str.substr(0,6))
  {alpha=true;}
  if((str.substr(0,12)=="Unrestricted"||str.substr(0,10)=="Restricted") && !alpha)
  {
   getline(input,str);
   copying=true;
  }
  if(copying && !alpha)
  {
   output<<str<<endl;
  }
 }
 input.close();
 output.close();
//Read from tmp88.txt the overlaps
 input.open("tmp88.txt");
 for(i=0;i<total;i++)
 {
  input>>sij[i];
 }
 input.close();
//Here we build SIJ matrix
 counter=0;
 k=0;
 do
 {
  for(i=0+k;i<nbasis;i++)
  {
   for(j=0+k;j<=i;j++)
   {
    if(j<k+5)
    {
     SIJ[i][j]=sij[counter];
     counter++;
     more=false;
    }
    else
    {
     j=i+1;
     more=true;
    }
   }
  }
  k=k+5;
 }while(more);
 output.open(("Transf_"+name).c_str());
 input.open((name).c_str());
 do
 {
  getline(input,str);
  output<<str<<endl;
 }while((str.substr(0,12)!="Unrestricted" && str.substr(0,10)!="Restricted"));
 input.close();
 output<<endl;
 output<<setprecision(10)<<fixed;
 for(i=0;i<nbasis;i++)
 {
  counter=0;
  for(j=0;j<=i;j++)
  {
   if(SIJ[i][j]<ZERO){output<<" ";}
   else{output<<"  ";}
   output<<SIJ[i][j];
   counter++;
   if(counter==8){output<<endl;counter=0;}
  }
  if(counter!=0)
  {output<<endl;}
 }
 output<<endl;
 input.open((name).c_str());
 while(getline(input,str))
 {
  if(" ALPHA"==str.substr(0,6))
  {
   output<<endl;
   output<<str<<endl;
   alpha=false;
  }
  if(!alpha && " ALPHA"!=str.substr(0,6))
  {output<<str<<endl;}
 }
 input.close();
 output.close();
 //Delete arrays and temporal files
 for(i=0;i<nbasis;i++)
 {delete[] SIJ[i];}
 delete[] sij;
 delete[] SIJ;
 system("rm tmp87.txt");
 system("rm tmp88.txt");
}
//Create a cube file
void cube_file(READ_FCHK_WFN &Read_fchk_wfn,string name_file,string op_cube,double cubex,double cubey, double cubez,double stepx,
double stepy,double stepz)
{
 int i,j,k,termsx,termsy,termsz;
 double eval,point[3],grad[3];
 double Density=ZERO,Density_alpha,Density_beta;
 double tauW_alpha,tauW_beta,tau_alpha,tau_beta,tcurr_alpha,tcurr_beta;
 double IND_alpha,IND_beta,ID_alpha,ID_beta;
 double coef_elf=(THREE/FIVE)*pow(SIX*PI*PI,TWO/THREE);
 bool not_elf,not_indic;
 ofstream cube;
 ifstream read_cube;
 string name_cube,aux;
 MO mo;
 NO no;
 wfn_fchk=Read_fchk_wfn.wfn;
 not_elf=false;
 not_indic=false;
 termsx=(int)cubex;termsy=(int)cubey;termsz=(int)cubez;
 FILE *pFile;
 pFile=fopen((name_file+".cub_tmp").c_str(),"w");
 fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f \n",termsx,stepx,ZERO,ZERO);
 fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f \n",termsy,ZERO,stepy,ZERO);
 fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f \n",termsz,ZERO,ZERO,stepz);
 if(op_cube!="densityp" && op_cube!="shannonp" && op_cube!="fisherp")
 {
  for(i=0;i<Read_fchk_wfn.natoms;i++)
  {
   point[0]=Read_fchk_wfn.Cartesian_Coor[i][0]+cubex*stepx/TWO;
   point[1]=Read_fchk_wfn.Cartesian_Coor[i][1]+cubey*stepy/TWO;
   point[2]=Read_fchk_wfn.Cartesian_Coor[i][2]+cubez*stepz/TWO;
   fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f    %8.6f \n",(int)Read_fchk_wfn.Nu_charge[i],ZERO,
   point[0],point[1],point[2]);
  }
 }
 else
 {
  point[0]=cubex*stepx/TWO;
  point[1]=cubey*stepy/TWO;
  point[2]=cubez*stepz/TWO;
  fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f    %8.6f \n",1,ZERO,point[0],point[1],point[2]);
 }
 point[0]=-cubex*stepx/TWO;
 for(i=0;i<termsx;i++)
 {
  point[1]=-cubey*stepy/TWO;
  for(j=0;j<termsy;j++)
  {
   point[2]=-cubez*stepz/TWO;
   for(k=0;k<termsz;k++)
   {
    if(op_cube=="density")
    {Read_fchk_wfn.rho_eval(point,eval);}
    if(op_cube=="laplacian")
    {Read_fchk_wfn.rho_lapl(point,eval);}
    if(op_cube=="neglaplacian")
    {
     Read_fchk_wfn.rho_lapl(point,eval);
     eval=-eval;
    }
    if(op_cube=="shannon")
    {
     Read_fchk_wfn.rho_eval(point,eval);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=-eval/(double)Read_fchk_wfn.nelectrons*log(eval/(double)Read_fchk_wfn.nelectrons);
    }
    if(op_cube=="fisher")
    {
     Read_fchk_wfn.rho_eval(point,eval);
     Read_fchk_wfn.rho_grad(point,grad);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=pow(norm3D(grad),TWO)/(eval*(double)Read_fchk_wfn.nelectrons);
    }
    if(op_cube=="elf" || op_cube=="elfa" || op_cube=="elfb")
    {
     if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
     {
      not_elf=true;
      cout<<"Unable to compute the ELF from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
      cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
      cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
     }
     if(!not_elf)
     {
      pre_elf(Read_fchk_wfn,point,Density_alpha,Density_beta,tauW_alpha,tauW_beta,tau_alpha,tau_beta,tcurr_alpha,tcurr_beta);
      Density=Density_alpha+Density_beta;
     }
     if(op_cube=="elf")
     {
      if(pow(Density,FIVE/THREE)>pow(TEN,-TEN))
      {
       eval=TWO*(tau_alpha+tau_beta-tauW_alpha-tauW_beta-tcurr_alpha-tcurr_beta)/(coef_elf*pow(Density,FIVE/THREE));
       eval=ONE/(ONE+eval*eval);
      }
      else
      {
       eval=ZERO;
      }
     }
     else if(op_cube=="elfa")
     {
      if(pow(Density_alpha,FIVE/THREE)>pow(TEN,-TEN))
      {
       eval=TWO*(tau_alpha-tauW_alpha-tcurr_alpha)/(coef_elf*pow(Density_alpha,FIVE/THREE));
       eval=ONE/(ONE+eval*eval);
      }
      else
      {
       eval=ZERO;
      }
     }
     else
     {
      if(pow(Density_beta,FIVE/THREE)>pow(TEN,-TEN))
      {
       eval=TWO*(tau_beta-tauW_beta-tcurr_beta)/(coef_elf*pow(Density_beta,FIVE/THREE));
       eval=ONE/(ONE+eval*eval);
      }
      else
      {
       eval=ZERO;
      }
     }
    }
    if(op_cube=="id" || op_cube=="ind" || op_cube=="ida" ||  op_cube=="idb" ||  op_cube=="inda" ||  op_cube=="indb" || op_cube=="idnd")
    {
     if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
     {
      not_indic=true;
      cout<<"Unable to compute the INDICATORS from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
      cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
      cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
     }
     if(!not_indic)
     {
      ID_IND_local(Read_fchk_wfn,point,ID_alpha,ID_beta,IND_alpha,IND_beta);
     }
     if(op_cube=="id")
     {
      eval=ID_alpha+ID_beta;
     }
     else if(op_cube=="ind")
     {
      eval=IND_alpha+IND_beta;
     }
     else if(op_cube=="inda")
     {
      eval=IND_alpha;
     }
     else if(op_cube=="indb")
     {
      eval=IND_beta;
     }
     else if(op_cube=="ida")
     {
      eval=ID_alpha;
     }
     else if(op_cube=="idb")
     {
      eval=ID_beta;
     }
     else
     {
      eval=ID_alpha+ID_beta+IND_alpha+IND_beta;
     }
    }
    if(op_cube=="mo")
    {
     mo=MO(Read_fchk_wfn,point,Read_fchk_wfn.Pair[0]);
     eval=mo.evaluation;
    }
    if(op_cube=="no")
    {
     no=NO(Read_fchk_wfn,point,Read_fchk_wfn.Pair[0]);
     eval=no.evaluation;
    }
    if(op_cube=="densityp")
    {Read_fchk_wfn.rho_p_eval(point,eval);}
    if(op_cube=="shannonp")
    {
     Read_fchk_wfn.rho_p_eval(point,eval);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=-eval/(double)Read_fchk_wfn.nelectrons*log(eval/(double)Read_fchk_wfn.nelectrons);
    }
    if(op_cube=="fisherp")
    {
     Read_fchk_wfn.rho_p_eval(point,eval);
     Read_fchk_wfn.rho_p_grad(point,grad);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=pow(norm3D(grad),TWO)/(eval*(double)Read_fchk_wfn.nelectrons);
    }
    if((k%6==0))
    {
     if(eval>=ZERO)
     {fprintf(pFile," %11.5E ",eval);}
     else
     {fprintf(pFile,"%11.5E ",eval);}
    }
    else
    {
     if(eval>=ZERO)
     {fprintf(pFile," %11.5E ",eval);}
     else
     {fprintf(pFile,"%11.5E ",eval);}
    }
    if((k+1)%6==0)
    {fprintf(pFile,"\n");}
    point[2]=point[2]+stepz;
   }
   fprintf(pFile,"\n");
   point[1]=point[1]+stepy;
  }
  point[0]=point[0]+stepx;
 }
 fclose(pFile);
 read_cube.open((name_file+".cub_tmp").c_str());
 if(Read_fchk_wfn.wfn)
 {
  name_cube=name_file.substr(0,(name_file.length()-4));
 }
 else
 {
  name_cube=name_file.substr(0,(name_file.length()-5));
 }
 cube.open((name_cube+".cube").c_str());
 cube<<"RHO_OPS CUBE FILE."<<endl;
 cube<<"OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"<<endl;
 if(op_cube!="densityp" && op_cube!="shannonp" && op_cube!="fisherp")
 {
  //Only up to 99 atoms. I think is enough for this code =P
  cube<<"  ";
  if(Read_fchk_wfn.natoms<10){cube<<" ";}
  cube<<Read_fchk_wfn.natoms;
  cube<<"    0.000000    0.000000    0.000000"<<endl;
 }
 else
 {
  cube<<"  ";
  if(Read_fchk_wfn.natoms<10){cube<<" ";}
  cube<<1;
  cube<<"    0.000000    0.000000    0.000000"<<endl;
 }
 while(getline(read_cube,aux))
 {
  cube<<aux<<endl;
 }
 cube.close();
 read_cube.close();
 system(("rm "+name_file+".cub_tmp").c_str());
}
//Calculate deviation from idempotency
double Deviation_idemp(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double I2=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {I2=ZERO;}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i]/TWO,TWO);}
     I2=ONE-I2/((double)Read_fchk_wfn.nelectrons*HALF);
    }
    else
    {
     cout<<"Warning! Unable to compute I2 open-shell correlated wfn file!"<<endl;
     I2=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i]/TWO,TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons*HALF);
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i],TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons);
   }
  }
  if(I2<pow(TEN,-FIVE)){I2=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i],TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons);
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i]/TWO,TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons*HALF);
   }
  }
  else
  {
   cout<<"Warning! Unable to compute I2 without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(I2 will be 0.0e0)."<<endl;
   I2=ZERO;
  }
  if(I2<pow(TEN,-FIVE)){I2=ZERO;}
 }
 return I2;
}
//Calculate ID
double ID_ni(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double ID=ZERO,occ;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {ID=ZERO;}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
      ID=ID-HALF*abs(occ)*(ONE-abs(occ));
     }
     ID=TWO*ID;
    }
    else
    {
     cout<<"Warning! Unable to compute ID open-shell correlated wfn file!"<<endl;
     ID=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
    ID=TWO*ID;
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i])<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     if(abs(Read_fchk_wfn.Ocupation[i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
    ID=ID;
   }
  }
  if(ID<pow(TEN,-FIVE)){ID=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i])<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     if(abs(Read_fchk_wfn.Ocupation[i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
    ID=ID*TWO;
   }
   if(ID<pow(TEN,-FIVE)){ID=ZERO;}
  }
  else
  {
   cout<<"Warning!  Unable to compute ID without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(ID will be 0.0e0)."<<endl;
   ID=ZERO;
  }
 }
 return ID;
}
//Calculate IND
double IND_ni(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double IND=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {IND=ZERO;}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      IND=IND+HALF*(Read_fchk_wfn.Ocupation[i]/TWO)*(ONE-Read_fchk_wfn.Ocupation[i]/TWO);
     }
     IND=IND*TWO;
    }
    else
    {
     cout<<"Warning! Unable to compute IND open-shell correlated wfn file!"<<endl;
     IND=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*(Read_fchk_wfn.Ocupation[i]/TWO)*(ONE-Read_fchk_wfn.Ocupation[i]/TWO);
    }
    IND=IND*TWO;
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*(Read_fchk_wfn.Ocupation[i])*(ONE-Read_fchk_wfn.Ocupation[i]);
    }
    IND=IND;
   }
  }
  if(IND<pow(TEN,-FIVE)){IND=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]);
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*(Read_fchk_wfn.Ocupation[i]/TWO)*(ONE-Read_fchk_wfn.Ocupation[i]/TWO);
    }
    IND=IND*TWO;
   }
   if(IND<pow(TEN,-FIVE)){IND=ZERO;}
  }
  else
  {
   cout<<"Warning! Unable to compute IND without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(IND will be 0.0e0)."<<endl;
   IND=ZERO;
  }
 }
 return IND;
}
//Calculate Shannon of ocupations
double Shannon_ni(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double S=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {S=(double)Read_fchk_wfn.nelectrons*log((double)Read_fchk_wfn.nelectrons);}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SIX))
      {
       S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons*TWO));
      }
      else
      {}
     }
    }
    else
    {
     cout<<"Warning! Unable to compute S(ni) open-shell correlated wfn file!"<<endl;
     S=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SIX))
     {
      S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons*TWO));
     }
     else
     {}
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SIX))
     {
      S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons));
     }
     else
     {}
    }
   }
  }
  if(S<ZERO)
  {S=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SEVEN))
     {
      S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/(double)Read_fchk_wfn.nelectrons);
     }
     else
     {}
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SEVEN))
     {S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons*TWO));}
     else
     {}
    }
   }
  }
  else
  {
   cout<<"Warning! Unable to compute S(ni) without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(S(ni) will be 0.0e0)."<<endl;
   S=ZERO;
  }
  if(S<ZERO)
  {S=ZERO;}
 }
 return S/(double)Read_fchk_wfn.nelectrons;
}
//Calculate deviation from idempotency from dmn
double Deviation_idemp_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double I2=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  I2=I2+pow(DMN.rho_matrixa[i][i],TWO)+pow(DMN.rho_matrixb[i][i],TWO);
 }
 I2=ONE-I2/N;
 return I2;
}
//Calculate ID
double ID_ni_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double occ,ID=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  if(abs(DMN.rho_matrixa[i][i])<pow(TEN,-SIX)){occ=ZERO;}
  else{occ=DMN.rho_matrixa[i][i];}
  if(abs(DMN.rho_matrixa[i][i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
  else{occ=DMN.rho_matrixa[i][i];}
  ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
  ID=ID-HALF*abs(occ)*(ONE-abs(occ));
  if(abs(DMN.rho_matrixb[i][i])<pow(TEN,-SIX)){occ=ZERO;}
  else{occ=DMN.rho_matrixb[i][i];}
  if(abs(DMN.rho_matrixb[i][i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
  else{occ=DMN.rho_matrixb[i][i];}
  ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
  ID=ID-HALF*abs(occ)*(ONE-abs(occ));
 }
 return ID;
}
//Calculate IND
double IND_ni_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double IND=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  IND=IND+HALF*DMN.rho_matrixa[i][i]*(ONE-DMN.rho_matrixa[i][i]);
  IND=IND+HALF*DMN.rho_matrixb[i][i]*(ONE-DMN.rho_matrixb[i][i]);
 }
 return IND;
}
//Calculate Shannon of ocupations
double Shannon_ni_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double S=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  if(DMN.rho_matrixa[i][i]>=pow(TEN,-SEVEN))
  {
   S=S-DMN.rho_matrixa[i][i]*log(DMN.rho_matrixa[i][i]/N);
  }
  else
  {}
  if(DMN.rho_matrixb[i][i]>=pow(TEN,-SEVEN))
  {
   S=S-DMN.rho_matrixb[i][i]*log(DMN.rho_matrixb[i][i]/N);
  }
  else
  {}
 }
 return S/N;
}

//Local version for Nondynamic and Dynamic corr indcators
void ID_IND_local(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &ID_alpha,double &ID_beta,double &IND_alpha,double &IND_beta)
{
 int nbasis,i,j;
 double **NO_orb_grad,occ;
 nbasis=Read_fchk_wfn.nbasis();
 NO_orb_grad=new double*[4];
 for(i=0;i<4;i++)
 {NO_orb_grad[i]=new double[nbasis];}
 for(i=0;i<4;i++)
 {
  for(j=0;j<nbasis;j++)
  {NO_orb_grad[i][j]=ZERO;}
 }
////////////////////////////////////////////
 if(Read_fchk_wfn.overlap || Read_fchk_wfn.wfn){Read_fchk_wfn.orb_grad(Point,NO_orb_grad);}//Get NOs and NO gradients
//////////////////////////////////////////
 if((Read_fchk_wfn.correlated && Read_fchk_wfn.open_shell) && (Read_fchk_wfn.wfn && firstcall))
 {cout<<"Warning open-shell correlated wfn file!"<<endl;}
 ID_alpha=ZERO;
 ID_beta=ZERO;
 IND_alpha=ZERO;
 IND_beta=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {
    ID_alpha=ZERO;
    ID_beta=ZERO;
    IND_alpha=ZERO;
    IND_beta=ZERO;
   }
   else
   {
    if(Read_fchk_wfn.open_shell)
    {
     if(firstcall)
     {
      cout<<"Warning! Unable to compute local indicators with open-shell"<<endl;
      cout<<"correlated wfn file!"<<endl;
      firstcall=false;
     }
     ID_alpha=ZERO;
     ID_beta=ZERO;
     IND_alpha=ZERO;
     IND_beta=ZERO;
    }
    else
    {
     for(i=0;i<nbasis;i++)
     {
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      ID_alpha=ID_alpha+(HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF)-HALF*abs(occ)*(ONE-abs(occ)))
              *NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]/TWO)
               *NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
     }
     ID_beta=ID_alpha;
     IND_beta=IND_alpha;
    }
   }
  }
  else
  {
   if(Read_fchk_wfn.open_shell)
   {
    for(i=0;i<nbasis;i++)
    {
     if(i%2==0)
     {
      ID_alpha=ID_alpha+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
              -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
               *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
     else
     {
      ID_beta=ID_beta+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
             -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_beta=IND_beta+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
              *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
    }
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID_alpha=ID_alpha+(HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF)-HALF*abs(occ)*(ONE-abs(occ)))
             *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]/TWO)
              *NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
    }
    ID_beta=ID_alpha;
    IND_beta=IND_alpha;
   }
  }
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<nbasis;i++)
    {
     if(i%2==0)
     {
      ID_alpha=ID_alpha+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
              -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
               *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
     else
     {
      ID_beta=ID_beta+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
             -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_beta=IND_beta+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
              *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
    }
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID_alpha=ID_alpha+(HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF)-HALF*abs(occ)*(ONE-abs(occ)))
             *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]/TWO)
              *NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
    }
    ID_beta=ID_alpha;
    IND_beta=IND_alpha;
   }
  }
  else if(firstcall)
  {
   cout<<"Warning! Unable to compute local indicators without the S^1/2 P S^1/2 matrix"<<endl;
   firstcall=false;
   ID_alpha=ZERO;
   ID_beta=ZERO;
  }
  else
  {
   firstcall=false;
   ID_alpha=ZERO;
   ID_beta=ZERO;
   IND_alpha=ZERO;
   IND_beta=ZERO;
  }
 }
 //delete dynamic arrays
 for(i=0;i<4;i++)
 {delete[] NO_orb_grad[i];}
 delete[] NO_orb_grad;
}
//Calculate the time that the program took
void calc_time(double DATE[2][4])
{
 long double timeI,timeF,timeTOTAL;
 timeI=DATE[0][3]+DATE[0][2]*SIXTY+DATE[0][1]*pow(SIXTY,TWO)
      +DATE[0][0]*pow(SIXTY,TWO)*THREE*EIGHT;
 timeF=DATE[1][3]+DATE[1][2]*SIXTY+DATE[1][1]*pow(SIXTY,TWO)
      +DATE[1][0]*pow(SIXTY,TWO)*THREE*EIGHT;
 timeTOTAL=timeF-timeI;
 DATE[0][0]=(int)(timeTOTAL/(pow(SIXTY,TWO)*THREE*EIGHT)) ;
 timeTOTAL=timeTOTAL-DATE[0][0]*pow(SIXTY,TWO)*THREE*EIGHT;
 DATE[0][1]=(int)(timeTOTAL/(pow(SIXTY,TWO)));
 timeTOTAL=timeTOTAL-DATE[0][1]*pow(SIXTY,TWO);
 DATE[0][2]=(int)(timeTOTAL/SIXTY);
 timeTOTAL=timeTOTAL-DATE[0][2]*SIXTY;
 DATE[0][3]=timeTOTAL;
}
