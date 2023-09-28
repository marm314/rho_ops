#include"mescal.h"

//Public functions.
MESCAL::MESCAL(){cout<<"Not allowed default constructor in MESCAL"<<endl;}
MESCAL::MESCAL(string name_output,string name_pdb)
{
 int ifrag,iatom,icoord,jcoord,Sum_Val_elect;
 double pos[3],**Im,**Urot;
 Urot=new double*[3];Im=new double*[3];
 for(icoord=0;icoord<3;icoord++)
 {
  Urot[icoord]=new double[3];
  Im[icoord]=new double[3];
 }
 // Init output file
 init_output(name_output);
 ofstream write_out(name_output,std::ios_base::app);
 write_out<<setprecision(8)<<fixed<<scientific;
 // Read PDB file to store Fragment and atomic information
 nfragments=0;
 read_pdb_file(name_pdb);
 write_out<<endl;
 write_out<<" Fragments read from the PDB file (distances in Bohr) "<<endl;
 write_out<<endl;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  Sum_Val_elect=0;
  write_out<<" Fragment "<<ifrag+1<<endl;
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   write_out<<"  Atom ";
   write_out<<setw(5)<<fragments[ifrag].atoms[iatom].Z;
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[0];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[1];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[2]<<endl;
   Sum_Val_elect+=Z2val_electrons(fragments[ifrag].atoms[iatom].Z);
  }
  Frag_T_inertia(ifrag,pos,Im,Urot);
  write_out<<"  Center of Mass "<<setw(20)<<pos[0]<<setw(20)<<pos[1]<<setw(20)<<pos[2]<<endl;
  write_out<<"  Normalized and diagonalized inertia tensor"<<endl;
  if(abs(Im[0][1])<tol8){Im[0][1]=0.0e0;}
  if(abs(Im[0][2])<tol8){Im[0][2]=0.0e0;}
  if(abs(Im[1][2])<tol8){Im[1][2]=0.0e0;}
  Im[1][0]=Im[0][1];
  Im[2][0]=Im[0][2];
  Im[2][1]=Im[1][2];
  for(icoord=0;icoord<3;icoord++)
  {
   for(jcoord=0;jcoord<3;jcoord++){write_out<<setw(20)<<Im[icoord][jcoord];}write_out<<endl;
  }
  write_out<<"  Rotation matrix (columns)"<<endl;
  for(icoord=0;icoord<3;icoord++)
  {
   for(jcoord=0;jcoord<3;jcoord++){write_out<<setw(20)<<Urot[icoord][order[jcoord]];}write_out<<endl;
  }
  write_out<<"  Total number of valence electrons of this fragment "<<setw(8)<<Sum_Val_elect<<endl;
  // Read fragment.dat files to store charges and asign alpha' atomic contributions using the number of val. electrons  
  // Note: We compute alpha' = U alpha U^T for each fragment 
  read_fragment_file((fragments[ifrag].name+".dat").c_str(),Im,Urot,ifrag,Sum_Val_elect);
  write_out<<endl;
 }
 // Delete Urot and Im because they are no longer needed
 for(icoord=0;icoord<3;icoord++){delete[] Urot[icoord];Urot[icoord]=NULL; delete[] Im[icoord];Im[icoord]=NULL;}
 delete[] Urot; Urot=NULL; delete[] Im; Im=NULL;
 write_out.close();
}

// Do self-consistent solution to find induced dipoles (mu) and fields (F_mu)
void MESCAL::mescal_scs(string name)
{
 int iter=0;
 bool tmp_false=false,conver;
 double E_old,E_new;
 // For permanent charges this is done only once because it does not change
 if(perm_q)
 {
  set_FV_q_inter_frag(tmp_false);
 }
 // Enter SC procedure 
 conver=false;
 do
 {
  if(induced_q)
  {
   set_FV_q_inter_frag(induced_q);
  }
  // mu = alpha F
  update_mu_ind();
  // F_mu_ind = D mu_ind
  set_F_mu_ind();
  // Check conver and compute energy
  // Print iter info 
  iter++;
 }while(iter<maxiter && !conver);
}

// Set F_ext and V_ext (due to a/many point charge(s))
void MESCAL::set_FV_ext_punct(double &q_ext,double Point_mescal[3])
{
 int ifrag,iatom,icoord;
 double r,r3,diff_xyz[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   r=0.0e0;
   for(icoord=0;icoord<3;icoord++)
   {
    diff_xyz[icoord]=fragments[ifrag].atoms[iatom].pos[icoord]-Point_mescal[icoord];
    r+=diff_xyz[icoord]*diff_xyz[icoord];
   }
   r=pow(r,0.5e0);
   r3=pow(r,3.0e0);
   for(icoord=0;icoord<3;icoord++)
   {
    fragments[ifrag].atoms[iatom].F_ext[icoord]+=q_ext*diff_xyz[icoord]/r3;
   }
   fragments[ifrag].atoms[iatom].V_ext+=q_ext/r;
  }
 }
}

MESCAL::~MESCAL()
{
 // Nth to be deleted manually
}

// Compute Inertia tensor for Fragments
void MESCAL::Frag_T_inertia(int &ifrag,double Rcm[3],double **Im,double **Urot)
{
 int iatom,icoord,jcoord,pivot;
 double mass,mass_tot,Norm,Norm_saved,pivot_doub;
 for(icoord=0;icoord<3;icoord++){Rcm[icoord]=0.0e0;}
 mass_tot=0.0e0;
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  mass=Z2mass(fragments[ifrag].atoms[iatom].Z);
  for(icoord=0;icoord<3;icoord++){Rcm[icoord]+=mass*fragments[ifrag].atoms[iatom].pos[icoord];}
  mass_tot+=mass;
 }
 for(icoord=0;icoord<3;icoord++){Rcm[icoord]=Rcm[icoord]/(mass_tot+tol8);}
 Norm_saved=-1.0e0;  
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  Norm=0.0e0;
  for(icoord=0;icoord<3;icoord++)
  {
   fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord]=fragments[ifrag].atoms[iatom].pos[icoord]-Rcm[icoord];
   Norm+=fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord]*fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord];
  }
  Norm=pow(Norm,0.5e0);
  if(abs(Norm)>Norm_saved){Norm_saved=Norm;}
 }
 for(icoord=0;icoord<3;icoord++)
 {
  for(jcoord=0;jcoord<3;jcoord++){Im[icoord][jcoord]=0.0e0;Urot[icoord][jcoord]=0.0e0;}
 }
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  for(icoord=0;icoord<3;icoord++){fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord]=fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord]/(Norm_saved+tol8);}
  mass=Z2mass(fragments[ifrag].atoms[iatom].Z);
  Im[0][0]+=mass*(pow(fragments[ifrag].atoms[iatom].pos_wrt_cm[1],2.0e0)+pow(fragments[ifrag].atoms[iatom].pos_wrt_cm[2],2.0e0));
  Im[1][1]+=mass*(pow(fragments[ifrag].atoms[iatom].pos_wrt_cm[0],2.0e0)+pow(fragments[ifrag].atoms[iatom].pos_wrt_cm[2],2.0e0));
  Im[2][2]+=mass*(pow(fragments[ifrag].atoms[iatom].pos_wrt_cm[0],2.0e0)+pow(fragments[ifrag].atoms[iatom].pos_wrt_cm[1],2.0e0));
  Im[0][1]+=-mass*fragments[ifrag].atoms[iatom].pos_wrt_cm[0]*fragments[ifrag].atoms[iatom].pos_wrt_cm[1];
  Im[0][2]+=-mass*fragments[ifrag].atoms[iatom].pos_wrt_cm[0]*fragments[ifrag].atoms[iatom].pos_wrt_cm[2];
  Im[1][2]+=-mass*fragments[ifrag].atoms[iatom].pos_wrt_cm[1]*fragments[ifrag].atoms[iatom].pos_wrt_cm[2];
 }
 if(abs(Im[0][1])<tol8){Im[0][1]=0.0e0;}
 if(abs(Im[0][2])<tol8){Im[0][2]=0.0e0;}
 if(abs(Im[1][2])<tol8){Im[1][2]=0.0e0;}
 Im[1][0]=Im[0][1];
 Im[2][0]=Im[0][2];
 Im[2][1]=Im[1][2];
 if(abs(Im[0][1])<tol8 && abs(Im[0][2])<tol8 && abs(Im[1][2])<tol8)
 {
  for(icoord=0;icoord<3;icoord++){Urot[icoord][icoord]=1.0e0;}   
 }
 else
 {
  jacobi(3,Im,Urot);
 }
 for(icoord=0;icoord<3;icoord++){order[icoord]=icoord;}
 for(icoord=0;icoord<3;icoord++)
 {
  Norm=abs(Im[icoord][icoord]);
  for(jcoord=icoord+1;jcoord<3;jcoord++)
  {
   if(abs(Im[jcoord][jcoord])>Norm)
   {
    Norm=abs(Im[jcoord][jcoord]);
    pivot_doub=Im[jcoord][jcoord];
    Im[jcoord][jcoord]=Im[icoord][icoord];
    Im[icoord][icoord]=pivot_doub;
    pivot=order[icoord];
    order[icoord]=order[jcoord];
    order[jcoord]=pivot;
   } 
  }
 }
}

// Update mu_ind
void MESCAL::update_mu_ind()
{
 int ifrag,iatom,icoord;
 double Field[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    Field[icoord]=fragments[ifrag].atoms[iatom].F_ext[icoord]
                 +fragments[ifrag].atoms[iatom].F_q_ind[icoord]
                 +fragments[ifrag].atoms[iatom].F_q_perm[icoord]
                 +fragments[ifrag].atoms[iatom].F_mu_ind[icoord];
   }
   alphaF2mu(ifrag,iatom,Field);
  }
 }
}

// Set F_q_inter_fragment and V_q_inter_fragment (due to INDUCED/PERMANENT point charge(s) of the other fragments)
void MESCAL::set_FV_q_inter_frag(bool &induced)
{
 int ifrag,jfrag,iatom,jatom,icoord;
 double r,r3,diff_xyz[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   if(induced)
   {
    fragments[ifrag].atoms[iatom].V_q_ind=0.0e0;
    for(icoord=0;icoord<3;icoord++)
    {
     fragments[ifrag].atoms[iatom].F_q_ind[icoord]=0.0e0;
    }
   }
   else
   {
    fragments[ifrag].atoms[iatom].V_q_perm=0.0e0;
    for(icoord=0;icoord<3;icoord++)
    {
     fragments[ifrag].atoms[iatom].F_q_perm[icoord]=0.0e0;
    }
   }
   for(jfrag=0;jfrag<nfragments;jfrag++)
   {
    if(ifrag!=jfrag) // Only inter-fragment contributions
    {
     for(jatom=0;jatom<fragments[jfrag].natoms;jatom++)
     {
      r=0.0e0;
      for(icoord=0;icoord<3;icoord++)
      {
       diff_xyz[icoord]=fragments[ifrag].atoms[iatom].pos[icoord]-fragments[jfrag].atoms[jatom].pos[icoord];
       r+=diff_xyz[icoord]*diff_xyz[icoord];
      }
      r=pow(r,0.5e0);
      r3=pow(r,3.0e0);
      if(induced)
      {
       fragments[ifrag].atoms[iatom].V_q_ind+=fragments[ifrag].atoms[iatom].charge_ind/r;
       for(icoord=0;icoord<3;icoord++)
       {
        fragments[ifrag].atoms[iatom].F_q_ind[icoord]+=fragments[ifrag].atoms[iatom].charge_ind*diff_xyz[icoord]/r3;
       }
      }
      else
      {
       fragments[ifrag].atoms[iatom].V_q_perm+=fragments[ifrag].atoms[iatom].charge/r;
       for(icoord=0;icoord<3;icoord++)
       {
        fragments[ifrag].atoms[iatom].F_q_perm[icoord]+=fragments[ifrag].atoms[iatom].charge*diff_xyz[icoord]/r3;
       }
      }
     }
    }
   }
  }
 }
}

// Set F_mu_ind (due to induced dipoles)
void MESCAL::set_F_mu_ind()
{
 int ifrag,jfrag,iatom,jatom,icoord;
 double r,r2,r5,diff_xyz[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    fragments[ifrag].atoms[iatom].F_mu_ind[icoord]=0.0e0;
   }
   for(jfrag=0;jfrag<nfragments;jfrag++)
   {
    if(ifrag!=jfrag)
    {
     for(jatom=0;jatom<fragments[jfrag].natoms;jatom++)
     {
      r=0.0e0;
      for(icoord=0;icoord<3;icoord++)
      {
       diff_xyz[icoord]=fragments[ifrag].atoms[iatom].pos[icoord]-fragments[jfrag].atoms[jatom].pos[icoord];
       r+=diff_xyz[icoord]*diff_xyz[icoord];
      }
      r=pow(r,0.5e0);
      r2=pow(r,2.0e0);
      r5=pow(r,5.0e0);
      fragments[ifrag].atoms[iatom].F_mu_ind[0]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[0]
                                                +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[0]*diff_xyz[1]
                                                +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[0]*diff_xyz[2])
                                                -fragments[jfrag].atoms[jatom].mu_ind[0]*r2)/r5;
      fragments[ifrag].atoms[iatom].F_mu_ind[1]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[1]
                                                +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]*diff_xyz[1]
                                                +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[1]*diff_xyz[2])
                                                -fragments[jfrag].atoms[jatom].mu_ind[1]*r2)/r5;
      fragments[ifrag].atoms[iatom].F_mu_ind[0]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[2]
                                                +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]*diff_xyz[2]
                                                +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[2]*diff_xyz[2])
                                                -fragments[jfrag].atoms[jatom].mu_ind[2]*r2)/r5;
     }
    }
   }
  } 
 }
}

// Do mu = alpha F
void MESCAL::alphaF2mu(int &ifrag, int &iatom, double Field[3])
{
 int icoord,jcoord;
 for(icoord=0;icoord<3;icoord++)
 {
  fragments[ifrag].atoms[iatom].mu_ind[icoord]=0.0e0;
  for(jcoord=0;jcoord<3;jcoord++)
  {
   fragments[ifrag].atoms[iatom].mu_ind[icoord]+=fragments[ifrag].atoms[iatom].alpha[icoord][jcoord]*Field[jcoord];
  }
 }
}
