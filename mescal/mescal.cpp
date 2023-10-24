#include"mescal.h"

// Public functions.
MESCAL::MESCAL(){cout<<"Not allowed default constructor in MESCAL"<<endl;}
MESCAL::MESCAL(string name_output,string name_pdb,bool &part_val_e_in, bool &induce_q)
{
 int ifrag,iatom,icoord,jcoord,Sum_Val_elect;
 double pos[3],**Im,**Urot,Sum_atomic_pol;
 part_val_e=part_val_e_in; 
 ind_q=induce_q;
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
  Sum_atomic_pol=0.0e0;
  write_out<<" Fragment "<<ifrag+1<<endl;
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   write_out<<"  Atom ";
   write_out<<setw(5)<<fragments[ifrag].atoms[iatom].Z;
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[0];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[1];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[2]<<endl;
   Sum_Val_elect+=Z2val_electrons(fragments[ifrag].atoms[iatom].Z);
   Sum_atomic_pol+=Z2atomic_pol(fragments[ifrag].atoms[iatom].Z);
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
  write_out<<"  Total number of valence electrons in this fragment "<<setw(8)<<Sum_Val_elect<<endl;
  write_out<<"    Sum of atomic polarizabilities for this fragment "<<setw(8)<<Sum_atomic_pol<<endl;
  // Read fragment.dat files to store charges, Pi, and asign alpha' atomic contributions using the number of val. electrons or tabulated atomic pols.
  // Note: We compute alpha' = U alpha U^T for each fragment 
  read_fragment_file((fragments[ifrag].name+".dat").c_str(),Im,Urot,ifrag,Sum_Val_elect,Sum_atomic_pol);
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
 bool tmp_false=false;
 // For permanent charges this is done only once because it does not change
 if(perm_q)
 {
  set_FV_q_inter_frag(tmp_false);
 }
 // Enter SC procedure
 if(!(w_damp>=0.0e0 && w_damp<1.0e0)){w_damp=0.0e0;cout<<" Warning! The weight_mu is recommended to be in [0,1)."<<endl;}
 if(r0<tol4){r0=0.0e0;cout<<" Warning! Screening r0 < 10^-4 found, setting r0 = 0.0e0."<<endl;}
 if(!ind_q){conver_q=true;}
 iter=0;
 print_init_sc(name); 
 do
 {
  // Set mu = alpha F  and check convergence on mu
  //  also set q_ind = - sum_j Pi_ij (if needed)
  mu_diff_max=-1.0e0;
  q_diff_max=-1.0e0;
  update_mu_q_ind();
  // Set F_mu_ind = D mu_ind and V_mu_ind = mu R if mu is not converged or iter=0
  //  also set q_ind contrib. to F_q_ind and V_q_ind (if needed)
  if(iter==0)
  {
   set_FV_mu_ind();
   if(ind_q){set_FV_q_inter_frag(ind_q);}
  }
  else
  {
   if(ind_q)
   {
    if(mu_diff_max>threshold_mu || q_diff_max>threshold_q)
    {
     set_FV_mu_ind();
     set_FV_q_inter_frag(ind_q);
    }
    else
    {
     conver_mu=true;
     conver_q=true;
    }
   }
   else
   {
    if(mu_diff_max>threshold_mu)
    {
     set_FV_mu_ind();
    }
    else
    {
     conver_mu=true;
    }
   }
  }
  // Check conver Energy (computig energy) and print iter info
  calc_E(name); 
  iter++;
 }while(iter<=maxiter && !(conver_E && conver_mu && conver_q));
 print_end_sc(name); 
 close_output(name);
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

// Send the number atoms in the PDB file
int MESCAL::natoms_tot()
{
 int ifrag,natoms=0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  natoms+=fragments[ifrag].natoms;
 }
 return natoms;
}

// Send coordinates on the PDB
void MESCAL::get_coords(double **Coords)
{
 int ifrag,iatom,icoord,jatom=0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    Coords[jatom][icoord]=fragments[ifrag].atoms[iatom].pos[icoord];
   }
   jatom++;
  }
 }
}

// Set F_ext and V_ext (due to a/many point charge(s))
void MESCAL::set_FV_ext_qm(double **F_ext,double *V_ext)
{
 int ifrag,iatom,icoord,jatom=0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    fragments[ifrag].atoms[iatom].F_ext[icoord]=F_ext[jatom][icoord];
   }
   fragments[ifrag].atoms[iatom].V_ext=V_ext[jatom];
   jatom++;
  }
 }
}

// Calc E and print iter info
void MESCAL::calc_E(string name)
{
 int ifrag,iatom,icoord;
 double E_mu=0.0e0,E_q=0.0e0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    E_mu+=(fragments[ifrag].atoms[iatom].F_ext[icoord]
          +fragments[ifrag].atoms[iatom].F_q_perm[icoord])*fragments[ifrag].atoms[iatom].mu_ind[icoord];
   }
   E_q+=(fragments[ifrag].atoms[iatom].V_ext
        +fragments[ifrag].atoms[iatom].V_q_perm)*fragments[ifrag].atoms[iatom].q_ind;
  }
 }
 Energy=0.5e0*(E_q-E_mu);
 if(iter>0)
 {
  E_diff=Energy-Energy_old;
  if(abs(E_diff)<threshold_E){conver_E=true;}
 }
 else
 {
  E_diff=0.0e0;
 }
 Energy_old=Energy;
 print_iter_info(name);
}

MESCAL::~MESCAL()
{
 int iatom,ifrag;
 if(ind_q)
 {
  for(ifrag=0;ifrag<nfragments;ifrag++)
  {
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    delete[] fragments[ifrag].Pi[iatom];fragments[ifrag].Pi[iatom]=NULL;
   }
   delete[] fragments[ifrag].Pi;fragments[ifrag].Pi=NULL;
  }
 } 
}

// Private functions.

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
   if(abs(fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord])<tol4){fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord]=0.0e0;}
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

// Compute Inertia tensor for reading Fragments (to compare)
void MESCAL::Frag_T_inertia_compare(int &ifrag, double **Cartes_coord, int *Zfrag, double Rcm[3],double **Im,double **Urot)
{
 int iatom,icoord,jcoord,pivot;
 double mass,mass_tot,Norm,Norm_saved,pivot_doub;
 for(icoord=0;icoord<3;icoord++){Rcm[icoord]=0.0e0;}
 mass_tot=0.0e0;
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  mass=Z2mass(Zfrag[iatom]);
  for(icoord=0;icoord<3;icoord++){Rcm[icoord]+=mass*Cartes_coord[iatom][icoord];}
  mass_tot+=mass;
 }
 for(icoord=0;icoord<3;icoord++){Rcm[icoord]=Rcm[icoord]/(mass_tot+tol8);if(abs(Rcm[icoord])>tol8){cout<<"Warning! The fragment was not computed at the CM"<<endl;}}
 Norm_saved=-1.0e0;  
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  Norm=0.0e0;
  for(icoord=0;icoord<3;icoord++)
  {
   Cartes_coord[iatom][icoord]=Cartes_coord[iatom][icoord]-Rcm[icoord];
   if(abs(Cartes_coord[iatom][icoord])<tol4){Cartes_coord[iatom][icoord]=0.0e0;}
   Norm+=Cartes_coord[iatom][icoord]*Cartes_coord[iatom][icoord];
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
  for(icoord=0;icoord<3;icoord++){Cartes_coord[iatom][icoord]=Cartes_coord[iatom][icoord]/(Norm_saved+tol8);}
  mass=Z2mass(Zfrag[iatom]);
  Im[0][0]+=mass*(pow(Cartes_coord[iatom][1],2.0e0)+pow(Cartes_coord[iatom][2],2.0e0));
  Im[1][1]+=mass*(pow(Cartes_coord[iatom][0],2.0e0)+pow(Cartes_coord[iatom][2],2.0e0));
  Im[2][2]+=mass*(pow(Cartes_coord[iatom][0],2.0e0)+pow(Cartes_coord[iatom][1],2.0e0));
  Im[0][1]+=-mass*Cartes_coord[iatom][0]*Cartes_coord[iatom][1];
  Im[0][2]+=-mass*Cartes_coord[iatom][0]*Cartes_coord[iatom][2];
  Im[1][2]+=-mass*Cartes_coord[iatom][1]*Cartes_coord[iatom][2];
  for(icoord=0;icoord<3;icoord++){Cartes_coord[iatom][icoord]=Cartes_coord[iatom][icoord]*Norm_saved;}
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

// Update mu_ind and q_ind
void MESCAL::update_mu_q_ind()
{
 int ifrag,iatom,jatom,icoord;
 double Field[3],V_atom,old_q_ind=0.0e0,q_diff=0.0e0,sum_q;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  sum_q=0.0e0;
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   // Field
   for(icoord=0;icoord<3;icoord++)
   {
    Field[icoord]=fragments[ifrag].atoms[iatom].F_ext[icoord]
                 +fragments[ifrag].atoms[iatom].F_q_perm[icoord]
                 +fragments[ifrag].atoms[iatom].F_q_ind[icoord]
                 +fragments[ifrag].atoms[iatom].F_mu_ind[icoord];
   }
   alphaF2mu(ifrag,iatom,Field);
   // Potential 
   if(ind_q)
   {
    if(iter>0){old_q_ind=fragments[ifrag].atoms[iatom].q_ind;}
    fragments[ifrag].atoms[iatom].q_ind=0.0e0;
    for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
    {
     V_atom=fragments[ifrag].atoms[jatom].V_ext
           +fragments[ifrag].atoms[jatom].V_q_perm
           +fragments[ifrag].atoms[jatom].V_q_ind
           +fragments[ifrag].atoms[jatom].V_mu_ind;
     fragments[ifrag].atoms[iatom].q_ind-=fragments[ifrag].Pi[iatom][jatom]*V_atom; 
    }
    if(iter>0)
    {
     fragments[ifrag].atoms[iatom].q_ind=(1.0e0-w_damp)*fragments[ifrag].atoms[iatom].q_ind+w_damp*old_q_ind;
     q_diff=abs(fragments[ifrag].atoms[iatom].q_ind-old_q_ind);
     if(q_diff>q_diff_max){q_diff_max=q_diff;}
    }
    sum_q+=fragments[ifrag].atoms[iatom].q_ind;
   }
  }
  if(abs(sum_q)>tol5)
  {
   cout<<"Warning! The total ind. charge in fragment "<<ifrag+1<<" does not conserve the num. of electrons "<<setprecision(4)<<fixed<<scientific<<sum_q<<endl;
   cout<<" iter "<<iter+1<<endl;
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
       fragments[ifrag].atoms[iatom].V_q_ind+=fragments[jfrag].atoms[jatom].q_ind/r;
       for(icoord=0;icoord<3;icoord++)
       {
        fragments[ifrag].atoms[iatom].F_q_ind[icoord]+=fragments[jfrag].atoms[jatom].q_ind*diff_xyz[icoord]/r3;
       }
      }
      else
      {
       fragments[ifrag].atoms[iatom].V_q_perm+=fragments[jfrag].atoms[jatom].q_perm/r;
       for(icoord=0;icoord<3;icoord++)
       {
        fragments[ifrag].atoms[iatom].F_q_perm[icoord]+=fragments[jfrag].atoms[jatom].q_perm*diff_xyz[icoord]/r3;
       }
      }
     }
    }
   }
  }
 }
}

// Set FV_mu_ind (due to induced dipoles)
void MESCAL::set_FV_mu_ind()
{
 int ifrag,jfrag,iatom,jatom,icoord;
 double r,r2,r3,r5,fr,diff_xyz[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   fragments[ifrag].atoms[iatom].V_mu_ind=0.0e0;
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
      r3=r*r2;
      r5=r3*r2;
      if(r0>=tol4) // Screen only for r0>=0.0001
      {
       fr=1.0e0-exp(-pow(r/r0,3.0e0));
      }
      else
      { 
       fr=1.0e0;
      }
      fragments[ifrag].atoms[iatom].F_mu_ind[0]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[0]
                                                +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[0]*diff_xyz[1]
                                                +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[0]*diff_xyz[2])
                                                -fragments[jfrag].atoms[jatom].mu_ind[0]*r2)*fr/r5;
      fragments[ifrag].atoms[iatom].F_mu_ind[1]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[1]
                                                +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]*diff_xyz[1]
                                                +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[1]*diff_xyz[2])
                                                -fragments[jfrag].atoms[jatom].mu_ind[1]*r2)*fr/r5;
      fragments[ifrag].atoms[iatom].F_mu_ind[2]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[2]
                                                +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]*diff_xyz[2]
                                                +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[2]*diff_xyz[2])
                                                -fragments[jfrag].atoms[jatom].mu_ind[2]*r2)*fr/r5;
      fragments[ifrag].atoms[iatom].V_mu_ind+=(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]
                                            +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]
                                            +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[2])/r3;
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
 double old_mu=0.0e0,mu_diff=0.0e0;
 for(icoord=0;icoord<3;icoord++)
 {
  if(iter>0){old_mu=fragments[ifrag].atoms[iatom].mu_ind[icoord];}
  fragments[ifrag].atoms[iatom].mu_ind[icoord]=0.0e0;
  for(jcoord=0;jcoord<3;jcoord++)
  {
   fragments[ifrag].atoms[iatom].mu_ind[icoord]+=fragments[ifrag].atoms[iatom].alpha[icoord][jcoord]*Field[jcoord];
  }
  if(iter>0)
  {
   fragments[ifrag].atoms[iatom].mu_ind[icoord]=(1.0e0-w_damp)*fragments[ifrag].atoms[iatom].mu_ind[icoord]+w_damp*old_mu;
   mu_diff+=pow(old_mu-fragments[ifrag].atoms[iatom].mu_ind[icoord],2.0e0);
  }
 }
 if(iter>0)
 {
  mu_diff=pow(mu_diff,0.5e0);
  if(mu_diff>mu_diff_max){mu_diff_max=mu_diff;}
 }
}
