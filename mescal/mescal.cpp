#include"Mescal.h"

// Public functions.
Mescal::Mescal()
{
 nfragments=0;
}

Mescal::Mescal(string name_output,string name_pdb,bool &part_val_e_in, bool &induce_q, bool &mute_in)
{
 part_val_e=part_val_e_in; 
 ind_q=induce_q;
 mute=mute_in;
 mescal_ofile=name_output;
 // Init output file
 if(!mute){init_output();}
 // Read PDB file to store Fragment arrays
 nfragments=0;
 read_pdb_file(name_pdb);
}

Mescal::Mescal(const Mescal&Mescal_obj)
{
 int ifrag,iatom,jatom,icoord,jcoord;
 sha=Mescal_obj.sha;
 perm_q=Mescal_obj.perm_q;
 ind_q=Mescal_obj.ind_q;
 part_val_e=Mescal_obj.part_val_e;
 mute=Mescal_obj.mute;
 nfragments=Mescal_obj.nfragments;
 maxiter=Mescal_obj.maxiter;
 if(Mescal_obj.ifrac_deact.size()>0)
 {
  for(ifrag=0;ifrag<(int)Mescal_obj.ifrac_deact.size();ifrag++)
  {
   ifrac_deact.push_back(Mescal_obj.ifrac_deact[ifrag]);
  }
 }
 iter=Mescal_obj.iter;
 r0=Mescal_obj.r0;
 w_damp=Mescal_obj.w_damp;
 mu_diff_max=Mescal_obj.mu_diff_max;
 q_diff_max=Mescal_obj.q_diff_max;
 E_diff=Mescal_obj.E_diff;
 threshold_mu=Mescal_obj.threshold_mu;
 threshold_E=Mescal_obj.threshold_E;
 threshold_q=Mescal_obj.threshold_q;
 Energy=Mescal_obj.Energy;
 deact_rad=Mescal_obj.deact_rad;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  fragments.push_back({Mescal_obj.fragments[ifrag].name,Mescal_obj.fragments[ifrag].natoms}); 
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   fragments[ifrag].atoms.push_back({Mescal_obj.fragments[ifrag].atoms[iatom].Z});
   for(icoord=0;icoord<3;icoord++)
   {
    fragments[ifrag].atoms[iatom].pos[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].pos[icoord];
    fragments[ifrag].atoms[iatom].mu_ind[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].mu_ind[icoord];
    fragments[ifrag].atoms[iatom].F_ext[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].F_ext[icoord];
    fragments[ifrag].atoms[iatom].F_q_ind[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].F_q_ind[icoord];
    fragments[ifrag].atoms[iatom].F_q_perm[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].F_q_perm[icoord];
    fragments[ifrag].atoms[iatom].F_mu_ind[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].F_mu_ind[icoord];
    fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].pos_wrt_cm[icoord];
    fragments[ifrag].atoms[iatom].dipole_ind[icoord]=Mescal_obj.fragments[ifrag].atoms[iatom].dipole_ind[icoord];
    for(jcoord=0;jcoord<3;jcoord++)
    {
     fragments[ifrag].atoms[iatom].alpha[icoord][jcoord]=Mescal_obj.fragments[ifrag].atoms[iatom].alpha[icoord][jcoord];
    }
   }
   fragments[ifrag].atoms[iatom].q_perm=Mescal_obj.fragments[ifrag].atoms[iatom].q_perm;
   fragments[ifrag].atoms[iatom].q_ind=Mescal_obj.fragments[ifrag].atoms[iatom].q_ind;
   fragments[ifrag].atoms[iatom].V_ext=Mescal_obj.fragments[ifrag].atoms[iatom].V_ext;
   fragments[ifrag].atoms[iatom].V_q_ind=Mescal_obj.fragments[ifrag].atoms[iatom].V_q_ind;
   fragments[ifrag].atoms[iatom].V_q_perm=Mescal_obj.fragments[ifrag].atoms[iatom].V_q_perm;
   fragments[ifrag].atoms[iatom].V_mu_ind=Mescal_obj.fragments[ifrag].atoms[iatom].V_mu_ind;
  }
  for(icoord=0;icoord<3;icoord++)
  {
   fragments[ifrag].Rcm[icoord]=Mescal_obj.fragments[ifrag].Rcm[icoord];
  }
  fragments[ifrag].dist_RcmO=Mescal_obj.fragments[ifrag].dist_RcmO;
  fragments[ifrag].active=Mescal_obj.fragments[ifrag].active;
  fragments[ifrag].i_was_active=Mescal_obj.fragments[ifrag].i_was_active;
  if(ind_q)
  {
   for(ifrag=0;ifrag<nfragments;ifrag++)
   {
    fragments[ifrag].Pi=vector<double>(fragments[ifrag].natoms*fragments[ifrag].natoms,0.0e0);
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
     {
      fragments[ifrag].Pi[iatom+jatom*fragments[ifrag].natoms]=Mescal_obj.fragments[ifrag].Pi[iatom+jatom*fragments[ifrag].natoms];
     }
    }
   }
  } 
 }
}

// Prepare fragments info 
void Mescal::mescal_get_frag_info()
{
 int ifrag,iatom,icoord,jcoord,Sum_Val_elect;
 double pos[3],**Im,**Urot,Urot_order[3],Sum_atomic_pol,norm_cm;
 Urot=new double*[3];Im=new double*[3];
 for(icoord=0;icoord<3;icoord++)
 {
  Urot[icoord]=new double[3];
  Im[icoord]=new double[3];
 }
 // Init output file
 ofstream write_out;
 if(!mute)
 {
  write_out.open(mescal_ofile,std::ios_base::app);
  write_out<<setprecision(8)<<fixed<<scientific;
  // Read atomic information for each fragment
  write_out<<endl;
  write_out<<" Fragments read from the PDB file (distances in Bohr) "<<endl;
  write_out<<endl;
 }
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  Sum_Val_elect=0;
  Sum_atomic_pol=0.0e0;
  if(!mute)
  {
   write_out<<" Fragment "<<ifrag+1<<endl;
  }
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   if(!mute)
   {
    write_out<<"  Atom ";
    write_out<<setw(5)<<fragments[ifrag].atoms[iatom].Z;
    write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[0];
    write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[1];
    write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[2]<<endl;
   }
   Sum_Val_elect+=Z2val_electrons(fragments[ifrag].atoms[iatom].Z);
   Sum_atomic_pol+=Z2atomic_pol(fragments[ifrag].atoms[iatom].Z);
  }
  Frag_T_inertia(ifrag,pos,Im,Urot);
  norm_cm=0.0e0;
  for(icoord=0;icoord<3;icoord++){fragments[ifrag].Rcm[icoord]=pos[icoord];norm_cm+=pow(pos[icoord],2.0e0);}
  if(!mute)
  {
   write_out<<"  Center of Mass ";
   for(icoord=0;icoord<3;icoord++){write_out<<setw(20)<<fragments[ifrag].Rcm[icoord];}
   write_out<<endl;
  }
  norm_cm=pow(norm_cm,0.5e0);
  fragments[ifrag].dist_RcmO=norm_cm;
  if(!mute)
  {
   write_out<<"  Distance of the center of mass to the origin "<<setw(20)<<norm_cm<<endl;
   write_out<<"  Normalized and diagonalized inertia tensor"<<endl;
  }
  if(abs(Im[0][1])<tol8){Im[0][1]=0.0e0;}
  if(abs(Im[0][2])<tol8){Im[0][2]=0.0e0;}
  if(abs(Im[1][2])<tol8){Im[1][2]=0.0e0;}
  Im[1][0]=Im[0][1];
  Im[2][0]=Im[0][2];
  Im[2][1]=Im[1][2];
  for(icoord=0;icoord<3;icoord++)
  {
   for(jcoord=0;jcoord<3;jcoord++)
   {
    Urot_order[jcoord]=Urot[icoord][order[jcoord]];
   }
   for(jcoord=0;jcoord<3;jcoord++)
   {
    Urot[icoord][jcoord]=Urot_order[jcoord];
   }
  }
  if(!mute)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    for(jcoord=0;jcoord<3;jcoord++){write_out<<setw(20)<<Im[icoord][jcoord];}write_out<<endl;
   }
   write_out<<"  Rotation matrix (columns)"<<endl;
   for(icoord=0;icoord<3;icoord++)
   {
    for(jcoord=0;jcoord<3;jcoord++){write_out<<setw(20)<<Urot[icoord][jcoord];}write_out<<endl;
   }
   write_out<<"  Total number of valence electrons in this fragment "<<setw(8)<<Sum_Val_elect<<endl;
   write_out<<"    Sum of atomic polarizabilities for this fragment "<<setw(8)<<Sum_atomic_pol<<endl;
  }
  // Read fragment.dat files to store charges, Pi, and asign alpha' atomic contributions using the number of val. electrons or tabulated atomic pols.
  // Note: We compute alpha' = U alpha U^T for each fragment 
  read_fragment_file((fragments[ifrag].name+".dat").c_str(),Im,Urot,ifrag,Sum_Val_elect,Sum_atomic_pol);
  if(!mute)
  {
   write_out<<endl;
  }
 }
 // Delete Urot and Im because they are no longer needed
 for(icoord=0;icoord<3;icoord++){delete[] Urot[icoord];Urot[icoord]=NULL; delete[] Im[icoord];Im[icoord]=NULL;}
 delete[] Urot; Urot=NULL; delete[] Im; Im=NULL;
 if(!mute)
 {
  write_out.close();
 }
}

// Do self-consistent solution to find induced dipoles (mu) and fields (F_mu)
void Mescal::mescal_scs()
{
 bool tmp_false=false,tmp_true=true;
 // For permanent charges this is done only once because it does not change
 if(perm_q)
 {
  set_FV_inter_frag(tmp_false,tmp_true); // induced_q = false  & permanent_q = true
 }
 // Enter SC procedure
 if(!(w_damp>=0.0e0 && w_damp<1.0e0)){w_damp=0.0e0;cout<<" Warning! The weight_mu is recommended to be in [0,1)."<<endl;}
 if(r0<tol4){r0=0.0e0;cout<<" Warning! Screening r0 < 10^-4 found, setting r0 = 0.0e0."<<endl;}
 if(!ind_q){conver_q=true;}
 iter=0;
 if(!mute){print_init_sc();}
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
   if(ind_q)
   {
    set_FV_inter_frag(tmp_true,tmp_false); // induced_q = true  & permanent_q = false
   }
   else
   {
    set_FV_inter_frag(tmp_false,tmp_false); // induced_q = false  & permanent_q = false
   }
  }
  else
  {
   if(ind_q)
   {
    if(mu_diff_max>threshold_mu || q_diff_max>threshold_q)
    {
     set_FV_inter_frag(tmp_true,tmp_false); // induced_q = true  & permanent_q = false
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
     set_FV_inter_frag(tmp_false,tmp_false); // induced_q = false  & permanent_q = false
    }
    else
    {
     conver_mu=true;
    }
   }
  }
  // Check conver Energy (computig energy) and print iter info
  calc_E(); 
  iter++;
 }while(iter<=maxiter && !(conver_E && conver_mu && conver_q));
 if(!mute)
 {
  print_charges_file();
  print_end_sc();
  close_output();
 }
}

// Set F_ext and V_ext (due to a/many point charge(s))
void Mescal::set_FV_ext_punct(double &q_ext,double Point_mescal[3])
{
 int ifrag,iatom,icoord;
 double r,r3,diff_xyz[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
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
}

// Clean data from previous run (also store info about active fragments) 
void Mescal::clean_converg()
{
 conver_E=false,conver_mu=false,conver_q=false;
}

// Clean data from previous run (also store info about active fragments) 
void Mescal::clean()
{
 int ifrag,iatom,icoord;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {	 
   fragments[ifrag].i_was_active=true;
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    for(icoord=0;icoord<3;icoord++)
    {
     fragments[ifrag].atoms[iatom].F_ext[icoord]=0.0e0;
    }
    fragments[ifrag].atoms[iatom].V_ext=0.0e0;
   }
  }
 }
 nactive=-1;
 conver_E=false,conver_mu=false,conver_q=false;
}

// Clean data from previous run (do not store info about active fragments) 
void Mescal::clean_all()
{
 int ifrag,iatom,icoord;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  fragments[ifrag].i_was_active=false;
  if(fragments[ifrag].active)
  {     
   fragments[ifrag].active=false;
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    for(icoord=0;icoord<3;icoord++)
    {
     fragments[ifrag].atoms[iatom].F_ext[icoord]=0.0e0;
     fragments[ifrag].atoms[iatom].F_q_ind[icoord]=0.0e0;
     fragments[ifrag].atoms[iatom].F_mu_ind[icoord]=0.0e0;
    }
    fragments[ifrag].atoms[iatom].V_ext=0.0e0;
    fragments[ifrag].atoms[iatom].V_q_ind=0.0e0;
    fragments[ifrag].atoms[iatom].V_mu_ind=0.0e0;
   }
  }
 }
 nactive=-1;
 conver_E=false,conver_mu=false,conver_q=false;
}

// Activate a particular fragment given the atomic number Z and its coordinates (a.u. are assumed)
void Mescal::activate_fragment(int &natoms_in, int *Z,double **coord)
{
 bool located=false;
 int ifrag,iatom;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(natoms_in==fragments[ifrag].natoms && !located)
  {
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    if( ( Z[iatom]=fragments[ifrag].atoms[iatom].Z && 
        abs(coord[iatom][0]-fragments[ifrag].atoms[iatom].pos[0])<tol4 ) &&
        ( abs(coord[iatom][1]-fragments[ifrag].atoms[iatom].pos[1])<tol4 &&
        abs(coord[iatom][2]-fragments[ifrag].atoms[iatom].pos[2])<tol4 )
      ){located=true;}
    else{located=false;iatom=fragments[ifrag].natoms;}
   }
  }
  if(located)
  {
   fragments[ifrag].active=true;nactive++;ifrac_deact.push_back(ifrag);ifrag=nfragments;
  }
 }
}

// Deactivate a particular fragment given the atomic number Z and its coordinates (a.u. are assumed)
void Mescal::deactivate_fragment(int &natoms_in, int *Z,double **coord)
{
 bool located=false;
 int ifrag,iatom;
 double rad=1e99;
 if(deact_rad)
 {
  for(ifrag=0;ifrag<nfragments;ifrag++)
  {
   if(natoms_in==fragments[ifrag].natoms && !located)
   {
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     if( ( Z[iatom]=fragments[ifrag].atoms[iatom].Z && 
         abs(coord[iatom][0]-fragments[ifrag].atoms[iatom].pos[0])<tol4 ) &&
         ( abs(coord[iatom][1]-fragments[ifrag].atoms[iatom].pos[1])<tol4 &&
         abs(coord[iatom][2]-fragments[ifrag].atoms[iatom].pos[2])<tol4 )
       ){located=true;}
     else{located=false;iatom=fragments[ifrag].natoms;}
    }
   }
   if(located)
   {
    fragments[ifrag].active=false;nactive--;ifrac_deact.push_back(ifrag);ifrag=nfragments;
   }
  }
 }
 else
 {
  deactivate_fragments(rad);
  for(ifrag=0;ifrag<nfragments;ifrag++)
  {
   if(natoms_in==fragments[ifrag].natoms && !located)
   {
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     if( ( Z[iatom]=fragments[ifrag].atoms[iatom].Z && 
         abs(coord[iatom][0]-fragments[ifrag].atoms[iatom].pos[0])<tol4 ) &&
         ( abs(coord[iatom][1]-fragments[ifrag].atoms[iatom].pos[1])<tol4 &&
         abs(coord[iatom][2]-fragments[ifrag].atoms[iatom].pos[2])<tol4 )
       ){located=true;}
     else{located=false;iatom=fragments[ifrag].natoms;}
    }
   }
   if(located)
   {
    fragments[ifrag].active=false;nactive--;ifrac_deact.push_back(ifrag);ifrag=nfragments;
   }
  }
 }
}

// Deactivate fragments at distance higher than rad
void Mescal::deactivate_fragments(double &rad)
{
 int ifrag;
 nactive=0;
 radius=rad;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].dist_RcmO<rad){fragments[ifrag].active=true;nactive++;}
  else{fragments[ifrag].active=false;}
 }
 deact_rad=true;
}

// Send the number atoms in the PDB file
int Mescal::natoms_tot()
{
 int ifrag,natoms=0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {
   natoms+=fragments[ifrag].natoms;
  }
 }
 return natoms;
}

// Send coordinates on the PDB
void Mescal::get_coords(double **Coords)
{
 int ifrag,iatom,icoord,jatom=0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
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
}

// Set F_ext and V_ext for a given atom of a fragment
void Mescal::set_FV_ext_one(int &ifrag,int &iatom,double F_ext[3],double &V_ext)
{
 int icoord;
 if(fragments[ifrag].active) 
 {
  for(icoord=0;icoord<3;icoord++)
  {
   fragments[ifrag].atoms[iatom].F_ext[icoord]=F_ext[icoord];
  }
  fragments[ifrag].atoms[iatom].V_ext=V_ext;
 }
}

// Set F_ext and V_ext (due to a/many point charge(s))
void Mescal::set_FV_ext_qm(double **F_ext,double *V_ext)
{
 int ifrag,iatom,icoord,jatom=0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active) 
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
}

// Calc E and print iter info
void Mescal::calc_E()
{
 int ifrag,iatom,icoord;
 double E_mu=0.0e0,E_q=0.0e0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
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
 if(!mute){print_iter_info();}
}

// Get V at Point_r (due to induced dipoles and charges)
void Mescal::get_V_punct(double &V_r,double Point_r[3])
{
 int ifrag,iatom,icoord;
 double r,r3,diff_xyz[3];
 V_r=0.0e0;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {	  
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    r=0.0e0;
    for(icoord=0;icoord<3;icoord++)
    {
     diff_xyz[icoord]=Point_r[icoord]-fragments[ifrag].atoms[iatom].pos[icoord];
     r+=diff_xyz[icoord]*diff_xyz[icoord];
    }
    r=pow(r,0.5e0);
    r3=pow(r,3.0e0);
    for(icoord=0;icoord<3;icoord++) // Pot. due to mu_ind 
    {
     V_r+=fragments[ifrag].atoms[iatom].mu_ind[icoord]*diff_xyz[icoord]/r3;
    }                               // Pot. due to q_ind
    V_r+=fragments[ifrag].atoms[iatom].q_ind/r;
   }
  }
 }
}

Mescal::~Mescal()
{

}

// Private functions.

// Compute Inertia tensor for Fragments
void Mescal::Frag_T_inertia(int &ifrag,double Rcm[3],double **Im,double **Urot)
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
  jacobi_mescal(3,Im,Urot);
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
void Mescal::Frag_T_inertia_compare(int &ifrag, double **Cartes_coord, int *Zfrag, double Rcm[3],double **Im,double **Urot)
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
  jacobi_mescal(3,Im,Urot);
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
void Mescal::update_mu_q_ind()
{
 int ifrag,iatom,jatom,icoord;
 double Field[3],V_atom,old_q_ind=0.0e0,q_diff=0.0e0,sum_q;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {
   if(!(iter==0 && fragments[ifrag].i_was_active)) // Using as guess the previous mu_ind and q_ind for the fragments that were active
   {	                                           // in the previous run of the SCS.
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
       fragments[ifrag].atoms[iatom].q_ind-=fragments[ifrag].Pi[iatom+jatom*fragments[ifrag].natoms]*V_atom; 
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
     cout<<"Warning! The total ind. charge in fragment "<<ifrag+1<<" does not conserve the num. of electrons ";
     cout<<setprecision(4)<<fixed<<scientific<<sum_q<<endl;
     cout<<" iter "<<iter+1<<endl;
    }
   }
  } 
 }
}

// Set F_inter_fragment and V_inter_fragment (due to INDUCED/PERMANENT point charge(s) and INDUCED dipoles of other fragments)
void Mescal::set_FV_inter_frag(bool &induced_q,bool &permanent_q)
{
 int ifrag,jfrag,iatom,jatom,icoord;
 double V_mu_ind,V_q_ind,V_q_perm;
 double pos[3],F_mu_ind[3],F_q_ind[3],F_q_perm[3];
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {	  
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    V_mu_ind=0.0e0;V_q_ind=0.0e0;V_q_perm=0.0e0;
    for(icoord=0;icoord<3;icoord++)
    {
     F_mu_ind[icoord]=0.0e0;
     F_q_ind[icoord]=0.0e0;
     F_q_perm[icoord]=0.0e0;
     pos[icoord]=fragments[ifrag].atoms[iatom].pos[icoord];
    }

    #pragma omp parallel num_threads(nthread) \
     private(jfrag,jatom,icoord) firstprivate(pos) \
     shared(fragments,ifrag,nfragments,induced_q,permanent_q) \
     reduction(+:V_mu_ind,V_q_ind,V_q_perm,F_mu_ind,F_q_ind,F_q_perm) 
    {
     //Variables that belong to each thread and are only needed here
     int nth=omp_get_num_threads();
     int ith=omp_get_thread_num();
     double r,r2,r3,r5,fr,diff_xyz[3];
     for(jfrag=ith;jfrag<nfragments;jfrag=jfrag+nth)
     {
      if(ifrag!=jfrag && fragments[jfrag].active) // Only inter-fragment contributions
      {
       for(jatom=0;jatom<fragments[jfrag].natoms;jatom++)
       {
        r=0.0e0;
        for(icoord=0;icoord<3;icoord++)
        {
         diff_xyz[icoord]=pos[icoord]-fragments[jfrag].atoms[jatom].pos[icoord];
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
        if(permanent_q)
        {
         V_q_perm+=fragments[jfrag].atoms[jatom].q_perm/r;
         for(icoord=0;icoord<3;icoord++)
         {
          F_q_perm[icoord]+=fragments[jfrag].atoms[jatom].q_perm*diff_xyz[icoord]/r3;
         }
        }
        else
        {
         // Induced charges? 
         if(induced_q)
         {
          V_q_ind+=fragments[jfrag].atoms[jatom].q_ind/r;
          for(icoord=0;icoord<3;icoord++)
          {
           F_q_ind[icoord]+=fragments[jfrag].atoms[jatom].q_ind*diff_xyz[icoord]/r3;
          }
         }
         // Induced dipoles are always updated 
         F_mu_ind[0]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[0]
                     +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[0]*diff_xyz[1]
                     +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[0]*diff_xyz[2])
                    -fragments[jfrag].atoms[jatom].mu_ind[0]*r2)*fr/r5;
         F_mu_ind[1]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[1]
                     +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]*diff_xyz[1]
                     +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[1]*diff_xyz[2])
                     -fragments[jfrag].atoms[jatom].mu_ind[1]*r2)*fr/r5;
         F_mu_ind[2]+=(3.0e0*(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]*diff_xyz[2]
                     +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]*diff_xyz[2]
                     +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[2]*diff_xyz[2])
                     -fragments[jfrag].atoms[jatom].mu_ind[2]*r2)*fr/r5;
         V_mu_ind+=(fragments[jfrag].atoms[jatom].mu_ind[0]*diff_xyz[0]
                  +fragments[jfrag].atoms[jatom].mu_ind[1]*diff_xyz[1]
                  +fragments[jfrag].atoms[jatom].mu_ind[2]*diff_xyz[2])/r3;
     
        }
       }
      }
     }
    }

    if(!permanent_q)
    {	   
     fragments[ifrag].atoms[iatom].V_mu_ind=V_mu_ind;
     for(icoord=0;icoord<3;icoord++)
     {
      fragments[ifrag].atoms[iatom].F_mu_ind[icoord]=F_mu_ind[icoord];
     }
    }
    if(induced_q)
    {
     fragments[ifrag].atoms[iatom].V_q_ind=V_q_ind;
     for(icoord=0;icoord<3;icoord++)
     {
      fragments[ifrag].atoms[iatom].F_q_ind[icoord]=F_q_ind[icoord];
     }
    }
    if(permanent_q)
    {
     fragments[ifrag].atoms[iatom].V_q_perm=V_q_perm;
     for(icoord=0;icoord<3;icoord++)
     {
      fragments[ifrag].atoms[iatom].F_q_perm[icoord]=F_q_perm[icoord];
     }
    }

   }
  }
 }
}

// Do mu = alpha F
void Mescal::alphaF2mu(int &ifrag, int &iatom, double Field[3])
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

// Get q_charge and their coords that reproduce q_ind and mu_ind (in AU) and consistent with Gabriele's code
void Mescal::get_ind_q_frag_atom(int &ifrag,int &iatom, double q_charge[2],double coord_q[2][3])
{
 int icoord;
 double distO,distAU,u_vec,norm_vec;
 distO=5.0e-4*Angs2au;distAU=2.0e0*distO;
 norm_vec=0.0e0;
 for(icoord=0;icoord<3;icoord++)
 {
  norm_vec=norm_vec+pow(fragments[ifrag].atoms[iatom].mu_ind[icoord],2.0e0);
 }
 norm_vec=pow(norm_vec,0.5e0);
 for(icoord=0;icoord<3;icoord++)
 {
  u_vec=fragments[ifrag].atoms[iatom].mu_ind[icoord]/norm_vec;
  coord_q[0][icoord]= u_vec*distO+fragments[ifrag].atoms[iatom].pos[icoord];
  coord_q[1][icoord]=-u_vec*distO+fragments[ifrag].atoms[iatom].pos[icoord];
 }
 q_charge[0]= norm_vec/distAU+0.5e0*fragments[ifrag].atoms[iatom].q_ind;
 q_charge[1]=-norm_vec/distAU+0.5e0*fragments[ifrag].atoms[iatom].q_ind;
}

// Functions to be used by MPI for read .pdb and .dat files stored in master

// Using read info (for parallel MPI coding)
void Mescal::use_pdb_info(int &natoms,string *pdb_file)
{
 int Z=1,ifrag,iatom,count_fragments=-1,old_fragment=-1,new_fragment;
 double pos[3];
 string line,line_aux;
 for(iatom=0;iatom<natoms;iatom++)
 {
  line=pdb_file[iatom];
  stringstream ss(line.substr(20,6));
  ss>>new_fragment;
  if(new_fragment!=old_fragment)
  {
   count_fragments++;
   line_aux=line.substr(17,3);
   line_aux.erase(std::remove_if(line_aux.begin(),line_aux.end(),::isspace),line_aux.end());
   old_fragment=new_fragment;
   fragments.push_back({line_aux,1});                                         // name fragment, natoms (init.)
   // Here add atomic info to push_back
   line_aux=line.substr(11,5);
   line_aux.erase(std::remove_if(line_aux.begin(),line_aux.end(),::isspace),line_aux.end());
   if(line_aux.length()>1)
   {
    if(islower(line_aux[1]))
    {
     line_aux=line_aux.substr(0,2);
    }
    else
    {
     line_aux=line_aux.substr(0,1);
    }
   }
   else
   {
    line_aux=line_aux.substr(0,1);
   }
   Asymbol2Z(Z,line_aux);
   // Det Atomic positions
   line_aux=line.substr(30,8);
   stringstream ss1(line_aux);
   ss1>>pos[0];
   pos[0]=pos[0]*Angs2au;
   line_aux=line.substr(38,8);
   stringstream ss2(line_aux);
   ss2>>pos[1];
   pos[1]=pos[1]*Angs2au;
   line_aux=line.substr(46,8);
   stringstream ss3(line_aux);
   ss3>>pos[2];
   pos[2]=pos[2]*Angs2au;
   fragments[count_fragments].atoms.push_back({Z,pos[0],pos[1],pos[2]});// Z, position
   nfragments++;
  }
  else
  {
   fragments[count_fragments].natoms++;                                        // Add an atom to the fragment
   // Here add atomic info to push_back
   // Det Z
   line_aux=line.substr(11,5);
   line_aux.erase(std::remove_if(line_aux.begin(),line_aux.end(),::isspace),line_aux.end());
   if(line_aux.length()>1)
   {
    if(islower(line_aux[1]))
    {
     line_aux=line_aux.substr(0,2);
    }
    else
    {
     line_aux=line_aux.substr(0,1);
    }
   }
   else
   {
    line_aux=line_aux.substr(0,1);
   }
   Asymbol2Z(Z,line_aux);
   // Det Atomic positions
   line_aux=line.substr(30,8);
   stringstream ss1(line_aux);
   ss1>>pos[0];
   pos[0]=pos[0]*Angs2au;
   line_aux=line.substr(38,8);
   stringstream ss2(line_aux);
   ss2>>pos[1];
   pos[1]=pos[1]*Angs2au;
   line_aux=line.substr(46,8);
   stringstream ss3(line_aux);
   ss3>>pos[2];
   pos[2]=pos[2]*Angs2au;
   fragments[count_fragments].atoms.push_back({Z,pos[0],pos[1],pos[2]});// Z, position  
  }
 }
 if(ind_q)
 {
  for(ifrag=0;ifrag<nfragments;ifrag++)
  {
   fragments[ifrag].Pi=vector<double>(fragments[ifrag].natoms*fragments[ifrag].natoms,0.0e0);
  }
 }
}
