#include"mescal.h"

//Public functions.
MESCAL::MESCAL(){cout<<"Not allowed default constructor in MESCAL"<<endl;}
MESCAL::MESCAL(string name_output,string name_pdb)
{
 int ifrag,iatom,icoord,jcoord;
 double pos[3],**Im,**Urot;
 Urot=new double*[3];Im=new double*[3];
 for(icoord=0;icoord<3;icoord++)
 {
  Urot[icoord]=new double[3];
  Im[icoord]=new double[3];
 }
 ofstream write_out(name_output);
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<"------                 MESCAL ( in C++ )                   -----"<<endl; 
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<"--          Developed by: Dr. M. Rodriguez-Mayorga            --"<<endl; 
 write_out<<"----               email: marm3.14@gmail.com               -----"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<endl;  
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<"--  MicroElectroStatic Calculations (MESCAL)               -----"<<endl;
 write_out<<"--  Based on MESCAL (in Fortran) code by                   -----"<<endl;
 write_out<<"--  Dr. Gabriele D'Avino (2013-2015)                       -----"<<endl;
 write_out<<"--  email: gabriele.davino@gmail.com                       -----"<<endl;
 write_out<<"--  download: https://gitlab.com/taphino/mescal            -----"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<endl;  
 write_out<<setprecision(8)<<fixed<<scientific;
 // Read PDB file to store Fragment and atomic information
 nfragments=0;
 read_pdb_file(name_pdb);

 // Read fragment.dat files to store alpha_ij and rot matrices 

 write_out<<endl;
 write_out<<" Fragments read from the PDB file (distances in Bohr) "<<endl;
 write_out<<endl;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  write_out<<" Fragment "<<ifrag+1<<endl;
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   write_out<<"  Atom ";
   write_out<<setw(5)<<fragments[ifrag].atoms[iatom].Z;
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[0];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[1];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[2]<<endl;
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
  write_out<<endl;
  // Compute alpha' = U alpha U^T for each fragment 
  
  // Set alpha atomic contributions 
 }
 // Delete Urot and Im because they are no longer needed
 for(icoord=0;icoord<3;icoord++){delete[] Urot[icoord];Urot[icoord]=NULL; delete[] Im[icoord];Im[icoord]=NULL;}
 delete[] Urot; Urot=NULL; delete[] Im; Im=NULL;

 // Do self-consistent sol. to get induced dipoles

 write_out<<endl;  
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<" /  \\     /  |/        | /      \\  /      \\           /  |"<<endl;
 write_out<<" $$  \\   /$$ |$$$$$$$$/ /$$$$$$  |/$$$$$$  |  ______  $$ |  "<<endl;
 write_out<<" $$$  \\ /$$$ |$$ |__    $$ \\__$$/ $$ |  $$/  /      \\ $$ | "<<endl;
 write_out<<" $$$$  /$$$$ |$$    |   $$      \\ $$ |       $$$$$$  |$$ |  "<<endl;
 write_out<<" $$ $$ $$/$$ |$$$$$/     $$$$$$  |$$ |   __  /    $$ |$$ |   "<<endl;
 write_out<<" $$ |$$$/ $$ |$$ |_____ /  \\__$$ |$$ \\__/  |/$$$$$$$ |$$ | "<<endl;
 write_out<<" $$ | $/  $$ |$$       |$$    $$/ $$    $$/ $$    $$ |$$ |   "<<endl;
 write_out<<" $$/      $$/ $$$$$$$$/  $$$$$$/   $$$$$$/   $$$$$$$/ $$/    "<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out<<endl;
 write_out<<"  Normal termination of MESCAL code          "<<endl;
 write_out<<endl;  
 write_out<<"----------------------------------------------------------------"<<endl;  
 write_out.close();
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
