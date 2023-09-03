#include"mescal.h"

//Public functions.
MESCAL::MESCAL(){cout<<"Not allowed default constructor in MESCAL"<<endl;}
MESCAL::MESCAL(string name_output,string name_pdb)
{
 string line,line_aux;
 bool space,space2;
 int ifrag,iatom,ichar,ichar1,icoord,jcoord,pivot,Z=1,old_fragment=-1,new_fragment,blank_spaces[14],order[3];
 double mass,mass_tot,Norm,Norm_saved,pivot_doub,pos[3],Rcm[3],**Im,**Urot;
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
 nfragments=0;
 ifstream read_pdb(name_pdb);
 while(getline(read_pdb,line))
 {
  // Find blank spaces in the lines of the PDB file
  ichar1=0;space=false;
  for(ichar=0;ichar<(int)line.length();ichar++)
  {
   if(line[ichar]==' ')
   {
    if(!space)
    {
     blank_spaces[ichar1]=ichar;
     ichar1++;
    }
    space=true;
   }
   else
   {
    if(space)
    {
     space=false;
     blank_spaces[ichar1]=ichar-1;
     ichar1++;
    }
   }
   if(ichar1==14){break;}
  }
  // Store information
  if(line.length()>27 && line.substr(0,4)=="ATOM")
  {
   stringstream ss(line.substr(blank_spaces[6],blank_spaces[9]-blank_spaces[6]+1));
   ss>>new_fragment;
   if(new_fragment!=old_fragment)
   {
    line_aux=line.substr(blank_spaces[4],blank_spaces[7]-blank_spaces[4]+1);
    line_aux.erase(std::remove_if(line_aux.begin(),line_aux.end(),::isspace),line_aux.end());
    old_fragment=new_fragment;
    fragments.push_back({line_aux,1});                                         // name fragment, natoms (init.)
    // Here add atomic info to push_back
    line_aux=line.substr(blank_spaces[2],blank_spaces[5]-blank_spaces[2]+1);
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
    line_aux=line.substr(blank_spaces[8],blank_spaces[11]-blank_spaces[8]+1);
    stringstream ss1(line_aux);
    ss1>>pos[0];
    pos[0]=pos[0]*Angs2au;
    line_aux=line.substr(blank_spaces[10],blank_spaces[13]-blank_spaces[10]+1);
    stringstream ss2(line_aux);
    ss2>>pos[1];
    pos[1]=pos[1]*Angs2au;
    line_aux=line.substr(blank_spaces[12],line.length()-blank_spaces[12]+1);
    space=false;space2=false;ichar1=0;
    for(ichar=0;ichar<(int)line_aux.length();ichar++)
    {
     if(line_aux[ichar]==' ')
     {
      space=true;
      if(space2)
      {
       ichar1=ichar;
       break;
      }
     }
     else
     {
      space=false;
      space2=true;
     }
    }
    if(ichar1!=0){line_aux=line_aux.substr(0,ichar1);}
    stringstream ss3(line_aux);
    ss3>>pos[2];
    pos[2]=pos[2]*Angs2au;
    fragments[old_fragment-1].atoms.push_back({Z,0.0,pos[0],pos[1],pos[2],0.0,0.0,0.0});// Z, charge, pos, dipole  
    nfragments++;
   }
   else
   {
    fragments[old_fragment-1].natoms++;                                        // Add an atom to the fragment
    // Here add atomic info to push_back
    // Det Z
    line_aux=line.substr(blank_spaces[2],blank_spaces[5]-blank_spaces[2]+1);
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
    line_aux=line.substr(blank_spaces[8],blank_spaces[11]-blank_spaces[8]+1);
    stringstream ss1(line_aux);
    ss1>>pos[0];
    pos[0]=pos[0]*Angs2au;
    line_aux=line.substr(blank_spaces[10],blank_spaces[13]-blank_spaces[10]+1);
    stringstream ss2(line_aux);
    ss2>>pos[1];
    pos[1]=pos[1]*Angs2au;
    line_aux=line.substr(blank_spaces[12],line.length()-blank_spaces[12]+1);
    space=false;space2=false;ichar1=0;
    for(ichar=0;ichar<(int)line_aux.length();ichar++)
    {
     if(line_aux[ichar]==' ')
     {
      space=true;
      if(space2)
      {
       ichar1=ichar;
       break;
      }
     }
     else
     {
      space=false;
      space2=true;
     }
    }
    if(ichar1!=0){line_aux=line_aux.substr(0,ichar1);}
    stringstream ss3(line_aux);
    ss3>>pos[2];
    pos[2]=pos[2]*Angs2au;
    fragments[old_fragment-1].atoms.push_back({Z,0.0,pos[0],pos[1],pos[2],0.0,0.0,0.0});// Z, charge, pos, dipole  
   }
  }
 } 
 read_pdb.close();

// Read here fragment.dat files to store alpha and rot matrices before computing U alpha U^T for each fragment

 write_out<<endl;
 write_out<<" Fragments read from the PDB file (distances in Bohr) "<<endl;
 write_out<<endl;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  write_out<<" Fragment "<<ifrag+1<<endl;
  for(icoord=0;icoord<3;icoord++){Rcm[icoord]=0.0e0;}
  mass_tot=0.0e0;
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   write_out<<"  Atom ";
   write_out<<setw(5)<<fragments[ifrag].atoms[iatom].Z;
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[0];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[1];
   write_out<<setw(17)<<fragments[ifrag].atoms[iatom].pos[2]<<endl;
   mass=Z2mass(fragments[ifrag].atoms[iatom].Z);
   for(icoord=0;icoord<3;icoord++){Rcm[icoord]+=mass*fragments[ifrag].atoms[iatom].pos[icoord];}
   mass_tot+=mass;
  }
  for(icoord=0;icoord<3;icoord++){Rcm[icoord]=Rcm[icoord]/(mass_tot+tol8);}
  write_out<<"  Center of Mass "<<setw(20)<<Rcm[0]<<setw(20)<<Rcm[1]<<setw(20)<<Rcm[2]<<endl;
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
 }
 for(icoord=0;icoord<3;icoord++){delete[] Urot[icoord];Urot[icoord]=NULL; delete[] Im[icoord];Im[icoord]=NULL;}
 delete[] Urot; Urot=NULL; delete[] Im; Im=NULL;

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
