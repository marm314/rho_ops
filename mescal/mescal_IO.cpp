#include"mescal.h"

// Read PDB file
void MESCAL::read_pdb_file(string name_pdb)
{
 int ichar,ichar1,Z=1,count_fragments=-1,old_fragment=-1,new_fragment,blank_spaces[14];
 double pos[3];
 string line,line_aux;
 bool space,space2;
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
    count_fragments++;
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
    fragments[count_fragments].atoms.push_back({Z,pos[0],pos[1],pos[2]});// Z, position
    nfragments++;
   }
   else
   {
    fragments[count_fragments].natoms++;                                        // Add an atom to the fragment
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
    fragments[count_fragments].atoms.push_back({Z,pos[0],pos[1],pos[2]});// Z, position  
   }
  }
 } 
 read_pdb.close();
}

// Read fragment file
void MESCAL::read_fragment_file(string name_frag,double **Im_frag,double **Urot,int &ifrag,int &Sum_Val_elect, double &Sum_atomic_pol)
{
 bool devItens=false,frag_file_good=true;
 int iindex,jindex,kindex,iatom,jatom,ialpha,jalpha;
 double tol2=pow(10.0e0,-2.0e0),fact_weight,Im_ref[3][3],alpha[3][3]={0.0e0},Temp_mat[3][3],*q_read;
 string line;
 if(ind_q)
 {
  fragments[ifrag].Pi=new double*[fragments[ifrag].natoms];
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   fragments[ifrag].Pi[iatom]=new double[fragments[ifrag].natoms];
   for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
   {fragments[ifrag].Pi[iatom][jatom]=0.0e0;}
  }
 }
 q_read=new double [fragments[ifrag].natoms];
 for(iindex=0;iindex<fragments[ifrag].natoms;iindex++){q_read[iindex]=0.0e0;}
 ifstream read_frag(name_frag);
 if(!read_frag.good()){cout<<"Warning! Unable to find the .dat file for fragment "<<setw(5)<<ifrag+1<<" "<<fragments[ifrag].name<<endl;frag_file_good=false;}
 while(getline(read_frag,line))
 {
  if(line.length()<25){line+="                         ";}
  if(line.substr(0,14)=="Polarizability")
  {
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<=iindex;jindex++){read_frag>>alpha[iindex][jindex];if(iindex!=jindex){alpha[jindex][iindex]=alpha[iindex][jindex];}}
   }
  }
  if(line.substr(0,25)=="Normalized inertia tensor")
  {
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++){read_frag>>Im_ref[iindex][jindex];}
   }
   devItens=false; 
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++)
    {
     Temp_mat[iindex][jindex]=0.0e0;
     if(abs(Im_ref[iindex][jindex]-Im_frag[iindex][jindex])>tol2){devItens=true;}
    }
   }
   if(devItens){cout<<"Comment: The Inert. tensor. of fragment "<<setw(5)<<ifrag+1<<" presents deviations >10^-2 w.r.t. reference."<<endl;}
  }
  if(line.substr(0,16)=="Mulliken Charges")
  {
   for(iindex=0;iindex<fragments[ifrag].natoms;iindex++){read_frag>>q_read[iindex];}
  }
  if(line.substr(0,14)=="Susceptibility")
  {
   if(ind_q)
   {
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<fragments[ifrag].natoms;jatom++){read_frag>>fragments[ifrag].Pi[iatom][jatom];}
    }
   }
  }
 }
 read_frag.close();
 if(frag_file_good)
 {
  // Transform alpha_read -> alpha_rot = U^T alpha U (first check Inertia tensor)
  for(iindex=0;iindex<3;iindex++)
  {
   for(jindex=0;jindex<3;jindex++)
   {
    for(kindex=0;kindex<3;kindex++)
    {
     Temp_mat[iindex][jindex]+=Urot[iindex][kindex]*alpha[kindex][jindex];
    }
   }
  }
  for(iindex=0;iindex<3;iindex++)
  {
   for(jindex=0;jindex<3;jindex++)
   {
    alpha[iindex][jindex]=0.0e0;
    for(kindex=0;kindex<3;kindex++)
    {
     alpha[iindex][jindex]+=Temp_mat[iindex][kindex]*Urot[jindex][kindex];
    }
    if(abs(alpha[iindex][jindex])<tol8){alpha[iindex][jindex]=0.0e0;}
   }
  }
 }
 // Transform alpha_rot  -> alpha_atomic (number of valence electrons or atomic polarizabilities used as weighting method)
 // Note: We are currently assuming that atoms are ordered in the same way in the PDB and in the fragment.dat file
 //       It is possible to improve this by using the rot matrix and the distances. 
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  if(part_val_e) // Use partition based on number of valence electrons
  {
   fact_weight=((double)Z2val_electrons(fragments[ifrag].atoms[iatom].Z)/(double)Sum_Val_elect);
  }
  else           // Use partition based on atomic polarizabilities
  {
   fact_weight=Z2atomic_pol(fragments[ifrag].atoms[iatom].Z)/Sum_atomic_pol;
  }
  for(ialpha=0;ialpha<3;ialpha++)
  {
   for(jalpha=0;jalpha<3;jalpha++)
   {
    fragments[ifrag].atoms[iatom].alpha[ialpha][jalpha]=fact_weight*alpha[ialpha][jalpha];
   }
  }
  fragments[ifrag].atoms[iatom].q_perm=q_read[iatom]; // Asign permanent charges to atoms (not used in SCF calc. and currently are the Mulliken ones) 
 }
 delete[] q_read; q_read=NULL;
}

// Print header ouput file
void MESCAL::init_output(string name_output)
{
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
 write_out.close();
}

// Print init SC procedure
void MESCAL::print_init_sc(string name_output)
{
 ofstream write_out(name_output,std::ios_base::app);
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"---------- Performing self-consistent procedure ----------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<endl;
 write_out<<setprecision(4)<<fixed<<scientific;
 write_out<<" Maxiter     "<<setw(20)<<maxiter<<endl;
 write_out<<" Threshold mu"<<setw(20)<<threshold_mu<<endl;
 write_out<<" Threshold  E"<<setw(20)<<threshold_E<<endl;
 if(ind_q)
 {
  write_out<<" Threshold  q"<<setw(20)<<threshold_q<<endl;
 }
 write_out<<" Screening r0 (au) "<<setw(14)<<r0<<endl;
 write_out<<" Damping weight mu "<<setw(14)<<w_mu<<endl;
 if(perm_q){write_out<<" Q_permanent option is ON"<<endl;}
 if(ind_q ){write_out<<" Q_induced option is ON"<<endl;}
 if(part_val_e){write_out<<" Partition of alpha using num. valence electrons is ON"<<endl;}
 else{write_out<<" Partition of alpha using atomic polarizabilities is ON"<<endl;}
 write_out<<endl;
 if(ind_q)
 {
  write_out<<"#   iter        Energy(au)          max(mu_diff)         max(q_diff)         E_diff(au)"<<endl;
 }
 else
 {
  write_out<<"#   iter        Energy(au)          max(mu_diff)         E_diff(au)"<<endl;
 }
 write_out.close();
}

// Print end SC procedure
void MESCAL::print_end_sc(string name_output)
{
 ofstream write_out(name_output,std::ios_base::app);
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"---------- Self-consistent procedure completed  ----------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<endl;
 write_out.close();
}

// Print iter information
void MESCAL::print_iter_info(string name_output)
{
 ofstream write_out(name_output,std::ios_base::app);
 write_out<<setprecision(8)<<fixed;
 if(ind_q)
 {
  if(iter!=0)
  {
   write_out<<setw(7)<<iter<<setw(20)<<Energy<<setw(20)<<mu_diff_max<<setw(20)<<q_diff_max<<setw(20)<<E_diff<<endl;
  }
  else
  {
   write_out<<setw(7)<<iter<<setw(20)<<Energy<<setw(20)<<0.0e0<<setw(20)<<0.0e0<<setw(20)<<E_diff<<endl;
  }
 }
 else
 {
  if(iter!=0)
  {
   write_out<<setw(7)<<iter<<setw(20)<<Energy<<setw(20)<<mu_diff_max<<setw(20)<<E_diff<<endl;
  }
  else
  {
   write_out<<setw(7)<<iter<<setw(20)<<Energy<<setw(20)<<0.0e0<<setw(20)<<E_diff<<endl;
  }
 }
 write_out.close();
}

// Print footer ouput file
void MESCAL::close_output(string name_output,string sha)
{
 int ifrag,iatom;
 ofstream write_out(name_output,std::ios_base::app);
 write_out<<setprecision(8)<<fixed;
 write_out<<endl;
 write_out<<"Final induced charges and dipoled"<<endl;
 write_out<<endl;
 write_out<<"#     Z                               R(au)                                     q_ind                                 mu_ind(au)"<<endl;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
  {
   write_out<<setw(7)<<fragments[ifrag].atoms[iatom].Z;
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].pos[0];
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].pos[1];
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].pos[2];
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].q_ind;
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].mu_ind[0];
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].mu_ind[1];
   write_out<<setw(20)<<fragments[ifrag].atoms[iatom].mu_ind[2]<<endl;
  }
 }
 write_out<<endl;
 write_out<<endl;
 write_out<<"Git sha: "<<sha<<endl;
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
