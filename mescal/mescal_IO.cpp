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
 bool devItens=false,frag_file_good=true,pimatrix_good=false,all_int=true;
 int iindex,jindex,kindex,iatom,jatom,ialpha,jalpha,imo,amo,nbasis=0,nocc=0;
 double tol2=pow(10.0e0,-2.0e0),fact_weight,Im_ref[3][3],alpha[3][3]={0.0e0},Temp_mat[3][3],*q_read,**Cartes_coord,***S_mat,*orb_ene,val,val2;
 string line;
 orb_ene=new double[1];orb_ene[0]=0.0e0;
 Cartes_coord=new double*[fragments[ifrag].natoms];
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  Cartes_coord[iatom]=new double[3];
  for(iindex=0;iindex<3;iindex++)
  {Cartes_coord[iatom][iindex]=0.0e0;}
 }
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
  if(line.length()<29){line+="                             ";}
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
   // Check Inertia tensor similarity
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
  if(line.substr(0,29)=="Current cartesian coordinates")
  {
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    for(iindex=0;iindex<3;iindex++){read_frag>>Cartes_coord[iatom][iindex];}
   }
  }
  if(line.substr(0,16)=="Mulliken Charges")
  {
   for(iindex=0;iindex<fragments[ifrag].natoms;iindex++){read_frag>>q_read[iindex];}
  }
  if(line.substr(0,25)=="Number of basis functions" && ind_q)
  {
   line=line.substr(44,line.length()-44);
   stringstream ss(line);
   ss>>nbasis;
   delete[] orb_ene;orb_ene=NULL;
   orb_ene=new double[nbasis];
  }
  if(line.substr(0,19)=="Number of electrons" && ind_q)
  {
   line=line.substr(44,line.length()-44);
   stringstream ss(line);
   ss>>nocc;
   nocc=nocc/2; 
  }
  if(line.substr(0,22)=="Alpha Orbital Energies" && ind_q && nbasis!=0)
  {
   for(imo=0;imo<nbasis;imo++)
   {
    read_frag>>orb_ene[imo];
   }
  }
  if(line.substr(0,14)=="Susceptibility")
  {
   pimatrix_good=true;
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
  // Allocate S_mat[natom][n_mo][n_mo]
  if(ind_q && !pimatrix_good)
  {
   S_mat=new double**[fragments[ifrag].natoms]; 
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    S_mat[iatom]=new double*[nbasis];
    for(imo=0;imo<nbasis;imo++)
    {
     S_mat[iatom][imo]=new double[nbasis];
     for(amo=0;amo<nbasis;amo++)
     {
      S_mat[iatom][imo][amo]=0.0e0;
     }
    }
   }
   // Read S_mat (int files) and check it
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    Z2label(fragments[ifrag].atoms[iatom].Z);
    jatom=iatom+1;
    string label_cc= static_cast<ostringstream*>( &(ostringstream() << jatom) )->str();
    label=label+label_cc+".int";
    ifstream read_int(label.c_str());
    if(!read_int.good() && all_int){cout<<"Warning! Unable to find the .int file "<<label+".int"<<endl;all_int=false;}
    while(getline(read_int,line))
    {
     if(line==" The Atomic Overlap Matrix:")
     {
      getline(read_int,line);
      getline(read_int,line);
      getline(read_int,line);
      for(imo=0;imo<nbasis;imo++)
      {
       for(amo=0;amo<=imo;amo++)
       {
        read_int>>S_mat[iatom][imo][amo];
        if(imo!=amo){S_mat[iatom][amo][imo]=S_mat[iatom][imo][amo];}
       }
      }
     }
    }  
    read_int.close();
   }
   for(imo=0;imo<nbasis;imo++)
   {
    for(amo=0;amo<=imo;amo++)
    {
     if(imo==amo){val2=1.0e0;}
     else{val2=0.0e0;}
     val=0.0e0; 
     for(iatom=0;iatom<fragments[ifrag].natoms-1;iatom++)
     {
      val+=S_mat[iatom][amo][imo];
     }
     S_mat[fragments[ifrag].natoms-1][amo][imo]=val2-val;
     S_mat[fragments[ifrag].natoms-1][imo][amo]=S_mat[fragments[ifrag].natoms-1][amo][imo];
    }
   }
   // Compute PI[iatom][jatom] = 4 sum _{ia} S_mat[iatom][i][a] S_mat[jatom][a][i] / ( e_i - e_a ) with i occ and a virtual
   if(all_int)
   {
    pimatrix_good=true;
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
     {
      fragments[ifrag].Pi[iatom][jatom]=0.0e0;
      for(imo=0;imo<nocc;imo++)
      {
       for(amo=nocc;amo<nbasis;amo++)
       {
        fragments[ifrag].Pi[iatom][jatom]+=S_mat[iatom][imo][amo]*S_mat[jatom][amo][imo]/(orb_ene[imo]-orb_ene[amo]); 
       }
      }
      fragments[ifrag].Pi[iatom][jatom]=4.0e0*fragments[ifrag].Pi[iatom][jatom];
     }
    }
    ofstream print_pi_mat((name_frag.substr(0,name_frag.length()-4)+".pi").c_str());
    print_pi_mat<<setprecision(10)<<fixed<<scientific;
    iindex=0;
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
     {
      print_pi_mat<<setw(25)<<fragments[ifrag].Pi[iatom][jatom];
      iindex++;
      if(iindex==5){iindex=0;print_pi_mat<<endl;}
     }
    } 
    print_pi_mat.close();
   }
   // Delete S_mat
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    for(imo=0;imo<nbasis;imo++)
    {
     delete[] S_mat[iatom][imo];S_mat[iatom][imo]=NULL;
    }
    delete[] S_mat[iatom];S_mat[iatom]=NULL;
   }
   delete[] S_mat;S_mat=NULL;
   if(!pimatrix_good){cout<<"Warning! Pi = 0; thus, q_ind = 0 (mu^CR _ind == mu^dipole _ind because alpha(PI) = 0)"<<endl;}
  }
  // If it we want charge redistribution model (ind_q = True), we must substract the alpha_c = alpha(Pi) contrib.
  if(ind_q && pimatrix_good)
  {
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++)
    {
     for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
     {
      for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
      {
       // alpha ^xy  = alpha ^xy - \sum ij Pi_ij x_i y_j
       alpha[iindex][jindex]-=fragments[ifrag].Pi[iatom][jatom]
                             *Cartes_coord[iatom][iindex]
                             *Cartes_coord[jatom][jindex];
      }
     }
    }
   }
  }
  // Transform alpha_read -> alpha_rot = U^T alpha U 
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
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  delete[] Cartes_coord[iatom];Cartes_coord[iatom]=NULL;
 }
 delete[] Cartes_coord;Cartes_coord=NULL;
 delete[] q_read; q_read=NULL;
 delete[] orb_ene;orb_ene=NULL;
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
 write_out<<" Damping weight    "<<setw(14)<<w_damp<<endl;
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
void MESCAL::close_output(string name_output)
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
 if(sha!="")
 {
  write_out<<endl;
  write_out<<"Git sha: "<<sha<<endl;
  write_out<<endl;
 }
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
