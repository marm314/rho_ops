#include"Mescal.h"

// Read PDB file
void Mescal::read_pdb_file(string name_pdb)
{
 int Z=1,count_fragments=-1,old_fragment=-1,new_fragment;
 double pos[3];
 string line,line_aux;
 ifstream read_pdb(name_pdb);
 while(getline(read_pdb,line))
 {
  if(line.length()>4)
  {
   if(line.substr(0,4)=="ATOM")
   {
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
  }
 } 
 read_pdb.close();
}

// Read fragment file
void Mescal::read_fragment_file(string name_frag,double **Im_frag,double **Urot2align,int &ifrag,int &Sum_Val_elect, double &Sum_atomic_pol)
{
 bool devItens=false,frag_file_good=true,pimatrix_good=false,all_int=true,Itensor=false;
 int iindex,jindex,kindex,iatom,jatom,ialpha,jalpha,imo,jmo,amo,bmo,ipair,jpair,nbasis=0,nocc=0,nvir=0,npair=1,npair_read;
 int *Zfrag;
 double *q_read,**Cartes_coord,***S_mat,**inv_ApB_mat,**U_cphf_cpks;
 double val,tol2=pow(10.0e0,-2.0e0),fact_weight,Im_ref[3][3],alpha[3][3]={0.0e0},Temp_mat[3][3]={0.0e0};
 string line;
 inv_ApB_mat=new double*[npair];
 inv_ApB_mat[0]=new double[npair];
 //orb_ene=new double[1];orb_ene[0]=0.0e0;
 U_cphf_cpks=new double*[fragments[ifrag].natoms];
 Cartes_coord=new double*[fragments[ifrag].natoms];
 Zfrag=new int[fragments[ifrag].natoms];
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  U_cphf_cpks[iatom]=new double[npair];
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
  if(line.substr(0,14)=="Atomic numbers")
  {
   ofstream tmp("tmp");
   getline(read_frag,line);
   do
   {
    tmp<<line<<endl;
    getline(read_frag,line);
   }while(line[0]==' ');
   tmp.close();
   ifstream tmp2("tmp");
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    tmp2>>Zfrag[iatom];
   } 
   tmp2.close();
   system("/bin/rm -rf tmp");
  }
  if(line.substr(0,29)=="Current cartesian coordinates")
  {
   ofstream tmp("tmp");
   getline(read_frag,line);
   do
   {
    tmp<<line<<endl;
    getline(read_frag,line);
   }while(line[0]==' ');
   tmp.close();
   ifstream tmp2("tmp");
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    for(iindex=0;iindex<3;iindex++){tmp2>>Cartes_coord[iatom][iindex];}
   }
   tmp2.close();
   system("/bin/rm -rf tmp");
  }
  if(line.substr(0,19)=="Number of electrons" && ind_q)
  {
   line=line.substr(44,line.length()-44);
   stringstream ss(line);
   ss>>nocc;
   nocc=nocc/2; 
  }
  if(line.substr(0,14)=="Polarizability")
  {
   ofstream tmp("tmp");
   getline(read_frag,line);
   do
   {
    tmp<<line<<endl;
    getline(read_frag,line);
   }while(line[0]==' ');
   tmp.close();
   ifstream tmp2("tmp");
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<=iindex;jindex++){tmp2>>alpha[iindex][jindex];if(iindex!=jindex){alpha[jindex][iindex]=alpha[iindex][jindex];}}
   }
   tmp2.close();
   system("/bin/rm -rf tmp"); 
  }
  if(line.substr(0,16)=="Mulliken Charges")
  {
   ofstream tmp("tmp");
   getline(read_frag,line);
   do
   {
    tmp<<line<<endl;
    getline(read_frag,line);
   }while(line[0]==' ');
   tmp.close();
   ifstream tmp2("tmp");
   for(iindex=0;iindex<fragments[ifrag].natoms;iindex++){tmp2>>q_read[iindex];}
   tmp2.close();
   system("/bin/rm -rf tmp");
  }
  if(line.substr(0,25)=="Normalized inertia tensor")
  {
   Itensor=true;
   ofstream tmp("tmp");
   getline(read_frag,line);
   do
   {
    tmp<<line<<endl;
    getline(read_frag,line);
   }while(line[0]==' ');
   tmp.close();
   ifstream tmp2("tmp");
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++){tmp2>>Im_ref[iindex][jindex];}
   }
   tmp2.close();
   system("/bin/rm -rf tmp"); 
   // Check Inertia tensor similarity
   devItens=false; 
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++)
    {
     if(abs(Im_ref[iindex][jindex]-Im_frag[iindex][jindex])>tol2){devItens=true;}
    }
   }
   if(devItens){cout<<"Comment: The read Inertia tensor of fragment "<<setw(5)<<ifrag+1<<" presents deviations >10^-2 w.r.t. reference."<<endl;}
  }
  if(line.substr(0,14)=="Susceptibility") // It is the last quantity, we can read it like this
  {
   pimatrix_good=true;
   if(ind_q)
   {
    ofstream tmp("tmp");
    getline(read_frag,line);
    do
    {
     tmp<<line<<endl;
     getline(read_frag,line);
    }while(line[0]==' ');
    tmp.close();
    ifstream tmp2("tmp");
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<fragments[ifrag].natoms;jatom++){tmp2>>fragments[ifrag].Pi[iatom][jatom];}
    }
    tmp2.close();
    system("/bin/rm -rf tmp"); 
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<=iatom;jatom++)
     {
      fragments[ifrag].Pi[iatom][jatom]=0.5e0*(fragments[ifrag].Pi[iatom][jatom]+fragments[ifrag].Pi[jatom][iatom]);
      if(iatom!=jatom)
      {
       fragments[ifrag].Pi[jatom][iatom]=fragments[ifrag].Pi[iatom][jatom];
      }
     }
    }
    for(iatom=0;iatom<fragments[ifrag].natoms-1;iatom++)
    {
     val=0.0e0;
     for(jatom=0;jatom<fragments[ifrag].natoms-1;jatom++)
     {
      val+=fragments[ifrag].Pi[jatom][iatom];
     }
     fragments[ifrag].Pi[fragments[ifrag].natoms-1][iatom]=-val;
     fragments[ifrag].Pi[iatom][fragments[ifrag].natoms-1]=fragments[ifrag].Pi[fragments[ifrag].natoms-1][iatom];
    }
   }
  }
 }
 read_frag.close();
 if(frag_file_good)
 {
  // Build and check the inertia tensor if it was not provided
  if(!Itensor)
  {
   double Rcm[3],**Im,**Urot;
   Im=new double*[3];
   Urot=new double*[3];
   for(iindex=0;iindex<3;iindex++)
   { 
    Im[iindex]=new double[3];
    Urot[iindex]=new double[3];
    for(jindex=0;jindex<3;jindex++)
    {Im[iindex][jindex]=0.0e0;Urot[iindex][jindex]=0.0e0;}
   }
   Frag_T_inertia_compare(ifrag,Cartes_coord,Zfrag,Rcm,Im,Urot);
   // Check Inertia tensor similarity
   devItens=false; 
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++)
    {
     if(abs(Im[iindex][jindex]-Im_frag[iindex][jindex])>tol2){devItens=true;}
    }
   }
   if(devItens){cout<<"Comment: The computed Inertia tensor of fragment "<<setw(5)<<ifrag+1<<" presents deviations >10^-2 w.r.t. reference."<<endl;}
   for(iindex=0;iindex<3;iindex++)
   { 
    delete[] Im[iindex];Im[iindex]=NULL;
    delete[] Urot[iindex];Urot[iindex]=NULL;
   }
   delete[] Im;Im=NULL;
   delete[] Urot;Urot=NULL;
  }
  // Allocate S_mat[natom][n_mo][n_mo]
  if(ind_q && !pimatrix_good)
  {
   ifstream read_invApB("inv_apb_mat");
   if(!read_invApB.good()){cout<<"Warning! Unable to find the inv_apb_mat file"<<endl;}
   else
   {	    
    cout<<"Reading the inv_apb_mat file"<<endl;
    delete[] inv_ApB_mat[0];inv_ApB_mat[0]=NULL;
    delete[] inv_ApB_mat;inv_ApB_mat=NULL;
    read_invApB>>npair_read>>jindex;
    read_invApB.close();
    inv_ApB_mat=new double*[npair_read];
    read_invApB.open("inv_apb_mat");
    getline(read_invApB,line); 
    for(iindex=0;iindex<npair_read;iindex++)
    {
     inv_ApB_mat[iindex]=new double[npair_read];
     for(jindex=0;jindex<npair_read;jindex++)
     {
      getline(read_invApB,line);
      for(kindex=0;kindex<(int)line.length();kindex++)
      {
       if(line[kindex]=='d' || line[kindex]=='D'){line[kindex]='E';}
      }
      stringstream ss(line);
      ss>>val;
      inv_ApB_mat[iindex][jindex]=val;
     }
    }
    read_invApB>>nbasis;
    read_invApB.close();
    nvir=nbasis-nocc;npair=nocc*nvir;
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     delete[] U_cphf_cpks[iatom];U_cphf_cpks[iatom]=NULL;
    }
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     U_cphf_cpks[iatom]=new double[npair];
     for(iindex=0;iindex<npair;iindex++){U_cphf_cpks[iatom][iindex]=0.0e0;}
    }
   }  
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
   // Read S_mat (int files) and make the sum produce the identity matrix
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    Z2label(fragments[ifrag].atoms[iatom].Z);
    jatom=iatom+1;
    ostringstream convert;
    convert<<jatom;
    string label_cc=convert.str();
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
     val=0.0e0; 
     for(iatom=0;iatom<fragments[ifrag].natoms-1;iatom++)
     {
      val+=S_mat[iatom][amo][imo];
     }
     if(imo==amo)
     {
      S_mat[fragments[ifrag].natoms-1][amo][imo]=1.0e0-val;
     }
     else
     {
      S_mat[fragments[ifrag].natoms-1][amo][imo]=-val;
     }
     S_mat[fragments[ifrag].natoms-1][imo][amo]=S_mat[fragments[ifrag].natoms-1][amo][imo];
    }
   }
   // We do not use the non-interacting approximation:
   // This is crap -> PI[iatom][jatom] = 4 sum _{ia} S_mat[iatom][i][a] S_mat[jatom][a][i] / ( e_i - e_a ) with i occ and a virtual
   // We use the CPHF/CPKS Pi:
   // Pi = 4 sum_{ia} S_mat[iatom][i][a] U_cphf_cpks[jatom][i][a] = 4 sum_{ia} S_mat[iatom][i][a] U_cphf_cpks[jatom][ ipair = (i-1)*nvir + a ]   
   if(all_int)
   {
    pimatrix_good=true;
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     ipair=0;
     for(imo=0;imo<nocc;imo++)
     {
      for(amo=nocc;amo<nbasis;amo++)
      {
       jpair=0;
       for(jmo=0;jmo<nocc;jmo++)
       {
        for(bmo=nocc;bmo<nbasis;bmo++)
        {
         U_cphf_cpks[iatom][ipair]+=inv_ApB_mat[ipair][jpair]*S_mat[iatom][jmo][bmo];
	 jpair++;
	}
       }
       ipair++;	
      }
     }
    }
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
     {
      fragments[ifrag].Pi[iatom][jatom]=0.0e0;
      ipair=0;
      for(imo=0;imo<nocc;imo++)
      {
       for(amo=nocc;amo<nbasis;amo++)
       {
        fragments[ifrag].Pi[iatom][jatom]+=S_mat[iatom][imo][amo]*U_cphf_cpks[jatom][ipair];
        ipair++;	
        //fragments[ifrag].Pi[iatom][jatom]+=S_mat[iatom][imo][amo]*S_mat[jatom][amo][imo]/(orb_ene[imo]-orb_ene[amo]); 
       }
      }
      fragments[ifrag].Pi[iatom][jatom]=4.0e0*fragments[ifrag].Pi[iatom][jatom];
     }
    }
    for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
    {
     for(jatom=0;jatom<=iatom;jatom++)
     {
      fragments[ifrag].Pi[iatom][jatom]=0.5e0*(fragments[ifrag].Pi[iatom][jatom]+fragments[ifrag].Pi[jatom][iatom]);
      if(iatom!=jatom)
      {
       fragments[ifrag].Pi[jatom][iatom]=fragments[ifrag].Pi[iatom][jatom];
      }
     }
    }
    for(iatom=0;iatom<fragments[ifrag].natoms-1;iatom++)
    {
     val=0.0e0;
     for(jatom=0;jatom<fragments[ifrag].natoms-1;jatom++)
     {
      val+=fragments[ifrag].Pi[jatom][iatom];
     }
     fragments[ifrag].Pi[fragments[ifrag].natoms-1][iatom]=-val;
     fragments[ifrag].Pi[iatom][fragments[ifrag].natoms-1]=fragments[ifrag].Pi[fragments[ifrag].natoms-1][iatom];
    }
    if(!mute)
    { 
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
   ofstream print_alpha_mat;
   if(!mute)
   {
    print_alpha_mat.open((name_frag.substr(0,name_frag.length()-4)+".alpha").c_str());
    print_alpha_mat<<setprecision(10)<<fixed<<scientific;
   }
   for(iindex=0;iindex<3;iindex++)
   {
    for(jindex=0;jindex<3;jindex++)
    {
     for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
     {
      for(jatom=0;jatom<fragments[ifrag].natoms;jatom++)
      {
       // alpha ^xy (Pi)  = \sum AB Pi_AB x_A y_B
       Temp_mat[iindex][jindex]+=fragments[ifrag].Pi[iatom][jatom]
	                          *Cartes_coord[iatom][iindex]
				  *Cartes_coord[jatom][jindex];
      }
     }
     // alpha ^xy = alpha ^xy,QM  - alpha ^xy (Pi)
     alpha[iindex][jindex]-=Temp_mat[iindex][jindex];
     if(abs(alpha[iindex][jindex])<tol8){alpha[iindex][jindex]=0.0e0;}
     if(abs(Temp_mat[iindex][jindex])<tol8){Temp_mat[iindex][jindex]=0.0e0;}
     if(!mute)
     {
      print_alpha_mat<<setw(25)<<Temp_mat[iindex][jindex];
     }
    }
    if(!mute)
    {
     print_alpha_mat<<endl;
    }
   }
   if(!mute)
   {
    print_alpha_mat.close();
   }
  }
  // Transform alpha_read -> alpha_rot = U^T alpha U 
  for(iindex=0;iindex<3;iindex++)
  {
   for(jindex=0;jindex<3;jindex++)
   {
    Temp_mat[iindex][jindex]=0.0e0;
    for(kindex=0;kindex<3;kindex++)
    {
     Temp_mat[iindex][jindex]+=Urot2align[iindex][kindex]*alpha[kindex][jindex];
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
     alpha[iindex][jindex]+=Temp_mat[iindex][kindex]*Urot2align[jindex][kindex];
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
 delete[] Zfrag;Zfrag=NULL;
 delete[] Cartes_coord;Cartes_coord=NULL;
 delete[] q_read; q_read=NULL;
 for(iindex=0;iindex<npair;iindex++)
 {delete[] inv_ApB_mat[iindex];inv_ApB_mat[iindex]=NULL;}
 delete[] inv_ApB_mat;inv_ApB_mat=NULL;
 for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
 {
  delete[] U_cphf_cpks[iatom];U_cphf_cpks[iatom]=NULL;
 }
 delete[] U_cphf_cpks;U_cphf_cpks=NULL;
}

// Print header ouput file
void Mescal::init_output()
{
 ofstream write_out(mescal_ofile);
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"------                 Mescal ( in C++ )                   -----"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"--          Developed by: Dr. M. Rodriguez-Mayorga            --"<<endl;
 write_out<<"----               email: marm3.14@gmail.com               -----"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"--  MicroElectroStatic Calculations (Mescal)               -----"<<endl;
 write_out<<"--  Based on Mescal (in Fortran) code by                   -----"<<endl;
 write_out<<"--  Dr. Gabriele D'Avino (2013-2015)                       -----"<<endl;
 write_out<<"--  email: gabriele.davino@gmail.com                       -----"<<endl;
 write_out<<"--  download: https://gitlab.com/taphino/mescal            -----"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<endl;
 write_out.close();
}

// Print init SC procedure
void Mescal::print_init_sc()
{
 ofstream write_out(mescal_ofile,std::ios_base::app);
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"---------- Performing self-consistent procedure ----------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<endl;
 write_out<<setprecision(4)<<fixed<<scientific;
 if(nactive!=-1)
 {
  write_out<<" Nactive Frag"<<setw(20)<<nactive<<endl;
  write_out<<" Radius Frag "<<setw(20)<<radius<<endl;
 }
 if(ifrac_deact!=-1)
 {
  write_out<<" Deact Frag  "<<setw(20)<<ifrac_deact+1<<endl;
 }
 write_out<<" Maxiter     "<<setw(20)<<maxiter<<endl;
 write_out<<" Threshold mu"<<setw(20)<<threshold_mu<<endl;
 write_out<<" Threshold  E"<<setw(20)<<threshold_E<<endl;
 if(ind_q)
 {
  write_out<<" Threshold  q"<<setw(20)<<threshold_q<<endl;
 }
 write_out<<" Screening r0 (au) "<<setw(14)<<r0<<endl;
 write_out<<" Damping weight    "<<setw(14)<<w_damp<<endl;
 write_out<<" OMP threads       "<<setw(14)<<nthread<<endl;
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
void Mescal::print_end_sc()
{
 ofstream write_out(mescal_ofile,std::ios_base::app);
 write_out<<endl;
 write_out<<setprecision(8)<<fixed;
 write_out<<" Total SCS energy  "<<setw(14)<<Energy_old<<" a.u."<<endl;
 write_out<<"                   "<<setw(14)<<Energy_old*au2eV<<" eV"<<endl;
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<"---------- Self-consistent procedure completed  ----------------"<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out<<endl;
 write_out.close();
}

// Print iter information
void Mescal::print_iter_info()
{
 ofstream write_out(mescal_ofile,std::ios_base::app);
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
void Mescal::close_output()
{
 int ifrag,iatom;
 ofstream write_out(mescal_ofile,std::ios_base::app);
 write_out<<setprecision(8)<<fixed;
 write_out<<endl;
 write_out<<"Final induced charges and dipoled"<<endl;
 write_out<<endl;
 write_out<<"#     Z                               R(au)                                     q_ind                                mu_ind(au)"<<endl;
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
 write_out<<"  Normal termination of Mescal code          "<<endl;
 write_out<<endl;
 write_out<<"----------------------------------------------------------------"<<endl;
 write_out.close();
}

// Print charges ouput file (in Angstrom)
void Mescal::print_charges_file()
{
 int ifrag,iatom,icoord,n_charges=0;
 double distO,distAU,u_vec,norm_vec,q_charge[2]={0.0e0},coord_q[2][3]={0.0e0}; 
 distO=5.0e-4;distAU=2.0e0*distO*Angs2au;
 string mescal_charges=mescal_ofile.substr(0,mescal_ofile.length()-3)+"charges";
 ofstream write_charges(mescal_charges);
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    n_charges++;
   }
  }
 }
 write_charges<<setw(9)<<n_charges*2<<endl;
 write_charges<<setprecision(8)<<fixed;
 for(ifrag=0;ifrag<nfragments;ifrag++)
 {
  if(fragments[ifrag].active)
  {
   for(iatom=0;iatom<fragments[ifrag].natoms;iatom++)
   {
    norm_vec=0.0e0;
    for(icoord=0;icoord<3;icoord++)
    {
     norm_vec=norm_vec+pow(fragments[ifrag].atoms[iatom].mu_ind[icoord],2.0e0);
    }
    norm_vec=pow(norm_vec,0.5e0);
    for(icoord=0;icoord<3;icoord++)
    {
     u_vec=fragments[ifrag].atoms[iatom].mu_ind[icoord]/norm_vec;
     coord_q[0][icoord]= u_vec*distO+fragments[ifrag].atoms[iatom].pos[icoord]/Angs2au;
     coord_q[1][icoord]=-u_vec*distO+fragments[ifrag].atoms[iatom].pos[icoord]/Angs2au;   
    }
    q_charge[0]= norm_vec/distAU+0.5e0*fragments[ifrag].atoms[iatom].q_ind;
    q_charge[1]=-norm_vec/distAU+0.5e0*fragments[ifrag].atoms[iatom].q_ind;
    write_charges<<setprecision(8)<<fixed;
    write_charges<<setw(16)<<q_charge[0];
    for(icoord=0;icoord<3;icoord++)
    {
     write_charges<<setw(16)<<coord_q[0][icoord];
    }
    write_charges<<endl;
    write_charges<<setw(16)<<q_charge[1];
    for(icoord=0;icoord<3;icoord++)
    {
     write_charges<<setw(16)<<coord_q[1][icoord];
    }
    write_charges<<endl;
   }
  }
 }
 write_charges.close();
}
