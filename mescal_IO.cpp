#include"mescal.h"

// Read PDB file
void MESCAL::read_pdb_file(string name_pdb)
{
 int ichar,ichar1,Z=1,old_fragment=-1,new_fragment,blank_spaces[14];
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
}
