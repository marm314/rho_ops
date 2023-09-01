#include"mescal.h"

//Public functions.
MESCAL::MESCAL(){cout<<"Not allowed default constructor in MESCAL"<<endl;}
MESCAL::MESCAL(string name_output,string name_pdb)
{
 string line,line_aux;
 bool space,space2;
 int ifrag,iatom,ichar,ichar1,Z=1,old_fragment=-1,new_fragment,blank_spaces[14];
 double pos[3];
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
 }



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


// Bellow wes point we have functions that only work as libraries to asing values. 
// Atomic symbol -> Z 
void MESCAL::Asymbol2Z(int &Z, string symbol)
{
 if(symbol.length()==2)
 {
  if(symbol=="He"){Z=2;}
  else if(symbol=="Li"){Z=3;}
  else if(symbol=="Be"){Z=4;}
  else if(symbol=="Ne"){Z=10;}
  else if(symbol=="Na"){Z=11;}
  else if(symbol=="Mg"){Z=12;}
  else if(symbol=="Al"){Z=13;}
  else if(symbol=="Si"){Z=14;}
  else if(symbol=="Cl"){Z=17;}
  else if(symbol=="Ar"){Z=18;}
  else if(symbol=="Ca"){Z=20;}
  else if(symbol=="Sc"){Z=21;}
  else if(symbol=="Ti"){Z=22;}
  else if(symbol=="Cr"){Z=24;}
  else if(symbol=="Mn"){Z=25;}
  else if(symbol=="Fe"){Z=26;}
  else if(symbol=="Co"){Z=27;}
  else if(symbol=="Ni"){Z=28;}
  else if(symbol=="Cu"){Z=29;}
  else if(symbol=="Zn"){Z=30;}
  else if(symbol=="Ga"){Z=31;}
  else if(symbol=="Ge"){Z=32;}
  else if(symbol=="As"){Z=33;}
  else if(symbol=="Se"){Z=34;}
  else if(symbol=="Br"){Z=35;}
  else if(symbol=="Kr"){Z=36;}
  else if(symbol=="Rb"){Z=37;}
  else if(symbol=="Sr"){Z=38;}
  else if(symbol=="Zr"){Z=40;}
  else if(symbol=="Nb"){Z=41;}
  else if(symbol=="Mo"){Z=42;}
  else if(symbol=="Tc"){Z=43;}
  else if(symbol=="Ru"){Z=44;}
  else if(symbol=="Rh"){Z=45;}
  else if(symbol=="Pd"){Z=46;}
  else if(symbol=="Ag"){Z=47;}
  else if(symbol=="Cd"){Z=48;}
  else if(symbol=="In"){Z=49;}
  else if(symbol=="Sn"){Z=50;}
  else if(symbol=="Sb"){Z=51;}
  else if(symbol=="Te"){Z=52;}
  else if(symbol=="Xe"){Z=54;}
  else if(symbol=="Cs"){Z=55;}
  else if(symbol=="Ba"){Z=56;}
  else if(symbol=="La"){Z=57;}
  else if(symbol=="Ce"){Z=58;}
  else if(symbol=="Pr"){Z=59;}
  else if(symbol=="Nd"){Z=60;}
  else if(symbol=="Pm"){Z=61;}
  else if(symbol=="Sm"){Z=62;}
  else if(symbol=="Eu"){Z=63;}
  else if(symbol=="Gd"){Z=64;}
  else if(symbol=="Tb"){Z=65;}
  else if(symbol=="Dy"){Z=66;}
  else if(symbol=="Ho"){Z=67;}
  else if(symbol=="Er"){Z=68;}
  else if(symbol=="Tm"){Z=69;}
  else if(symbol=="Yb"){Z=70;}
  else if(symbol=="Lu"){Z=71;}
  else if(symbol=="Hf"){Z=72;}
  else if(symbol=="Ta"){Z=73;}
  else if(symbol=="Re"){Z=75;}
  else if(symbol=="Os"){Z=76;}
  else if(symbol=="Ir"){Z=77;}
  else if(symbol=="Pt"){Z=78;}
  else if(symbol=="Au"){Z=79;}
  else if(symbol=="Hg"){Z=80;}
  else if(symbol=="Tl"){Z=81;}
  else if(symbol=="Pb"){Z=82;}
  else if(symbol=="Bi"){Z=83;}
  else if(symbol=="Po"){Z=84;}
  else if(symbol=="At"){Z=85;}
  else if(symbol=="Rn"){Z=86;}
  else if(symbol=="Fr"){Z=87;}
  else if(symbol=="Ra"){Z=88;}
  else if(symbol=="Ac"){Z=89;}
  else if(symbol=="Th"){Z=90;}
  else if(symbol=="Pa"){Z=91;}
  else if(symbol=="Np"){Z=93;}
  else if(symbol=="Pu"){Z=94;}
  else if(symbol=="Am"){Z=95;}
  else if(symbol=="Cm"){Z=96;}
  else if(symbol=="Bk"){Z=97;}
  else if(symbol=="Cf"){Z=98;}
  else if(symbol=="Es"){Z=99;}
  else if(symbol=="Fm"){Z=100;}
  else if(symbol=="Md"){Z=101;}
  else if(symbol=="No"){Z=102;}
  else if(symbol=="Lr"){Z=103;}
  else if(symbol=="Rf"){Z=104;}
  else if(symbol=="Db"){Z=105;}
  else if(symbol=="Sg"){Z=106;}
  else if(symbol=="Bh"){Z=107;}
  else if(symbol=="Hs"){Z=108;}
  else if(symbol=="Mt"){Z=109;}
  else if(symbol=="Ds"){Z=110;}
  else if(symbol=="Rg"){Z=111;}
  else if(symbol=="Cn"){Z=112;}
  else if(symbol=="Nh"){Z=113;}
  else if(symbol=="Fl"){Z=114;}
  else if(symbol=="Mc"){Z=115;}
  else if(symbol=="Lv"){Z=116;}
  else if(symbol=="Ts"){Z=117;}
  else if(symbol=="Og"){Z=118;}
  else{cout<<"Warning! Atomic number not found for Symbol "<<symbol<<endl;}
 }
 else
 {
  if(symbol=="H"){Z=1;}
  else if(symbol=="B"){Z=5;}
  else if(symbol=="C"){Z=6;}
  else if(symbol=="N"){Z=7;}
  else if(symbol=="O"){Z=8;}
  else if(symbol=="F"){Z=9;}
  else if(symbol=="P"){Z=15;}
  else if(symbol=="S"){Z=16;}
  else if(symbol=="K"){Z=19;}
  else if(symbol=="V"){Z=23;}
  else if(symbol=="Y"){Z=39;}
  else if(symbol=="I"){Z=53;}
  else if(symbol=="W"){Z=74;}
  else if(symbol=="U"){Z=92;}
  else{cout<<"Warning! Atomic number not found for Symbol "<<symbol<<endl;}
 }
} 
