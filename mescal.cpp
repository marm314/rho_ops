#include"mescal.h"

//Public functions.
MESCAL::MESCAL(){cout<<"Not allowed default constructor in MESCAL"<<endl;}
MESCAL::MESCAL(string name_output,string name_pdb)
{
 string line;
 int ifrag,iatom,old_fragment=-1,new_fragment;
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
 write_out<<setprecision(10)<<fixed<<scientific;
 nfragments=0;
 ifstream read_pdb(name_pdb);
 while(getline(read_pdb,line))
 {
  if(line.length()>27)
  {
   stringstream ss(line.substr(20,7)); // This is not standard take care and first check spaces for each line to do this.
   ss>>new_fragment;
   if(new_fragment!=old_fragment)
   {
    old_fragment=new_fragment;
    fragments.push_back({"a",1});                                              // name fragment, natoms 
    fragments[old_fragment-1].atoms.push_back({0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});// Z, charge, pos, dipole  
    nfragments++;
   }
   else
   {
    fragments[old_fragment-1].natoms++;
    fragments[old_fragment-1].atoms.push_back({0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});// Z, charge, pos, dipole  
   }
  }
 } 
 read_pdb.close();
 write_out<<endl;
 write_out<<" Fragments read from the PDB file"<<endl;
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
