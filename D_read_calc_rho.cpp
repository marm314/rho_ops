#include"D_read_calc_rho.h"
   ///////////////////////////
   //Declare Class Functions//
   ///////////////////////////

/////////////////////////////////////////
/////////////////////////////////////////
// Creation and anihilation of objects //
// ----------------------------------- //
/////////////////////////////////////////
/////////////////////////////////////////

//Public functions.
READ_FCHK_WFN::READ_FCHK_WFN(){cout<<"Not allowed default constructor in READ_FCHK_WFN"<<endl;}
//Build Density object from fchk or wfn (read all required variables)
READ_FCHK_WFN::READ_FCHK_WFN(string name,string name_log,bool WFN,bool log_file,bool cas,bool CM_in,int mult_in)
{
  int i,j,k;
  wfn=false;
  wfn=WFN;
  name_file=name;
  nelectrons=0;multiplicity=0;nprimitv=0;MO_coef=0;MO_beta_coef=0;
  prim_exp=0;smap=0;natoms=0;counter=0;counter2=0;
  cart_coor=false;nprimsh=false;satmap=false;exponents=false;mocoef=false;
  mo_beta_coef=false;open_shell=true;rhf=false;uhf=false;CI=false;CAS=false;
  correlated=false;virtuals=false;relaxed=false;error_opens_wfn=false;PS_bool=false;overlap=false;
  BETA_MOS=false;no_beta_wfn=false;scf_dens_found=false;wfx=false;im_wfn_wfx=false;
  if(!wfn)
  {
   string aux;
   cartes_coord=0;nalphael=0; nbetael=0;nbasisf=0;
   nindepbasisf=0;nshells=0; largcontr=0;
   npure_d=0;npure_f=0; Hi_ang_moment=0;nu_charge=0;
   stype=0;nprim_shell=0;contr_coef=0;spcontr_coef=0;
   nshells=0;SCF_rho=0;
   spin_SCF_rho=0;rho_CI=0;
   nu_ch=false;shtype=false;contract=false;contractSP=false;scfrho=false;
   spinscfrho=false;extra0=false;extra1=false;rho_CI_bool=false;
   activeSP=false;spinCIrho=false;
  //go is use to break if pure d or f shells are found
  //extras activate the destruction and control in extra0 case the use of SP
  //coefficients also extra1 controls the use of coefficients
 /////////////////////////////////////////////////////////////////////////////////
 //Read all information required
   input_fchk.open((name_file).c_str());
   if(input_fchk.is_open())
    { //Read information and castings
      while(getline(input_fchk,line))
      {
      //look for special phrases and activate keywords
       if(line.substr(0,15)=="Number of atoms")
       {line_fill(natoms,line);}
       else if(line.substr(0,12)=="Multiplicity")
       {line_fill(multiplicity,line);
        if(multiplicity!=1){open_shell=true;}
        else{open_shell=false;}
       }
       else if(line.substr(0,19)=="Number of electrons")
       {line_fill(nelectrons,line);}
       else if(line.substr(0,25)=="Number of alpha electrons")
       {line_fill(nalphael,line);}
       else if(line.substr(0,24)=="Number of beta electrons")
       {line_fill(nbetael,line);}
       else if(line.substr(0,25)=="Number of basis functions")
       {line_fill(nbasisf,line);}
       else if(line.substr(0,31)=="Number of independent functions" ||
               line.substr(0,31)=="Number of independant functions")
       {line_fill(nindepbasisf,line);}
       else if(line.substr(0,15)=="Nuclear charges")
       {line_fill(nu_charge,line);nu_ch=true;}
       else if(line.substr(0,29)=="Current cartesian coordinates")
       {line_fill(cartes_coord,line);cart_coor=true;}
       else if(line.substr(0,27)=="Number of contracted shells")
       {line_fill(nshells,line);}
       else if(line.substr(0,26)=="Number of primitive shells")
       {line_fill(nprimitv,line);}
       else if(line.substr(0,23)=="Pure/Cartesian d shells")
       {line_fill(npure_d,line);}
       else if(line.substr(0,23)=="Pure/Cartesian f shells")
       {line_fill(npure_f,line);}
       else if(line.substr(0,24)=="Highest angular momentum")
       {line_fill(Hi_ang_moment,line);}
       else if(line.substr(0,29)=="Largest degree of contraction")
       {line_fill(largcontr,line);}
       else if(line.substr(0,11)=="Shell types")
       {line_fill(stype,line);shtype=true;}
       else if(line.substr(0,30)=="Number of primitives per shell")
       {line_fill(nprim_shell,line);nprimsh=true;}
       else if(line.substr(0,17)=="Shell to atom map")
       {line_fill(smap,line);satmap=true;}
       else if(line.substr(0,19)=="Primitive exponents")
       {line_fill(prim_exp,line);exponents=true;}
       else if(line.substr(0,24)=="Contraction coefficients")
       {line_fill(contr_coef,line);contract=true;}
       else if(line.substr(0,31)=="P(S=P) Contraction coefficients")
       {line_fill(spcontr_coef,line);contractSP=true;extra0=true;}
       else if(line.substr(0,21)=="Alpha MO coefficients")
       {line_fill(MO_coef,line);mocoef=true;}
       else if(line.substr(0,20)=="Beta MO coefficients" || line.substr(0,20)=="B3TA MO coefficients")
       {line_fill(MO_beta_coef,line);mo_beta_coef=true;extra1=true;uhf=true;}
       else if(line.substr(0,17)=="Total SCF Density")
       {line_fill(SCF_rho,line);scfrho=true;scf_dens_found=true;}
       else if(line.substr(0,16)=="Spin SCF Density")
       {line_fill(spin_SCF_rho,line);spinscfrho=true;}
       else if(line.substr(0,16)=="Total CI Density" || line.substr(0,16)== "Total CC Density" ||
               line.substr(0,17)=="Total MP2 Density")
       {line_fill(rho_CI,line);rho_CI_bool=true;correlated=true;}
       else if(line.substr(0,15)=="Spin CI Density" || line.substr(0,15)=="Spin CC Density" ||
               line.substr(0,16)=="Spin MP2 Density")
       {line_fill(spin_CI_rho,line);spinCIrho=true;}
       else
       {//Casting strings from getline to double vectors and matrixes
       if(nu_ch)
       {counter=0;
        Nu_charge=new double[nu_charge];
        cast_doub(Nu_charge,line,counter,nu_charge);
        while(counter<nu_charge)
        {
         getline(input_fchk,line);
         cast_doub(Nu_charge,line,counter,nu_charge);
        }
        nu_ch=false;
       }
       else if(cart_coor)
       {counter=0;
        Cartesian_Coor=new double*[natoms];
        for(i=0;i<natoms;i++)
        {Cartesian_Coor[i]=new double[3];}
        AUX=new double[cartes_coord];
        cast_doub(AUX,line,counter,cartes_coord);
        while(counter<cartes_coord)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,cartes_coord);
        }
        counter=0;
        for(i=0;i<natoms;i++)
        {for(j=0;j<3;j++)
         {
          Cartesian_Coor[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        cart_coor=false;
       }
       else if(shtype)
       {counter=0;
        shell_type=new int[stype];
        cast_int(shell_type,line,counter,stype);
        while(counter<stype)
        {
         getline(input_fchk,line);
         cast_int(shell_type,line,counter,stype);
        }
        shtype=false;
       }
       else if(nprimsh)
       {counter=0;
        n_prim_per_shell=new int[nprim_shell];
        cast_int(n_prim_per_shell,line,counter,nprim_shell);
        while(counter<stype)
        {
         getline(input_fchk,line);
         cast_int(n_prim_per_shell,line,counter,nprim_shell);
        }
        nprimsh=false;
       }
       else if(satmap)
       {counter=0;
        shell_map=new int[smap];
        cast_int(shell_map,line,counter,smap);
        while(counter<stype)
        {
         getline(input_fchk,line);
         cast_int(shell_map,line,counter,smap);
        }
        satmap=false;
       }
       else if(exponents)
       {counter=0;
        Prim_exp=new double[prim_exp];
        cast_doub(Prim_exp,line,counter,prim_exp);
        while(counter<prim_exp)
        {
         getline(input_fchk,line);
         cast_doub(Prim_exp,line,counter,prim_exp);
        }
        exponents=false;
       }
       else if(contract)
       {counter=0;
        Contr_Coef=new double[contr_coef];
        cast_doub(Contr_Coef,line,counter,contr_coef);
        while(counter<contr_coef)
        {
         getline(input_fchk,line);
         cast_doub(Contr_Coef,line,counter,contr_coef);
        }
        contract=false;
       }
       else if(contractSP)
       {counter=0;
        SP_Contr_Coef=new double[spcontr_coef];
        cast_doub(SP_Contr_Coef,line,counter,spcontr_coef);
        while(counter<spcontr_coef)
        {
         getline(input_fchk,line);
         cast_doub(SP_Contr_Coef,line,counter,spcontr_coef);
        }
        contractSP=false;
       }
       else if(mocoef)
       {counter=0;
        MOcoefA=new double*[nbasisf];
        SPIN=new bool[nbasisf];
        for(i=0;i<nbasisf;i++)
        {
         MOcoefA[i]=new double[nbasisf];
         SPIN[i]=true;
        }
        AUX=new double[MO_coef];
        cast_doub(AUX,line,counter,MO_coef);
        while(counter<MO_coef)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,MO_coef);
        }
        counter=0;
        for(i=0;i<nbasisf;i++)
        {for(j=0;j<nbasisf;j++)
         {
          MOcoefA[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        mocoef=false;
       }
       else if(mo_beta_coef)
       {counter=0;
        MOcoefB=new double*[nbasisf];
        delete[] SPIN;
        SPIN=new bool[nbasisf+nbasisf];
        for(i=0;i<nbasisf+nbasisf;i++)
        {
         SPIN[i]=true;
         if(i%2 !=0)
         {SPIN[i]=false;}
        }
        for(i=0;i<nbasisf;i++)
        {MOcoefB[i]=new double[nbasisf];}
        AUX=new double[MO_beta_coef];
        cast_doub(AUX,line,counter,MO_beta_coef);
        while(counter<MO_beta_coef)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,MO_beta_coef);
        }
        counter=0;
        for(i=0;i<nbasisf;i++)
        {for(j=0;j<nbasisf;j++)
         {
          MOcoefB[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        mo_beta_coef=false;
       }
       else if(scfrho)
       {counter=1;
        Total_rho=new double*[nbasisf];
        for(i=0;i<nbasisf;i++)
        {
         Total_rho[i]=new double[counter];
         counter++;
        }
        counter=0;
        AUX=new double[SCF_rho];
        cast_doub(AUX,line,counter,SCF_rho);
        while(counter<SCF_rho)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,SCF_rho);
        }
        counter=0;
        for(i=0;i<nbasisf;i++)
        {for(j=0;j<i+1;j++)
         {
          Total_rho[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        scfrho=false;
       }
       else if(spinscfrho)
       {counter=1;
        Spin_rho=new double*[nbasisf];
        for(i=0;i<nbasisf;i++)
        {
         Spin_rho[i]=new double[counter];
         counter++;
        }
        counter=0;
        AUX=new double[spin_SCF_rho];
        cast_doub(AUX,line,counter,spin_SCF_rho);
        while(counter<spin_SCF_rho)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,spin_SCF_rho);
        }
        counter=0;
        for(i=0;i<nbasisf;i++)
        {for(j=0;j<i+1;j++)
         {
          Spin_rho[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        spinscfrho=false;
       }
       else if(rho_CI_bool)
       {
        if(scf_dens_found)
        {
         for(i=0;i<nbasisf;i++)
         {
          delete[] Total_rho[i];
          Total_rho[i]=NULL;
         }
         delete[] Total_rho;Total_rho=NULL;
        }
        counter=1;
        Total_rho=new double*[nbasisf];
        for(i=0;i<nbasisf;i++)
        {
         Total_rho[i]=new double[counter];
         counter++;
        }
        counter=0;
        AUX=new double[rho_CI];
        cast_doub(AUX,line,counter,rho_CI);
        while(counter<rho_CI)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,rho_CI);
        }
        counter=0;
        for(i=0;i<nbasisf;i++)
        {for(j=0;j<i+1;j++)
         {
          Total_rho[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        rho_CI_bool=false;
       }
       else if(spinCIrho)
       {
        counter=0;
        AUX=new double[spin_CI_rho];
        cast_doub(AUX,line,counter,spin_CI_rho);
        while(counter<spin_CI_rho)
        {
         getline(input_fchk,line);
         cast_doub(AUX,line,counter,spin_CI_rho);
        }
        counter=0;
        for(i=0;i<nbasisf;i++)
        {for(j=0;j<i+1;j++)
         {
          Spin_rho[i][j]=AUX[counter];
          counter++;
         }
        }
        delete[] AUX;
        spinCIrho=false;
       }
       else{}
      }
     }
     input_fchk.close();
   }
  else{cout<<"Unable to open file "<<name_file<<endl;}
  if(nbasisf!= nindepbasisf)
  {
   cout<<"Warning the number of independent function \n";
   cout<<"is different from the total number of functions \n";
  }
  //Build S^1/2 P S^1/2 diagonalization
  if(log_file)
  {
   string aux;
   double **Eigenvec,**Temp,**Temp2;
   bool check=true;
   PS_bool=true;
   P=new double*[nbasisf];
   if(extra1)
   {Ocupation=new double[nbasisf+nbasisf];}
   else
   {Ocupation=new double[nbasisf];}
   if(extra1)
   {
    Pbeta=new double*[nbasisf];
    Sbeta=new double*[nbasisf];
   }
   Sao=new double*[nbasisf];
   S=new double*[nbasisf];
   Temp=new double*[nbasisf];
   Temp2=new double*[nbasisf];
   Eigenvec=new double*[nbasisf];
   for(i=0;i<nbasisf;i++)
   {
    P[i]=new double[nbasisf];
    S[i]=new double[nbasisf];
    Sao[i]=new double[nbasisf];
    Eigenvec[i]=new double[nbasisf];
    Temp[i]=new double[nbasisf];
    Temp2[i]=new double[nbasisf];
    if(extra1)
    {
     Pbeta[i]=new double[nbasisf];
     Sbeta[i]=new double[nbasisf];
    }
   }
   for(i=0;i<nbasisf;i++)
   {
    for(j=0;j<nbasisf;j++)
    {
     P[i][j]=ZERO;
     if(extra1)
     {
      Pbeta[i][j]=ZERO;
      Sbeta[i][j]=ZERO;
     }
     S[i][j]=ZERO;
     Eigenvec[i][j]=ZERO;
     Temp[i][j]=ZERO;
     Temp2[i][j]=ZERO;
    }
    Ocupation[i]=ZERO;
    if(extra1)
    {Ocupation[i+nbasisf]=ZERO;}
   }
   ifstream log;
   log.open((name_log).c_str());
   counter=0;counter2=0;counter3=0;
   if(log.good())
   {
    while(getline(log,line))
    {
     if(line.length()>12)
     {
      if(line.substr(5,7)=="Overlap")
      {
       getline(log,line);
       counter=nbasisf;
       overlap=true;
       counter2=0;
       counter3=0;
       do
       {
        getline(log,line);
        aux=line.substr(6,5);
        nan(aux,check);
        for(i=0;i<(int)line.length();i++)
        {
         if(line[i]=='D')
         {
          line[i]='E';
         }
        }
        if(check)
        {
         line=line.substr(7,line.length());
         if((line.length()==14||line.length()==28)||
         (line.length()==42||line.length()==56)|| (line.length()==70))
         {
          if(line.length()==14)
          {
           stringstream ss(line);
           ss>>S[counter2][counter3];
           if(S[counter2][counter3]==ZERO){S[counter2][counter3]=9999;}
           counter2++;
          }
          else if(line.length()==28)
          {
           stringstream ss(line.substr(0,14));
           ss>>S[counter2][counter3];
           if(S[counter2][counter3]==ZERO){S[counter2][counter3]=9999;}
           stringstream ss1(line.substr(14,14));
           ss1>>S[counter2][counter3+1];
           if(S[counter2][counter3+1]==ZERO){S[counter2][counter3+1]=9999;}
           counter2++;
          }
          else if(line.length()==42)
          {
           stringstream ss(line.substr(0,14));
           ss>>S[counter2][counter3];
           if(S[counter2][counter3]==ZERO){S[counter2][counter3]=9999;}
           stringstream ss1(line.substr(14,14));
           ss1>>S[counter2][counter3+1];
           if(S[counter2][counter3+1]==ZERO){S[counter2][counter3+1]=9999;}
           stringstream ss2(line.substr(28,14));
           ss2>>S[counter2][counter3+2];
           if(S[counter2][counter3+2]==ZERO){S[counter2][counter3+2]=9999;}
           counter2++;
          }
          else if(line.length()==56)
          {
           stringstream ss(line.substr(0,14));
           ss>>S[counter2][counter3];
           if(S[counter2][counter3]==ZERO){S[counter2][counter3]=9999;}
           stringstream ss1(line.substr(14,14));
           ss1>>S[counter2][counter3+1];
           if(S[counter2][counter3+1]==ZERO){S[counter2][counter3+1]=9999;}
           stringstream ss2(line.substr(28,14));
           ss2>>S[counter2][counter3+2];
           if(S[counter2][counter3+2]==ZERO){S[counter2][counter3+2]=9999;}
           stringstream ss3(line.substr(42,14));
           ss3>>S[counter2][counter3+3];
           if(S[counter2][counter3+3]==ZERO){S[counter2][counter3+3]=9999;}
           counter2++;
          }
          else
          {
           stringstream ss(line.substr(0,14));
           ss>>S[counter2][counter3];
           if(S[counter2][counter3]==ZERO){S[counter2][counter3]=9999;}
           stringstream ss1(line.substr(14,14));
           ss1>>S[counter2][counter3+1];
           if(S[counter2][counter3+1]==ZERO){S[counter2][counter3+1]=9999;}
           stringstream ss2(line.substr(28,14));
           ss2>>S[counter2][counter3+2];
           if(S[counter2][counter3+2]==ZERO){S[counter2][counter3+2]=9999;}
           stringstream ss3(line.substr(42,14));
           ss3>>S[counter2][counter3+3];
           if(S[counter2][counter3+3]==ZERO){S[counter2][counter3+3]=9999;}
           stringstream ss4(line.substr(56,14));
           ss4>>S[counter2][counter3+4];
           if(S[counter2][counter3+4]==ZERO){S[counter2][counter3+4]=9999;}
           counter2++;
           if(counter==counter2)
           {
            counter=counter-5;
            counter2=0;
            counter3=counter3+5;
           }
          }
         }
         else
         {}
        }
       }while(check);
       check=true;
      }
     }
    }
    if(nbasisf>5)
    {
     for(counter=5;counter<nbasisf;counter++)
     {
      for(i=0;i<nbasisf;i++)
      {
       for(j=i+1;j<nbasisf;j++)
       {
        if(S[nbasisf-1-j][counter]!=ZERO)
        {
          S[nbasisf-1-i][counter]=S[nbasisf-1-j][counter];
          S[nbasisf-1-j][counter]=ZERO;
          j=nbasisf;
        }
        else
        {}
       }
      }
     }
    }
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      if(S[i][j]==9999){S[i][j]=ZERO;}
     }
    }
    for(i=0;i<nbasisf;i++)
    {
     for(j=i;j<nbasisf;j++)
     {
      S[i][j]=S[j][i];
      Sao[i][j]=S[j][i];
      Sao[j][i]=S[j][i];
     }
    }
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<i+1;j++)
     {
      P[i][j]=Total_rho[i][j];
      if(extra1)
      {
       Pbeta[i][j]=(Total_rho[i][j]-Spin_rho[i][j])/TWO;
       P[i][j]=(Total_rho[i][j]+Spin_rho[i][j])/TWO;
       Pbeta[j][i]=Pbeta[i][j];
      }
      P[j][i]=P[i][j];
     }
    }
    jacobi(nbasisf,S,Eigenvec);
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      if(i==j)
      {
       S[i][j]=sqrt(S[i][j]);
       Temp2[i][j]=ONE/(sqrt(S[i][j])+pow(TEN,-TEN));
      }
      else
      {
       S[i][j]=ZERO;
       Temp2[i][j]=ZERO;
      }
     }
    }
  //  mat_inverse2(nbasisf,S,Temp2);
    matmul(nbasisf,Eigenvec,S,Temp);//Save in Temp the result of the product
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      S[i][j]=Temp[i][j];
     }
    }
    matmul(nbasisf,Eigenvec,Temp2,Temp);//Save in Temp the result of the product
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      Temp2[i][j]=Temp[i][j]; //Save in Temp2 = U s^-1/2
     }
    }
    //Save in Temp the inverse of the Eigenvectors Matrix. U^-1 =U^Tranpose
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      Temp[i][j]=Eigenvec[j][i];
     }
    }
    //Save in Eigenvectors U^-1
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      Eigenvec[i][j]=Temp[i][j];
     }
    }
    matmul(nbasisf,S,Eigenvec,Temp); //So I have S^1/2 in S (Save in Temp the result of the product)
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      S[i][j]=Temp[i][j];
     }
    }
    matmul(nbasisf,Temp2,Eigenvec,Temp); //So I have S^-1/2 in Temp2 (Save in Temp the result of the product)
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      Temp2[i][j]=Temp[i][j];
     }
    }
    matmul(nbasisf,P,S,Temp); //Save in Temp the result of the product
    for(i=0;i<nbasisf;i++) //Save in P=PS^1/2
    {
     for(j=0;j<nbasisf;j++)
     {
      P[i][j]=Temp[i][j];
     }
    }
    matmul(nbasisf,S,P,Temp); //Save in Temp the result of the product
    for(i=0;i<nbasisf;i++) //Save in P=S^1/2PS^1/2
    {
     for(j=0;j<nbasisf;j++)
     {
      P[i][j]=Temp[i][j]; //Save P' = S^1/2 P S^1/2
     }
    }
    if(extra1)
    {
     matmul(nbasisf,Pbeta,S,Temp); //Save in Temp the result of the product
     for(i=0;i<nbasisf;i++) //Save in Pbeta=PbetaS^1/2
     {
      for(j=0;j<nbasisf;j++)
      {
       Pbeta[i][j]=Temp[i][j];
      }
     }
     matmul(nbasisf,S,Pbeta,Temp); //Save in Temp the result of the product
     for(i=0;i<nbasisf;i++) //Save in Pbeta=S^1/2PbetaS^1/2
     {
      for(j=0;j<nbasisf;j++)
      {
       Pbeta[i][j]=Temp[i][j]; //Save Pbeta' = S^1/2 Pbeta S^1/2
      }
     }
    }
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      S[i][j]=Temp2[i][j]; //Save S^-1/2
      if(extra1){Sbeta[i][j]=Temp2[i][j];}
     }
    }
    jacobi(nbasisf,P,Eigenvec); //Diagonalize P=> p=U^-1 S^1/2PS^1/2 U
    matmul(nbasisf,S,Eigenvec,Temp); //S^-1/2 Eigenv=Coefs AOs (columns) => NO
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nbasisf;j++)
     {
      S[i][j]=Temp[i][j]; //Saved in S the coefs of NO alpha
     }
    }
    if(extra1)
    {
     jacobi(nbasisf,Pbeta,Eigenvec); //Diagonalize P=> p=U^-1 S^1/2PS^1/2 U
     matmul(nbasisf,Sbeta,Eigenvec,Temp); //S^-1/2 Eigenv=Coefs AOs (columns) => NObeta
     for(i=0;i<nbasisf;i++)
     {
      for(j=0;j<nbasisf;j++)
      {
       Sbeta[i][j]=Temp[i][j]; //Saved in Sbeta the coefs of NO beta
      }
     }
    }
    //pseudo wfn and wfx files are printed by default
    pseudowfn.open(("pseudo_"+name_file.substr(0,(name_file.length()-5))+".wfn").c_str());
    pseudowfx.open(("pseudo_"+name_file.substr(0,(name_file.length()-5))+".wfx").c_str());
    pseudowfx<<"<Number of Occupied Molecular Orbitals>"<<endl;
    if(extra1)
    {
     counter=0;
     counter2=0;
     pseudowfx<<2*nbasisf<<endl;
     pseudowfx<<"</Number of Occupied Molecular Orbitals>"<<endl;
     pseudowfx<<"<Molecular Orbital Occupation Numbers>"<<endl;
     pseudowfx<<setprecision(10)<<fixed<<scientific;
     for(i=0;i<nbasisf+nbasisf;i++)
     {
      if(i%2==0)
      {Ocupation[i]=P[counter][counter];counter++;}
      else
      {Ocupation[i]=Pbeta[counter2][counter2];counter2++;}
      if(Ocupation[i]<ZERO)
      {
       pseudowfn<<setprecision(0)<<fixed;
       if(i+1 <10)
       {
        pseudowfn<<"MO     "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       else if(i+1<100)
       {
        pseudowfn<<"MO    "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       else if(i+1<1000)
       {
        pseudowfn<<"MO   "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       else
       {
        pseudowfn<<"MO  "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       pseudowfn<<setprecision(7)<<fixed;
       pseudowfn<<Ocupation[i]<<"  ORB. ENERGY =    0.000000"<<endl;
      }
      else
      {
       pseudowfn<<setprecision(0)<<fixed;
       if(i+1 <10)
       {
        pseudowfn<<"MO     "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       else if(i+1<100)
       {
        pseudowfn<<"MO    "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       else if(i+1<1000)
       {
        pseudowfn<<"MO   "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       else
       {
        pseudowfn<<"MO  "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       pseudowfn<<setprecision(7)<<fixed;
       pseudowfn<<Ocupation[i]<<"  ORB. ENERGY =    0.000000"<<endl;
      }
      pseudowfx<<setw(17)<<Ocupation[i]<<endl;
     }
     counter=0;
     counter2=0;
    }
    else
    {
     pseudowfx<<nbasisf<<endl;
     pseudowfx<<"</Number of Occupied Molecular Orbitals>"<<endl;
     pseudowfx<<"<Molecular Orbital Occupation Numbers>"<<endl;
     pseudowfx<<setprecision(10)<<fixed<<scientific;
     for(i=0;i<nbasisf;i++)
     {
      Ocupation[i]=P[i][i];
      if(Ocupation[i]<ZERO)
      {
       pseudowfn<<setprecision(0)<<fixed;
       if(i+1 <10)
       {
        pseudowfn<<"MO     "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       else if(i+1<100)
       {
        pseudowfn<<"MO    "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       else if(i+1<1000)
       {
        pseudowfn<<"MO   "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       else
       {
        pseudowfn<<"MO  "<<i+1<<"    MO 0.0        OCC NO =   ";
       }
       pseudowfn<<setprecision(7)<<fixed;
       pseudowfn<<Ocupation[i]<<"  ORB. ENERGY =    0.000000"<<endl;
      }
      else
      {
       pseudowfn<<setprecision(0)<<fixed;
       if(i+1 <10)
       {
        pseudowfn<<"MO     "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       else if(i+1<100)
       {
        pseudowfn<<"MO    "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       else if(i+1<1000)
       {
        pseudowfn<<"MO   "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       else
       {
        pseudowfn<<"MO  "<<i+1<<"    MO 0.0        OCC NO =    ";
       }
       pseudowfn<<setprecision(7)<<fixed;
       pseudowfn<<Ocupation[i]<<"  ORB. ENERGY =    0.000000"<<endl;
      }
      pseudowfx<<setw(17)<<Ocupation[i]<<endl;
     }
    }
    pseudowfn.close();
    pseudowfx<<"</Molecular Orbital Occupation Numbers>"<<endl;
    pseudowfx.close();
    log.close();
    if(!overlap)
    {
     cout<<"(Warning, the AO overlap matrix was not found! Include the keyword iop(3/33=4) and"<<endl;
     cout<<"repeat the gaussian calculation. Then, use the new .log file for RHO.OPS)"<<endl;
    }
   }
   else
   {cout<<"Unable to open the .log file"<<endl;}
   for(i=0;i<nbasisf;i++)
   {
    delete[] Eigenvec[i];
    delete[] Temp[i];
    delete[] Temp2[i];
   }
   delete[] Eigenvec;
   delete[] Temp;
   delete[] Temp2;
   counter=0;
   counter2=0;
   counter3=0;
  }
  else
  {}
 }
 else
 {
  string aux;
  ifstream logical;
  bool gamess=false;
  double sum;
  //The information available in wfx file is almost the same as wfn but with extra info!
  if(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X')
  {wfx=true;}
  /////////////////////////////////////////////////////////////////////////////////
  //Read all information required
  input_wfn.open((name_file).c_str());
  if(input_wfn.is_open())
  { //Read information and castings
   while(getline(input_wfn,line))
   {
     if(!wfx)
     {
      //save out of range of last-1 line!!
      if(line.length()<100)
      {line=line+"                                                                                                             ";}
      //look for special phrases and activate keywords
      if(line.substr(44,10)=="PRIMITIVES")
      {
       aux=line.substr(16,24);
       stringstream ss(aux);
       //Mo_coef number of NOs
       ss>>MO_coef;
       nbasisf=MO_coef; //nbasisf=MO_coef number of NO
       aux=line.substr(37,44);
       stringstream ss1(aux);
       //nprimitiv number of primitives
       ss1>>nprimitv;
       aux=line.substr(55,64);
       stringstream ss2(aux);
       //number of atoms
       ss2>>natoms;
       Cartesian_Coor=new double*[natoms];
       for(i=0;i<natoms;i++)
       {Cartesian_Coor[i]=new double[3];}
       Nu_charge=new double[natoms];
       Prim_exp=new double[nprimitv];
       //shell_map is primitive assignment in this context!!!
       //shell_type is the type of the primitive
       shell_map=new int[nprimitv];
       shell_type=new int[nprimitv];
       //MOcoefA are the coefficients which go with the primitive
       //and form the MO. MO= sum coef*primitive
       MOcoefA=new double*[MO_coef];
       for(i=0;i<MO_coef;i++)
       {MOcoefA[i]=new double[nprimitv];}
       for(i=0;i<MO_coef;i++)
       {
        for(j=0;j<nprimitv;j++)
        {
         MOcoefA[i][j]=ZERO;
        }
       }
       //Ocupation of the MO
       Ocupation=new double[MO_coef];
      }
      else if(line.substr(13,6)=="CENTRE")
      {
       for(i=0;i<natoms-1;i++)
       {
        cast_wfn_coord(line.substr(25,34),Cartesian_Coor,counter);
        cast_wfn_nuch(line.substr(71,5),Nu_charge,counter2);
        getline(input_wfn,line);
       }
       cast_wfn_coord(line.substr(25,34),Cartesian_Coor,counter);
       cast_wfn_nuch(line.substr(71,5),Nu_charge,counter2);
      }
      else if(line.substr(0,6)=="CENTRE")
      { counter=0;
        cast_wfn_type_asig(line.substr(20,line.length()),shell_map,counter);
        while(counter<nprimitv)
        {
         getline(input_wfn,line);
         cast_wfn_type_asig(line.substr(20,line.length()),shell_map,counter);
        }
      }
      else if(line.substr(0,4)=="TYPE")
      {counter=0;
       cast_wfn_type_asig(line.substr(20,line.length()),shell_type,counter);
       while(counter<nprimitv)
       {
        getline(input_wfn,line);
        cast_wfn_type_asig(line.substr(20,line.length()),shell_type,counter);
       }
      }
      else if(line.substr(0,9)=="EXPONENTS")
      {counter=0;
       cast_wfn_expon(line.substr(10,line.length()),Prim_exp,counter);
       while(counter<nprimitv)
       {
        getline(input_wfn,line);
        cast_wfn_expon(line.substr(10,line.length()),Prim_exp,counter);
       }
       counter=0;
      }
      else if(line.substr(0,2)=="MO" && line.substr(0,6)!="MOLCAS")
      {
       counter2=0;
       aux=line.substr(35,12);
       stringstream ss(aux);
       ss>>Ocupation[counter];
       do
       {
        getline(input_wfn,line);
        cast_wfn_MO(line,MOcoefA,counter,counter2);
       }while(counter2<nprimitv);
       counter++;
      }
      else
      {}
    }
    else
    {
      if(line=="<Number of Nuclei>")
      {
       input_wfn>>natoms;
       Cartesian_Coor=new double*[natoms];
       for(i=0;i<natoms;i++)
       {Cartesian_Coor[i]=new double[3];}
       Nu_charge=new double[natoms];
      }
      if(line=="<Number of Occupied Molecular Orbitals>")
      {
       input_wfn>>MO_coef;
       nbasisf=MO_coef;
       Ocupation=new double[MO_coef];
       SPIN=new bool[MO_coef];
      }
      if(line=="<Electronic Spin Multiplicity>")
      {
       input_wfn>>multiplicity;
      }
      if(line=="<Nuclear Charges>")
      {
       for(i=0;i<natoms;i++)
       {
        input_wfn>>Nu_charge[i];
       }
      }
      if(line=="<Nuclear Cartesian Coordinates>")
      {
       for(i=0;i<natoms;i++)
       {
        input_wfn>>Cartesian_Coor[i][0]>>Cartesian_Coor[i][1]>>Cartesian_Coor[i][2];
       }
      }
      if(line=="<Number of Primitives>")
      {
       input_wfn>>nprimitv;
       Prim_exp=new double[nprimitv];
       //shell_map is primitive assignment in this context!!!
       //shell_type is the type of the primitive
       shell_map=new int[nprimitv];
       shell_type=new int[nprimitv];
       MOcoefA=new double*[MO_coef];
       for(i=0;i<MO_coef;i++)
       {MOcoefA[i]=new double[nprimitv];}
      }
      if(line=="<Primitive Centers>")
      {
       for(i=0;i<nprimitv;i++)
       {
        input_wfn>>shell_map[i];
       }
      }
      if(line=="<Primitive Types>")
      {
       for(i=0;i<nprimitv;i++)
       {
        input_wfn>>shell_type[i];
       }
      }
      if(line=="<Primitive Exponents>")
      {
       for(i=0;i<nprimitv;i++)
       {
        input_wfn>>Prim_exp[i];
       }
      }
      if(line=="<Molecular Orbital Occupation Numbers>")
      {
       for(i=0;i<MO_coef;i++)
       {
        input_wfn>>Ocupation[i];
       }
      }
      if(line=="<Molecular Orbital Spin Types>")
      {
       for(i=0;i<MO_coef;i++)
       {
        SPIN[i]=true;
        getline(input_wfn,line);
        if(line.substr(0,5)!=" Alph")
        {
         SPIN[i]=false;
        }
       }
      }
      if(line=="<Molecular Orbital Primitive Coefficients>")
      {
       i=0;
       do
       {
        do
        {
         getline(input_wfn,line);
        }while(line!="</MO Number>");
        line="clean";
        for(j=0;j<nprimitv;j++)
        {
         input_wfn>>MOcoefA[i][j];
        }
        i++;
       }while(i<MO_coef);
      }
    }
   }
   input_wfn.close();
     //Take the multiplicity as parameter from Gaussian or Gamess .log
     if(!wfx)
     {
      multiplicity=0;
      multiplicity=mult_in;
     }
     if(multiplicity==0)
     {
      logical.open((name_log).c_str());
      if(logical.good())
      {
       while(getline(logical,line))
       {
        if(line.substr(0,26)=="          *         GAMESS" || gamess)
        {
         gamess=true;
         if(line.substr(0,18)==" SPIN MULTIPLICITY")
         {
          aux=line.substr(50,2);
          stringstream ss(aux);
          ss>>multiplicity;
         }
         else if(line.substr(0,8)==" MULT  =")
         {
          aux=line.substr(14,2);
          stringstream ss(aux);
          ss>>multiplicity;
         }
        }
       else if(!gamess)
        {
         if(line.substr(0,9)==" Charge =")
         {
          aux=line.substr(27,29);
          stringstream ss(aux);
          ss>>multiplicity;
         }
        }
        else
        {}
       }
      }
      logical.close();
     }
  ///////////////////////////////
  //Define type of calculation///
  ///////////////////////////////
     sum=ZERO;
     for(i=0;i<MO_coef;i++)
     {
      sum=sum+Ocupation[i];
      if(Ocupation[i]==ZERO){virtuals=true;}
     }
     nelectrons=(int)round(sum);
     for(i=0;i<MO_coef;i++)
     {
      k=abs((int) Ocupation[i]);
      if(abs(Ocupation[i]-k)==ZERO)
      {
       identity="1-det.";
      }
      else
      {
       identity="N-det.";
       correlated=true;
       break;
      }
     }
     if(cas && identity=="1-det.")
     {
      correlated=true;
      identity="N-det.";
      virtuals=false;
     }
     if(!correlated)
     {
      for(i=0;i<MO_coef;i++)
      {
       k=abs((int) Ocupation[i]);
       if(k==2 || k==0)
       {
        open_shell=false;
        if(k==0){virtuals=true;}
       }
       else
       {
        open_shell=true;
        break;
       }
      }
      if(open_shell)
      {
       uhf=true;
       if(nelectrons%2==0 && multiplicity==1)
       {
        identity=identity+"closed-shell.unrestricted.";
        //but treated as open-shell MOs for calculating
       }
       else{identity=identity+"open-shell.";}
      }
      else{identity=identity+"closed-shell.";rhf=true;uhf=false;}
     }
     else
     {
      if((nelectrons%2==0)&&(multiplicity==1))
      {
       open_shell=false;
       identity=identity+"closed-shell.";
       for(i=0;i<MO_coef;i++)
       {
        if(Ocupation[i]>=ZERO && Ocupation[i]<=TWO)
        {
         relaxed=false;
        }
        else
        {
         relaxed=true;
         identity=identity+"relaxed.CI";
         cout<<"Warning density is "<<identity<<"!"<<endl;
         CI=true;
         break;
        }
       }
      }
      else
      {
       open_shell=true;
       identity=identity+"open-shell.";
       for(i=0;i<MO_coef;i++)
       {
        if(Ocupation[i]>=ZERO && Ocupation[i]<=ONE)
        {
         relaxed=false;
        }
        else
        {
         relaxed=true;
         identity=identity+"relaxed.CI";
         cout<<"Warning density is "<<identity<<" !"<<endl;
         CI=true;
         break;
        }
       }
      }
      if(!relaxed)
      {
       if(cas)
       {
        relaxed=true; //For CAS relaxed=unrelaxed
        identity=identity+"relaxed.CAS";
        CAS=true;
       }
       else
       {
        relaxed=false;
        identity=identity+"unrelaxed.CI";
        CI=true;
       }
      }
     }
    if(!wfx)
    {
     SPIN=new bool[nbasisf];
     if(!open_shell)
     {for(i=0;i<nbasisf;i++){SPIN[i]=true;}}
     else
     {
      if(!correlated)
      {
       for(i=0;i<nbasisf;i++)
       {
        SPIN[i]=false;
        if(i%2==0){SPIN[i]=true;}
       }
       if(virtuals)
       {
        counter=0;
        counter2=0;
        for(i=0;i<nbasisf;i++)
        {
         counter++;
         if((int)Ocupation[i]!=0)
         {
          if(counter2!=0)
          {
           break; //counter is 1 > number of alpha NO
          }
          else{}
         }
         else
         {counter2=i;}
        }
        if(counter!=nbasisf)
        {
         counter=counter-1;
         double **TEMP;
         TEMP=new double*[nbasisf];
         AUX=new double[nbasisf];
         for(i=0;i<nbasisf;i++)
         {TEMP[i]=new double[nprimitv];}
         counter2=0;
         for(i=0;i<nbasisf;i++)
         {
          if(i%2==0)
          {
           for(j=0;j<nprimitv;j++)
           {
            TEMP[i][j]=MOcoefA[counter2][j];
            AUX[i]=Ocupation[counter2];
           }
           counter2++;
          }
          else
          {
           for(j=0;j<nprimitv;j++)
           {
            TEMP[i][j]=MOcoefA[counter][j];
            AUX[i]=Ocupation[counter];
           }
           counter++;
          }
         }
         for(i=0;i<nbasisf;i++)
         {
          for(j=0;j<nprimitv;j++)
          {
           MOcoefA[i][j]=TEMP[i][j];
          }
          Ocupation[i]=AUX[i];
         }
         for(i=0;i<nbasisf;i++)
         {delete[] TEMP[i];}
         delete[] TEMP;
         delete[] AUX;
         counter=0;
         counter2=0;
        }
        else
        {
         //If we only have alpha NOs in the wfn file, correct spin component
         //and do not asume that the first orbital and second orbitals have diff spins.
         for(i=0;i<nbasisf;i++)
         {
          SPIN[i]=true;
          no_beta_wfn=true;
         }
        }
       }
       else
       {
        int unpaired=multiplicity-1;
        if(unpaired==0)
        {
         double **TEMP;
         TEMP=new double*[nbasisf];
         for(i=0;i<nbasisf;i++)
         {TEMP[i]=new double[nprimitv];}
         AUX=new double[nbasisf];
         counter=0;
         counter2=nbasisf/2;
         for(i=0;i<nbasisf;i++)
         {
          for(j=0;j<nprimitv;j++)
          {
           if(i%2==0)
           {
            TEMP[i][j]=MOcoefA[counter][j];
            AUX[i]=Ocupation[counter];
           }
           else
           {
            TEMP[i][j]=MOcoefA[counter2][j];
            AUX[i]=Ocupation[counter2];
           }
          }
          if(i%2==0)
          {counter++;}
          else
          {counter2++;}
         }
         for(i=0;i<nbasisf;i++)
         {
          for(j=0;j<nprimitv;j++)
          {
           MOcoefA[i][j]=TEMP[i][j];
          }
          Ocupation[i]=AUX[i];
         }
         for(i=0;i<nbasisf;i++)
         {delete[] TEMP[i];}
         delete[] TEMP;
         delete[] AUX;
         counter=0;
         counter2=0;
        }
        else
        {
         int paired,alpha,beta,holes;
         paired=(nelectrons-unpaired)/2;
         alpha=paired+unpaired;
         beta=nelectrons-alpha;
         holes=alpha-1;
         holes=holes-beta;
         if(holes!=0)
         {cout<<"Warning trivial orbitals have been built!"<<endl;}
         double **TEMP;
         TEMP=new double*[nbasisf+holes];
         AUX=new double[nbasisf+holes];
         for(i=0;i<nbasisf+holes;i++)
         {TEMP[i]=new double[nprimitv];}
         for(i=0;i<nbasisf+holes;i++)
         {
          for(j=0;j<nprimitv;j++)
          {
           if(i<nbasisf)
           {
            TEMP[i][j]=MOcoefA[i][j];
            AUX[i]=Ocupation[i];
           }
           else
           {
            TEMP[i][j]=ZERO; //Trivial virtual orbitals
            AUX[i]=ZERO;
           }
          }
         }
         for(i=0;i<nbasisf;i++)
         {delete[] MOcoefA[i];}
         delete[] MOcoefA;
         delete[] Ocupation;
         delete[] SPIN;
         nbasisf=nbasisf+holes;
         MOcoefA=new double*[nbasisf];
         Ocupation=new double[nbasisf];
         SPIN=new bool[nbasisf];
         for(i=0;i<nbasisf;i++)
         {
          SPIN[i]=false;
          if(i%2==0){SPIN[i]=true;}
         }
         for(i=0;i<nbasisf;i++)
         {MOcoefA[i]=new double[nprimitv];}
         for(i=0;i<alpha;i++)
         {
          for(j=0;j<nprimitv;j++)
          {
           MOcoefA[2*i][j]=TEMP[i][j];
          }
          Ocupation[2*i]=AUX[i];
          beta=i+1;
         }
         for(i=0;i<alpha-1;i++)
         {
          for(j=0;j<nprimitv;j++)
          {
           MOcoefA[2*i+1][j]=TEMP[beta][j];
          }
          Ocupation[2*i+1]=AUX[beta];
          beta++;
         }
         for(i=0;i<nbasisf;i++)
         {delete[] TEMP[i];}
         delete[] TEMP;
         delete[] AUX;
        }
       }
      }
      else
      {for(i=0;i<nbasisf;i++){SPIN[i]=true;}}
     }
    }
    else
    {
     //We turn off correlated variable problem since wfx file contains occupancies divide by spin (the 1-RDM was diagonalized with spins separated)
     //this will not affect any density dependent quantity but only any property that was not accesible from correlated open-shell wfn files (properties that require correlated open shell spins)
     //are now available with wfx files! =)
     correlated=false;
     bool orderNOs=false;
     counter=0;
     for(i=0;i<MO_coef;i++)
     {
      if(!SPIN[i]){orderNOs=true;break;}
      else{counter++;}
     }
     if(orderNOs)
     {
      counter2=MO_coef-counter;
      double **TEMP;
      if(counter2==counter)
      {
       //Same number of alpha and beta NOs (probably some unrestricted formalism)
       TEMP=new double*[MO_coef];
       AUX=new double[MO_coef];
       for(i=0;i<MO_coef;i++)
       {
        TEMP[i]=new double[nprimitv];
       }
       //We only need to order the NOs and Occ.: Alpha, Beta, Alpha, Beta...
       for(i=0;i<counter;i++)
       {
        for(j=0;j<nprimitv;j++)
        {
         TEMP[2*i][j]=MOcoefA[i][j];
         TEMP[2*i+1][j]=MOcoefA[i+counter][j];
        }
        AUX[2*i]=Ocupation[i];
        AUX[2*i+1]=Ocupation[i+counter];
       }
       for(i=0;i<MO_coef;i++)
       {
        for(j=0;j<nprimitv;j++)
        {
         MOcoefA[i][j]=TEMP[i][j];
        }
        Ocupation[i]=AUX[i];
       }
       for(i=0;i<MO_coef;i++)
       {
        delete[] TEMP[i];TEMP[i]=NULL;
       }
       delete[] TEMP; TEMP=NULL;
       delete[] AUX; AUX=NULL;
      }
      else
      {
       cout<<"Warning "<<counter-counter2<<" trivial orbital(s) was/were built"<<endl;
       TEMP=new double*[MO_coef];
       AUX=new double[MO_coef];
       for(i=0;i<MO_coef;i++)
       {
        TEMP[i]=new double[nprimitv];
       }
       for(i=0;i<MO_coef;i++)
       {
        for(j=0;j<nprimitv;j++)
        {
         TEMP[i][j]=MOcoefA[i][j];
        }
        AUX[i]=Ocupation[i];
       }
       for(i=0;i<MO_coef;i++)
       {
        delete[] MOcoefA[i];MOcoefA[i]=NULL;
       }
       delete[] MOcoefA;MOcoefA=NULL;
       delete[] Ocupation,Ocupation=NULL;
       delete[] SPIN;SPIN=NULL;
       counter3=MO_coef;
       if(counter2>counter)
       {
        cout<<"Normal excess of beta NOs found in wfx file"<<endl;
        MO_coef=2*counter2;
        nbasisf=MO_coef;
        SPIN=new bool[MO_coef];
        Ocupation=new double[MO_coef];
        MOcoefA=new double*[MO_coef];
        for(i=0;i<MO_coef;i++)
        {
         MOcoefA[i]=new double[nprimitv];
         for(j=0;j<nprimitv;j++)
         {
          MOcoefA[i][j]=ZERO;
          Ocupation[i]=ZERO;
          SPIN[i]=true;
         }
        }
        for(i=0;i<counter;i++)
        {
         for(j=0;j<nprimitv;j++)
         {
          MOcoefA[2*i][j]=TEMP[i][j];
         }
         Ocupation[2*i]=AUX[i];
        }
        k=0;
        for(i=counter;i<counter3;i++)
        {
         for(j=0;j<nprimitv;j++)
         {
          MOcoefA[2*k+1][j]=TEMP[i][j];
         }
         Ocupation[2*k+1]=AUX[i];
         SPIN[2*k+1]=false;
         k++;
        }
       }
       else
       {
        cout<<"Normal excess of alpha NOs found in wfx file"<<endl;
        MO_coef=2*counter;
        nbasisf=MO_coef;
        SPIN=new bool[MO_coef];
        Ocupation=new double[MO_coef];
        MOcoefA=new double*[MO_coef];
        for(i=0;i<MO_coef;i++)
        {
         MOcoefA[i]=new double[nprimitv];
         for(j=0;j<nprimitv;j++)
         {
          MOcoefA[i][j]=ZERO;
          Ocupation[i]=ZERO;
          SPIN[i]=false;
         }
        }
        for(i=0;i<counter;i++)
        {
         for(j=0;j<nprimitv;j++)
         {
          MOcoefA[2*i][j]=TEMP[i][j];
         }
         Ocupation[2*i]=AUX[i];
         SPIN[2*i]=true;
        }
        k=0;
        for(i=counter;i<counter3;i++)
        {
         for(j=0;j<nprimitv;j++)
         {
          MOcoefA[2*k+1][j]=TEMP[i][j];
         }
         Ocupation[2*k+1]=AUX[i];
         k++;
        }
       }
       for(i=0;i<counter3;i++)
       {
        delete[] TEMP[i];TEMP[i]=NULL;
       }
       delete[] TEMP; TEMP=NULL;
       delete[] AUX; AUX=NULL;
      }
     }
     counter=0;
     counter2=0;
     counter3=0;
    }
    //Once we have define the type of wfn of wfx file, we proceed
  }
  else  cout<<"Unable to open file xxx.wfn or xxx.wfx"<<endl;
 }
 if(CM_in){Center_of_mass();}
}
//Copy density objects
//WARNING! WARNING! WARNING!
//Beta MOs info is not copied neither alpha MO coefs!
//They are define only TO MAKE THE DESTRUCTOR NOT GIVE SEGMENTATION FAULT!
READ_FCHK_WFN::READ_FCHK_WFN(const READ_FCHK_WFN&RHO)
{
 int i,j;
 if(!RHO.wfn)
 {
  wfn=RHO.wfn;
  extra0=RHO.extra0;
  prim_exp=RHO.prim_exp;
  counter=1;
  nbasisf=RHO.nbasisf;
  nshells=RHO.smap;
  PS_bool=false;
  //extra1 which is refered to MO beta is not used by copy
  //be aware that if you copy READ_FCHK_WFN objects the beta MO
  //info is not copied! (Is defined only for the constructor)
  extra1=false;
  natoms=RHO.natoms;
  Total_rho=new double*[nbasisf];
  MOcoefA=new double*[nbasisf];
  //Alpha MO coef inf. is NOT COPIED!!
  //WARNING! WARNING! WARNING!
  for(i=0;i<nbasisf;i++)
  {
   Total_rho[i]=new double[counter];
   MOcoefA[i]=new double[1];
   counter++;
  }
  counter=0;
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<i+1;j++)
   {
    Total_rho[i][j]=RHO.Total_rho[i][j];
   }
  }
  Cartesian_Coor=new double*[natoms];
  for(i=0;i<natoms;i++)
  {
   Cartesian_Coor[i]=new double[3];
  }
  for(i=0;i<natoms;i++)
  {
   for(j=0;j<3;j++)
   {Cartesian_Coor[i][j]=RHO.Cartesian_Coor[i][j];}
  }
  shell_type=new int[nshells];
  n_prim_per_shell=new int[nshells];
  shell_map=new int[nshells];
  Nu_charge=new double[nbasisf];
  for(i=0;i<nshells;i++)
  {
   shell_type[i]=RHO.shell_type[i];
   n_prim_per_shell[i]=RHO.n_prim_per_shell[i];
   shell_map[i]=RHO.shell_map[i];
   Nu_charge[i]=RHO.Nu_charge[i];
  }
  Prim_exp=new double[prim_exp];
  Contr_Coef=new double[prim_exp];
  if(extra0){SP_Contr_Coef=new double[prim_exp];}
  for(i=0;i<prim_exp;i++)
  {
   Prim_exp[i]=RHO.Prim_exp[i];
   Contr_Coef[i]=RHO.Contr_Coef[i];
   if(extra0){SP_Contr_Coef[i]=RHO.SP_Contr_Coef[i];}
  }
 }
 else
 {
   wfn=RHO.wfn;
   im_wfn_wfx=RHO.im_wfn_wfx;
   nbasisf=RHO.nbasisf;
   MO_coef=nbasisf;
   natoms=RHO.natoms;
   nprimitv=RHO.nprimitv;
   Nu_charge=new double[natoms];
   Cartesian_Coor=new double*[natoms];
   for(i=0;i<natoms;i++)
   {Cartesian_Coor[i]=new double[3];}
   for(i=0;i<natoms;i++)
   {
    for(j=0;j<3;j++)
    {Cartesian_Coor[i][j]=RHO.Cartesian_Coor[i][j];}
    Nu_charge[i]=RHO.Nu_charge[i];
   }
   Prim_exp=new double[nprimitv];
   shell_map=new int[nprimitv];
   shell_type=new int[nprimitv];
   for(i=0;i<nprimitv;i++)
   {
    Prim_exp[i]=RHO.Prim_exp[i];
    shell_map[i]=RHO.shell_map[i];
    shell_type[i]=RHO.shell_type[i];
   }
   MOcoefA=new double*[MO_coef];
   for(i=0;i<MO_coef;i++)
   {MOcoefA[i]=new double[nprimitv];}
   Ocupation=new double[MO_coef];
   for(i=0;i<nbasisf;i++)
   {
    for(j=0;j<nprimitv;j++)
    {MOcoefA[i][j]=RHO.MOcoefA[i][j];}
    Ocupation[i]=RHO.Ocupation[i];
   }
   if(RHO.im_wfn_wfx)
   {
    MOcoefA_im=new double*[MO_coef];
    for(i=0;i<MO_coef;i++)
    {MOcoefA_im[i]=new double[nprimitv];}
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<nprimitv;j++)
     {MOcoefA_im[i][j]=RHO.MOcoefA_im[i][j];}
    }
   }
   SPIN=new bool[1]; //Check for another use
 }
}
//Cast from strings to numbers
//Function for reading the dimension of the array need.
void READ_FCHK_WFN::line_fill(int &number, string line_read)
{
  line_read=line_read.substr(55,line_read.length());
  stringstream ss(line_read); //string to int in ss
  ss>>number;
}
//Casting from getline to double in the array
void READ_FCHK_WFN::cast_doub(double *filling, string line_read,int &counter,int control)
{
 string aux;
 aux=line_read.substr(0,17);
 stringstream ss(aux);
 ss>>filling[counter];
 counter++;
 if(counter<control)
 {
  aux=line_read.substr(17,33);
  stringstream ss(aux);
  ss>>filling[counter];
  counter++;
  if(counter<control)
  {
   aux=line_read.substr(33,49);
   stringstream ss(aux);
   ss>>filling[counter];
   counter++;
   if(counter<control)
   {
    aux=line_read.substr(49,65);
    stringstream ss(aux);
    ss>>filling[counter];
    counter++;
    if(counter<control)
    {
     aux=line_read.substr(65,line_read.length());
     stringstream ss(aux);
     ss>>filling[counter];
     counter++;
    }
   }
  }
 }
}
//Casting from getline to int in the array
void READ_FCHK_WFN::cast_int(int *filling, string line_read,int &counter,int control)
{
 string aux;
 aux=line_read.substr(0,17);
 stringstream ss(aux);
 ss>>filling[counter];
 counter++;
 if(counter<control)
 {
  aux=line_read.substr(17,33);
  stringstream ss(aux);
  ss>>filling[counter];
  counter++;
  if(counter<control)
  {
   aux=line_read.substr(33,44);
   stringstream ss(aux);
   ss>>filling[counter];
   counter++;
   if(counter<control)
   {
    aux=line_read.substr(44,55);
    stringstream ss(aux);
    ss>>filling[counter];
    counter++;
    if(counter<control)
    {
     aux=line_read.substr(55,66);
     stringstream ss(aux);
     ss>>filling[counter];
     counter++;
     if(counter<control)
     {
     aux=line_read.substr(66,line_read.length());
     stringstream ss(aux);
     ss>>filling[counter];
     counter++;
     }
    }
   }
  }
 }
}
//Casting nuclear coord. for wfn
void READ_FCHK_WFN::cast_wfn_coord(string in, double **Cartesian_Coor,int &counter)
{
 string aux;
 aux=in.substr(0,11);
 stringstream ss(aux);
 ss>>Cartesian_Coor[counter][0];
 aux=in.substr(11,12);
 stringstream ss1(aux);
 ss1>>Cartesian_Coor[counter][1];
 aux=in.substr(23,15);
 stringstream ss2(aux);
 ss2>>Cartesian_Coor[counter][2];
 counter++;
}
//Casting nuclear charges for wfn
void READ_FCHK_WFN::cast_wfn_nuch(string in, double *Nu_charge,int &counter2)
{
 string aux;
 aux=in.substr(0,6);
 stringstream ss(aux);
 ss>>Nu_charge[counter2];
 counter2++;
}
//casting type and center assignment wfn
void READ_FCHK_WFN::cast_wfn_type_asig(string in,int *filling,int &counter)
{
 int i,j;
 i=0;
 j=1;
 string aux;
 do
 {
  aux=in.substr(i,3);
  stringstream ss(aux);
  ss>>filling[counter];
  i=i+3;
  j++;
  counter++;
 }while((counter<nprimitv) && (j<21));
}
//casting exponents wfn
void READ_FCHK_WFN::cast_wfn_expon(string in,double *filling,int &counter)
{
 int i,j;
 unsigned k;
 for(k=0;k<(in.length());k++){if(in[k]=='D'){in[k]='e';}}
 i=0;
 j=1;
 string aux;
 do
 {
  aux=in.substr(i,14);
  stringstream ss(aux);
  ss>>filling[counter];
  i=i+14;
  j++;
  counter++;
 }while((counter<nprimitv) && (j<6));
}
//casting coefs wfn
void READ_FCHK_WFN::cast_wfn_MO(string in,double **filling,int &counter,int &counter2)
{
 int i,j;
 unsigned k;
 for(k=0;k<(in.length());k++){if(in[k]=='D'){in[k]='e';}}
 i=0;
 j=1;
 string aux;
 do
 {
  aux=in.substr(i,16);
  stringstream ss(aux);
  ss>>filling[counter][counter2];
  i=i+16;
  j++;
  counter2++;
 }while((counter2<nprimitv) && (j<6));
}
//Send fchk MOs coef
void READ_FCHK_WFN::bring_mo_coefs(double **coefs)
{
 int i,j;
 //Build MO
 if(extra1)
 {
  for(i=0;i<2*nbasisf;i++)
  {
   if(i%2==0)
   {
    for(j=0;j<nbasisf;j++)
    {
     coefs[i][j]=MOcoefA[i/2][j];
    }
   }
   else
   {
    for(j=0;j<nbasisf;j++)
    {
     if(BETA_MOS)
     {
      coefs[i][j]=MOcoefB[(i-1)/2][j];
     }
     else
     {
      coefs[i][j]=MOcoefA[(i-1)/2][j];
     }
    }
   }
  }
 }
 else
 {
  for(i=0;i<2*nbasisf;i++)
  {
   if(i%2==0)
   {
    for(j=0;j<nbasisf;j++)
    {
     coefs[i][j]=MOcoefA[i/2][j];
    }
   }
   else
   {
    for(j=0;j<nbasisf;j++)
    {
     coefs[i][j]=MOcoefA[(i-1)/2][j];
    }
   }
  }
 }
}
//Delete object
READ_FCHK_WFN::~READ_FCHK_WFN()
{
 int i;
 if(!wfn)
 {
 //Free memory allocated
  delete[] Nu_charge;
  delete[] shell_type;
  delete[] n_prim_per_shell;
  delete[] shell_map;
  if(PS_bool)
  {
   for(i=0;i<nbasisf;i++)
   {
    delete[] P[i];
    delete[] S[i];
    if(extra1)
    {
     delete[] Pbeta[i];
     delete[] Sbeta[i];
    }
   }
   delete[] P;
   delete[] S;
   delete[] Sao;
   delete[] Ocupation;
   if(extra1)
   {
    delete[] Pbeta;
    delete[] Sbeta;
   }
  }
  for(i=0;i<natoms;i++)
  {delete[] Cartesian_Coor[i];}
  delete[] Cartesian_Coor;
  for(i=0;i<nbasisf;i++)
  {delete[] Total_rho[i];}
  delete[] Total_rho;
  delete[] Prim_exp;   //Only the pointers must be deleted
  delete[] Contr_Coef;
  for(i=0;i<nbasisf;i++)
  {delete[] MOcoefA[i];}
  delete[] MOcoefA;
  if(extra0)
  {delete[] SP_Contr_Coef;}
  if(extra1)
  {
   for(i=0;i<nbasisf;i++)
   {delete[] MOcoefB[i];}
    delete[] MOcoefB;
   for(i=0;i<nbasisf;i++)
   {delete[] Spin_rho[i];}
    delete[] Spin_rho;
  }
 }
 else
 {
  if(im_wfn_wfx)
  {
   for(i=0;i<MO_coef;i++)
   {delete[] MOcoefA_im[i];}
   delete[] MOcoefA_im;
  }
  for(i=0;i<MO_coef;i++)
  {delete[] MOcoefA[i];}
  delete[] MOcoefA;
  for(i=0;i<natoms;i++)
  {delete[] Cartesian_Coor[i];}
  delete[] Cartesian_Coor;
  delete[] Nu_charge;
  delete[] shell_map;
  delete[] shell_type;
  delete[] Prim_exp;
  delete[] Ocupation;
  delete[] SPIN;
 }
}
//Set BETA_MOS
void READ_FCHK_WFN::set_BETA_MOS(bool &setmos)
{
 BETA_MOS=setmos;
}

//////////////////////////////////
//////////////////////////////////
//      Position Space          //
//      ---------------         //
//////////////////////////////////
//////////////////////////////////

////////////////////////////////////////////////////////
//Functions used by the class for external evaluations//
////////////////////////////////////////////////////////
void READ_FCHK_WFN::rho_eval(double point[3],double &result)
{
 int *i,*j;
 i=new int[1];
 j=new int[1];
 double *AO,**AO_grad;
 if(!wfn)
 {
  //point comes from outside
  AO_grad=new double *[3];
  for(i[0]=0;i[0]<3;i[0]++)
  {AO_grad[i[0]]=new double[nbasisf];}
  for(i[0]=0;i[0]<3;i[0]++)
  {for(j[0]=0;j[0]<nbasisf;j[0]++){AO_grad[i[0]][j[0]]=ZERO;}}
  AO=new double[nbasisf];
  for(i[0]=0;i[0]<nbasisf;i[0]++){AO[i[0]]=ZERO;}
  //Run for each shell calculating Point at each shell
  //construction of AOs and AOs gradients
  build_AO_AOgrad(AO,AO_grad,point);
  //Evaluate the density
   density=ZERO;
   for(i[0]=0;i[0]<nbasisf;i[0]++)
   {
    for(j[0]=i[0];j[0]<nbasisf;j[0]++)
    {
     density=density+AO[i[0]]*AO[j[0]]*Total_rho[j[0]][i[0]];
     if(i[0]!=j[0])
     {
     density=density+AO[j[0]]*AO[i[0]]*Total_rho[j[0]][i[0]];
     }
    }
   }
  delete[] AO;
  AO=NULL;
  for(i[0]=0;i[0]<3;i[0]++)
  {delete[] AO_grad[i[0]];AO_grad[i[0]]=NULL;}
  delete[] AO_grad;
  AO_grad=NULL;
 //Return result
 result=density;
 }
 else
 {
  if(im_wfn_wfx)
  {
   complex<double> *NOs,**NOs_grad;
   NOs=new complex<double>[MO_coef];
   NOs_grad=new complex<double>*[3];
   for(i[0]=0;i[0]<3;i[0]++)
   { 
    NOs_grad[i[0]]=new complex<double>[MO_coef];
   }
   build_NO_grad_wfn_allCC(NOs,NOs_grad,point);
   density=ZERO;
   result=ZERO;
   for(i[0]=0;i[0]<MO_coef;i[0]++)
   {
    if(Ocupation[i[0]]!=ZERO)
    {
     density=density+pow(abs(NOs[i[0]]),TWO)*Ocupation[i[0]];
    }
   }
   //Return result
   result=density;
   for(i[0]=0;i[0]<3;i[0]++)
   { 
    delete[] NOs_grad[i[0]]; NOs_grad[i[0]]=NULL;
   }
   delete[] NOs; NOs=NULL;
  }
  else
  {
   double *NOs,**NOs_grad;
   NOs=new double[MO_coef];
   NOs_grad=new double*[3];
   for(i[0]=0;i[0]<3;i[0]++)
   { 
    NOs_grad[i[0]]=new double[MO_coef];
   }
   build_NO_grad_wfn_all(NOs,NOs_grad,point);
   density=ZERO;
   result=ZERO;
   for(i[0]=0;i[0]<MO_coef;i[0]++)
   {
    if(Ocupation[i[0]]!=ZERO)
    {
     density=density+NOs[i[0]]*NOs[i[0]]*Ocupation[i[0]];
    }
   }
   //Return result
   result=density;
   for(i[0]=0;i[0]<3;i[0]++)
   { 
    delete[] NOs_grad[i[0]]; NOs_grad[i[0]]=NULL;
   }
   delete[] NOs; NOs=NULL;
   delete[] NOs_grad; NOs_grad=NULL;
  }
 }
 delete[] i;
 delete[] j;
 i=NULL;
 j=NULL;
}
void READ_FCHK_WFN::rho_eval_a_b(double point[3],double &rhoa,double &rhob)
{
 int *i,*j;
 i=new int[1];j=new int[1];
 double *AO,**AO_grad;
 if(!wfn)
 {
  if(extra1) //If there are alpha and beta MO
  {
   double *density2;
   density2=new double[1];density2[0]=ZERO;
   AO_grad=new double *[3];
   for(i[0]=0;i[0]<3;i[0]++)
   {AO_grad[i[0]]=new double[nbasisf];}
   for(i[0]=0;i[0]<3;i[0]++)
   {for(j[0]=0;j[0]<nbasisf;j[0]++){AO_grad[i[0]][j[0]]=ZERO;}}
   AO=new double[nbasisf];
   for(i[0]=0;i[0]<nbasisf;i[0]++){AO[i[0]]=ZERO;}
   //Run for each shell calculating Point at each shell
   //construction of AOs and AOs gradients
   build_AO_AOgrad(AO,AO_grad,point);
   //Evaluate the density
    density=ZERO;
    for(i[0]=0;i[0]<nbasisf;i[0]++)
    {
     for(j[0]=i[0];j[0]<nbasisf;j[0]++)
     {
      density=density+AO[i[0]]*AO[j[0]]*(Total_rho[j[0]][i[0]]+Spin_rho[j[0]][i[0]])/TWO;
      density2[0]=density2[0]+AO[i[0]]*AO[j[0]]*(Total_rho[j[0]][i[0]]-Spin_rho[j[0]][i[0]])/TWO;
      if(i[0]!=j[0])
      {
      density=density+AO[j[0]]*AO[i[0]]*(Total_rho[j[0]][i[0]]+Spin_rho[j[0]][i[0]])/TWO;
      density2[0]=density2[0]+AO[i[0]]*AO[j[0]]*(Total_rho[j[0]][i[0]]-Spin_rho[j[0]][i[0]])/TWO;
      }
     }
    }
   delete[] AO;
   AO=NULL;
   for(i[0]=0;i[0]<3;i[0]++)
   {delete[] AO_grad[i[0]];AO_grad[i[0]]=NULL;}
   delete[] AO_grad;
   AO_grad=NULL;
  //Return result
  rhoa=density;
  rhob=density2[0];
  delete[] density2;
  density2=NULL;
  }
  else //No alpha beta MO only alpha->Close-Shell
  {
   rho_eval(point,rhoa);
   rhoa=rhoa/TWO;
   rhob=rhoa;
  }
 }
 else
 {
  if(!open_shell)
  {
   rho_eval(point,rhoa);
   rhoa=rhoa/TWO;
   rhob=rhoa;
  }
  else
  {
   rhoa=ZERO;
   rhob=ZERO;
   if(!correlated)
   {
    if(im_wfn_wfx)
    {
     complex<double> *NOs,**NOs_grad;
     NOs=new complex<double>[MO_coef];
     NOs_grad=new complex<double>*[3];
     for(i[0]=0;i[0]<3;i[0]++)
     { 
      NOs_grad[i[0]]=new complex<double>[MO_coef];
     }
     build_NO_grad_wfn_allCC(NOs,NOs_grad,point);
     for(i[0]=0;i[0]<MO_coef;i[0]++)
     {
      if(Ocupation[i[0]]!=ZERO)
      {
       if(i[0]%2==0)
       {
        rhoa=rhoa+pow(abs(NOs[i[0]]),TWO)*Ocupation[i[0]]; //calculate rhoa
       }
       else
       {
        rhob=rhob+pow(abs(NOs[i[0]]),TWO)*Ocupation[i[0]]; //calculate rhob
       }
      }
     }
     for(i[0]=0;i[0]<3;i[0]++)
     { 
      delete[] NOs_grad[i[0]]; NOs_grad[i[0]]=NULL;
     }
     delete[] NOs; NOs=NULL;
     delete[] NOs_grad; NOs_grad=NULL;
    }
    else
    {
     double *NOs,**NOs_grad;
     NOs=new double[MO_coef];
     NOs_grad=new double*[3];
     for(i[0]=0;i[0]<3;i[0]++)
     { 
      NOs_grad[i[0]]=new double[MO_coef];
     }
     build_NO_grad_wfn_all(NOs,NOs_grad,point);
     for(i[0]=0;i[0]<MO_coef;i[0]++)
     {
      if(Ocupation[i[0]]!=ZERO)
      {
       if(i[0]%2==0)
       {
        rhoa=rhoa+NOs[i[0]]*NOs[i[0]]*Ocupation[i[0]]; //calculate rhoa
       }
       else
       {
        rhob=rhob+NOs[i[0]]*NOs[i[0]]*Ocupation[i[0]]; //calculate rhob
       }
      }
     }
     for(i[0]=0;i[0]<3;i[0]++)
     { 
      delete[] NOs_grad[i[0]]; NOs_grad[i[0]]=NULL;
     }
     delete[] NOs; NOs=NULL;
     delete[] NOs_grad; NOs_grad=NULL;
    }
   }
   else      //Not calculate for Open Shell correlated .wfn files
   {
    error_opens_wfn=true;
   }
  }
 }
 delete[] i;
 delete[] j;
 i=NULL;
 j=NULL;
}
void READ_FCHK_WFN::rho_grad(double point[3],double Grad[3])
{
 int *i,*j,*k;
 i=new int[1];j=new int[1];k=new int[1];
 double *AO,**AO_grad;
 if(!wfn)
 {
  AO_grad=new double *[3];
  for(i[0]=0;i[0]<3;i[0]++)
  {AO_grad[i[0]]=new double[nbasisf];}
  for(i[0]=0;i[0]<3;i[0]++)
  {for(j[0]=0;j[0]<nbasisf;j[0]++){AO_grad[i[0]][j[0]]=ZERO;}}
  AO=new double[nbasisf];
  for(i[0]=0;i[0]<nbasisf;i[0]++){AO[i[0]]=ZERO;}
  //Run for each shell calculating Point at each shell
  //construction of AOs and AOs gradients
  build_AO_AOgrad(AO,AO_grad,point);
  //Evaluate Grad[3]
  for(i[0]=0;i[0]<3;i[0]++){Grad[i[0]]=ZERO;}
   for(i[0]=0;i[0]<nbasisf;i[0]++)
   {
    for(j[0]=i[0];j[0]<nbasisf;j[0]++)
    {
     for(k[0]=0;k[0]<3;k[0]++)
     {
      Grad[k[0]]=Grad[k[0]]+(AO[i[0]]*AO_grad[k[0]][j[0]]+AO[j[0]]*AO_grad[k[0]][i[0]])*Total_rho[j[0]][i[0]];
      if(i[0]!=j[0])
      {
      Grad[k[0]]=Grad[k[0]]+(AO[i[0]]*AO_grad[k[0]][j[0]]+AO[j[0]]*AO_grad[k[0]][i[0]])*Total_rho[j[0]][i[0]];
      }
     }
    }
   }
  delete[] AO;
  AO=NULL;
  for(i[0]=0;i[0]<3;i[0]++)
  {delete[] AO_grad[i[0]];AO_grad[i[0]]=NULL;}
  delete[] AO_grad;
  AO_grad=NULL;
 }
 else
 { 
  if(im_wfn_wfx)
  {
   complex<double> *MOs,**MO_grad;
   MOs=new complex<double>[MO_coef];
   MO_grad=new complex<double>*[3];
   for(i[0]=0;i[0]<3;i[0]++)
   {
    MO_grad[i[0]]=new complex<double>[MO_coef];
   }
   build_NO_grad_wfn_allCC(MOs,MO_grad,point);
   for(i[0]=0;i[0]<3;i[0]++)
   {
    complex<double>res(ZERO,ZERO);
    for(j[0]=0;j[0]<MO_coef;j[0]++)
    {
     // OCC Grad (MO* x MO) 
     res=res+Ocupation[j[0]]*(conj(MOs[j[0]])*MO_grad[i[0]][j[0]]+conj(MO_grad[i[0]][j[0]])*MOs[j[0]]);
    }
    Grad[i[0]]=real(res);
   }
   for(i[0]=0;i[0]<3;i[0]++)
   {
    delete[] MO_grad[i[0]];MO_grad[i[0]]=NULL;
   }
   delete[] MO_grad;MO_grad=NULL;
   delete[] MOs;MOs=NULL;
  }
  else
  {
   double *MOs,**MO_grad;
   MO_grad=new double*[3];
   for(i[0]=0;i[0]<3;i[0]++)
   {MO_grad[i[0]]=new double[MO_coef];}
   MOs=new double[MO_coef];
   for(i[0]=0;i[0]<3;i[0]++)
   {
    for(j[0]=0;j[0]<MO_coef;j[0]++)
    {
     MO_grad[i[0]][j[0]]=ZERO;
    }
   }
   for(i[0]=0;i[0]<MO_coef;i[0]++)
   {
    MOs[i[0]]=ZERO;
   }
   build_NO_grad_wfn_all(MOs,MO_grad,point);
   for(i[0]=0;i[0]<3;i[0]++)
   {
    Grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<MO_coef;j[0]++)
    {
     Grad[i[0]]=Grad[i[0]]+Ocupation[j[0]]*TWO*MO_grad[i[0]][j[0]]*MOs[j[0]];
    }
   }
   for(i[0]=0;i[0]<3;i[0]++)
   {delete[] MO_grad[i[0]];MO_grad[i[0]]=NULL;}
   delete[] MO_grad;MO_grad=NULL;
   delete[] MOs;MOs=NULL;
  }
 }
 delete[] i;
 delete[] j;
 delete[] k;
 i=NULL;
 j=NULL;
 k=NULL;
}
void READ_FCHK_WFN::rho_grad_a_b(double point[3],double Grada[3],double Gradb[3])
{
 int i,j,k;
 double *AO,**AO_grad;
 for(i=0;i<3;i++){Grada[i]=ZERO;Gradb[i]=ZERO;}
 if(!wfn)
 {
  if(extra1) //extra1 alpha and beta MO
  {
   AO_grad=new double *[3];
   for(i=0;i<3;i++)
   {AO_grad[i]=new double[nbasisf];}
   for(i=0;i<3;i++)
   {for(j=0;j<nbasisf;j++){AO_grad[i][j]=ZERO;}}
   AO=new double[nbasisf];
   for(i=0;i<nbasisf;i++){AO[i]=ZERO;}
   //Run for each shell calculating Point at each shell
   //construction of AOs and AOs gradients
   build_AO_AOgrad(AO,AO_grad,point);
   //Evaluate Grad[3]
   for(i=0;i<3;i++)
   {
    Grada[i]=ZERO;
    Gradb[i]=ZERO;
   }
    for(i=0;i<nbasisf;i++)
    {
     for(j=i;j<nbasisf;j++)
     {
      for(k=0;k<3;k++)
      {
       Grada[k]=Grada[k]+(AO[i]*AO_grad[k][j]+AO[j]*AO_grad[k][i])*(Total_rho[j][i]
       +Spin_rho[j][i])/TWO;
       Gradb[k]=Gradb[k]+(AO[i]*AO_grad[k][j]+AO[j]*AO_grad[k][i])*(Total_rho[j][i]
       -Spin_rho[j][i])/TWO;
       if(i!=j)
       {
       Grada[k]=Grada[k]+(AO[i]*AO_grad[k][j]+AO[j]*AO_grad[k][i])*(Total_rho[j][i]
       +Spin_rho[j][i])/TWO;
       Gradb[k]=Gradb[k]+(AO[i]*AO_grad[k][j]+AO[j]*AO_grad[k][i])*(Total_rho[j][i]
       -Spin_rho[j][i])/TWO;
       }
      }
     }
    }
   delete[] AO;
   AO=NULL;
   for(i=0;i<3;i++)
   {delete[] AO_grad[i];AO_grad[i]=NULL;}
   delete[] AO_grad;
   AO_grad=NULL;
  }
  else //Only alpha NO -> Close-Shell
  {
   rho_grad(point,Grada);
   for(i=0;i<3;i++)
   {
    Grada[i]=Grada[i]/TWO;
    Gradb[i]=Grada[i];
   }
  }
 }
 else
 {
  for(i=0;i<3;i++){Grada[i]=ZERO;}
  if(!open_shell)
  {
   rho_grad(point,Gradb);
   for(i=0;i<3;i++)
   {
    Grada[i]=Gradb[i]/TWO;
    Gradb[i]=Grada[i];
   }
  }
  else
  {
   if(!correlated)
   {
    if(im_wfn_wfx)
    {
     complex<double> *MOs,**MO_grad;
     MOs=new complex<double>[MO_coef];
     MO_grad=new complex<double>*[3];
     for(i=0;i<3;i++)
     {
      MO_grad[i]=new complex<double>[MO_coef];
     }
     build_NO_grad_wfn_allCC(MOs,MO_grad,point);
     for(i=0;i<3;i++)
     {
      complex<double>resa(ZERO,ZERO);
      complex<double>resb(ZERO,ZERO);
      for(j=0;j<MO_coef;j++)
      {
       // OCC Grad (MO* x MO) 
       if(i%2==0)
       {
        resa=resa+Ocupation[j]*(conj(MOs[j])*MO_grad[i][j]+conj(MO_grad[i][j])*MOs[j]);
       }
       else
       {
        resb=resb+Ocupation[j]*(conj(MOs[j])*MO_grad[i][j]+conj(MO_grad[i][j])*MOs[j]);
       }
      }
      Grada[i]=real(resa);
      Gradb[i]=real(resb);
     }
     for(i=0;i<3;i++)
     {
      delete[] MO_grad[i];MO_grad[i]=NULL;
     }
     delete[] MO_grad;MO_grad=NULL;
     delete[] MOs;MOs=NULL;
    }
    else
    {
     double **MO_grad,*MOs;
     MO_grad=new double*[3];
     for(i=0;i<3;i++)
     {MO_grad[i]=new double[MO_coef];}
     MOs=new double[MO_coef];
     build_NO_grad_wfn_all(MOs,MO_grad,point);
     for(i=0;i<MO_coef;i++)
     {
      if(Ocupation[i]!=ZERO)
      {
       if(i%2==0)
       {
        for(k=0;k<3;k++)
        {
         Grada[k]=Grada[k]+TWO*MOs[i]*Ocupation[i]*MO_grad[i][k]; //calculate grada
        }
       }
       else
       {
        for(k=0;k<3;k++)
        {
         Gradb[k]=Gradb[k]+TWO*MOs[i]*Ocupation[i]*MO_grad[i][k]; //calculate grada
        }
       }
      }
     }
     //delete arrays
     for(i=0;i<3;i++)
     {delete[] MO_grad[i];MO_grad[i]=NULL;}
     delete[] MO_grad;MO_grad=NULL;
     delete[] MOs;MOs=NULL;
    }
   }
   else
   {}
  }
 }
}
void READ_FCHK_WFN::rho_lapl(double point[3],double &result)
{
//point comes from outside
 int l,m;//Do not use i or j as they exist only once!
 double delta=1.0e-4,temp;
 double AUX1[3];
 result=ZERO; 
 for(l=0;l<3;l++)
 {for(m=0;m<3;m++){AUX1[m]=point[m];}
    AUX1[l]=AUX1[l]+delta;
    rho_eval(AUX1,temp);
    result=result+temp;
    AUX1[l]=AUX1[l]-TWO*delta;
    rho_eval(AUX1,temp);
    result=result+temp;
 }
 rho_eval(point,temp);
 result=(result-SIX*temp)/pow(delta,TWO);
}
void READ_FCHK_WFN::rho_lapl_a_b(double point[3],double &result1,double &result2)
{
 if(!wfn)
 {
  if(extra1)
  {
   //point comes from outside
   int l,m;//Do not use i or j as they exist only once!
   double delta=1.0e-4,temp,temp2;
   double AUX1[3];
   result1=ZERO;
   result2=ZERO;
   for(l=0;l<3;l++)
   {for(m=0;m<3;m++){AUX1[m]=point[m];}
    AUX1[l]=AUX1[l]+delta;
    rho_eval_a_b(AUX1,temp,temp2);
    result1=result1+temp;
    result2=result2+temp2;
    AUX1[l]=AUX1[l]-TWO*delta;
    rho_eval_a_b(AUX1,temp,temp2);
    result1=result1+temp;
    result2=result2+temp2;
   }
  rho_eval_a_b(point,temp,temp2);
  result1=(result1-SIX*temp)/pow(delta,TWO);
  result2=(result2-SIX*temp2)/pow(delta,TWO);
  }
  else
  {
   rho_lapl(point,result1);
   result1=result1/TWO;
   result2=result1;
  }
 }
 else
 {
  if((!correlated) || (!open_shell))
  {      //point comes from outside
   int l,m;//Do not use i or j as they exist only once!
   double delta=1.0e-4,temp,temp2;
   double AUX1[3];
   result1=ZERO;
   result2=ZERO;
   for(l=0;l<3;l++)
   {for(m=0;m<3;m++){AUX1[m]=point[m];}
    AUX1[l]=AUX1[l]+delta;
    rho_eval_a_b(AUX1,temp,temp2);
    result1=result1+temp;
    result2=result2+temp2;
    AUX1[l]=AUX1[l]-TWO*delta;
    rho_eval_a_b(AUX1,temp,temp2);
    result1=result1+temp;
    result2=result2+temp2;
   }
   rho_eval_a_b(point,temp,temp2);
   result1=(result1-SIX*temp)/pow(delta,TWO);
   result2=(result2-SIX*temp2)/pow(delta,TWO);
  }
  else
  {}
 }
}
void READ_FCHK_WFN::rho_hessian(double point[3],double **hess,double *grad,double &density)
{
 int l;
 double h=1.0e-4,h2,temp[3];
 double f_x_y_z;
 double f_xph_y_z,f_x_yph_z,f_x_y_zph;
 double f_xmh_y_z,f_x_ymh_z,f_x_y_zmh;
 double f_xph_yph_z,f_xph_y_zph,f_x_yph_zph;
 double f_xmh_ymh_z,f_xmh_y_zmh,f_x_ymh_zmh;
 // h^2
 h2=h*h;
 // f_x_y_z
 rho_eval(point,f_x_y_z);
 // f_xph_y_z
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[0]+=h;
 rho_eval(temp,f_xph_y_z);
 // f_x_yph_z
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[1]+=h;
 rho_eval(temp,f_x_yph_z);
 // f_x_y_zph
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[2]+=h;
 rho_eval(temp,f_x_y_zph);
 // f_xph_yph_z
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[0]+=h;
 temp[1]+=h;
 rho_eval(temp,f_xph_yph_z);
 // f_xph_y_zph
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[0]+=h;
 temp[2]+=h;
 rho_eval(temp,f_xph_y_zph);
 // f_x_yph_zph
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[1]+=h;
 temp[2]+=h;
 rho_eval(temp,f_x_yph_zph);
 // f_xmh_y_z
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[0]-=h;
 rho_eval(temp,f_xmh_y_z);
 // f_x_ymh_z
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[1]-=h;
 rho_eval(temp,f_x_ymh_z);
 // f_x_y_zmh
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[2]-=h;
 rho_eval(temp,f_x_y_zmh);
 // f_xmh_ymh_z
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[0]-=h;
 temp[1]-=h;
 rho_eval(temp,f_xmh_ymh_z);
 // f_xmh_y_zmh
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[0]-=h;
 temp[2]-=h;
 rho_eval(temp,f_xmh_y_zmh);
 // f_x_ymh_zmh
 for(l=0;l<3;l++){temp[l]=point[l];}
 temp[1]-=h;
 temp[2]-=h;
 rho_eval(temp,f_x_ymh_zmh);
 // Hessian matrix
 hess[0][0]=(f_xph_y_z+f_xmh_y_z-TWO*f_x_y_z)/h2;
 hess[1][1]=(f_x_yph_z+f_x_ymh_z-TWO*f_x_y_z)/h2;
 hess[2][2]=(f_x_y_zph+f_x_y_zmh-TWO*f_x_y_z)/h2;
 hess[0][1]=(f_xph_yph_z-f_xph_y_z-f_x_yph_z+TWO*f_x_y_z-f_xmh_y_z-f_x_ymh_z+f_xmh_ymh_z)/(TWO*h2);
 hess[0][2]=(f_xph_y_zph-f_xph_y_z-f_x_y_zph+TWO*f_x_y_z-f_xmh_y_z-f_x_y_zmh+f_xmh_y_zmh)/(TWO*h2);
 hess[1][2]=(f_x_yph_zph-f_x_yph_z-f_x_y_zph+TWO*f_x_y_z-f_x_ymh_z-f_x_y_zmh+f_x_ymh_zmh)/(TWO*h2);
 hess[1][0]=hess[0][1];
 hess[2][0]=hess[0][2];
 hess[2][1]=hess[1][2];
 // Gradient cart. components
 grad[0]=(f_xph_y_z-f_xmh_y_z)/(TWO*h);
 grad[1]=(f_x_yph_z-f_x_ymh_z)/(TWO*h);
 grad[2]=(f_x_y_zph-f_x_y_zmh)/(TWO*h);
 // Density
 density=f_x_y_z;
}
void READ_FCHK_WFN::orb_grad(double point[3],double **res)
{
 int i,j,numbasis;
 double evaluation,AUX[3];
 if(!wfn)
 {
  if(extra1) //extra1 true for alpha and beta
  {numbasis=(MO_coef+MO_beta_coef)/nbasisf;} //numbasis all posible MO
  else                                     //nbasisf number of AO
  {numbasis=MO_coef/nbasisf;}
  for(i=0;i<4;i++)
  {for(j=0;j<numbasis;j++)
   {res[i][j]=ZERO;}
  }
  for(i=0;i<numbasis;i++) //numbasis all posible MOs
  {
   build_NO_grad_fchk(evaluation,AUX,point,i);
   res[0][i]=evaluation;
   for(j=0;j<3;j++)
   {res[j+1][i]=AUX[j];}
  }
 }
 else
 {
  double **MO_grad,*MOs;
  //dynamic arrays
  MO_grad=new double*[3];
  for(i=0;i<3;i++)
  {MO_grad[i]=new double[MO_coef];}
  MOs=new double[MO_coef];
  build_NO_grad_wfn_all(MOs,MO_grad,point);
  for(i=0;i<nbasisf;i++)
  {
   res[0][i]=MOs[i];
   for(j=0;j<3;j++)
   {res[j+1][i]=MO_grad[j][i];}
  }
  //delete arrays
  for(i=0;i<3;i++)
  {delete[] MO_grad[i];MO_grad[i]=NULL;}
  delete[] MO_grad;MO_grad=NULL;
  delete[] MOs;MOs=NULL;
 }
}

void READ_FCHK_WFN::orb_gradCC(double point[3],complex<double> **res)
{
 int i,j;
 if(!wfn)
 {
  cout<<"Warning! Not implemented for FCHK files!"<<endl;
 }
 else
 {
  complex<double> **MO_grad,*MOs;
  //dynamic arrays
  MO_grad=new complex<double>*[3];
  MOs=new complex<double>[MO_coef];
  for(i=0;i<3;i++)
  {MO_grad[i]=new complex<double>[MO_coef];}
  build_NO_grad_wfn_allCC(MOs,MO_grad,point);
  for(i=0;i<nbasisf;i++)
  {
   res[0][i]=MOs[i];
   for(j=0;j<3;j++)
   {res[j+1][i]=MO_grad[j][i];}
  }
  //delete arrays
  for(i=0;i<3;i++)
  {delete[] MO_grad[i];MO_grad[i]=NULL;}
  delete[] MO_grad;MO_grad=NULL;
  delete[] MOs;MOs=NULL;
 }
}

void READ_FCHK_WFN::build_NO_grad_wfn(double &NO,double grad[3],double Point[3], int &numMO)
{
  int *i,*j,*k,*nlm;
  double *pos_nuclei,*exponent;
  int **Quant;
  i=new int[1];j=new int[1];k=new int[1];
  nlm=new int[3];
  pos_nuclei=new double[3];exponent=new double[1];
  ///////////////////////////////////////
  //dynamic arrays
  Quant=new int*[35];
  for(i[0]=0;i[0]<35;i[0]++)
  {Quant[i[0]]=new int[3];}
  ///////////////////////////////////////
  //Initialize
  Quant_fill(Quant,0); //0 because no type is needed
  for(j[0]=0;j[0]<3;j[0]++)
  {grad[j[0]]=ZERO;}
  NO=ZERO;
  //////////////////////////////////////////
  //build gradients for the NOs and the NOs
  for(i[0]=0;i[0]<3;i[0]++)
  {
   for(j[0]=0;j[0]<nprimitv;j[0]++)
   {
    for(k[0]=0;k[0]<3;k[0]++)
    {
     pos_nuclei[k[0]]=Cartesian_Coor[shell_map[j[0]]-1][k[0]];
     nlm[k[0]]=Quant[shell_type[j[0]]-1][k[0]];
    }
    exponent[0]=Prim_exp[j[0]];
    grad[i[0]]=grad[i[0]]+MOcoefA[numMO][j[0]]*grad_Primitive_wfn(Point,i[0],
    pos_nuclei,nlm,exponent[0]);
    if(i[0]==0)
    {
     NO=NO+MOcoefA[numMO][j[0]]*eval_Primitive_wfn(Point,pos_nuclei,nlm,exponent[0]);
    }
   }
  }
  for(i[0]=0;i[0]<35;i[0]++)
  {delete[] Quant[i[0]];Quant[i[0]]=NULL;}
  delete[] Quant;
  Quant=NULL;
  delete[] i; delete[] j; delete[] k; delete[] nlm;
  delete[] pos_nuclei; delete[] exponent;
  i=NULL; j=NULL; k=NULL; nlm=NULL;
  pos_nuclei=NULL; exponent=NULL;
}

void READ_FCHK_WFN::build_NO_grad_wfn_all(double *NO,double **NO_grad,double Point[3])
{
 int i,iprim,ishell,idir,iNO,nlm[3];
 double pos_nuclei[3],exponent,grad_prim,val_prim;
 int **Quant;
 ///////////////////////////////////////
 //Dynamic array
 Quant=new int*[35];
 for(i=0;i<35;i++)
 {Quant[i]=new int[3];}
 ///////////////////////////////////////
 //Initialize
 Quant_fill(Quant,0); //0 because no type is needed
 for(idir=0;idir<3;idir++)
 {
  for(iNO=0;iNO<MO_coef;iNO++)
  {
   NO_grad[idir][iNO]=ZERO;
   if(idir==0) 
   {
    NO[iNO]=ZERO;
   }
  }
 }
 //////////////////////////////////////////
 //Build NOs and their gradients 
 for(iprim=0;iprim<nprimitv;iprim++)
 {
  for(ishell=0;ishell<3;ishell++)
  {
   pos_nuclei[ishell]=Cartesian_Coor[shell_map[iprim]-1][ishell];
   nlm[ishell]=Quant[shell_type[iprim]-1][ishell];
  }
  exponent=Prim_exp[iprim];
  val_prim=eval_Primitive_wfn(Point,pos_nuclei,nlm,exponent);
  for(idir=0;idir<3;idir++)
  {
   grad_prim=grad_Primitive_wfn(Point,idir,pos_nuclei,nlm,exponent);
   for(iNO=0;iNO<MO_coef;iNO++)
   {
    NO_grad[idir][iNO]=NO_grad[idir][iNO]+MOcoefA[iNO][iprim]*grad_prim;
    if(idir==0)
    {
     NO[iNO]=NO[iNO]+MOcoefA[iNO][iprim]*val_prim;
    }
   }
  }
 }
 //Delete allocated arrays
 for(i=0;i<35;i++)
 {delete[] Quant[i];Quant[i]=NULL;}
 delete[] Quant;Quant=NULL;
}

void READ_FCHK_WFN::build_NO_grad_wfn_allCC(complex<double> *NO,complex<double> **NO_grad,double Point[3])
{
 int i,iprim,ishell,idir,iNO,nlm[3];
 double pos_nuclei[3],exponent,grad_prim,val_prim;
 int **Quant;
 complex<double>ztmp0(ZERO,ZERO);
 ///////////////////////////////////////
 //Dynamic array
 Quant=new int*[35];
 for(i=0;i<35;i++)
 {Quant[i]=new int[3];}
 ///////////////////////////////////////
 //Initialize
 Quant_fill(Quant,0); //0 because no type is needed
 for(idir=0;idir<3;idir++)
 {
  for(iNO=0;iNO<MO_coef;iNO++)
  {
   NO_grad[idir][iNO]=ztmp0;
   if(idir==0)
   {
    NO[iNO]=ztmp0;
   }
  }
 }
 //////////////////////////////////////////
 //Build NOs and their gradients 
 for(iprim=0;iprim<nprimitv;iprim++)
 {
  for(ishell=0;ishell<3;ishell++)
  {
   pos_nuclei[ishell]=Cartesian_Coor[shell_map[iprim]-1][ishell];
   nlm[ishell]=Quant[shell_type[iprim]-1][ishell];
  }
  exponent=Prim_exp[iprim];
  val_prim=eval_Primitive_wfn(Point,pos_nuclei,nlm,exponent);
  complex<double>eval_prim(val_prim,ZERO);
  for(idir=0;idir<3;idir++)
  {
   grad_prim=grad_Primitive_wfn(Point,idir,pos_nuclei,nlm,exponent);
   complex<double>eval_grad(grad_prim,ZERO);
   for(iNO=0;iNO<MO_coef;iNO++)
   {
    complex<double>mo_coef(MOcoefA[iNO][iprim],MOcoefA_im[iNO][iprim]);
    NO_grad[idir][iNO]=NO_grad[idir][iNO]+mo_coef*eval_grad;
    if(idir==0)
    {
     NO[iNO]=NO[iNO]+mo_coef*eval_prim;
    }
   }
  }
 }
 //Delete allocated arrays
 for(i=0;i<35;i++)
 {delete[] Quant[i];Quant[i]=NULL;}
 delete[] Quant;Quant=NULL;
}

void READ_FCHK_WFN::build_NO_grad_fchk(double &NO,double grad[3],double point[3],int &numMO)
{
  int *i,*j;
  i=new int [1];j=new int[1];
  double *AO,**AO_grad;
  //point comes from outside
  NO=ZERO;
  AO_grad=new double *[3];
  for(i[0]=0;i[0]<3;i[0]++)
  {AO_grad[i[0]]=new double[nbasisf];}
  for(i[0]=0;i[0]<3;i[0]++)
  {for(j[0]=0;j[0]<nbasisf;j[0]++){AO_grad[i[0]][j[0]]=ZERO;}}
  AO=new double[nbasisf];
  for(i[0]=0;i[0]<nbasisf;i[0]++)
  {
   AO[i[0]]=ZERO;
  }
  build_AO_AOgrad(AO,AO_grad,point);
 //Build NO
 if(extra1)
 {
  if(numMO%2==0)
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    NO=NO+S[j[0]][numMO/2]*AO[j[0]];
   } //S already has S^-1/2 Eigenv for NO alpha
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     grad[i[0]]=grad[i[0]]+S[j[0]][numMO/2]*AO_grad[i[0]][j[0]];
    } //S already has S^-1/2 Eigenv for NO alpha
   }
  }
  else
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    NO=NO+Sbeta[j[0]][(numMO-1)/2]*AO[j[0]];
   } //Sbeta already has S^-1/2 Eigenv for NO beta
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     grad[i[0]]=grad[i[0]]+Sbeta[j[0]][(numMO-1)/2]*AO_grad[i[0]][j[0]];
    } //Sbeta already has S^-1/2 Eigenv for NO beta
   }
  }
 }
 else
 {
  for(j[0]=0;j[0]<nbasisf;j[0]++)
  {
   NO=NO+S[j[0]][numMO]*AO[j[0]];
  } //S already has S^-1/2 Eigenv for NO alpha
  for(i[0]=0;i[0]<3;i[0]++)
  {
   grad[i[0]]=ZERO;
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    grad[i[0]]=grad[i[0]]+S[j[0]][numMO]*AO_grad[i[0]][j[0]];
   } //S already has S^-1/2 Eigenv for NO alpha
  }
 }
///////////////////////////////////////////////////////////////////////////////
  delete[] AO;
  AO=NULL;
  for(i[0]=0;i[0]<3;i[0]++)
  {delete[] AO_grad[i[0]];AO_grad[i[0]]=NULL;}
  delete[] AO_grad;
  AO_grad=NULL;
  delete[] i;
  i=NULL;
  delete[] j;
  j=NULL;
}
//Same as the previous one but reading AOs and AOs_grad from outside
void READ_FCHK_WFN::build_NO_grad_fchk2(double *AO,double **AO_grad,double &NO,double grad[3],int &numMO)
{
  int *i,*j;
  i=new int [1];j=new int[1];
  //point comes from outside
  NO=ZERO;
 //Build NO
 if(extra1)
 {
  if(numMO%2==0)
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    NO=NO+S[j[0]][numMO/2]*AO[j[0]];
   } //S already has S^-1/2 Eigenv for NO alpha
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     grad[i[0]]=grad[i[0]]+S[j[0]][numMO/2]*AO_grad[i[0]][j[0]];
    } //S already has S^-1/2 Eigenv for NO alpha
   }
  }
  else
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    NO=NO+Sbeta[j[0]][(numMO-1)/2]*AO[j[0]];
   } //Sbeta already has S^-1/2 Eigenv for NO beta
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     grad[i[0]]=grad[i[0]]+Sbeta[j[0]][(numMO-1)/2]*AO_grad[i[0]][j[0]];
    } //Sbeta already has S^-1/2 Eigenv for NO beta
   }
  }
 }
 else
 {
  for(j[0]=0;j[0]<nbasisf;j[0]++)
  {
   NO=NO+S[j[0]][numMO]*AO[j[0]];
  } //S already has S^-1/2 Eigenv for NO alpha
  for(i[0]=0;i[0]<3;i[0]++)
  {
   grad[i[0]]=ZERO;
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    grad[i[0]]=grad[i[0]]+S[j[0]][numMO]*AO_grad[i[0]][j[0]];
   } //S already has S^-1/2 Eigenv for NO alpha
  }
 }
///////////////////////////////////////////////////////////////////////////////
  delete[] i;
  i=NULL;
  delete[] j;
  j=NULL;
}
void READ_FCHK_WFN::build_MO_grad_fchk(double &MO,double grad[3],double point[3],int &numMO)
{
 int *i,*j;
 i=new int[1];j=new int[1];
 double *AO,**AO_grad;
 //point comes from outside
 MO=ZERO;
 AO_grad=new double *[3];
 for(i[0]=0;i[0]<3;i[0]++)
 {AO_grad[i[0]]=new double[nbasisf];}
 for(i[0]=0;i[0]<3;i[0]++)
 {for(j[0]=0;j[0]<nbasisf;j[0]++){AO_grad[i[0]][j[0]]=ZERO;}}
 AO=new double[nbasisf];
 for(i[0]=0;i[0]<nbasisf;i[0]++)
 {
  AO[i[0]]=ZERO;
 }
 build_AO_AOgrad(AO,AO_grad,point);
 //Build MO
 if(extra1)
 {
  if(numMO%2==0)
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    MO=MO+MOcoefA[numMO/2][j[0]]*AO[j[0]];
   }
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     grad[i[0]]=grad[i[0]]+MOcoefA[numMO/2][j[0]]*AO_grad[i[0]][j[0]];
    }
   }
  }
  else
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    if(BETA_MOS)
    {
     MO=MO+MOcoefB[(numMO-1)/2][j[0]]*AO[j[0]];
    }
    else
    {
     MO=MO+MOcoefA[(numMO-1)/2][j[0]]*AO[j[0]];
    }
   }
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     if(BETA_MOS)
     {
      grad[i[0]]=grad[i[0]]+MOcoefB[(numMO-1)/2][j[0]]*AO_grad[i[0]][j[0]];
     }
     else
     {
      grad[i[0]]=grad[i[0]]+MOcoefA[(numMO-1)/2][j[0]]*AO_grad[i[0]][j[0]];
     }
    }
   }
  }
 }
 else
 {
  for(j[0]=0;j[0]<nbasisf;j[0]++)
  {
   MO=MO+MOcoefA[numMO][j[0]]*AO[j[0]];
  }
  for(i[0]=0;i[0]<3;i[0]++)
  {
   grad[i[0]]=ZERO;
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    grad[i[0]]=grad[i[0]]+MOcoefA[numMO][j[0]]*AO_grad[i[0]][j[0]];
   }
  }
 }
///////////////////////////////////////////////////////////////////////////////
  delete[] AO;
  AO=NULL;
  for(i[0]=0;i[0]<3;i[0]++)
  {delete[] AO_grad[i[0]]; AO_grad[i[0]]=NULL;}
  delete[] AO_grad;
  AO_grad=NULL;
  delete[] i;
  i=NULL;
  delete[] j;
  j=NULL;
}
//Build with AOs comming from outside
void READ_FCHK_WFN::build_MO_grad_fchk2(double *AO,double **AO_grad,double &MO,double grad[3],int &numMO)
{
 int *i,*j;
 i=new int[1];j=new int[1];
 //point comes from outside
 MO=ZERO;
 //Build MO
 if(extra1)
 {
  if(numMO%2==0)
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    MO=MO+MOcoefA[numMO/2][j[0]]*AO[j[0]];
   }
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     grad[i[0]]=grad[i[0]]+MOcoefA[numMO/2][j[0]]*AO_grad[i[0]][j[0]];
    }
   }
  }
  else
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    if(BETA_MOS)
    {
     MO=MO+MOcoefB[(numMO-1)/2][j[0]]*AO[j[0]];
    }
    else
    {
     MO=MO+MOcoefA[(numMO-1)/2][j[0]]*AO[j[0]];
    }
   }
   for(i[0]=0;i[0]<3;i[0]++)
   {
    grad[i[0]]=ZERO;
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     if(BETA_MOS)
     {
      grad[i[0]]=grad[i[0]]+MOcoefB[(numMO-1)/2][j[0]]*AO_grad[i[0]][j[0]];
     }
     else
     {
      grad[i[0]]=grad[i[0]]+MOcoefA[(numMO-1)/2][j[0]]*AO_grad[i[0]][j[0]];
     }
    }
   }
  }
 }
 else
 {
  for(j[0]=0;j[0]<nbasisf;j[0]++)
  {
   MO=MO+MOcoefA[numMO][j[0]]*AO[j[0]];
  }
  for(i[0]=0;i[0]<3;i[0]++)
  {
   grad[i[0]]=ZERO;
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    grad[i[0]]=grad[i[0]]+MOcoefA[numMO][j[0]]*AO_grad[i[0]][j[0]];
   }
  }
 }
///////////////////////////////////////////////////////////////////////////////
  delete[] i;
  i=NULL;
  delete[] j;
  j=NULL;
}
//Give the total basis size,
//half size for close shell
int READ_FCHK_WFN::nbasis()
{
 int number_basis_fun;
 if(wfn)
 {number_basis_fun=nbasisf;}
 else
 {
  if(extra1) //extra1 true for alpha and beta
  {number_basis_fun=(MO_coef+MO_beta_coef)/nbasisf;}
  else
  {number_basis_fun=MO_coef/nbasisf;}
 }
 return number_basis_fun;
}
//Dipole moment contribution from nuclei
void READ_FCHK_WFN::muATOMS(double mu[3])
{
 int i;
 for(i=0;i<3;i++){mu[i]=ZERO;}
 for(i=0;i<natoms;i++)
 {
  mu[0]=mu[0]+(double)Nu_charge[i]*Cartesian_Coor[i][0];
  mu[1]=mu[1]+(double)Nu_charge[i]*Cartesian_Coor[i][1];
  mu[2]=mu[2]+(double)Nu_charge[i]*Cartesian_Coor[i][2];
 }
}
//Quadrupole moment contribution from nuclei
void READ_FCHK_WFN::quadrupoleATOMS(double **quadrupole)
{
 int i,j;
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   quadrupole[i][j]=ZERO;
  }
 }
 for(i=0;i<natoms;i++)
 {
  quadrupole[0][0]=quadrupole[0][0]+(double)Nu_charge[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][0];
  quadrupole[0][1]=quadrupole[0][1]+(double)Nu_charge[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][1];
  quadrupole[0][2]=quadrupole[0][2]+(double)Nu_charge[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][2];
  quadrupole[1][1]=quadrupole[1][1]+(double)Nu_charge[i]*Cartesian_Coor[i][1]*Cartesian_Coor[i][1];
  quadrupole[1][2]=quadrupole[1][2]+(double)Nu_charge[i]*Cartesian_Coor[i][1]*Cartesian_Coor[i][2];
  quadrupole[2][2]=quadrupole[2][2]+(double)Nu_charge[i]*Cartesian_Coor[i][2]*Cartesian_Coor[i][2];
 }
 quadrupole[1][0]=quadrupole[0][1];
 quadrupole[2][0]=quadrupole[0][2];
 quadrupole[2][1]=quadrupole[1][2];
}
//Nuclei contribution for the potential at point r
double READ_FCHK_WFN::Vnuclear(double Point_Vr[3])
{
 int i;
 double res=ZERO,XYZ[3]={ZERO};
 for(i=0;i<natoms;i++)
 {
  XYZ[0]=Point_Vr[0]-Cartesian_Coor[i][0];
  XYZ[1]=Point_Vr[1]-Cartesian_Coor[i][1];
  XYZ[2]=Point_Vr[2]-Cartesian_Coor[i][2];
  res=res+((double)Nu_charge[i])/(norm3D(XYZ)+pow(TEN,-SIX*TWO));
 }
 return res;
}
//Build AOs and AO_Grads outside the class
//Warning also AOs_Grads are needed
void READ_FCHK_WFN::build_AO_AOgrad2(double *AO,double **AO_grad,double Point[3])
{
 if(!wfn)
 {
  int *i,*j,*k,*counter,*counter2,*iprim,*iprim_tot;
  double *Coord_At;
  double *P_exp_send, *C_coef_send, *SPC_coef_send;
  i=new int[1];j=new int[1];k=new int[1];counter=new int[1];
  counter2=new int[1];iprim=new int[1];iprim_tot=new int[1];
  Coord_At=new double[3];
  counter[0]=0;
  counter2[0]=0;
  for(i[0]=0;i[0]<nshells;i[0]++)
  {//save atomic coordinates of the shell, number of primitives
    //send only coefficients and exponents needed for that shell for building AOs.
    for(j[0]=0;j[0]<3;j[0]++)
    {Coord_At[j[0]]=Cartesian_Coor[(shell_map[i[0]]-1)][j[0]];}
    iprim[0]=n_prim_per_shell[i[0]];
    iprim_tot[0]=0;
    for(j[0]=0;j[0]<i[0];j[0]++)
    {iprim_tot[0]=iprim_tot[0]+n_prim_per_shell[j[0]];}
    C_coef_send=new double[iprim[0]];
    if(extra0)
    {SPC_coef_send=new double[iprim[0]];}
    P_exp_send=new double[iprim[0]];
    k[0]=0;
    for(j[0]=iprim_tot[0];j[0]<iprim[0]+iprim_tot[0];j[0]++)
    {C_coef_send[k[0]]=Contr_Coef[j[0]];
     if(extra0)
     {SPC_coef_send[k[0]]=SP_Contr_Coef[j[0]];}
     P_exp_send[k[0]]=Prim_exp[j[0]];
     k[0]++;
    }
    build(iprim[0],shell_type[i[0]],counter[0],C_coef_send,SPC_coef_send,P_exp_send,Point,
    AO,Coord_At,activeSP);
    build_grad(iprim[0],shell_type[i[0]],counter2[0],C_coef_send,SPC_coef_send,P_exp_send,Point,
    AO_grad,Coord_At,activeSP);
    delete[] C_coef_send;
    C_coef_send=NULL;
    if(extra0)
    {delete[] SPC_coef_send;SPC_coef_send=NULL;}
    delete[] P_exp_send;
    P_exp_send=NULL;
  }
 delete[] i;delete[] j;delete[] k;delete[] counter;delete[] counter2;delete[] iprim;delete[] iprim_tot;
 delete[] Coord_At;
 i=NULL;j=NULL;k=NULL;counter=NULL;counter2=NULL;iprim=NULL;iprim_tot=NULL;Coord_At=NULL;
 }
}
////////////////////////////////////////////////////////
//Functions used by the class for internal evaluations//
////////////////////////////////////////////////////////
//build AO and AOgradients for fchk
void READ_FCHK_WFN::build_AO_AOgrad(double *AO,double **AO_grad,double Point[3])
{
  int *i,*j,*k,*counter,*counter2,*iprim,*iprim_tot;
  double *Coord_At;
  double *P_exp_send, *C_coef_send, *SPC_coef_send;
  i=new int[1];j=new int[1];k=new int[1];counter=new int[1];
  counter2=new int[1];iprim=new int[1];iprim_tot=new int[1];
  Coord_At=new double[3];
  counter[0]=0;
  counter2[0]=0;
  for(i[0]=0;i[0]<nshells;i[0]++)
  {//save atomic coordinates of the shell, number of primitives
    //send only coefficients and exponents needed for that shell for building AOs.
    for(j[0]=0;j[0]<3;j[0]++)
    {Coord_At[j[0]]=Cartesian_Coor[(shell_map[i[0]]-1)][j[0]];}
    iprim[0]=n_prim_per_shell[i[0]];
    iprim_tot[0]=0;
    for(j[0]=0;j[0]<i[0];j[0]++)
    {iprim_tot[0]=iprim_tot[0]+n_prim_per_shell[j[0]];}
    C_coef_send=new double[iprim[0]];
    if(extra0)
    {SPC_coef_send=new double[iprim[0]];}
    P_exp_send=new double[iprim[0]];
    k[0]=0;
    for(j[0]=iprim_tot[0];j[0]<iprim[0]+iprim_tot[0];j[0]++)
    {C_coef_send[k[0]]=Contr_Coef[j[0]];
     if(extra0)
     {SPC_coef_send[k[0]]=SP_Contr_Coef[j[0]];}
     P_exp_send[k[0]]=Prim_exp[j[0]];
     k[0]++;
    }
    build(iprim[0],shell_type[i[0]],counter[0],C_coef_send,SPC_coef_send,P_exp_send,Point,
    AO,Coord_At,activeSP);
    build_grad(iprim[0],shell_type[i[0]],counter2[0],C_coef_send,SPC_coef_send,P_exp_send,Point,
    AO_grad,Coord_At,activeSP);
    delete[] C_coef_send;
    C_coef_send=NULL;
    if(extra0)
    {delete[] SPC_coef_send;SPC_coef_send=NULL;}
    delete[] P_exp_send;
    P_exp_send=NULL;
  }
 delete[] i;delete[] j;delete[] k;delete[] counter;delete[] counter2;delete[] iprim;delete[] iprim_tot;
 delete[] Coord_At;
 i=NULL;j=NULL;k=NULL;counter=NULL;counter2=NULL;iprim=NULL;iprim_tot=NULL;Coord_At=NULL;
}
//Build builds the AOs. Recives the atomic coordinates, coeficients, exponents and
//the point where density must be evaluated it decides whether or not use SP and
//calls for a function(Quant_fill) to make the n,l,m.
void READ_FCHK_WFN::build(int &iprim,int &styp,int &counter,double *C_coef_send,
double *SPC_coef_send,double *P_exp_send,double Point[3],double *AO,double Coord_At[3],
bool &activeSP)
{ int *i,*j,*qdim;
  int **Quant;
  double *evaluation;
  i=new int[1];j=new int[1];qdim=new int[1];
  evaluation=new double[1];
  if((styp==0 || styp==1) || (styp==2 || styp==3) || styp==4)
  {
   qdim[0]=(styp+1)*(styp+2)/2;
   Quant=new int*[qdim[0]];
   for(i[0]=0;i[0]<qdim[0];i[0]++)
   {Quant[i[0]]=new int[3];}
   //Quantum array n,l,m
   Quant_fill(Quant,styp);
   for(i[0]=0;i[0]<qdim[0];i[0]++)
   {for(j[0]=0;j[0]<iprim;j[0]++)
    {//eval contains the value of primitive at point
     evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],Quant[i[0]][0],
     Quant[i[0]][1],Quant[i[0]][2]);
     if(!activeSP)
     {
      AO[counter]=AO[counter]+C_coef_send[j[0]]*evaluation[0];
     }
     else
     {
      AO[counter]=AO[counter]+SPC_coef_send[j[0]]*evaluation[0];
     }
    }
    counter++;
   }
   if(activeSP){activeSP=false;}
   for(i[0]=0;i[0]<styp+1;i[0]++)
   {delete[] Quant[i[0]];Quant[i[0]]=NULL;}
   delete Quant;
   Quant=NULL;
  }
  else if(styp==-1)
  {
    int *send_type;
    send_type=new int[1];
    send_type[0]=0;
    build(iprim,send_type[0],counter,C_coef_send,SPC_coef_send,P_exp_send,Point,AO,
    Coord_At,activeSP);
    activeSP=true;
    send_type[0]=1;
    build(iprim,send_type[0],counter,C_coef_send,SPC_coef_send,P_exp_send,Point,AO,
    Coord_At,activeSP);
    delete[] send_type;
    send_type=NULL;
  }
  else if(styp==-2)
  {
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
   int *n,*l,*m;
   n=new int[1];l=new int[1];m=new int[1];
   n[0]=0;l[0]=1;m[0]=2;
   for(j[0]=0;j[0]<iprim;j[0]++)
   {
   //(2,0)
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],n[0],n[0],m[0])
              -HALF*(eval(Coord_At,Point,P_exp_send[j[0]],m[0],n[0],n[0])
              +eval(Coord_At,Point,P_exp_send[j[0]],n[0],m[0],n[0]));
    AO[counter]=AO[counter]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(2,1)+(2,-1)}
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],l[0],n[0],l[0]);
    AO[counter+1]=AO[counter+1]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(2,1)-(2,-1)}
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],n[0],l[0],l[0]);
    AO[counter+2]=AO[counter+2]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(2,2)+(2,-2)}
    evaluation[0]=pow(THREE/FOUR,HALF)*(eval(Coord_At,Point,P_exp_send[j[0]],m[0],n[0],n[0])
              -eval(Coord_At,Point,P_exp_send[j[0]],n[0],m[0],n[0]));
    AO[counter+3]=AO[counter+3]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(2,2)-(2,-2)}
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],l[0],l[0],n[0]);
    AO[counter+4]=AO[counter+4]+C_coef_send[j[0]]*evaluation[0];
   }
   counter=counter+5;
   delete[] n;delete[] m; delete[] l;
   l=NULL;m=NULL;n=NULL;
  }
  else if(styp==-3)
  {
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
   int *n,*l,*m,*o;
   n=new int[1];l=new int[1];m=new int[1];o=new int[1];
   n[0]=0;l[0]=1;m[0]=2;o[0]=3;
   for(j[0]=0;j[0]<iprim;j[0]++)
   {
   //(3,0)
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],n[0],n[0],o[0])
              -(THREE/(TWO*pow(FIVE,HALF)))*(eval(Coord_At,Point,P_exp_send[j[0]],m[0],n[0],l[0])
              +eval(Coord_At,Point,P_exp_send[j[0]],n[0],m[0],l[0]));
    AO[counter]=AO[counter]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(3,1)+(3,-1)}
    evaluation[0]=pow(SIX/FIVE,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],n[0],m[0])
              -(pow(SIX,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],o[0],n[0],n[0])
              -(pow(SIX/FIVE,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],m[0],n[0]);
    AO[counter+1]=AO[counter+1]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(3,1)-(3,-1)}
    evaluation[0]=pow(SIX/FIVE,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],l[0],m[0])
              -(pow(SIX,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],o[0],n[0])
              -(pow(SIX/FIVE,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],m[0],l[0],n[0]);
    AO[counter+2]=AO[counter+2]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(3,2)+(3,-2)}
    evaluation[0]=pow(SIX/EIGHT,HALF)*(eval(Coord_At,Point,P_exp_send[j[0]],m[0],n[0],l[0])
              -eval(Coord_At,Point,P_exp_send[j[0]],n[0],m[0],l[0]));
    AO[counter+3]=AO[counter+3]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(3,2)-(3,-2)}
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],l[0],l[0],l[0]);
    AO[counter+4]=AO[counter+4]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(3,3)+(3,-3)}
    evaluation[0]=(pow(TEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],o[0],n[0],n[0])
              -(THREE/(TWO*pow(TWO,HALF)))*eval(Coord_At,Point,P_exp_send[j[0]],l[0],m[0],n[0]);
    AO[counter+5]=AO[counter+5]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(3,3)-(3,-3)}
    evaluation[0]=-(pow(TEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],o[0],n[0])
              +(THREE/(TWO*pow(TWO,HALF)))*eval(Coord_At,Point,P_exp_send[j[0]],m[0],l[0],n[0]);
    AO[counter+6]=AO[counter+6]+C_coef_send[j[0]]*evaluation[0];
   }
   counter=counter+7;
   delete[] n;delete[] m; delete[] l; delete[] o;
   l=NULL;m=NULL;n=NULL;o=NULL;
  }
  else if(styp==-4)
  {
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
   int *n,*l,*m,*o,*p;
   n=new int[1];l=new int[1];m=new int[1];o=new int[1];p=new int[1];
   n[0]=0;l[0]=1;m[0]=2;o[0]=3;p[0]=4;
   for(j[0]=0;j[0]<iprim;j[0]++)
   {
   //(4,0)
    evaluation[0]=eval(Coord_At,Point,P_exp_send[j[0]],n[0],n[0],p[0])
              +(THREE/EIGHT)*(eval(Coord_At,Point,P_exp_send[j[0]],p[0],n[0],n[0])
              +eval(Coord_At,Point,P_exp_send[j[0]],n[0],p[0],n[0]))
              -(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*(eval(Coord_At,Point,P_exp_send[j[0]],m[0],n[0],m[0])
              +eval(Coord_At,Point,P_exp_send[j[0]],n[0],m[0],m[0])
              -HALF*HALF*eval(Coord_At,Point,P_exp_send[j[0]],m[0],m[0],n[0]));
    AO[counter]=AO[counter]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,1)+(4,-1)}
    evaluation[0]=(pow(TEN/SEVEN,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],n[0],o[0])
                  -(THREE*pow(TEN/SEVEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],o[0],n[0],l[0])
                  -(THREE*pow(TWO/SEVEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],m[0],l[0]));
    AO[counter+1]=AO[counter+1]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,1)-(4,-1)}
    evaluation[0]=(pow(TEN/SEVEN,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],l[0],o[0])
                  -(THREE*pow(TEN/SEVEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],o[0],l[0])
                  -(THREE*pow(TWO/SEVEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],m[0],l[0],l[0]));
    AO[counter+2]=AO[counter+2]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,2)+(4,-2)}
    evaluation[0]=(THREE*pow(THREE/SEVEN,HALF)/TWO)*eval(Coord_At,Point,P_exp_send[j[0]],m[0],n[0],m[0])
                 -(THREE*pow(THREE/SEVEN,HALF)/TWO)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],m[0],m[0])
                 -HALF*HALF*pow(FIVE,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],p[0],n[0],n[0])
                 +HALF*HALF*pow(FIVE,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],p[0],n[0]);
    AO[counter+3]=AO[counter+3]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,2)-(4,-2)}
    evaluation[0]=(THREE/pow(SEVEN,HALF))*eval(Coord_At,Point,P_exp_send[j[0]],l[0],l[0],m[0])
                 -HALF*pow(FIVE/SEVEN,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],o[0],l[0],n[0])
                 -HALF*pow(FIVE/SEVEN,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],o[0],n[0]);
    AO[counter+4]=AO[counter+4]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,3)+(4,-3)}
    evaluation[0]=(pow(TEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],o[0],n[0],l[0])
                 -THREE*(pow(TWO,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],m[0],l[0]);
    AO[counter+5]=AO[counter+5]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,3)-(4,-3)}
    evaluation[0]=(-pow(TEN,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],o[0],l[0])
                 +THREE*(pow(TWO,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],m[0],l[0],l[0]);
    AO[counter+6]=AO[counter+6]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,4)+(4,-4)}
    evaluation[0]=(pow(SEVEN*FIVE,HALF)/EIGHT)*eval(Coord_At,Point,P_exp_send[j[0]],p[0],n[0],n[0])
                 +(pow(SEVEN*FIVE,HALF)/EIGHT)*eval(Coord_At,Point,P_exp_send[j[0]],n[0],p[0],n[0])
                 -(THREE*pow(THREE,HALF)/FOUR)*eval(Coord_At,Point,P_exp_send[j[0]],m[0],m[0],n[0]);
    AO[counter+7]=AO[counter+7]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,4)-(4,-4)}
    evaluation[0]=pow(FIVE/FOUR,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],o[0],l[0],n[0])
                 -pow(FIVE/FOUR,HALF)*eval(Coord_At,Point,P_exp_send[j[0]],l[0],o[0],n[0]);
    AO[counter+8]=AO[counter+8]+C_coef_send[j[0]]*evaluation[0];
   }
   counter=counter+9;
   delete[] n;delete[] m; delete[] l; delete[] o; delete[] p;
   l=NULL;m=NULL;n=NULL;o=NULL;p=NULL;
  }
  else
  {}
  delete[] i;delete[] j;delete[] qdim;
  delete[] evaluation;
 i=NULL;j=NULL;qdim=NULL;evaluation=NULL;
}
//Build builds the AO_gradients. Recives the atomic coordinates, coeficients,
// exponents and the point where density must be evaluated it decides whether
// or not use SP and calls for a function(Quant_fill) to make the n,l,m.
void READ_FCHK_WFN::build_grad(int &iprim,int &styp,int &counter2,double *C_coef_send,
double *SPC_coef_send,double *P_exp_send,double Point[3],double **AO_grad,
double Coord_At[3],bool &activeSP)
{ int *i,*j,*qdim,*dir;
  int **Quant;
  double *evaluation;
  i=new int[1];j=new int[1];qdim=new int[1];dir=new int[1];
  evaluation=new double[1];
  if((styp==0 || styp==1) || (styp==2 || styp==3) || styp==4)
  {
   qdim[0]=(styp+1)*(styp+2)/2;
   Quant=new int*[qdim[0]];
   for(i[0]=0;i[0]<qdim[0];i[0]++)
   {Quant[i[0]]=new int[3];}
   //Quantum array n,l,m
   Quant_fill(Quant,styp);
   for(i[0]=0;i[0]<qdim[0];i[0]++)
   {for(j[0]=0;j[0]<iprim;j[0]++)
    {if(!activeSP)
     {
      for(dir[0]=0;dir[0]<3;dir[0]++)
      {//dir 0 x direction, dir 1 y direction and dir 2 z direction
      evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],Quant[i[0]][0],
      Quant[i[0]][1],Quant[i[0]][2]);
      AO_grad[dir[0]][counter2]=AO_grad[dir[0]][counter2]+C_coef_send[j[0]]*evaluation[0];
      }
     }
     else
     {
      for(dir[0]=0;dir[0]<3;dir[0]++)
      {
      evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],Quant[i[0]][0],
      Quant[i[0]][1],Quant[i[0]][2]);
      AO_grad[dir[0]][counter2]=AO_grad[dir[0]][counter2]+SPC_coef_send[j[0]]*evaluation[0];
      }
     }
    }
    counter2++;
   }
   if(activeSP){activeSP=false;}
   for(i[0]=0;i[0]<styp+1;i[0]++)
   {delete[] Quant[i[0]];Quant[i[0]]=NULL;}
   delete Quant;
   Quant=NULL;
  }
  else if(styp==-1)
  {
    int *send_type;
    send_type=new int[1];
    send_type[0]=0;
    build_grad(iprim,send_type[0],counter2,C_coef_send,SPC_coef_send,P_exp_send,Point,AO_grad,
    Coord_At,activeSP);
    activeSP=true;
    send_type[0]=1;
    build_grad(iprim,send_type[0],counter2,C_coef_send,SPC_coef_send,P_exp_send,Point,AO_grad,
    Coord_At,activeSP);
    delete[] send_type;
    send_type=NULL;
  }
  else if(styp==-2)
  {
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
   int *n,*l,*m;
   n=new int[1];l=new int[1];m=new int[1];
   n[0]=0;l[0]=1;m[0]=2;
   for(j[0]=0;j[0]<iprim;j[0]++)
   {
    for(dir[0]=0;dir[0]<3;dir[0]++)
    {
    //(2,0)
     evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],n[0],m[0])
               -HALF*(eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],n[0],n[0])
               +eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],m[0],n[0]));
     AO_grad[dir[0]][counter2]=AO_grad[dir[0]][counter2]+C_coef_send[j[0]]*evaluation[0];
    //2^(-1/2){(2,1)+(2,-1)}
     evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],n[0],l[0]);
     AO_grad[dir[0]][counter2+1]=AO_grad[dir[0]][counter2+1]+C_coef_send[j[0]]*evaluation[0];
    //(-2)^(-1/2){(2,1)+(2,-1)}
     evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],l[0],l[0]);
     AO_grad[dir[0]][counter2+2]=AO_grad[dir[0]][counter2+2]+C_coef_send[j[0]]*evaluation[0];
    //2^(-1/2){(2,2)+(2,-2)}
     evaluation[0]=pow(THREE/FOUR,HALF)*(eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],n[0],n[0])-
     eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],m[0],n[0]));
     AO_grad[dir[0]][counter2+3]=AO_grad[dir[0]][counter2+3]+C_coef_send[j[0]]*evaluation[0];
    //(-2)^(-1/2){(2,2)-(2,-2)}
     evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],l[0],n[0]);
     AO_grad[dir[0]][counter2+4]=AO_grad[dir[0]][counter2+4]+C_coef_send[j[0]]*evaluation[0];
    }
   }
   counter2=counter2+5;
   delete[] n;delete[] m; delete[] l;
   l=NULL;m=NULL;n=NULL;
  }
  else if(styp==-3)
  {
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
   int *n,*l,*m,*o;
   n=new int[1];l=new int[1];m=new int[1];o=new int[1];
   n[0]=0;l[0]=1;m[0]=2;o[0]=3;
   for(j[0]=0;j[0]<iprim;j[0]++)
   {
    for(dir[0]=0;dir[0]<3;dir[0]++)
    {
    //(3,0)
     evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],n[0],o[0])
               -(THREE/(TWO*pow(FIVE,HALF)))*(eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],n[0],l[0])
               +eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],m[0],l[0]));
     AO_grad[dir[0]][counter2]=AO_grad[dir[0]][counter2]+C_coef_send[j[0]]*evaluation[0];
    //2^(-1/2){(3,1)+(3,-1)}
     evaluation[0]=pow(SIX/FIVE,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],n[0],m[0])
               -(pow(SIX,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],o[0],n[0],n[0])
               -(pow(SIX/FIVE,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],m[0],n[0]);
     AO_grad[dir[0]][counter2+1]=AO_grad[dir[0]][counter2+1]+C_coef_send[j[0]]*evaluation[0];
    //(-2)^(-1/2){(3,1)-(3,-1)}
     evaluation[0]=pow(SIX/FIVE,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],l[0],m[0])
               -(pow(SIX,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],o[0],n[0])
               -(pow(SIX/FIVE,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],l[0],n[0]);
     AO_grad[dir[0]][counter2+2]=AO_grad[dir[0]][counter2+2]+C_coef_send[j[0]]*evaluation[0];
    //2^(-1/2){(3,2)+(3,-2)}
     evaluation[0]=pow(SIX/EIGHT,HALF)*(eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],n[0],l[0])
               -eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],m[0],l[0]));
     AO_grad[dir[0]][counter2+3]=AO_grad[dir[0]][counter2+3]+C_coef_send[j[0]]*evaluation[0];
    //(-2)^(-1/2){(3,2)-(3,-2)}
     evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],l[0],l[0]);
     AO_grad[dir[0]][counter2+4]=AO_grad[dir[0]][counter2+4]+C_coef_send[j[0]]*evaluation[0];
    //2^(-1/2){(3,3)+(3,-3)}
     evaluation[0]=(pow(TEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],o[0],n[0],n[0])
               -(THREE/(TWO*pow(TWO,HALF)))*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],m[0],n[0]);
     AO_grad[dir[0]][counter2+5]=AO_grad[dir[0]][counter2+5]+C_coef_send[j[0]]*evaluation[0];
    //(-2)^(-1/2){(3,3)-(3,-3)}
     evaluation[0]=-(pow(TEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],o[0],n[0])
               +(THREE/(TWO*pow(TWO,HALF)))*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],l[0],n[0]);
     AO_grad[dir[0]][counter2+6]=AO_grad[dir[0]][counter2+6]+C_coef_send[j[0]]*evaluation[0];
    }
   }
  counter2=counter2+7;
  delete[] n;delete[] m; delete[] l; delete[] o;
  l=NULL;m=NULL;n=NULL;o=NULL;
  }
  else if(styp==-4)
  {
   //See H. Bernhard Schlegel and Michael J. Frisch
   //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
   //Int. J. Quant. Chem., 54, 83-87 (1995).
   int *n,*l,*m,*o,*p;
   n=new int[1];l=new int[1];m=new int[1];o=new int[1];p=new int[1];
   n[0]=0;l[0]=1;m[0]=2;o[0]=3;p[0]=4;
   for(j[0]=0;j[0]<iprim;j[0]++)
   {
    for(dir[0]=0;dir[0]<3;dir[0]++)
    {
   //(4,0)
    evaluation[0]=eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],n[0],p[0])
              +(THREE/EIGHT)*(eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],p[0],n[0],n[0])
              +eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],p[0],n[0]))
              -(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*(eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],n[0],m[0])
              +eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],m[0],m[0])
              -HALF*HALF*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],m[0],n[0]));
    AO_grad[dir[0]][counter2]=AO_grad[dir[0]][counter2]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,1)+(4,-1)}
    evaluation[0]=(pow(TEN/SEVEN,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],n[0],o[0])
                  -(THREE*pow(TEN/SEVEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],o[0],n[0],l[0])
                  -(THREE*pow(TWO/SEVEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],m[0],l[0]));
    AO_grad[dir[0]][counter2+1]=AO_grad[dir[0]][counter2+1]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,1)-(4,-1)}
    evaluation[0]=(pow(TEN/SEVEN,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],l[0],o[0])
                  -(THREE*pow(TEN/SEVEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],o[0],l[0])
                  -(THREE*pow(TWO/SEVEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],l[0],l[0]));
    AO_grad[dir[0]][counter2+2]=AO_grad[dir[0]][counter2+2]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,2)+(4,-2)}
    evaluation[0]=(THREE*pow(THREE/SEVEN,HALF)/TWO)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],n[0],m[0])
                 -(THREE*pow(THREE/SEVEN,HALF)/TWO)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],m[0],m[0])
                 -HALF*HALF*pow(FIVE,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],p[0],n[0],n[0])
                 +HALF*HALF*pow(FIVE,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],p[0],n[0]);
    AO_grad[dir[0]][counter2+3]=AO_grad[dir[0]][counter2+3]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,2)-(4,-2)}
    evaluation[0]=(THREE/pow(SEVEN,HALF))*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],l[0],m[0])
                 -HALF*pow(FIVE/SEVEN,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],o[0],l[0],n[0])
                 -HALF*pow(FIVE/SEVEN,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],o[0],n[0]);
    AO_grad[dir[0]][counter2+4]=AO_grad[dir[0]][counter2+4]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,3)+(4,-3)}
    evaluation[0]=(pow(TEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],o[0],n[0],l[0])
                 -THREE*(pow(TWO,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],m[0],l[0]);
    AO_grad[dir[0]][counter2+5]=AO_grad[dir[0]][counter2+5]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,3)-(4,-3)}
    evaluation[0]=(-pow(TEN,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],o[0],l[0])
                 +THREE*(pow(TWO,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],l[0],l[0]);
    AO_grad[dir[0]][counter2+6]=AO_grad[dir[0]][counter2+6]+C_coef_send[j[0]]*evaluation[0];
   //2^(-1/2){(4,4)+(4,-4)}
    evaluation[0]=(pow(SEVEN*FIVE,HALF)/EIGHT)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],p[0],n[0],n[0])
                 +(pow(SEVEN*FIVE,HALF)/EIGHT)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],n[0],p[0],n[0])
                 -(THREE*pow(THREE,HALF)/FOUR)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],m[0],m[0],n[0]);
    AO_grad[dir[0]][counter2+7]=AO_grad[dir[0]][counter2+7]+C_coef_send[j[0]]*evaluation[0];
   //(-2)^(-1/2){(4,4)-(4,-4)}
    evaluation[0]=pow(FIVE/FOUR,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],o[0],l[0],n[0])
                 -pow(FIVE/FOUR,HALF)*eval_grad(dir[0],Coord_At,Point,P_exp_send[j[0]],l[0],o[0],n[0]);
    AO_grad[dir[0]][counter2+8]=AO_grad[dir[0]][counter2+8]+C_coef_send[j[0]]*evaluation[0];
    }
   }
   counter2=counter2+9;
   delete[] n;delete[] m; delete[] l; delete[] o; delete[] p;
   l=NULL;m=NULL;n=NULL;o=NULL;p=NULL;
  }
  else
  {}
  delete[] i;delete[] j;delete[] qdim;delete[] dir;
  delete[] evaluation;
  i=NULL; j=NULL; qdim=NULL; dir=NULL;
  evaluation=NULL;
}
//Evaluate point in each primitive calls dfact() to do the double factorial
double READ_FCHK_WFN::eval(double Coord_At[3],double Point[3],double &expon, int &n, int &l,int &m)
{
 double evaluation,*RR,*Pnorm,*norm,*Gx,*Gy,*Gz;
 int *n_df,*m_df,*l_df;
 RR=new double[1];Pnorm=new double[1];
 norm=new double[1];Gx=new double[1];Gy=new double[1];Gz=new double[1];
 n_df=new int[1];l_df=new int[1];m_df=new int[1];
 norm[0]=(double)n+l+m;
 RR[0]=pow((Point[0]-Coord_At[0]),TWO)+pow((Point[1]-Coord_At[1]),TWO)+
 pow((Point[2]-Coord_At[2]),TWO);
 n_df[0]=dfact(2*n-1);
 l_df[0]=dfact(2*l-1);
 m_df[0]=dfact(2*m-1);
 Pnorm[0]=pow((TWO*expon/PI),(THREE/FOUR))*pow((FOUR*expon),(norm[0]/TWO));
 Pnorm[0]=Pnorm[0]/pow(n_df[0]*m_df[0]*l_df[0],HALF);
 Gx[0]=pow((Point[0]-Coord_At[0]),n);
 Gy[0]=pow((Point[1]-Coord_At[1]),l);
 Gz[0]=pow((Point[2]-Coord_At[2]),m);
 evaluation=Pnorm[0]*Gx[0]*Gy[0]*Gz[0]*exp(-expon*RR[0]);
 delete[] RR;delete[] Pnorm;delete[] norm; delete[] Gx; delete[] Gy;
 delete[] Gz; delete[] n_df;delete[] m_df;delete[] l_df;
 RR=NULL;Pnorm=NULL;norm=NULL;Gx=NULL;Gy=NULL;Gz=NULL;n_df=NULL;
 m_df=NULL;l_df=NULL;
 return evaluation;
}
//Evaluate point in each primitive calls dfact() to do the double factorial
double READ_FCHK_WFN::eval_grad(int &dir, double Coord_At[3],double Point[3],double &expon,
int &n, int &l,int &m)
{ //dir 0 x axis, 1 y axis and 2 is z axis
 double evaluation,*RR,*Pnorm,*norm,*Gx,*Gy,*Gz,*dervG;
 double *nd,*ld,*md;
 RR=new double[1];Pnorm=new double[1];dervG=new double[1];
 norm=new double[1];Gx=new double[1];Gy=new double[1];Gz=new double[1];
 nd=new double[1];ld=new double[1];md=new double[1];
 nd[0]=(double)n;
 ld[0]=(double)l;
 md[0]=(double)m;
 int *n_df,*m_df,*l_df;
 n_df=new int[1];l_df=new int[1];m_df=new int[1];
 norm[0]=(double)n+l+m;
 RR[0]=pow((Point[0]-Coord_At[0]),TWO)+pow((Point[1]-Coord_At[1]),TWO)+
 pow((Point[2]-Coord_At[2]),TWO);
 n_df[0]=dfact(2*n-1);
 l_df[0]=dfact(2*l-1);
 m_df[0]=dfact(2*m-1);
 Pnorm[0]=pow((TWO*expon/PI),(THREE/FOUR))*pow((FOUR*expon),(norm[0]/TWO));
 Pnorm[0]=Pnorm[0]/pow(n_df[0]*m_df[0]*l_df[0],HALF);
 if(dir==0)
 {
 Gx[0]=pow((Point[0]-Coord_At[0]),nd[0]);
 dervG[0]=pow((Point[0]-Coord_At[0]),(nd[0]-ONE));
 if((Point[0]-Coord_At[0])==ZERO){dervG[0]=pow((1.0e-30),(nd[0]-ONE));}
 Gy[0]=pow((Point[1]-Coord_At[1]),ld[0]);
 Gz[0]=pow((Point[2]-Coord_At[2]),md[0]);
 evaluation=Pnorm[0]*dervG[0]*Gy[0]*Gz[0]*exp(-expon*RR[0])*nd[0];
 evaluation=evaluation-TWO*expon*Pnorm[0]*Gx[0]*Gy[0]*Gz[0]*(Point[0]-Coord_At[0])*exp(-expon*RR[0]);
 }
 else if(dir==1)
 {
 Gx[0]=pow((Point[0]-Coord_At[0]),nd[0]);
 Gy[0]=pow((Point[1]-Coord_At[1]),ld[0]);
 dervG[0]=pow((Point[1]-Coord_At[1]),ld[0]-ONE);
 if((Point[1]-Coord_At[1])==ZERO){dervG[0]=pow((1.0e-30),(ld[0]-ONE));}
 Gz[0]=pow((Point[2]-Coord_At[2]),md[0]);
 evaluation=Pnorm[0]*dervG[0]*Gx[0]*Gz[0]*exp(-expon*RR[0])*ld[0];
 evaluation=evaluation-TWO*expon*Pnorm[0]*Gx[0]*Gy[0]*Gz[0]*(Point[1]-Coord_At[1])*exp(-expon*RR[0]);
 }
 else
 {
 Gx[0]=pow((Point[0]-Coord_At[0]),nd[0]);
 Gy[0]=pow((Point[1]-Coord_At[1]),ld[0]);
 Gz[0]=pow((Point[2]-Coord_At[2]),md[0]);
 dervG[0]=pow((Point[2]-Coord_At[2]),md[0]-ONE);
 if((Point[2]-Coord_At[2])==ZERO){dervG[0]=pow((1.0e-30),(md[0]-ONE));}
 evaluation=Pnorm[0]*dervG[0]*Gx[0]*Gy[0]*exp(-expon*RR[0])*md[0];
 evaluation=evaluation-TWO*expon*Pnorm[0]*Gx[0]*Gy[0]*Gz[0]*(Point[2]-Coord_At[2])*exp(-expon*RR[0]);
 }
 delete[] RR;delete[] Pnorm;delete[] norm; delete[] Gx; delete[] Gy;
 delete[] Gz; delete[] n_df;delete[] m_df;delete[] l_df;delete[] dervG;
 delete[] nd; delete[] ld; delete[] md;
 RR=NULL;Pnorm=NULL;norm=NULL;Gx=NULL;Gy=NULL;Gz=NULL;n_df=NULL;
 m_df=NULL;l_df=NULL;dervG=NULL; nd=NULL;md=NULL;ld=NULL;
 return evaluation;
}
//evaluate point in each primitive for wfn type of file to build NOs
//Note: The coefficients in wfn already have the normalization coef of the gaussian.
double READ_FCHK_WFN::eval_Primitive_wfn(double Point[3],double pos_nuclei[3],int nlm[3],double exponent)
{
 int i;
 double RR,result,n,l,m;
 n=(double) nlm[0];
 l=(double) nlm[1];
 m=(double) nlm[2];
 RR=ZERO;
 for(i=0;i<3;i++)
 {
 RR=RR+(pow((Point[i]-pos_nuclei[i]),TWO));
 }
 result=pow((Point[0]-pos_nuclei[0]),n)*pow((Point[1]-pos_nuclei[1]),l);
 result=result*pow((Point[2]-pos_nuclei[2]),m)*exp(-exponent*RR);
 return result;
}
//evaluate point in each primitive for wfn type of file to build MOs
double READ_FCHK_WFN::grad_Primitive_wfn(double Point[3],int dir,double pos_nuclei[3],int nlm[3],double exponent)
{
 int i;
 double RR,evaluation,nd,ld,md,Gx,Gy,Gz,dervG;
 double expon,Coord_At[3];
 nd=(double) nlm[0];
 ld=(double) nlm[1];
 md=(double) nlm[2];
 expon=exponent;
 RR=ZERO;
 for(i=0;i<3;i++)
 {
  RR=RR+(pow((Point[i]-pos_nuclei[i]),TWO));
  Coord_At[i]=pos_nuclei[i];
 }
 if(dir==0)
 {
 Gx=pow((Point[0]-Coord_At[0]),nd);
 dervG=pow((Point[0]-Coord_At[0]),(nd-ONE));
 if((Point[0]-Coord_At[0])==ZERO){dervG=pow((1.0e-30),(nd-ONE));}
 Gy=pow((Point[1]-Coord_At[1]),ld);
 Gz=pow((Point[2]-Coord_At[2]),md);
 evaluation=dervG*Gy*Gz*exp(-expon*RR)*nd;
 evaluation=evaluation-TWO*expon*Gx*Gy*Gz*(Point[0]-Coord_At[0])*exp(-expon*RR);
 }
 else if(dir==1)
 {
 Gx=pow((Point[0]-Coord_At[0]),nd);
 Gy=pow((Point[1]-Coord_At[1]),ld);
 dervG=pow((Point[1]-Coord_At[1]),ld-ONE);
 if((Point[1]-Coord_At[1])==ZERO){dervG=pow((1.0e-30),(ld-ONE));}
 Gz=pow((Point[2]-Coord_At[2]),md);
 evaluation=dervG*Gx*Gz*exp(-expon*RR)*ld;
 evaluation=evaluation-TWO*expon*Gx*Gy*Gz*(Point[1]-Coord_At[1])*exp(-expon*RR);
 }
 else
 {
 Gx=pow((Point[0]-Coord_At[0]),nd);
 Gy=pow((Point[1]-Coord_At[1]),ld);
 Gz=pow((Point[2]-Coord_At[2]),md);
 dervG=pow((Point[2]-Coord_At[2]),md-ONE);
 if((Point[2]-Coord_At[2])==ZERO){dervG=pow((1.0e-30),(md-ONE));}
 evaluation=dervG*Gx*Gy*exp(-expon*RR)*md;
 evaluation=evaluation-TWO*expon*Gx*Gy*Gz*(Point[2]-Coord_At[2])*exp(-expon*RR);
 }
 return evaluation;
}


//////////////////////////////////
//////////////////////////////////
//      Momentum Space          //
//      ---------------         //
//  (Fourier transformations)   //
//////////////////////////////////
//////////////////////////////////

////////////////////////////////////////////////////////
//Functions used by the class for external evaluations//
////////////////////////////////////////////////////////
void READ_FCHK_WFN::rho_p_eval(double point_p[3],double &result)
{
 int *i,*j;
 i=new int[1];j=new int[1];
 if(!wfn)
 {
  double **AOp;
  AOp=new double*[2];
  for(i[0]=0;i[0]<2;i[0]++)
  {AOp[i[0]]=new double[nbasisf];}
  for(i[0]=0;i[0]<nbasisf;i[0]++)
  {AOp[0][i[0]]=ZERO;AOp[1][i[0]]=ZERO;}
  build_AOp(AOp,point_p);
  //Evaluate the density 
   result=ZERO;
   for(i[0]=0;i[0]<nbasisf;i[0]++)
   {
    for(j[0]=i[0];j[0]<nbasisf;j[0]++)
    {
     result=result+(AOp[0][i[0]]*AOp[0][j[0]]+AOp[1][i[0]]*AOp[1][j[0]])*Total_rho[j[0]][i[0]]; 
     if(i[0]!=j[0])
     {
     result=result+(AOp[0][j[0]]*AOp[0][i[0]]+AOp[1][j[0]]*AOp[1][i[0]])*Total_rho[j[0]][i[0]];
     }
    }
   }
  for(i[0]=0;i[0]<2;i[0]++)
  {delete[] AOp[i[0]];AOp[i[0]]=NULL;}
  delete[] AOp;
  AOp=NULL;
 }
 else //Build for wfn
 {
  complex<double> *NOp;
  NOp=new complex<double>[MO_coef];
  result=ZERO;
  density=ZERO;
  build_NOp_wfn_all(NOp,point_p);
  for(i[0]=0;i[0]<MO_coef;i[0]++)
  {
   if(Ocupation[i[0]]!=ZERO)
   {
    density=density+Ocupation[i[0]]*pow(abs(NOp[i[0]]),TWO);
   }
  }
  //Return result
  result=density;
  delete[] NOp;NOp=NULL;
 }
 delete[] i; delete[] j;
 i=NULL;j=NULL;
}
void READ_FCHK_WFN::rho_p_eval_a_b(double point_p[3],double &rhoa,double &rhob)
{
 int i,j;
 if(!wfn)
 {
  if(extra1)
  {
   double **AOp;
   AOp=new double*[2];
   for(i=0;i<2;i++)
   {AOp[i]=new double[nbasisf];}
   for(i=0;i<nbasisf;i++)
   {AOp[0][i]=ZERO;AOp[1][i]=ZERO;}
   build_AOp(AOp,point_p);
   //Evaluate the density
    rhoa=ZERO;
    rhob=ZERO;
    for(i=0;i<nbasisf;i++)
    {
     for(j=i;j<nbasisf;j++)
     {
      rhoa=rhoa+(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*(Total_rho[j][i]+Spin_rho[j][i])/TWO;
      rhob=rhob+(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*(Total_rho[j][i]-Spin_rho[j][i])/TWO;
      if(i!=j)
      {
       rhoa=rhoa+(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*(Total_rho[j][i]+Spin_rho[j][i])/TWO;
       rhob=rhob+(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*(Total_rho[j][i]-Spin_rho[j][i])/TWO;
      }
     }
    }
   for(i=0;i<2;i++)
   {delete[] AOp[i];AOp[i]=NULL;}
   delete[] AOp;
   AOp=NULL;
   }
   else
   {
    rho_p_eval(point_p,rhoa);
    rhoa=rhoa/TWO;
    rhob=rhoa;
   }
 }
 else
 {
  if(!open_shell)
  {
   rho_p_eval(point_p,rhoa);
   rhoa=rhoa/TWO;
   rhob=rhoa;
  }
  else
  {
   complex<double> *NOp;
   NOp=new complex<double>[MO_coef];
   rhoa=ZERO;
   rhob=ZERO;
   build_NOp_wfn_all(NOp,point_p);
   if(!correlated)
   {
    for(i=0;i<MO_coef;i++)
    {
     if(Ocupation[i]!=ZERO)
     {
      if(i%2==0)
      {
       rhoa=rhoa+Ocupation[i]*pow(abs(NOp[i]),TWO); //calculate rhoa
      }
      else
      {
       rhob=rhob+Ocupation[i]*pow(abs(NOp[i]),TWO); //calculate rhob
      }
     }
    }
   }
   else      //Only not calculable in Open Shell Correlated .wfn
   {
    error_opens_wfn=true;
    rhoa=ZERO;
    rhob=ZERO;
   }
   delete[] NOp;NOp=NULL;
  }
 }
}
void READ_FCHK_WFN::rho_p_grad(double point_p[3], double grad_p[3])
{
 int i,j,k;
 double *aux1,*aux2,*h;
 aux1=new double[3]; aux2=new double[3];h=new double[1];
 h[0]=pow(TEN,-FIVE);
 for(i=0;i<3;i++){grad_p[i]=ZERO;}
 if(!wfn)
 {
  double **AOp;
  AOp=new double*[2];
  for(i=0;i<2;i++)
  {AOp[i]=new double[nbasisf];}
  //Density for finite diferences(px+-h,py+-h,pz+-h)
  for(k=0;k<3;k++)
  {
   if(k==0)
   {
    aux1[0]=point_p[0]+h[0];aux1[1]=point_p[1];aux1[2]=point_p[2];
    aux2[0]=point_p[0]-h[0];aux2[1]=point_p[1];aux2[2]=point_p[2];
   }
   else if(k==1)
   {
    aux1[0]=point_p[0];aux1[1]=point_p[1]+h[0];aux1[2]=point_p[2];
    aux2[0]=point_p[0];aux2[1]=point_p[1]-h[0];aux2[2]=point_p[2];
   }
   else
   {
    aux1[0]=point_p[0];aux1[1]=point_p[1];aux1[2]=point_p[2]+h[0];
    aux2[0]=point_p[0];aux2[1]=point_p[1];aux2[2]=point_p[2]-h[0];
   }
   for(i=0;i<nbasisf;i++)
   {AOp[0][i]=ZERO;AOp[1][i]=ZERO;}
   build_AOp(AOp,aux1);
   //Evaluate the density
   for(i=0;i<nbasisf;i++)
   {
    for(j=i;j<nbasisf;j++)
    {
     grad_p[k]=grad_p[k]+(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*Total_rho[j][i];
     if(i!=j)
     {
     grad_p[k]=grad_p[k]+(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*Total_rho[j][i];
     }
    }
   }
   for(i=0;i<nbasisf;i++)
   {AOp[0][i]=ZERO;AOp[1][i]=ZERO;}
   build_AOp(AOp,aux2);
   //Evaluate the density
   for(i=0;i<nbasisf;i++)
   {
    for(j=i;j<nbasisf;j++)
    {
     grad_p[k]=grad_p[k]-(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*Total_rho[j][i];
     if(i!=j)
     {
     grad_p[k]=grad_p[k]-(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*Total_rho[j][i];
     }
    }
   }
  }
  for(k=0;k<3;k++){grad_p[k]=grad_p[k]/(TWO*h[0]);}
  for(i=0;i<2;i++)
  {delete[] AOp[i];AOp[i]=NULL;}
  delete[] AOp;
  AOp=NULL;
 }
 else //Build for wfn
 {
  complex<double> *NOp1;
  NOp1=new complex<double>[MO_coef];
  complex<double> *NOp2;
  NOp2=new complex<double>[MO_coef];
  for(k=0;k<3;k++)
  {
   if(k==0)
   {
    aux1[0]=point_p[0]+h[0];aux1[1]=point_p[1];aux1[2]=point_p[2];
    aux2[0]=point_p[0]-h[0];aux2[1]=point_p[1];aux2[2]=point_p[2];
   }
   else if(k==1)
   {
    aux1[0]=point_p[0];aux1[1]=point_p[1]+h[0];aux1[2]=point_p[2];
    aux2[0]=point_p[0];aux2[1]=point_p[1]-h[0];aux2[2]=point_p[2];
   }
   else
   {
    aux1[0]=point_p[0];aux1[1]=point_p[1];aux1[2]=point_p[2]+h[0];
    aux2[0]=point_p[0];aux2[1]=point_p[1];aux2[2]=point_p[2]-h[0];
   }
   build_NOp_wfn_all(NOp1,aux1);
   build_NOp_wfn_all(NOp2,aux2);
   for(i=0;i<MO_coef;i++)
   {
    if(Ocupation[i]!=ZERO)
    {
     grad_p[k]=grad_p[k]+Ocupation[i]*(pow(abs(NOp1[i]),TWO)-pow(abs(NOp2[i]),TWO));
    }
   }
  }
  for(k=0;k<3;k++){grad_p[k]=grad_p[k]/(TWO*h[0]);}
  delete[] NOp1; NOp1=NULL;
  delete[] NOp2; NOp2=NULL;
 }
 delete[] aux1; delete[] aux2; delete[] h;
 aux1=NULL;aux2=NULL;h=NULL;
}
void READ_FCHK_WFN::rho_p_grad_a_b(double point_p[3], double grad_p_a[3], double grad_p_b[3])
{
  int i,j,k;
  double aux1[3],aux2[3],h=pow(TEN,-FIVE);
  for(i=0;i<3;i++)
  {
   grad_p_a[i]=ZERO;
   grad_p_b[i]=ZERO;
  }
  if(!wfn)
  {
   if(extra1)
   {
    double **AOp;
    AOp=new double*[2];
    for(i=0;i<2;i++)
    {AOp[i]=new double[nbasisf];}
    //Density for finite diferences(px+-h,py+-h,pz+-h)
    for(k=0;k<3;k++)
    {
     if(k==0)
     {
      aux1[0]=point_p[0]+h;aux1[1]=point_p[1];aux1[2]=point_p[2];
      aux2[0]=point_p[0]-h;aux2[1]=point_p[1];aux2[2]=point_p[2];
     }
     else if(k==1)
     {
      aux1[0]=point_p[0];aux1[1]=point_p[1]+h;aux1[2]=point_p[2];
      aux2[0]=point_p[0];aux2[1]=point_p[1]-h;aux2[2]=point_p[2];
     }
     else
     {
      aux1[0]=point_p[0];aux1[1]=point_p[1];aux1[2]=point_p[2]+h;
      aux2[0]=point_p[0];aux2[1]=point_p[1];aux2[2]=point_p[2]-h;
     }
     for(i=0;i<nbasisf;i++)
     {AOp[0][i]=ZERO;AOp[1][i]=ZERO;}
     build_AOp(AOp,aux1);
     //Evaluate the density
     for(i=0;i<nbasisf;i++)
     {
      for(j=i;j<nbasisf;j++)
      {
       grad_p_a[k]=grad_p_a[k]+(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*(Total_rho[j][i]+Spin_rho[j][i])/TWO;
       grad_p_b[k]=grad_p_b[k]+(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*(Total_rho[j][i]-Spin_rho[j][i])/TWO;
       if(i!=j)
       {
       grad_p_a[k]=grad_p_a[k]+(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*(Total_rho[j][i]+Spin_rho[j][i])/TWO;
       grad_p_b[k]=grad_p_b[k]+(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*(Total_rho[j][i]-Spin_rho[j][i])/TWO;
       }
      }
     }
     for(i=0;i<nbasisf;i++)
     {AOp[0][i]=ZERO;AOp[1][i]=ZERO;}
     build_AOp(AOp,aux2);
     //Evaluate the density
     for(i=0;i<nbasisf;i++)
     {
      for(j=i;j<nbasisf;j++)
      {
       grad_p_a[k]=grad_p_a[k]-(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*(Total_rho[j][i]+Spin_rho[j][i])/TWO;
       grad_p_b[k]=grad_p_b[k]-(AOp[0][i]*AOp[0][j]+AOp[1][i]*AOp[1][j])*(Total_rho[j][i]-Spin_rho[j][i])/TWO;
       if(i!=j)
       {
       grad_p_a[k]=grad_p_a[k]-(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*(Total_rho[j][i]+Spin_rho[j][i])/TWO;
       grad_p_b[k]=grad_p_b[k]-(AOp[0][j]*AOp[0][i]+AOp[1][j]*AOp[1][i])*(Total_rho[j][i]-Spin_rho[j][i])/TWO;
       }
      }
     }
    }
    for(k=0;k<3;k++)
    {
     grad_p_a[k]=grad_p_a[k]/(TWO*h);
     grad_p_b[k]=grad_p_b[k]/(TWO*h);
    }
    for(i=0;i<2;i++)
    {delete[] AOp[i];AOp[i]=NULL;}
    delete[] AOp;
    AOp=NULL;
   }
   else
   {
    rho_p_grad(point_p,grad_p_a);
    for(k=0;k<3;k++)
    {
     grad_p_a[k]=grad_p_a[k]/TWO;
     grad_p_b[k]=grad_p_a[k];
    }
   }
  }
  else
  {
   for(i=0;i<3;i++)
   {
    grad_p_a[i]=ZERO;
    grad_p_b[i]=ZERO;
   }
   if(!open_shell)
   {
    rho_p_grad(point_p,grad_p_a);
    for(i=0;i<3;i++)
    {
     grad_p_a[i]=grad_p_a[i]/TWO;
     grad_p_b[i]=grad_p_a[i];
    }
   }
   else
   {
    if(!correlated)
    {
     complex<double> NOp;
     for(k=0;k<3;k++)
     {
      if(k==0)
      {
       aux1[0]=point_p[0]+h;aux1[1]=point_p[1];aux1[2]=point_p[2];
       aux2[0]=point_p[0]-h;aux2[1]=point_p[1];aux2[2]=point_p[2];
      }
      else if(k==1)
      {
       aux1[0]=point_p[0];aux1[1]=point_p[1]+h;aux1[2]=point_p[2];
       aux2[0]=point_p[0];aux2[1]=point_p[1]-h;aux2[2]=point_p[2];
      }
      else
      {
       aux1[0]=point_p[0];aux1[1]=point_p[1];aux1[2]=point_p[2]+h;
       aux2[0]=point_p[0];aux2[1]=point_p[1];aux2[2]=point_p[2]-h;
      }
      for(i=0;i<MO_coef;i++)
      {
       if(Ocupation[i]!=ZERO)
       {
        build_NOp_wfn(NOp,aux1,i);
        if(i%2==0)
        {grad_p_a[k]=grad_p_a[k]+Ocupation[i]*pow(abs(NOp),TWO);}
        else
        {grad_p_b[k]=grad_p_b[k]+Ocupation[i]*pow(abs(NOp),TWO);}
        build_NOp_wfn(NOp,aux2,i);
        if(i%2==0)
        {grad_p_a[k]=grad_p_a[k]-Ocupation[i]*pow(abs(NOp),TWO);}
        else
        {grad_p_b[k]=grad_p_b[k]-Ocupation[i]*pow(abs(NOp),TWO);}
       }
      }
     }
     for(k=0;k<3;k++)
     {
      grad_p_a[k]=grad_p_a[k]/(TWO*h);
      grad_p_b[k]=grad_p_b[k]/(TWO*h);
     }
    }
    else
    {}
   }
  }
}
void READ_FCHK_WFN::rho_p_lapl(double point_p[3], double &laplacian)
{
 int l,m;
 double h=pow(TEN,-FIVE),temp;
 double AUX1[3];
 laplacian=ZERO;
 for(l=0;l<3;l++)
 {for(m=0;m<3;m++){AUX1[m]=point_p[m];}
    AUX1[l]=AUX1[l]+h;
    rho_p_eval(AUX1,temp);
    laplacian=laplacian+temp;
    AUX1[l]=AUX1[l]-TWO*h;
    rho_p_eval(AUX1,temp);
    laplacian=laplacian+temp;
 }
 rho_p_eval(point_p,temp);
 laplacian=(laplacian-SIX*temp)/pow(h,TWO);
}
void READ_FCHK_WFN::rho_p_lapl_a_b(double point_p[3], double &laplacian_a,double &laplacian_b)
{
 int l,m;
 double h=pow(TEN,-FIVE),temp,temp2;
 double AUX1[3];
 laplacian_a=ZERO;
 laplacian_b=ZERO;
 for(l=0;l<3;l++)
 {for(m=0;m<3;m++){AUX1[m]=point_p[m];}
    AUX1[l]=AUX1[l]+h;
    rho_p_eval_a_b(AUX1,temp,temp2);
    laplacian_a=laplacian_a+temp;
    laplacian_b=laplacian_b+temp2;
    AUX1[l]=AUX1[l]-TWO*h;
    rho_p_eval_a_b(AUX1,temp,temp2);
    laplacian_a=laplacian_a+temp;
    laplacian_b=laplacian_b+temp2;
 }
 rho_p_eval_a_b(point_p,temp,temp2);
 laplacian_a=(laplacian_a-SIX*temp)/pow(h,TWO);
 laplacian_b=(laplacian_b-SIX*temp2)/pow(h,TWO);
}
void READ_FCHK_WFN::build_NOp_wfn(complex<double> &NOp,double Point_p[3],int &numMO)
{
 int *i,*j,*nlm, **Quant;
 i=new int [1];j=new int [1];nlm=new int[3];
 double *pos_nuclei,*exponent,*evaluation_prim_real,*evaluation_prim_imag;
 pos_nuclei=new double[3];exponent=new double[1];
 evaluation_prim_real=new double[1];evaluation_prim_imag=new double[1];
 ///////////////////////////////////////
 //dynamic arrays
 Quant=new int*[35];
 for(i[0]=0;i[0]<35;i[0]++)
 {Quant[i[0]]=new int[3];}
 ///////////////////////////////////////
 //Initialize
 Quant_fill(Quant,0); //0 because no type is needed
 NOp=ZERO*NOp;
 //////////////////////////////////////////
 //build NOps
 for(i[0]=0;i[0]<nprimitv;i[0]++)
 {
  for(j[0]=0;j[0]<3;j[0]++)
  {
   pos_nuclei[j[0]]=Cartesian_Coor[shell_map[i[0]]-1][j[0]];
   nlm[j[0]]=Quant[shell_type[i[0]]-1][j[0]];
  }
  exponent[0]=Prim_exp[i[0]];
  eval_Primitive_p_wfn(Point_p,pos_nuclei,nlm,exponent[0],evaluation_prim_real[0],evaluation_prim_imag[0]);
  complex<double>ztmp(ZERO,evaluation_prim_imag[0]);
  NOp=NOp+(MOcoefA[numMO][i[0]]*evaluation_prim_real[0]+MOcoefA[numMO][i[0]]*ztmp);
 }
 for(i[0]=0;i[0]<35;i[0]++)
 {delete[] Quant[i[0]];Quant[i[0]]=NULL;}
 delete[] Quant;
 delete[] i; delete [] j;delete[] nlm;
 delete[] pos_nuclei; delete [] evaluation_prim_real;
 delete[] evaluation_prim_imag;
 i=NULL;j=NULL;nlm=NULL;Quant=NULL;pos_nuclei=NULL;
 evaluation_prim_real=NULL; evaluation_prim_imag=NULL;
}
//Build all NOp 
void READ_FCHK_WFN::build_NOp_wfn_all(complex<double> *NOp,double Point_p[3])
{
 int i,j,k,nlm[3],**Quant;
 double pos_nuclei[3],exponent,evaluation_prim_real,evaluation_prim_imag;
 ///////////////////////////////////////
 //dynamic arrays
 Quant=new int*[35];
 for(i=0;i<35;i++)
 {Quant[i]=new int[3];}
 ///////////////////////////////////////
 //Initialize
 Quant_fill(Quant,0); //0 because no type is needed
 complex<double>ztmp0(ZERO,ZERO);
 for(k=0;k<MO_coef;k++)
 {
  NOp[k]=ztmp0;
 }
 //////////////////////////////////////////
 //build NOps
 for(i=0;i<nprimitv;i++)
 {
  for(j=0;j<3;j++)
  {
   pos_nuclei[j]=Cartesian_Coor[shell_map[i]-1][j];
   nlm[j]=Quant[shell_type[i]-1][j];
  }
  exponent=Prim_exp[i];
  eval_Primitive_p_wfn(Point_p,pos_nuclei,nlm,exponent,evaluation_prim_real,evaluation_prim_imag);
  complex<double>ztmp(ZERO,evaluation_prim_imag);
  for(k=0;k<MO_coef;k++)
  {
   complex<double>mocoef_var(MOcoefA[k][i],ZERO);
   if(im_wfn_wfx){complex<double>mocoef_var2(ZERO,MOcoefA_im[k][i]);mocoef_var=mocoef_var+mocoef_var2;}
   NOp[k]=NOp[k]+mocoef_var*(evaluation_prim_real+ztmp);
  }
 }
 for(i=0;i<35;i++)
 {delete[] Quant[i];Quant[i]=NULL;}
 delete[] Quant;Quant=NULL;
}
//Build MOp
void READ_FCHK_WFN::build_MOp_fchk(complex<double> &MOp,double Point_p[3],int &numMO)
{
  int *i,*j;
  i=new int[1];j=new int[1];
  double  **AOp;
  complex<double>ztmp0(ZERO,ZERO);
  //Point comes from outside
  MOp=ztmp0;
  AOp=new double*[2];
  for(i[0]=0;i[0]<2;i[0]++)
  {AOp[i[0]]=new double[nbasisf];}
  for(i[0]=0;i[0]<nbasisf;i[0]++)
  {
   AOp[0][i[0]]=ZERO;
   AOp[1][i[0]]=ZERO;
  }
  build_AOp(AOp,Point_p);
  //Build MOp
  if(extra1)
  {
   if(numMO%2==0)
   {
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    { 
     complex<double>ztmp(ZERO,AOp[1][j[0]]);
     MOp=MOp+(MOcoefA[numMO/2][j[0]])*(AOp[0][j[0]]+ztmp);
    }
   }
   else
   {
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     complex<double>ztmp(ZERO,AOp[1][j[0]]);
     if(BETA_MOS)
     {
      MOp=MOp+(MOcoefB[(numMO-1)/2][j[0]])*(AOp[0][j[0]]+ztmp);
     }
     else
     {
      MOp=MOp+(MOcoefA[(numMO-1)/2][j[0]])*(AOp[0][j[0]]+ztmp);
     }
    }
   }
  }
  else
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    complex<double>ztmp2(ZERO,AOp[1][j[0]]);
    MOp=MOp+(MOcoefA[numMO][j[0]])*(AOp[0][j[0]]+ztmp2);
   }
  }
///////////////////////////////////////////////////////////////////////////////
  for(i[0]=0;i[0]<2;i[0]++)
  {delete[] AOp[i[0]];AOp[i[0]]=NULL;}
  delete[] AOp; delete[] i; delete[] j;
  AOp=NULL;i=NULL;j=NULL;
}
//Build MOp taking AOps created from outside
void READ_FCHK_WFN::build_MOp_fchk2(double  **AOp,complex<double> &MOp,int &numMO)
{
  int *j;
  j=new int[1];
  complex<double>ztmp0(ZERO,ZERO);
  MOp=ztmp0;
  //Build MOp
  if(extra1)
  {
   if(numMO%2==0)
   {
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     complex<double> ztmp(AOp[0][j[0]],AOp[1][j[0]]);
     MOp=MOp+(MOcoefA[numMO/2][j[0]])*ztmp;
    }
   }
   else
   {
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     if(BETA_MOS)
     {
      complex<double> ztmp(AOp[0][j[0]],AOp[1][j[0]]);
      MOp=MOp+(MOcoefB[(numMO-1)/2][j[0]])*ztmp;
     }
     else
     {
      complex<double> ztmp(AOp[0][j[0]],AOp[1][j[0]]);
      MOp=MOp+(MOcoefA[(numMO-1)/2][j[0]])*ztmp;
     }
    }
   }
  }
  else
  {
   for(j[0]=0;j[0]<nbasisf;j[0]++)
   {
    complex<double> ztmp(AOp[0][j[0]],AOp[1][j[0]]);
    MOp=MOp+(MOcoefA[numMO][j[0]])*ztmp;
   }
  }
///////////////////////////////////////////////////////////////////////////////
 delete[] j;
 j=NULL;
}
//Calculate the gradient of the MOp (numerical gradient!)
void READ_FCHK_WFN::grad_MOp_fchk(complex<double> Grad[3],double Point_p[3],int &numMO)
{
 int *i,*j;
 double *h,*AUX1,*AUX2;
 complex<double> *eval1,*eval2;
 complex<double>ztmp0(ZERO,ZERO);
 i=new int[1];j=new int[1];
 h=new double[1]; AUX1=new double[3];
 AUX2=new double[3]; eval1=new complex<double>[1];
 eval2=new complex<double>[1];
 h[0]=pow(TEN,-SIX);
 for(i[0]=0;i[0]<3;i[0]++)
 {Grad[i[0]]=ztmp0;}
 for(i[0]=0;i[0]<3;i[0]++)
 {
  for(j[0]=0;j[0]<3;j[0]++)
  {
   AUX1[j[0]]=Point_p[j[0]];
   AUX2[j[0]]=Point_p[j[0]];
  }
  AUX1[i[0]]=AUX1[i[0]]+h[0];
  AUX2[i[0]]=AUX2[i[0]]-h[0];
  build_MOp_fchk(eval1[0],AUX1,numMO);
  build_MOp_fchk(eval2[0],AUX2,numMO);
  Grad[i[0]]=(eval1[0]-eval2[0])/(TWO*h[0]);
 }
 delete[] i; delete[] j;
 delete[] h; delete[] AUX1;
 delete[] AUX2;delete eval1;
 delete[] eval2;
 i=NULL;j=NULL;h=NULL;AUX1=NULL;
 AUX2=NULL;eval1=NULL;eval2=NULL;
}
//Build grad MOps using grads from outside
void READ_FCHK_WFN::grad_MOp_fchk2(complex<double> **AOps_Grad,complex<double> Grad[3],int &numMO)
{
  int *i,*j;
  i=new int[1];j=new int[1];
  complex<double>ztmp0(ZERO,ZERO);
  for(i[0]=0;i[0]<3;i[0]++)
  {
   Grad[i[0]]=ztmp0;
   //Build MOp
   if(extra1)
   {
    if(numMO%2==0)
    {
     for(j[0]=0;j[0]<nbasisf;j[0]++)
     {
      Grad[i[0]]=Grad[i[0]]+(MOcoefA[numMO/2][j[0]])*(AOps_Grad[i[0]][j[0]]);
     }
    }
    else
    {
     for(j[0]=0;j[0]<nbasisf;j[0]++)
     {
      Grad[i[0]]=Grad[i[0]]+(MOcoefB[(numMO-1)/2][j[0]])*(AOps_Grad[i[0]][j[0]]);
     }
    }
   }
   else
   {
    for(j[0]=0;j[0]<nbasisf;j[0]++)
    {
     Grad[i[0]]=Grad[i[0]]+(MOcoefA[numMO][j[0]])*(AOps_Grad[i[0]][j[0]]);
    }
   }
  }
///////////////////////////////////////////////////////////////////////////////
 delete[] i; delete[] j;
 i=NULL;j=NULL;
}
//Same as build AOps but may be use outside the class!
void READ_FCHK_WFN::build_AOp2(double **AOp,double point_p[3])
{
 int *i,*j,*k,*counter,*iprim,*iprim_tot;
 double *Coord_At, *P_exp_send, *C_coef_send, *SPC_coef_send;
 i=new int[1];j=new int[1];k=new int[1];counter=new int[1];iprim=new int[1];iprim_tot=new int[1];
 Coord_At=new double[3];
 counter[0]=0;
 for(i[0]=0;i[0]<nshells;i[0]++)
 {
  //save atomic coordinates of the shell, number of primitives
  //send only coefficients and exponents needed for that shell for building AOps.
  for(j[0]=0;j[0]<3;j[0]++)
  {Coord_At[j[0]]=Cartesian_Coor[(shell_map[i[0]]-1)][j[0]];}
  iprim[0]=n_prim_per_shell[i[0]];
  iprim_tot[0]=0;
  for(j[0]=0;j[0]<i[0];j[0]++)
  {iprim_tot[0]=iprim_tot[0]+n_prim_per_shell[j[0]];}
  C_coef_send=new double[iprim[0]];
  if(extra0)
  {SPC_coef_send=new double[iprim[0]];}
  P_exp_send=new double[iprim[0]];
  k[0]=0;
  for(j[0]=iprim_tot[0];j[0]<iprim[0]+iprim_tot[0];j[0]++)
  {
   C_coef_send[k[0]]=Contr_Coef[j[0]];
   if(extra0)
   {SPC_coef_send[k[0]]=SP_Contr_Coef[j[0]];}
   P_exp_send[k[0]]=Prim_exp[j[0]];
   k[0]++;
  }
  build_prim_AOp(iprim[0],shell_type[i[0]],counter[0],C_coef_send,SPC_coef_send,P_exp_send,point_p,
  AOp,Coord_At,activeSP);
  delete[] C_coef_send;
  C_coef_send=NULL;
  if(extra0)
  {delete[] SPC_coef_send;SPC_coef_send=NULL;}
  delete[] P_exp_send;
  P_exp_send=NULL;
 }
 delete[] i;delete[] j;delete[] k;delete[] counter;delete[] iprim;delete[] iprim_tot;
 delete[] Coord_At;
 i=NULL;j=NULL;k=NULL;counter=NULL; iprim=NULL; iprim_tot=NULL;
}
void READ_FCHK_WFN::build_AOp_grad2(complex<double> **AOps_Grad,double point_p[3])
{
 int *i,*j,*k;
 i=new int[1];j=new int[1];k=new int[1];
 //Now we build the Gradients Numerically
 //Someday analytically... =$
 double *h,point_step[3],point_step2[3];
 h=new double[1];
 double **AOps_aux,**AOps_aux2;
 h[0]=pow(TEN,-SIX);
 AOps_aux=new double*[2];
 AOps_aux2=new double*[2];
 for(i[0]=0;i[0]<2;i[0]++)
 {
  AOps_aux[i[0]]=new double[nbasisf];
  AOps_aux2[i[0]]=new double[nbasisf];
 }
 for(i[0]=0;i[0]<3;i[0]++)
 {
  for(j[0]=0;j[0]<3;j[0]++)
  {
   point_step[j[0]]=point_p[j[0]];
   point_step2[j[0]]=point_p[j[0]];
  }
  for(j[0]=0;j[0]<nbasisf;j[0]++)
  {
   AOps_aux[0][j[0]]=ZERO;AOps_aux[1][j[0]]=ZERO;
   AOps_aux2[0][j[0]]=ZERO;AOps_aux2[1][j[0]]=ZERO;
  }
  point_step[i[0]]=point_step[i[0]]+h[0];
  point_step2[i[0]]=point_step2[i[0]]-h[0];
  build_AOp2(AOps_aux,point_step);
  build_AOp2(AOps_aux2,point_step2);
  for(k[0]=0;k[0]<nbasisf;k[0]++)
  {
   AOps_Grad[i[0]][k[0]]=((AOps_aux[0][k[0]]-AOps_aux2[0][k[0]]),(AOps_aux[1][k[0]]-AOps_aux2[1][k[0]]))/(TWO*h[0]);
  }
 }
 for(i[0]=0;i[0]<2;i[0]++)
 {
  delete[] AOps_aux[i[0]];AOps_aux[i[0]]=NULL;
  delete[] AOps_aux2[i[0]];AOps_aux2[i[0]]=NULL;
 }
 delete[] AOps_aux;
 AOps_aux=NULL;
 delete[] AOps_aux2;
 AOps_aux2=NULL;
 delete[] i;delete[] j;delete[] k;
 i=NULL;j=NULL;k=NULL;
}
////////////////////////////////////////////////////////
//Functions used by the class for internal evaluations//
////////////////////////////////////////////////////////
void READ_FCHK_WFN::build_AOp(double **AOp,double point_p[3])
{
 int *i,*j,*k,*counter,*iprim,*iprim_tot;
 double *Coord_At, *P_exp_send, *C_coef_send, *SPC_coef_send;
 i=new int[1];j=new int[1];k=new int[1];counter=new int[1];iprim=new int[1];iprim_tot=new int[1];
 Coord_At=new double[3];
 counter[0]=0;
 for(i[0]=0;i[0]<nshells;i[0]++)
 {
  //save atomic coordinates of the shell, number of primitives
  //send only coefficients and exponents needed for that shell for building AOps.
  for(j[0]=0;j[0]<3;j[0]++)
  {Coord_At[j[0]]=Cartesian_Coor[(shell_map[i[0]]-1)][j[0]];}
  iprim[0]=n_prim_per_shell[i[0]];
  iprim_tot[0]=0;
  for(j[0]=0;j[0]<i[0];j[0]++)
  {iprim_tot[0]=iprim_tot[0]+n_prim_per_shell[j[0]];}
  C_coef_send=new double[iprim[0]];
  if(extra0)
  {SPC_coef_send=new double[iprim[0]];}
  P_exp_send=new double[iprim[0]];
  k[0]=0;
  for(j[0]=iprim_tot[0];j[0]<iprim[0]+iprim_tot[0];j[0]++)
  {
   C_coef_send[k[0]]=Contr_Coef[j[0]];
   if(extra0)
   {SPC_coef_send[k[0]]=SP_Contr_Coef[j[0]];}
   P_exp_send[k[0]]=Prim_exp[j[0]];
   k[0]++;
  }
  build_prim_AOp(iprim[0],shell_type[i[0]],counter[0],C_coef_send,SPC_coef_send,P_exp_send,point_p,
  AOp,Coord_At,activeSP);
  delete[] C_coef_send;
  C_coef_send=NULL;
  if(extra0)
  {delete[] SPC_coef_send;SPC_coef_send=NULL;}
  delete[] P_exp_send;
  P_exp_send=NULL;
 }
 delete[] i;delete[] j;delete[] k;delete[] counter;delete[] iprim;delete[] iprim_tot;
 delete[] Coord_At;
 i=NULL;j=NULL;k=NULL;counter=NULL; iprim=NULL; iprim_tot=NULL;
}
void READ_FCHK_WFN::build_AOp_grad(complex<double> **AOps_Grad,double point_p[3])
{
 int *i,*j,*k;
 i=new int[1];j=new int[1];k=new int[1];
 //Now we build the Gradients Numerically
 //Someday analytically... =$
 double *h,point_step[3],point_step2[3];
 h=new double[1];
 double **AOps_aux,**AOps_aux2;
 h[0]=pow(TEN,-SIX);
 AOps_aux=new double*[2];
 for(i[0]=0;i[0]<2;i[0]++)
 {AOps_aux[i[0]]=new double[nbasisf];}
 for(i[0]=0;i[0]<nbasisf;i[0]++)
 {AOps_aux[0][i[0]]=ZERO;AOps_aux[1][i[0]]=ZERO;}
 AOps_aux2=new double*[2];
 for(i[0]=0;i[0]<2;i[0]++)
 {AOps_aux2[i[0]]=new double[nbasisf];}
 for(i[0]=0;i[0]<nbasisf;i[0]++)
 {AOps_aux2[0][i[0]]=ZERO;AOps_aux2[1][i[0]]=ZERO;}
 for(i[0]=0;i[0]<3;i[0]++)
 {
  for(j[0]=0;j[0]<3;j[0]++)
  {
   point_step[j[0]]=point_p[j[0]];
   point_step2[j[0]]=point_p[j[0]];
  }
  point_step[i[0]]=point_step[i[0]]+h[0];
  point_step2[i[0]]=point_step2[i[0]]-h[0];
  build_AOp2(AOps_aux,point_step);
  build_AOp2(AOps_aux2,point_step2);
  for(k[0]=0;k[0]<nbasisf;k[0]++)
  {
   AOps_Grad[k[0]][i[0]]=((AOps_aux[0][k[0]]-AOps_aux2[0][k[0]]),(AOps_aux[1][k[0]]-AOps_aux2[1][k[0]]))/(TWO*h[0]);
  }
 }
 for(i[0]=0;i[0]<2;i[0]++)
 {delete[] AOps_aux[i[0]];AOps_aux[i[0]]=NULL;}
 delete[] AOps_aux;
 AOps_aux=NULL;
 for(i[0]=0;i[0]<2;i[0]++)
 {delete[] AOps_aux2[i[0]];AOps_aux2[i[0]]=NULL;}
 delete[] AOps_aux2;
 AOps_aux2=NULL;
 delete[] i;delete[] j;delete[] k;
 i=NULL;j=NULL;k=NULL;
}
void READ_FCHK_WFN::build_prim_AOp(int iprim,int styp,int &counter,double *C_coef_send,
double *SPC_coef_send,double *P_exp_send,double point_p[3],double **AOp,double Coord_At[3],
bool &activeSP)
{
 int *i,*j,*qdim;
 int **Quant;
 double *evaluationRE,*evaluationIM;
 i=new int[1];j=new int[1];qdim=new int[1];
 evaluationRE=new double [1]; evaluationIM=new double[1];
 if((styp==0 || styp==1) || (styp==2 || styp==3) || styp==4)
 {
  qdim[0]=(styp+1)*(styp+2)/2;
  Quant=new int*[qdim[0]];
  for(i[0]=0;i[0]<qdim[0];i[0]++)
  {Quant[i[0]]=new int[3];}
  //Quantum array n,l,m
  Quant_fill(Quant,styp);
  for(i[0]=0;i[0]<qdim[0];i[0]++)
  {for(j[0]=0;j[0]<iprim;j[0]++)
   {//eval contains the value of primitive at point_p
    eval_p(Coord_At,point_p,P_exp_send[j[0]],Quant[i[0]][0],
    Quant[i[0]][1],Quant[i[0]][2],evaluationRE[0],evaluationIM[0]);
    if(!activeSP)
    {
     AOp[0][counter]=AOp[0][counter]+C_coef_send[j[0]]*evaluationRE[0];
     AOp[1][counter]=AOp[1][counter]+C_coef_send[j[0]]*evaluationIM[0];
    }
    else
    {
     AOp[0][counter]=AOp[0][counter]+SPC_coef_send[j[0]]*evaluationRE[0];
     AOp[1][counter]=AOp[1][counter]+SPC_coef_send[j[0]]*evaluationIM[0];
    }
   }
   counter++;
  }
  if(activeSP){activeSP=false;}
  for(i[0]=0;i[0]<styp+1;i[0]++)
  {delete[] Quant[i[0]];Quant[i[0]]=NULL;}
  delete Quant;
  Quant=NULL;
 }
 else if(styp==-1)
 {
  build_prim_AOp(iprim,0,counter,C_coef_send,SPC_coef_send,P_exp_send,point_p,AOp,
  Coord_At,activeSP);
  activeSP=true;
  build_prim_AOp(iprim,1,counter,C_coef_send,SPC_coef_send,P_exp_send,point_p,AOp,
  Coord_At,activeSP);
 }
 else if(styp==-2)
 {
  //See H. Bernhard Schlegel and Michael J. Frisch
  //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
  //Int. J. Quant. Chem., 54, 83-87 (1995).
  int *n,*l,*m;
  n=new int[1];l=new int[1];m=new int[1];
  n[0]=0;l[0]=1;m[0]=2;
  for(j[0]=0;j[0]<iprim;j[0]++)
  {
   //(2,0)
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],n[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]+C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]-HALF*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]-HALF*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],m[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]-HALF*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]-HALF*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(2,1)+(2,-1)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],n[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]+C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(2,1)-(2,-1)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],l[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]+C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(2,2)+(2,-2)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]+pow(THREE/FOUR,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]+pow(THREE/FOUR,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],m[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]-C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]-C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(2,2)-(2,-2)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],l[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+4]=AOp[0][counter+4]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+4]=AOp[1][counter+4]+C_coef_send[j[0]]*evaluationIM[0];
  }
  counter=counter+5;
  delete[] n; delete[] l; delete[] m;
  n=NULL; l=NULL; m=NULL;
 }
 else if(styp==-3)
 {
  //See H. Bernhard Schlegel and Michael J. Frisch
  //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
  //Int. J. Quant. Chem., 54, 83-87 (1995).
  int *n,*l,*m,*o;
  n=new int[1];l=new int[1];m=new int[1];o=new int[1];
  n[0]=0;l[0]=1;m[0]=2;o[0]=3;
  for(j[0]=0;j[0]<iprim;j[0]++)
  {
   //(3,0)
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],n[0],o[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]+C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],n[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]-(THREE/(TWO*pow(FIVE,HALF)))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]-(THREE/(TWO*pow(FIVE,HALF)))*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],m[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]-(THREE/(TWO*pow(FIVE,HALF)))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]-(THREE/(TWO*pow(FIVE,HALF)))*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(3,1)+(3,-1)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],n[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]+pow(SIX/FIVE,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]+pow(SIX/FIVE,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],o[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]-(pow(SIX,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]-(pow(SIX,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],m[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]-(pow(SIX/FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]-(pow(SIX/FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(3,1)-(3,-1)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],l[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]+pow(SIX/FIVE,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]+pow(SIX/FIVE,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],o[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]-(pow(SIX,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]-(pow(SIX,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],l[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]-(pow(SIX/FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]-(pow(SIX/FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(3,2)+(3,-2)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],n[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]+pow(SIX/EIGHT,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]+pow(SIX/EIGHT,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],m[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]-pow(SIX/EIGHT,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]-pow(SIX/EIGHT,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(3,2)-(3,-2)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],l[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+4]=AOp[0][counter+4]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+4]=AOp[1][counter+4]+C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(3,3)+(3,-3)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],o[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+5]=AOp[0][counter+5]+(pow(TEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+5]=AOp[1][counter+5]+(pow(TEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],m[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+5]=AOp[0][counter+5]-(THREE/(TWO*pow(TWO,HALF)))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+5]=AOp[1][counter+5]-(THREE/(TWO*pow(TWO,HALF)))*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(3,3)-(3,-3)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],o[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+6]=AOp[0][counter+6]-(pow(TEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+6]=AOp[1][counter+6]-(pow(TEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],l[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+6]=AOp[0][counter+6]+(THREE/(TWO*pow(TWO,HALF)))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+6]=AOp[1][counter+6]+(THREE/(TWO*pow(TWO,HALF)))*C_coef_send[j[0]]*evaluationIM[0];
  }
  counter=counter+7;
  delete[] n; delete[] l; delete[] m; delete[] o;
  n=NULL; l=NULL; m=NULL; o=NULL;
 }
 else if(styp==-4)
 {
  //See H. Bernhard Schlegel and Michael J. Frisch
  //"Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
  //Int. J. Quant. Chem., 54, 83-87 (1995).
  int *n,*l,*m,*o,*p;
  n=new int[1];l=new int[1];m=new int[1];o=new int[1];p=new int[1];
  n[0]=0;l[0]=1;m[0]=2;o[0]=3;p[0]=4;
  for(j[0]=0;j[0]<iprim;j[0]++)
  {
   //(4,0)
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],n[0],p[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]+C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]+C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],p[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]+(THREE/EIGHT)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]+(THREE/EIGHT)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],p[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]+(THREE/EIGHT)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]+(THREE/EIGHT)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],n[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]-(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]-(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],m[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]-(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]-(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],m[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter]=AOp[0][counter]+(ONE/FOUR)*(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter]=AOp[1][counter]+(ONE/FOUR)*(THREE*pow(THREE,HALF)/pow(SEVEN*FIVE,HALF))*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(4,1)+(4,-1)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],n[0],o[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]+pow(TEN/SEVEN,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]+pow(TEN/SEVEN,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],o[0],n[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]-(THREE*pow(TEN/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]-(THREE*pow(TEN/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],m[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+1]=AOp[0][counter+1]-(THREE*pow(TWO/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+1]=AOp[1][counter+1]-(THREE*pow(TWO/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(4,1)-(4,-1)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],l[0],o[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]+pow(TEN/SEVEN,HALF)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]+pow(TEN/SEVEN,HALF)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],o[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]-(THREE*pow(TEN/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]-(THREE*pow(TEN/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],l[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+2]=AOp[0][counter+2]-(THREE*pow(TWO/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+2]=AOp[1][counter+2]-(THREE*pow(TWO/SEVEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(4,2)+(4,-2)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],n[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]+(THREE*pow(SIX/(TWO*SEVEN),HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]+(THREE*pow(SIX/(TWO*SEVEN),HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],m[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]-(THREE*pow(SIX/(TWO*SEVEN),HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]-(THREE*pow(SIX/(TWO*SEVEN),HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],p[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]-(pow(FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]-(pow(FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],p[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+3]=AOp[0][counter+3]+(pow(FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+3]=AOp[1][counter+3]+(pow(FIVE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(4,2)-(4,-2)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],l[0],m[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+4]=AOp[0][counter+4]+(THREE/pow(SEVEN,HALF))*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+4]=AOp[1][counter+4]+(THREE/pow(SEVEN,HALF))*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],o[0],l[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+4]=AOp[0][counter+4]-(pow(FIVE/SEVEN,HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+4]=AOp[1][counter+4]-(pow(FIVE/SEVEN,HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],o[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+4]=AOp[0][counter+4]-(pow(FIVE/SEVEN,HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+4]=AOp[1][counter+4]-(pow(FIVE/SEVEN,HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(4,3)+(4,-3)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],o[0],n[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+5]=AOp[0][counter+5]+(pow(TEN,HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+5]=AOp[1][counter+5]+(pow(TEN,HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],m[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+5]=AOp[0][counter+5]-(THREE*pow(TWO,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+5]=AOp[1][counter+5]-(THREE*pow(TWO,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(4,3)-(4,-3)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],o[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+6]=AOp[0][counter+6]-(pow(TEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+6]=AOp[1][counter+6]-(pow(TEN,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],l[0],l[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+6]=AOp[0][counter+6]+(THREE*pow(TWO,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+6]=AOp[1][counter+6]+(THREE*pow(TWO,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //2^(-1/2){(4,4)+(4,-4)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],p[0],n[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+7]=AOp[0][counter+7]+(pow(SEVEN*FIVE,HALF)/EIGHT)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+7]=AOp[1][counter+7]+(pow(SEVEN*FIVE,HALF)/EIGHT)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],n[0],p[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+7]=AOp[0][counter+7]+(pow(SEVEN*FIVE,HALF)/EIGHT)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+7]=AOp[1][counter+7]+(pow(SEVEN*FIVE,HALF)/EIGHT)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],m[0],m[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+7]=AOp[0][counter+7]-(THREE*pow(THREE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+7]=AOp[1][counter+7]-(THREE*pow(THREE,HALF)/FOUR)*C_coef_send[j[0]]*evaluationIM[0];
   //(-2)^(-1/2){(4,4)-(4,-4)}
   eval_p(Coord_At,point_p,P_exp_send[j[0]],o[0],l[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+8]=AOp[0][counter+8]+(pow(FIVE,HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+8]=AOp[1][counter+8]+(pow(FIVE,HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
   eval_p(Coord_At,point_p,P_exp_send[j[0]],l[0],o[0],n[0],evaluationRE[0],evaluationIM[0]);
   AOp[0][counter+8]=AOp[0][counter+8]-(pow(FIVE,HALF)/TWO)*C_coef_send[j[0]]*evaluationRE[0];
   AOp[1][counter+8]=AOp[1][counter+8]-(pow(FIVE,HALF)/TWO)*C_coef_send[j[0]]*evaluationIM[0];
  }
  counter=counter+9;
  delete[] n; delete[] l; delete[] m; delete[] o; delete[] p;
  n=NULL; l=NULL; m=NULL; o=NULL;p=NULL;
 }
 else
 {}
 delete[] i; delete[] j; delete[] qdim;
 delete[] evaluationRE; delete[] evaluationIM;
 i=NULL;j=NULL;qdim=NULL;
 evaluationRE=NULL;evaluationIM=NULL;
}
void READ_FCHK_WFN::eval_p(double Coord_At[3],double point_p[3],double &expon, int &n, int &l,int &m,double &re,double &im)
{
 int *i,*n_df,*m_df,*l_df,*nn;
 double *PP,*Pnorm,*norm;
 complex<double> *z,*z1,*z2,*z3,*fpx,*fpy,*fpz;
 double *coord_at,*point;
 i=new int[1];n_df=new int[1];m_df=new int[1];l_df=new int[1];nn=new int[1];
 PP=new double[1];Pnorm=new double[1];norm=new double[1];
 z=new complex<double>[1];z1=new complex<double>[1];z2=new complex<double>[1];
 z3=new complex<double>[1];fpx=new complex<double>[1];fpy=new complex<double>[1];
 fpz=new complex<double>[1];
 norm[0]=(double)(n+l+m);
 PP[0]=pow(norm3D(point_p),TWO);
 //Normalization of gaussian primitive
 n_df[0]=dfact(2*n-1);
 l_df[0]=dfact(2*l-1);
 m_df[0]=dfact(2*m-1);
 Pnorm[0]=pow((TWO*expon/PI),(THREE/FOUR))*pow((FOUR*expon),(norm[0]/TWO));
 Pnorm[0]=Pnorm[0]/pow(n_df[0]*m_df[0]*l_df[0],HALF);
 //f(ps) see the notes. Is the function that is multiplied by the 3D gaussian and Euler term
 fps(fpx[0],point_p[0],n,expon);
 fps(fpy[0],point_p[1],l,expon);
 fps(fpz[0],point_p[2],m,expon);
 //Euler term
 nn[0]=3;
 coord_at=new double[nn[0]];
 point=new double[nn[0]];
 for(i[0]=0;i[0]<3;i[0]++)
 {
  coord_at[i[0]]=Coord_At[i[0]];
  point[i[0]]=point_p[i[0]];
 }
 complex<double> ztmp(cos(dot(nn[0],coord_at,point)),sin(dot(nn[0],coord_at,point)));
 z1[0]=ztmp;
 //Multiplication
 z2[0]=Pnorm[0];
 z3[0]=exp(-PP[0]/(FOUR*expon));
 z[0]=z1[0]*z2[0]*z3[0]*fpx[0]*fpy[0]*fpz[0];
 //Send only RE and IM
 re=real(z[0]);
 im=imag(z[0]);
 delete[] i;delete[] n_df;delete[] m_df;delete[] l_df;delete[] nn;
 delete[] PP;delete[] Pnorm;delete[] norm;
 delete[] coord_at; delete[] point;
 delete[] z;delete[] z1;delete[] z2;delete[] z3;delete[] fpx;delete[] fpy;delete[] fpz;
 coord_at=NULL;
 point=NULL;
}
//evaluate point_p in each primitive for wfn type of file and build NOps
void READ_FCHK_WFN::eval_Primitive_p_wfn(double point_p[3],double pos_nuclei[3],int nlm[3],double &expon,double &re,double &im)
{
 int *i,*nn;
 double *PP;
 complex<double> *z,*z1,*z2,*fpx,*fpy,*fpz;
 double *coord_at,*point;
 i=new int[1];nn=new int[1];
 PP=new double[1];
 z=new complex<double>[1];z1=new complex<double>[1];z2=new complex<double>[1];
 fpx=new complex<double>[1];fpy=new complex<double>[1];
 fpz=new complex<double>[1];
 PP[0]=pow(norm3D(point_p),TWO);
 //f(ps) see the notes. Is the function that is multiplied by the 3D gaussian and Euler term
 fps(fpx[0],point_p[0],nlm[0],expon);
 fps(fpy[0],point_p[1],nlm[1],expon);
 fps(fpz[0],point_p[2],nlm[2],expon);
 //Euler term
 nn[0]=3;
 coord_at=new double[nn[0]];
 point=new double[nn[0]];
 for(i[0]=0;i[0]<3;i[0]++)
 {
  coord_at[i[0]]=pos_nuclei[i[0]];
  point[i[0]]=point_p[i[0]];
 }
 complex<double>ztmp(cos(dot(nn[0],coord_at,point)),sin(dot(nn[0],coord_at,point)));
 z1[0]=ztmp;
 //Multiplication
 z2[0]=exp(-PP[0]/(FOUR*expon));
 z[0]=z1[0]*z2[0]*fpx[0]*fpy[0]*fpz[0];
 //Send only RE and IM
 re=real(z[0]);
 im=imag(z[0]);
 delete[] coord_at;delete[] point;
 delete[] i; delete[] nn; delete[] z;
 delete[] z1; delete[] z2; delete[] fpx;
 delete[] fpy; delete[] fpz;delete[] PP;
 coord_at=NULL;point=NULL;
 i=NULL; nn=NULL;z=NULL;
 z1=NULL; z2=NULL; fpx=NULL;
 fpy=NULL; fpz=NULL;PP=NULL;
}
//For f(ps) see the notes about Momentum Space. Is the function that is multiplied by the 3D gaussian and Euler term
//used by fchk and wfn
void READ_FCHK_WFN::fps(complex<double> &f,double &p,int &quant_num,double &expon)
{
 if(quant_num==0)
 {
  complex<double> f1((ONE/pow(TWO*expon,HALF)),ZERO);
  f=f1;
 }
 else if(quant_num==1)
 { 
  complex<double> f1(ZERO,(p/(TWO*pow(TWO,HALF)*pow(expon,THREE/TWO))));
  f=f1;
 }
 else if(quant_num==2)
 {
  complex<double> f1((TWO*expon-pow(p,TWO))/(FOUR*pow(TWO,HALF)*pow(expon,FIVE/TWO)),ZERO);
  f=f1;
 }
 else if(quant_num==3)
 {
  complex<double> f1(ZERO,((SIX*p/pow(expon,FIVE/TWO)-pow(p,THREE)/pow(expon,SEVEN/TWO))/(EIGHT*pow(TWO,HALF))));
  f=f1;
 }
 else if(quant_num==4)
 {
  complex<double> f1(((TWO*SIX*pow(expon,TWO)-TWO*SIX*pow(p,TWO)*expon+pow(p,FOUR))/(pow(expon,NINE/TWO)*TWO*EIGHT*pow(TWO,HALF))),ZERO);
  f=f1;
 }
 else
 {cout<<"Warning! Basis set not available!"<<endl;}
}

////////////////////////////////////////
////////////////////////////////////////
//Quantum numbers for cart. gaussians //
//----------------------------------- //
////////////////////////////////////////
////////////////////////////////////////

//Generate n, l and m
void READ_FCHK_WFN::Quant_fill(int **Quant,int styp)
{
 if(!wfn)
 {
  int i,j,k,qfil,permute[3];
  qfil=0;
  for(i=styp;i>-1;i=i-1)
  {for(j=i;j>-1;j=j-1)
    {for(k=j;k>-1;k=k-1)
     {
      if((i+j+k)==styp)
      {
        Quant[qfil][0]=i;
        Quant[qfil][1]=j;
        Quant[qfil][2]=k;
        qfil++;
        if(j!=k)
        {
         Quant[qfil][0]=i;
         Quant[qfil][1]=k;
         Quant[qfil][2]=j;
         qfil++;
        }
        if(i!=j)
        {
         Quant[qfil][0]=j;
         Quant[qfil][1]=i;
         Quant[qfil][2]=k;
         qfil++;
        }
        if((i!=j) && (j!=k))
        {
         Quant[qfil][0]=j;
         Quant[qfil][1]=k;
         Quant[qfil][2]=i;
         qfil++;
         Quant[qfil][0]=k;
         Quant[qfil][1]=i;
         Quant[qfil][2]=j;
         qfil++;
        }
        if(k!=i)
        {
         Quant[qfil][0]=k;
         Quant[qfil][1]=j;
         Quant[qfil][2]=i;
         qfil++;
        }
      }
    }
   }
  }
  if(styp==3)
  {
   permute[0]=Quant[5][0];permute[1]=Quant[5][1];permute[2]=Quant[5][2];
   Quant[5][0]=Quant[4][0];Quant[5][1]=Quant[4][1];Quant[5][2]=Quant[4][2];
   Quant[4][0]=Quant[3][0];Quant[4][1]=Quant[3][1];Quant[4][2]=Quant[3][2];
   Quant[3][0]=permute[0];Quant[3][1]=permute[1];Quant[3][2]=permute[2];
   permute[0]=Quant[7][0];permute[1]=Quant[7][1];permute[2]=Quant[7][2];
   Quant[7][0]=Quant[8][0];Quant[7][1]=Quant[8][1];Quant[7][2]=Quant[8][2];
   Quant[8][0]=permute[0];Quant[8][1]=permute[1];Quant[8][2]=permute[2];
  }
  if(styp==4)
  {
   permute[0]=Quant[6][0];permute[1]=Quant[6][1];permute[2]=Quant[6][2];
   Quant[6][0]=Quant[7][0];Quant[6][1]=Quant[7][1];Quant[6][2]=Quant[7][2];
   Quant[7][0]=permute[0];Quant[7][1]=permute[1];Quant[7][2]=permute[2];
  }
 }
 else
 {
  Quant[0][0]=0;Quant[0][1]=0;Quant[0][2]=0;

  Quant[1][0]=1;Quant[1][1]=0;Quant[1][2]=0;
  Quant[2][0]=0;Quant[2][1]=1;Quant[2][2]=0;
  Quant[3][0]=0;Quant[3][1]=0;Quant[3][2]=1;

  Quant[4][0]=2;Quant[4][1]=0;Quant[4][2]=0;
  Quant[5][0]=0;Quant[5][1]=2;Quant[5][2]=0;
  Quant[6][0]=0;Quant[6][1]=0;Quant[6][2]=2;
  Quant[7][0]=1;Quant[7][1]=1;Quant[7][2]=0;
  Quant[8][0]=1;Quant[8][1]=0;Quant[8][2]=1;
  Quant[9][0]=0;Quant[9][1]=1;Quant[9][2]=1;

  Quant[10][0]=3;Quant[10][1]=0;Quant[10][2]=0;
  Quant[11][0]=0;Quant[11][1]=3;Quant[11][2]=0;
  Quant[12][0]=0;Quant[12][1]=0;Quant[12][2]=3;
  Quant[13][0]=2;Quant[13][1]=1;Quant[13][2]=0;
  Quant[14][0]=2;Quant[14][1]=0;Quant[14][2]=1;
  Quant[15][0]=0;Quant[15][1]=2;Quant[15][2]=1;
  Quant[16][0]=1;Quant[16][1]=2;Quant[16][2]=0;
  Quant[17][0]=1;Quant[17][1]=0;Quant[17][2]=2;
  Quant[18][0]=0;Quant[18][1]=1;Quant[18][2]=2;
  Quant[19][0]=1;Quant[19][1]=1;Quant[19][2]=1;

  Quant[20][0]=4;Quant[20][1]=0;Quant[20][2]=0;
  Quant[21][0]=0;Quant[21][1]=4;Quant[21][2]=0;
  Quant[22][0]=0;Quant[22][1]=0;Quant[22][2]=4;
  Quant[23][0]=3;Quant[23][1]=1;Quant[23][2]=0;
  Quant[24][0]=3;Quant[24][1]=0;Quant[24][2]=1;
  Quant[25][0]=1;Quant[25][1]=3;Quant[25][2]=0;
  Quant[26][0]=0;Quant[26][1]=3;Quant[26][2]=1;
  Quant[27][0]=1;Quant[27][1]=0;Quant[27][2]=3;
  Quant[28][0]=0;Quant[28][1]=1;Quant[28][2]=3;
  Quant[29][0]=2;Quant[29][1]=2;Quant[29][2]=0;
  Quant[30][0]=2;Quant[30][1]=0;Quant[30][2]=2;
  Quant[31][0]=0;Quant[31][1]=2;Quant[31][2]=2;
  Quant[32][0]=2;Quant[32][1]=1;Quant[32][2]=1;
  Quant[33][0]=1;Quant[33][1]=2;Quant[33][2]=1;
  Quant[34][0]=1;Quant[34][1]=1;Quant[34][2]=2;
 }
}

//Store im part of wfn coefs
void READ_FCHK_WFN::Init_MOim(double **MOcoef_in)
{
 int i,j;
 MOcoefA_im=new double*[MO_coef];
 for(i=0;i<MO_coef;i++)
 {
  MOcoefA_im[i]=new double[nprimitv];
  for(j=0;j<nprimitv;j++)
  {
   MOcoefA_im[i][j]=MOcoef_in[i][j];
  }
 }
}

//Move the origin to the center of mass
void READ_FCHK_WFN::Center_of_mass()
{
 int i,j,k,pivot,order[3];
 double **Im,**EigenV,*mass,Rcm[3]={ZERO},Mtot=ZERO,Norm,Norm_cm,Auxiliar,Auxiliar2;
 ofstream results_inert;
 if((name_file[name_file.length()-1]=='n' || name_file[name_file.length()-1]=='N')||(name_file[name_file.length()-1]=='x' || name_file[name_file.length()-1]=='X'))
 {
  results_inert.open((name_file.substr(0,(name_file.length()-4))+"_inert.tmp").c_str());
  Im=new double*[3];
  EigenV=new double*[3];
  mass=new double[natoms];
  Auxiliar2=ZERO;
  Auxiliar2=Auxiliar2;
  k=0;
  k=k;
  for(i=0;i<3;i++)
  {
   Im[i]=new double[3];
   EigenV[i]=new double[3];
   for(j=0;j<3;j++)
   {
    Im[i][j]=ZERO;
    EigenV[i][j]=ZERO;
   }
  }
  for(i=0;i<natoms;i++)
  {
   mass[i]=ZERO;
   //Silly if to check that the initialization was correct!
 if(mass[i]==ZERO)
 {
  if(Nu_charge[i]==ZERO)
  {
   mass[i]=ZERO;
  }
  else if(Nu_charge[i]==ONE)
  {
   mass[i]=1.00079;
  }
  else if(Nu_charge[i]==TWO)
  {
   mass[i]=4.0026;
  }
  else if(Nu_charge[i]==THREE)
  {
   mass[i]=6.941;
  }
  else if(Nu_charge[i]==FOUR)
  {
   mass[i]=9.0122;
  }
  else if(Nu_charge[i]==FIVE)
  {
   mass[i]=10.811;
  }
  else if(Nu_charge[i]==SIX)
  {
   mass[i]=12.011;
  }
  else if(Nu_charge[i]==SEVEN)
  {
   mass[i]=14.007;
  }
  else if(Nu_charge[i]==EIGHT)
  {
   mass[i]=15.999;
  }
  else if(Nu_charge[i]==NINE)
  {
   mass[i]=18.998;
  }
  else if(Nu_charge[i]==TEN)
  {
   mass[i]=20.180;
  }
  else if(Nu_charge[i]==(TEN+ONE))
  {
   mass[i]=22.990;
  }
  else if(Nu_charge[i]==(TEN+TWO))
  {
   mass[i]=24.305;
  }
  else if(Nu_charge[i]==(TEN+THREE))
  {
   mass[i]=26.982;
  }
  else if(Nu_charge[i]==(TEN+FOUR))
  {
   mass[i]=28.086;
  }
  else if(Nu_charge[i]==(TEN+FIVE))
  {
   mass[i]=30.974;
  }
  else if(Nu_charge[i]==(TEN+SIX))
  {
   mass[i]=32.065;
  }
  else if(Nu_charge[i]==(TEN+SEVEN))
  {
   mass[i]=35.453;
  }
  else if(Nu_charge[i]==(TEN+EIGHT))
  {
   mass[i]=39.948;
  }
  else if(Nu_charge[i]==(TEN+NINE))
  {
   mass[i]=39.098;
  }
  else if(Nu_charge[i]==(TWO*TEN))
  {
   mass[i]=40.078;
  }
  else if(Nu_charge[i]==(TWO*TEN+ONE))
  {
   mass[i]=44.956;
  }
  else if(Nu_charge[i]==(TWO*TEN+TWO))
  {
   mass[i]=47.867;
  }
  else if(Nu_charge[i]==(TWO*TEN+THREE))
  {
   mass[i]=50.942;
  }
  else if(Nu_charge[i]==(TWO*TEN+FOUR))
  {
   mass[i]=51.996;
  }
  else if(Nu_charge[i]==(TWO*TEN+FIVE))
  {
   mass[i]=54.938;
  }
  else if(Nu_charge[i]==(TWO*TEN+SIX))
  {
   mass[i]=55.845;
  }
  else if(Nu_charge[i]==(TWO*TEN+SEVEN))
  {
   mass[i]=58.933;
  }
  else if(Nu_charge[i]==(TWO*TEN+EIGHT))
  {
   mass[i]=58.693;
  }
  else if(Nu_charge[i]==(TWO*TEN+NINE))
  {
   mass[i]=63.546;
  }
  else if(Nu_charge[i]==(THREE*TEN))
  {
   mass[i]=65.39;
  }
  else if(Nu_charge[i]==(THREE*TEN+ONE))
  {
   mass[i]=69.723;
  }
  else if(Nu_charge[i]==(THREE*TEN+TWO))
  {
   mass[i]=72.61;
  }
  else if(Nu_charge[i]==(THREE*TEN+THREE))
  {
   mass[i]=74.922;
  }
  else if(Nu_charge[i]==(THREE*TEN+FOUR))
  {
   mass[i]=78.96;
  }
  else if(Nu_charge[i]==(THREE*TEN+FIVE))
  {
   mass[i]=79.904;
  }
  else if(Nu_charge[i]==(THREE*TEN+SIX))
  {
   mass[i]=83.80;
  }
  else if(Nu_charge[i]==(THREE*TEN+SEVEN))
  {
   mass[i]=85.468;
  }
  else if(Nu_charge[i]==(THREE*TEN+EIGHT))
  {
   mass[i]=87.62;
  }
  else if(Nu_charge[i]==(THREE*TEN+NINE))
  {
   mass[i]=89.906;
  }
  else if(Nu_charge[i]==(FOUR*TEN))
  {
   mass[i]=91.224;
  }
  else if(Nu_charge[i]==(FOUR*TEN+ONE))
  {
   mass[i]=92.906;
  }
  else if(Nu_charge[i]==(FOUR*TEN+TWO))
  {
   mass[i]=95.94;
  }
  else if(Nu_charge[i]==(FOUR*TEN+THREE))
  {
   mass[i]=98;
  }
  else if(Nu_charge[i]==(FOUR*TEN+FOUR))
  {
   mass[i]=101.07;
  }
  else if(Nu_charge[i]==(FOUR*TEN+FIVE))
  {
   mass[i]=102.91;
  }
  else if(Nu_charge[i]==(FOUR*TEN+SIX))
  {
   mass[i]=106.42;
  }
  else if(Nu_charge[i]==(FOUR*TEN+SEVEN))
  {
   mass[i]=107.87;
  }
  else if(Nu_charge[i]==(FOUR*TEN+EIGHT))
  {
   mass[i]=112.41;
  }
  else if(Nu_charge[i]==(FOUR*TEN+NINE))
  {
   mass[i]=114.82;
  }
  else if(Nu_charge[i]==(FIVE*TEN))
  {
   mass[i]=118.71;
  }
  else if(Nu_charge[i]==(FIVE*TEN+ONE))
  {
   mass[i]=121.76;
  }
  else if(Nu_charge[i]==(FIVE*TEN+TWO))
  {
   mass[i]=127.60;
  }
  else if(Nu_charge[i]==(FIVE*TEN+THREE))
  {
   mass[i]=126.90;
  }
  else if(Nu_charge[i]==(FIVE*TEN+FOUR))
  {
   mass[i]=131.29;
  }
  else if(Nu_charge[i]==(FIVE*TEN+FIVE))
  {
   mass[i]=132.91;
  }
  else if(Nu_charge[i]==(FIVE*TEN+SIX))
  {
   mass[i]=137.33;
  }
  else if(Nu_charge[i]==(FIVE*TEN+SEVEN))
  {
   mass[i]=138.91;
  }
  else if(Nu_charge[i]==(FIVE*TEN+EIGHT))
  {
   mass[i]=140.12;
  }
  else if(Nu_charge[i]==(FIVE*TEN+NINE))
  {
   mass[i]=140.91;
  }
  else if(Nu_charge[i]==(SIX*TEN))
  {
   mass[i]=144.24;
  }
  else if(Nu_charge[i]==(SIX*TEN+ONE))
  {
   mass[i]=145;
  }
  else if(Nu_charge[i]==(SIX*TEN+TWO))
  {
   mass[i]=150.36;
  }
  else if(Nu_charge[i]==(SIX*TEN+THREE))
  {
   mass[i]=151.96;
  }
  else if(Nu_charge[i]==(SIX*TEN+FOUR))
  {
   mass[i]=157.25;
  }
  else if(Nu_charge[i]==(SIX*TEN+FIVE))
  {
   mass[i]=158.93;
  }
  else if(Nu_charge[i]==(SIX*TEN+SIX))
  {
   mass[i]=162.50;
  }
  else if(Nu_charge[i]==(SIX*TEN+SEVEN))
  {
   mass[i]=164.93;
  }
  else if(Nu_charge[i]==(SIX*TEN+EIGHT))
  {
   mass[i]=167.26;
  }
  else if(Nu_charge[i]==(SIX*TEN+NINE))
  {
   mass[i]=168.93;
  }
  else if(Nu_charge[i]==(SEVEN*TEN))
  {
   mass[i]=173.04;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+ONE))
  {
   mass[i]=174.97;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+TWO))
  {
   mass[i]=178.49;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+THREE))
  {
   mass[i]=180.95;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+FOUR))
  {
   mass[i]=183.84;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+FIVE))
  {
   mass[i]=186.21;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+SIX))
  {
   mass[i]=190.23;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+SEVEN))
  {
   mass[i]=192.22;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+EIGHT))
  {
   mass[i]=195.08;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+NINE))
  {
   mass[i]=196.97;
  }
  else if(Nu_charge[i]==(EIGHT*TEN))
  {
   mass[i]=200.59;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+ONE))
  {
   mass[i]=204.38;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+TWO))
  {
   mass[i]=207.2;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+THREE))
  {
   mass[i]=208.98;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+FOUR))
  {
   mass[i]=209;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+FIVE))
  {
   mass[i]=210;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+SIX))
  {
   mass[i]=222;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+SEVEN))
  {
   mass[i]=223;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+EIGHT))
  {
   mass[i]=226;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+NINE))
  {
   mass[i]=227;
  }
  else if(Nu_charge[i]==(NINE*TEN))
  {
   mass[i]=232.04;
  }
  else if(Nu_charge[i]==(NINE*TEN+ONE))
  {
   mass[i]=231.04;
  }
  else if(Nu_charge[i]==(NINE*TEN+TWO))
  {
   mass[i]=238.03;
  }
  else if(Nu_charge[i]==(NINE*TEN+THREE))
  {
   mass[i]=237;
  }
  else if(Nu_charge[i]==(NINE*TEN+FOUR))
  {
   mass[i]=244;
  }
  else if(Nu_charge[i]==(NINE*TEN+FIVE))
  {
   mass[i]=243;
  }
  else if(Nu_charge[i]==(NINE*TEN+SIX))
  {
   mass[i]=247;
  }
  else if(Nu_charge[i]==(NINE*TEN+SEVEN))
  {
   mass[i]=247;
  }
  else if(Nu_charge[i]==(NINE*TEN+EIGHT))
  {
   mass[i]=251;
  }
  else if(Nu_charge[i]==(NINE*TEN+NINE))
  {
   mass[i]=252;
  }
  else if(Nu_charge[i]==(TEN*TEN))
  {
   mass[i]=257;
  }
  else if(Nu_charge[i]==(TEN*TEN+ONE))
  {
   mass[i]=258;
  }
  else if(Nu_charge[i]==(TEN*TEN+TWO))
  {
   mass[i]=259;
  }
  else if(Nu_charge[i]==(TEN*TEN+THREE))
  {
   mass[i]=262;
  }
  else if(Nu_charge[i]==(TEN*TEN+FOUR))
  {
   mass[i]=261;
  }
  else if(Nu_charge[i]==(TEN*TEN+FIVE))
  {
   mass[i]=262;
  }
  else if(Nu_charge[i]==(TEN*TEN+SIX))
  {
   mass[i]=266;
  }
  else if(Nu_charge[i]==(TEN*TEN+SEVEN))
  {
   mass[i]=264;
  }
  else if(Nu_charge[i]==(TEN*TEN+EIGHT))
  {
   mass[i]=269;
  }
  else if(Nu_charge[i]==(TEN*TEN+NINE))
  {
   mass[i]=268;
  }
  else
  {
   mass[i]=ZERO;
  }
  }
  for(j=0;j<3;j++)
  {
   Rcm[j]=Rcm[j]+mass[i]*Cartesian_Coor[i][j];
  }
  Mtot=Mtot+mass[i];
 }
 for(i=0;i<3;i++)
 {
  Rcm[i]=Rcm[i]/(Mtot+pow(TEN,-TEN));
 }
 results_inert<<"Coordinates of the center of mass:"<<endl;
 results_inert<<setprecision(10)<<fixed;
 for(i=0;i<3;i++)
 {
  results_inert<<setw(15)<<Rcm[i]<<" ";
 }
 results_inert<<endl;
 results_inert<<endl;
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<3;j++)
  {
   Cartesian_Coor[i][j]=Cartesian_Coor[i][j]-Rcm[j];
  }
 }
 Norm_cm=norm3D(Cartesian_Coor[0]);
 for(i=0;i<natoms;i++)
 {
  Auxiliar=norm3D(Cartesian_Coor[i]);
  if(Auxiliar>Norm_cm)
  {
   Norm_cm=Auxiliar;
  }
 }
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<3;j++)
  {
   Cartesian_Coor[i][j]=Cartesian_Coor[i][j]/(Norm_cm+pow(TEN,-TEN));
  }
 }
 for(i=0;i<natoms;i++)
 {
  Im[0][0]=Im[0][0]+mass[i]*(pow(Cartesian_Coor[i][1],TWO)+pow(Cartesian_Coor[i][2],TWO));
  Im[1][1]=Im[1][1]+mass[i]*(pow(Cartesian_Coor[i][0],TWO)+pow(Cartesian_Coor[i][2],TWO));
  Im[2][2]=Im[2][2]+mass[i]*(pow(Cartesian_Coor[i][0],TWO)+pow(Cartesian_Coor[i][1],TWO));
  Im[0][1]=Im[0][1]-mass[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][1];
  Im[0][2]=Im[0][2]-mass[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][2];
  Im[1][2]=Im[1][2]-mass[i]*Cartesian_Coor[i][1]*Cartesian_Coor[i][2];
 }
 if(abs(Im[0][1])<pow(TEN,-NINE)){Im[0][1]=ZERO;}
 if(abs(Im[0][2])<pow(TEN,-NINE)){Im[0][2]=ZERO;}
 if(abs(Im[1][2])<pow(TEN,-NINE)){Im[1][2]=ZERO;}
 Im[1][0]=Im[0][1];
 Im[2][0]=Im[0][2];
 Im[2][1]=Im[1][2];
 results_inert<<"Normalized inertia tensor for the atomic masses (pre-diagonalization):"<<endl;
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Im[i][j]<<" ";
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 if(abs(Im[0][1])==ZERO && abs(Im[0][2])==ZERO && abs(Im[1][2])==ZERO)
 {
  results_inert<<"The inertia tensor is alredy diagonal"<<endl;
  results_inert<<endl;
  for(i=0;i<3;i++)
  {
   EigenV[i][i]=ONE;
  }
 }
 else
 {
  results_inert<<"Using Jacobi diagonalization for the inertia tensor"<<endl;
  results_inert<<endl;
  jacobi(3,Im,EigenV);
 }
 results_inert<<"Normalized inertia tensor for the atomic masses (post-diagonalization and ordering):"<<endl;
 for(i=0;i<3;i++)
 {
  order[i]=i;
 }
 for(i=0;i<3;i++)
 {
  Norm=abs(Im[i][i]);
  for(j=i+1;j<3;j++)
  {
   if(abs(Im[j][j])>Norm)
   {
    Norm=abs(Im[j][j]);
    Auxiliar2=Im[j][j];
    Im[j][j]=Im[i][i];
    Im[i][i]=Auxiliar2;
    pivot=order[i];
    order[i]=order[j];
    order[j]=pivot;
   }
  }
 }
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Im[i][j]<<" ";
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 results_inert<<"Eigenvectors for the diagonalized inertia tensor (columns):"<<endl;
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<EigenV[i][order[j]]<<" ";
   //Also store the rotation matrix but in rows!
   Rot_ICM[j][i]=EigenV[i][order[j]];
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 // Recover real distances
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<3;j++)
  {
   Cartesian_Coor[i][j]=Norm_cm*Cartesian_Coor[i][j];
  }
 }
 results_inert<<"Final coordinates (the order of atoms is the same as in the FCHK/WFN/WFX file):"<<endl;
 results_inert<<"       Charge                     Coordinates"<<endl;
 for(i=0;i<natoms;i++)
 {
  results_inert<<setw(15)<<Nu_charge[i];
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Cartesian_Coor[i][j];
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 results_inert<<"Final rotated coordinates (the order of atoms is the same as in the FCHK/WFN/WFX file):"<<endl;
 results_inert<<"       Charge                     Coordinates"<<endl;
 for(i=0;i<natoms;i++)
 {
  results_inert<<setw(15)<<Nu_charge[i];
  for(j=0;j<3;j++){Rcm[j]=ZERO;}
  for(k=0;k<3;k++)
  {
   for(j=0;j<3;j++)
   {
    Rcm[k]=Rcm[k]+EigenV[j][order[k]]*Cartesian_Coor[i][j];
   }
  }
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Rcm[j];
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 results_inert.close();
 }
 else
 {
  results_inert.open((name_file.substr(0,(name_file.length()-5))+"_inert.tmp").c_str());
  Im=new double*[3];
  EigenV=new double*[3];
  mass=new double[natoms];
  Auxiliar2=ZERO;
  Auxiliar2=Auxiliar2;
  k=0;
  k=k;
  for(i=0;i<3;i++)
  {
   Im[i]=new double[3];
   EigenV[i]=new double[3];
   for(j=0;j<3;j++)
   {
    Im[i][j]=ZERO;
    EigenV[i][j]=ZERO;
   }
  }
  for(i=0;i<natoms;i++)
  {
   mass[i]=ZERO;
  //Silly if to cheack that initialization was correct!
  if(mass[i]==ZERO)
  {
  if(Nu_charge[i]==ZERO)
  {
   mass[i]=ZERO;
  }
  else if(Nu_charge[i]==ONE)
  {
   mass[i]=1.00079;
  }
  else if(Nu_charge[i]==TWO)
  {
   mass[i]=4.0026;
  }
  else if(Nu_charge[i]==THREE)
  {
   mass[i]=6.941;
  }
  else if(Nu_charge[i]==FOUR)
  {
   mass[i]=9.0122;
  }
  else if(Nu_charge[i]==FIVE)
  {
   mass[i]=10.811;
  }
  else if(Nu_charge[i]==SIX)
  {
   mass[i]=12.011;
  }
  else if(Nu_charge[i]==SEVEN)
  {
   mass[i]=14.007;
  }
  else if(Nu_charge[i]==EIGHT)
  {
   mass[i]=15.999;
  }
  else if(Nu_charge[i]==NINE)
  {
   mass[i]=18.998;
  }
  else if(Nu_charge[i]==TEN)
  {
   mass[i]=20.180;
  }
  else if(Nu_charge[i]==(TEN+ONE))
  {
   mass[i]=22.990;
  }
  else if(Nu_charge[i]==(TEN+TWO))
  {
   mass[i]=24.305;
  }
  else if(Nu_charge[i]==(TEN+THREE))
  {
   mass[i]=26.982;
  }
  else if(Nu_charge[i]==(TEN+FOUR))
  {
   mass[i]=28.086;
  }
  else if(Nu_charge[i]==(TEN+FIVE))
  {
   mass[i]=30.974;
  }
  else if(Nu_charge[i]==(TEN+SIX))
  {
   mass[i]=32.065;
  }
  else if(Nu_charge[i]==(TEN+SEVEN))
  {
   mass[i]=35.453;
  }
  else if(Nu_charge[i]==(TEN+EIGHT))
  {
   mass[i]=39.948;
  }
  else if(Nu_charge[i]==(TEN+NINE))
  {
   mass[i]=39.098;
  }
  else if(Nu_charge[i]==(TWO*TEN))
  {
   mass[i]=40.078;
  }
  else if(Nu_charge[i]==(TWO*TEN+ONE))
  {
   mass[i]=44.956;
  }
  else if(Nu_charge[i]==(TWO*TEN+TWO))
  {
   mass[i]=47.867;
  }
  else if(Nu_charge[i]==(TWO*TEN+THREE))
  {
   mass[i]=50.942;
  }
  else if(Nu_charge[i]==(TWO*TEN+FOUR))
  {
   mass[i]=51.996;
  }
  else if(Nu_charge[i]==(TWO*TEN+FIVE))
  {
   mass[i]=54.938;
  }
  else if(Nu_charge[i]==(TWO*TEN+SIX))
  {
   mass[i]=55.845;
  }
  else if(Nu_charge[i]==(TWO*TEN+SEVEN))
  {
   mass[i]=58.933;
  }
  else if(Nu_charge[i]==(TWO*TEN+EIGHT))
  {
   mass[i]=58.693;
  }
  else if(Nu_charge[i]==(TWO*TEN+NINE))
  {
   mass[i]=63.546;
  }
  else if(Nu_charge[i]==(THREE*TEN))
  {
   mass[i]=65.39;
  }
  else if(Nu_charge[i]==(THREE*TEN+ONE))
  {
   mass[i]=69.723;
  }
  else if(Nu_charge[i]==(THREE*TEN+TWO))
  {
   mass[i]=72.61;
  }
  else if(Nu_charge[i]==(THREE*TEN+THREE))
  {
   mass[i]=74.922;
  }
  else if(Nu_charge[i]==(THREE*TEN+FOUR))
  {
   mass[i]=78.96;
  }
  else if(Nu_charge[i]==(THREE*TEN+FIVE))
  {
   mass[i]=79.904;
  }
  else if(Nu_charge[i]==(THREE*TEN+SIX))
  {
   mass[i]=83.80;
  }
  else if(Nu_charge[i]==(THREE*TEN+SEVEN))
  {
   mass[i]=85.468;
  }
  else if(Nu_charge[i]==(THREE*TEN+EIGHT))
  {
   mass[i]=87.62;
  }
  else if(Nu_charge[i]==(THREE*TEN+NINE))
  {
   mass[i]=89.906;
  }
  else if(Nu_charge[i]==(FOUR*TEN))
  {
   mass[i]=91.224;
  }
  else if(Nu_charge[i]==(FOUR*TEN+ONE))
  {
   mass[i]=92.906;
  }
  else if(Nu_charge[i]==(FOUR*TEN+TWO))
  {
   mass[i]=95.94;
  }
  else if(Nu_charge[i]==(FOUR*TEN+THREE))
  {
   mass[i]=98;
  }
  else if(Nu_charge[i]==(FOUR*TEN+FOUR))
  {
   mass[i]=101.07;
  }
  else if(Nu_charge[i]==(FOUR*TEN+FIVE))
  {
   mass[i]=102.91;
  }
  else if(Nu_charge[i]==(FOUR*TEN+SIX))
  {
   mass[i]=106.42;
  }
  else if(Nu_charge[i]==(FOUR*TEN+SEVEN))
  {
   mass[i]=107.87;
  }
  else if(Nu_charge[i]==(FOUR*TEN+EIGHT))
  {
   mass[i]=112.41;
  }
  else if(Nu_charge[i]==(FOUR*TEN+NINE))
  {
   mass[i]=114.82;
  }
  else if(Nu_charge[i]==(FIVE*TEN))
  {
   mass[i]=118.71;
  }
  else if(Nu_charge[i]==(FIVE*TEN+ONE))
  {
   mass[i]=121.76;
  }
  else if(Nu_charge[i]==(FIVE*TEN+TWO))
  {
   mass[i]=127.60;
  }
  else if(Nu_charge[i]==(FIVE*TEN+THREE))
  {
   mass[i]=126.90;
  }
  else if(Nu_charge[i]==(FIVE*TEN+FOUR))
  {
   mass[i]=131.29;
  }
  else if(Nu_charge[i]==(FIVE*TEN+FIVE))
  {
   mass[i]=132.91;
  }
  else if(Nu_charge[i]==(FIVE*TEN+SIX))
  {
   mass[i]=137.33;
  }
  else if(Nu_charge[i]==(FIVE*TEN+SEVEN))
  {
   mass[i]=138.91;
  }
  else if(Nu_charge[i]==(FIVE*TEN+EIGHT))
  {
   mass[i]=140.12;
  }
  else if(Nu_charge[i]==(FIVE*TEN+NINE))
  {
   mass[i]=140.91;
  }
  else if(Nu_charge[i]==(SIX*TEN))
  {
   mass[i]=144.24;
  }
  else if(Nu_charge[i]==(SIX*TEN+ONE))
  {
   mass[i]=145;
  }
  else if(Nu_charge[i]==(SIX*TEN+TWO))
  {
   mass[i]=150.36;
  }
  else if(Nu_charge[i]==(SIX*TEN+THREE))
  {
   mass[i]=151.96;
  }
  else if(Nu_charge[i]==(SIX*TEN+FOUR))
  {
   mass[i]=157.25;
  }
  else if(Nu_charge[i]==(SIX*TEN+FIVE))
  {
   mass[i]=158.93;
  }
  else if(Nu_charge[i]==(SIX*TEN+SIX))
  {
   mass[i]=162.50;
  }
  else if(Nu_charge[i]==(SIX*TEN+SEVEN))
  {
   mass[i]=164.93;
  }
  else if(Nu_charge[i]==(SIX*TEN+EIGHT))
  {
   mass[i]=167.26;
  }
  else if(Nu_charge[i]==(SIX*TEN+NINE))
  {
   mass[i]=168.93;
  }
  else if(Nu_charge[i]==(SEVEN*TEN))
  {
   mass[i]=173.04;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+ONE))
  {
   mass[i]=174.97;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+TWO))
  {
   mass[i]=178.49;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+THREE))
  {
   mass[i]=180.95;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+FOUR))
  {
   mass[i]=183.84;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+FIVE))
  {
   mass[i]=186.21;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+SIX))
  {
   mass[i]=190.23;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+SEVEN))
  {
   mass[i]=192.22;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+EIGHT))
  {
   mass[i]=195.08;
  }
  else if(Nu_charge[i]==(SEVEN*TEN+NINE))
  {
   mass[i]=196.97;
  }
  else if(Nu_charge[i]==(EIGHT*TEN))
  {
   mass[i]=200.59;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+ONE))
  {
   mass[i]=204.38;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+TWO))
  {
   mass[i]=207.2;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+THREE))
  {
   mass[i]=208.98;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+FOUR))
  {
   mass[i]=209;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+FIVE))
  {
   mass[i]=210;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+SIX))
  {
   mass[i]=222;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+SEVEN))
  {
   mass[i]=223;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+EIGHT))
  {
   mass[i]=226;
  }
  else if(Nu_charge[i]==(EIGHT*TEN+NINE))
  {
   mass[i]=227;
  }
  else if(Nu_charge[i]==(NINE*TEN))
  {
   mass[i]=232.04;
  }
  else if(Nu_charge[i]==(NINE*TEN+ONE))
  {
   mass[i]=231.04;
  }
  else if(Nu_charge[i]==(NINE*TEN+TWO))
  {
   mass[i]=238.03;
  }
  else if(Nu_charge[i]==(NINE*TEN+THREE))
  {
   mass[i]=237;
  }
  else if(Nu_charge[i]==(NINE*TEN+FOUR))
  {
   mass[i]=244;
  }
  else if(Nu_charge[i]==(NINE*TEN+FIVE))
  {
   mass[i]=243;
  }
  else if(Nu_charge[i]==(NINE*TEN+SIX))
  {
   mass[i]=247;
  }
  else if(Nu_charge[i]==(NINE*TEN+SEVEN))
  {
   mass[i]=247;
  }
  else if(Nu_charge[i]==(NINE*TEN+EIGHT))
  {
   mass[i]=251;
  }
  else if(Nu_charge[i]==(NINE*TEN+NINE))
  {
   mass[i]=252;
  }
  else if(Nu_charge[i]==(TEN*TEN))
  {
   mass[i]=257;
  }
  else if(Nu_charge[i]==(TEN*TEN+ONE))
  {
   mass[i]=258;
  }
  else if(Nu_charge[i]==(TEN*TEN+TWO))
  {
   mass[i]=259;
  }
  else if(Nu_charge[i]==(TEN*TEN+THREE))
  {
   mass[i]=262;
  }
  else if(Nu_charge[i]==(TEN*TEN+FOUR))
  {
   mass[i]=261;
  }
  else if(Nu_charge[i]==(TEN*TEN+FIVE))
  {
   mass[i]=262;
  }
  else if(Nu_charge[i]==(TEN*TEN+SIX))
  {
   mass[i]=266;
  }
  else if(Nu_charge[i]==(TEN*TEN+SEVEN))
  {
   mass[i]=264;
  }
  else if(Nu_charge[i]==(TEN*TEN+EIGHT))
  {
   mass[i]=269;
  }
  else if(Nu_charge[i]==(TEN*TEN+NINE))
  {
   mass[i]=268;
  }
  else
  {
   mass[i]=ZERO;
  }
  }
  for(j=0;j<3;j++)
  {
   Rcm[j]=Rcm[j]+mass[i]*Cartesian_Coor[i][j];
  }
  Mtot=Mtot+mass[i];
 }
 for(i=0;i<3;i++)
 {
  Rcm[i]=Rcm[i]/(Mtot+pow(TEN,-TEN));
 }
 results_inert<<"Coordinates of the center of mass:"<<endl;
 results_inert<<setprecision(10)<<fixed;
 for(i=0;i<3;i++)
 {
  results_inert<<setw(15)<<Rcm[i]<<" ";
 }
 results_inert<<endl;
 results_inert<<endl;
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<3;j++)
  {
   Cartesian_Coor[i][j]=Cartesian_Coor[i][j]-Rcm[j];
  }
 }
 Norm_cm=norm3D(Cartesian_Coor[0]);
 for(i=0;i<natoms;i++)
 {
  Auxiliar=norm3D(Cartesian_Coor[i]);
  if(Auxiliar>Norm_cm)
  {
   Norm_cm=Auxiliar;
  }
 }
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<3;j++)
  {
   Cartesian_Coor[i][j]=Cartesian_Coor[i][j]/(Norm_cm+pow(TEN,-TEN));
  }
 }
 for(i=0;i<natoms;i++)
 {
  Im[0][0]=Im[0][0]+mass[i]*(pow(Cartesian_Coor[i][1],TWO)+pow(Cartesian_Coor[i][2],TWO));
  Im[1][1]=Im[1][1]+mass[i]*(pow(Cartesian_Coor[i][0],TWO)+pow(Cartesian_Coor[i][2],TWO));
  Im[2][2]=Im[2][2]+mass[i]*(pow(Cartesian_Coor[i][0],TWO)+pow(Cartesian_Coor[i][1],TWO));
  Im[0][1]=Im[0][1]-mass[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][1];
  Im[0][2]=Im[0][2]-mass[i]*Cartesian_Coor[i][0]*Cartesian_Coor[i][2];
  Im[1][2]=Im[1][2]-mass[i]*Cartesian_Coor[i][1]*Cartesian_Coor[i][2];
 }
 if(abs(Im[0][1])<pow(TEN,-NINE)){Im[0][1]=ZERO;}
 if(abs(Im[0][2])<pow(TEN,-NINE)){Im[0][2]=ZERO;}
 if(abs(Im[1][2])<pow(TEN,-NINE)){Im[1][2]=ZERO;}
 Im[1][0]=Im[0][1];
 Im[2][0]=Im[0][2];
 Im[2][1]=Im[1][2];
 results_inert<<"Normalized inertia tensor for the atomic masses (pre-diagonalization):"<<endl;
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Im[i][j]<<" ";
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 if(abs(Im[0][1])==ZERO && abs(Im[0][2])==ZERO && abs(Im[1][2])==ZERO)
 {
  results_inert<<"The inertia tensor is alredy diagonal"<<endl;
  results_inert<<endl;
  for(i=0;i<3;i++)
  {
   EigenV[i][i]=ONE;
  }
 }
 else
 {
  results_inert<<"Using Jacobi diagonalization for the inertia tensor"<<endl;
  results_inert<<endl;
  jacobi(3,Im,EigenV);
 }
 results_inert<<"Normalized inertia tensor for the atomic masses (post-diagonalization and ordering):"<<endl;
 for(i=0;i<3;i++)
 {
  order[i]=i;
 }
 for(i=0;i<3;i++)
 {
  Norm=abs(Im[i][i]);
  for(j=i+1;j<3;j++)
  {
   if(abs(Im[j][j])>Norm)
   {
    Norm=abs(Im[j][j]);
    Auxiliar2=Im[j][j];
    Im[j][j]=Im[i][i];
    Im[i][i]=Auxiliar2;
    pivot=order[i];
    order[i]=order[j];
    order[j]=pivot;
   }
  }
 }
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Im[i][j]<<" ";
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 results_inert<<"Eigenvectors for the diagonalized inertia tensor (columns):"<<endl;
 for(i=0;i<3;i++)
 {
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<EigenV[i][order[j]]<<" ";
   //Also store the rotation matrix but in rows!
   Rot_ICM[j][i]=EigenV[i][order[j]];
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 // Recover real distances
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<3;j++)
  {
   Cartesian_Coor[i][j]=Norm_cm*Cartesian_Coor[i][j];
  }
 }
 results_inert<<"Final coordinates (the order of atoms is the same as in the FCHK/WFN/WFX file):"<<endl;
 results_inert<<"       Charge                     Coordinates"<<endl;
 for(i=0;i<natoms;i++)
 {
  results_inert<<setw(15)<<Nu_charge[i];
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Cartesian_Coor[i][j];
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 results_inert<<"Final rotated coordinates (the order of atoms is the same as in the FCHK/WFN/WFX file):"<<endl;
 results_inert<<"       Charge                     Coordinates"<<endl;
 for(i=0;i<natoms;i++)
 {
  results_inert<<setw(15)<<Nu_charge[i];
  for(j=0;j<3;j++){Rcm[j]=ZERO;}
  for(k=0;k<3;k++)
  {
   for(j=0;j<3;j++)
   {
    Rcm[k]=Rcm[k]+EigenV[j][order[k]]*Cartesian_Coor[i][j];
   }
  }
  for(j=0;j<3;j++)
  {
   results_inert<<setw(15)<<Rcm[j];
  }
  results_inert<<endl;
 }
 results_inert<<endl;
 results_inert.close();
 }
 counter=0;
 for(i=0;i<3;i++)
 {
  delete[] Im[i];Im[i]=NULL;
  delete[] EigenV[i];EigenV[i]=NULL;
 }
 delete[] Im;Im=NULL;
 delete[] EigenV;EigenV=NULL;
 delete[] mass;mass=NULL;
}
