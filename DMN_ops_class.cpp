#include"DMN_ops_class.h"


DMN_OPS::DMN_OPS(){}
DMN_OPS::DMN_OPS(string name_dmn_file,int order)
{
 dmn_file=false;
 ifstream check;
 check.open((name_dmn_file).c_str());
 if(check.good())
 {
  fchk=false;
  dm1=false;dm2=false;dm3=false;dm4=false;
  evaluate=false;trace=false;grad_density=false;
  name_file=name_dmn_file;diagonalize_dm1=false;
  diagonalize_dm1_ab=false;diag_done=false;
  diag_ab_done=false;dm2store=false;DMNinCORE=false;
  dmn_file=true;
  RECORD_DELIMITER_LENGTH=4;
  if(order==1)
  {dm1=true;}
  else if(order==2)
  {dm2=true;}
  else if(order==3)
  {dm3=true;}
  else if(order==4)
  {dm4=true;}
  threshold=pow(TEN,-EIGHT);
  send_threshold=threshold;
  check.close();
 }
 else
 {cout<<"The file "<<name_dmn_file<<" was not found!"<<endl;}
}
//Copy object
DMN_OPS::DMN_OPS(const DMN_OPS &DMN)
{
 //Not copying neither diagonalization info nor fchk info
 dm1=DMN.dm1;
 dm2=DMN.dm2;
 dm3=DMN.dm3;
 dm4=DMN.dm4;
 threshold=DMN.send_threshold;
 name_file=DMN.name_file;
 RECORD_DELIMITER_LENGTH=4;
}
//Delete object
DMN_OPS::~DMN_OPS()
{
 int i;
 if(fchk)
 {nbasisf=FCHK_for_DMN[0].nbasisf;}
 if(diag_done && fchk)
 {
  for(i=0;i<nbasisf;i++)
  {
   delete[] rho_matrix[i];rho_matrix[i]=NULL;
   delete[] NO_coef[i];NO_coef[i]=NULL;
  }
  delete[] rho_matrix;
  delete[] NO_coef;
  rho_matrix=NULL;
  NO_coef=NULL;
 }
 if(diag_ab_done && fchk)
 {
  for(i=0;i<nbasisf;i++)
  {
   delete[] rho_matrixa[i];rho_matrixa[i]=NULL;
   delete[] NO_coefa[i];NO_coefa[i]=NULL;
   delete[] rho_matrixb[i];rho_matrixb[i]=NULL;
   delete[] NO_coefb[i];NO_coefb[i]=NULL;
  }
  delete[] rho_matrixa;
  delete[] NO_coefa;
  delete[] rho_matrixb;
  delete[] NO_coefb;
  rho_matrixa=NULL;
  NO_coefa=NULL;
  rho_matrixb=NULL;
  NO_coefb=NULL;
 }
 if(fchk)
 {FCHK_for_DMN[0].READ_FCHK_WFN::~READ_FCHK_WFN();}
 if(DMNinCORE && fchk)
 {
  if(dm1)
  {
   for(i=0;i<2*nbasisf;i++)
   {
    delete[] DM1[i];
    DM1[i]=NULL;
   }
   delete[] DM1;
   DM1=NULL;
  }
 }
}
//Make a copy of the fchk file for DMN (get the MOs)
void DMN_OPS::set_fchk(string name_fchk,string name_log,bool wfn_fchk,bool log,bool cas,bool cm,int multiplicity)
{
 ifstream test;
 test.open((name_fchk).c_str());
 if(test.good())
 {
  test.close();
  fchk=true;
  FCHK_for_DMN=new READ_FCHK_WFN(name_fchk,name_log,wfn_fchk,log,cas,cm,multiplicity);
 }
}
//Set a threshold
void DMN_OPS::set_thershold(double THRESHOLD)
{
 threshold=THRESHOLD;
 send_threshold=threshold;
}

//Size of the basis
int DMN_OPS::nbasis()
{
 int n=0;
 if(fchk)
 {n=FCHK_for_DMN[0].nbasisf;}
 return n;
}
//Set store or not the DMN
void DMN_OPS::dmnincore(bool store_or_not,double mem)
{
 int total_terms=0,storable_terms=0,i,j;
 if(dmn_file)
 {
  if(fchk)
  {
   if(store_or_not)
   {
    if(dm1)
    {
     total_terms=2*FCHK_for_DMN[0].nbasisf*2*FCHK_for_DMN[0].nbasisf;
     storable_terms=(int)(mem*BILLION/(EIGHT));
     cout<<"Number of DM1 terms "<<total_terms<<endl;
    }
    if(dm2)
    {
     Trace=calc_tr();
     total_terms=nterms;
     nterms=0;
     storable_terms=(int)(mem*BILLION/(EIGHT+FOUR*FOUR));//Gb to one double and four integers
    }
    if(storable_terms>=total_terms)
    {
     DMNinCORE=true;
     if(dm1)
     {
      //We always have alpha and beta separated from Edu's DMN program!!
      DM1=new double*[2*FCHK_for_DMN[0].nbasisf];
      for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
      {
       DM1[i]=new double[2*FCHK_for_DMN[0].nbasisf];
       for(j=0;j<2*FCHK_for_DMN[0].nbasisf;j++)
       {
        DM1[i][j]=ZERO;
       }
      }
      cout<<"DM1 array stored"<<endl;
      Read_use_DMN();
      cout<<"DM1 stored"<<endl;
     }
     if(dm2)
     {
      dm2store=true;
      Read_use_DMN();
     }
    }
    else
    {
     DMNinCORE=false;
     cout<<"Warning! Unable to store the DMN matrix due to the small ram memory given."<<endl;
     cout<<"[Calculations will be preformed without storing for the DM1 matrix (on the fly)]"<<endl;
    }
   }
  }
  else
  {
   cout<<"Warning! Include the fchk file to get the size of the basis set before storing the the DMN matrix!"<<endl;
  }
 }
 else
 {
  cout<<"Could not store the DMN file since it is missing"<<endl;
 }
}
//Trace
double DMN_OPS::calc_tr()
{
 Trace=ZERO;
 nterms=0;
 if(dmn_file)
 {
  trace=true;
  Read_use_DMN();
 }
 return Trace;
}
//Use the DM1 file to calc the Total Density for fchk files
void DMN_OPS::tot_dens_fchk(double *tot_density,int &aos)
{
 int i,j,k,counter=0;
 double **coef_MOs,**auxiliar;
 if(DMNinCORE)
 {
  //Succeed if the fchk and dm1 files are correct and the DM1 is already stored.
  coef_MOs=new double *[2*FCHK_for_DMN[0].nbasisf];
  auxiliar=new double *[2*FCHK_for_DMN[0].nbasisf];
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   coef_MOs[i]=new double[FCHK_for_DMN[0].nbasisf];
   auxiliar[i]=new double[FCHK_for_DMN[0].nbasisf];
  }
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<FCHK_for_DMN[0].nbasisf;j++)
   {
    coef_MOs[i][j]=ZERO;
    auxiliar[i][j]=ZERO;
   }
  }
  FCHK_for_DMN[0].bring_mo_coefs(coef_MOs);
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<FCHK_for_DMN[0].nbasisf;j++)
   {
    for(k=0;k<2*FCHK_for_DMN[0].nbasisf;k++)
    {auxiliar[i][j]=auxiliar[i][j]+DM1[i][k]*coef_MOs[k][j];}
   }
  }
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<=i;j++)
   {
    for(k=0;k<2*FCHK_for_DMN[0].nbasisf;k++)
    {
     tot_density[counter]=tot_density[counter]+coef_MOs[k][i]*auxiliar[k][j];
    }
    counter++;
   }
  }
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   delete[] coef_MOs[i];coef_MOs[i]=NULL;
   delete[] auxiliar[i];auxiliar[i]=NULL;
  }
  delete[] coef_MOs;delete[] auxiliar;
  coef_MOs=NULL;auxiliar=NULL;
 }
 else
 {cout<<"Could not calculate de 'Total SCF Density' since the fchk and/or dm1 file(s) is/are missing"<<endl;}
}
//Use the DM1 file to calc the Spin Density for fchk files
void DMN_OPS::spin_dens_fchk(double *tot_density,int &aos)
{
 int i,j,k,counter=0;
 double **coef_MOs,**auxiliar,**DM1a;
 if(DMNinCORE)
 {
  //Succeed if the fchk and dm1 files are correct and the DM1 is already stored.
  coef_MOs=new double *[2*FCHK_for_DMN[0].nbasisf];
  auxiliar=new double *[2*FCHK_for_DMN[0].nbasisf];
  DM1a=new double *[2*FCHK_for_DMN[0].nbasisf];
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   coef_MOs[i]=new double[FCHK_for_DMN[0].nbasisf];
   auxiliar[i]=new double[FCHK_for_DMN[0].nbasisf];
   DM1a[i]=new double[2*FCHK_for_DMN[0].nbasisf];
  }
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<FCHK_for_DMN[0].nbasisf;j++)
   {
    coef_MOs[i][j]=ZERO;
    auxiliar[i][j]=ZERO;
   }
  }
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<2*FCHK_for_DMN[0].nbasisf;j++)
   {
    DM1a[i][j]=-DM1[i][j];
    if(i%2==0 && j%2==0)
    {DM1a[i][j]=DM1a[i][j]+TWO*DM1[i][j];}
   }
  }
  FCHK_for_DMN[0].bring_mo_coefs(coef_MOs);
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<FCHK_for_DMN[0].nbasisf;j++)
   {
    for(k=0;k<2*FCHK_for_DMN[0].nbasisf;k++)
    {auxiliar[i][j]=auxiliar[i][j]+DM1a[i][k]*coef_MOs[k][j];}
   }
  }
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
  {
   for(j=0;j<=i;j++)
   {
    for(k=0;k<2*FCHK_for_DMN[0].nbasisf;k++)
    {
     tot_density[counter]=tot_density[counter]+coef_MOs[k][i]*auxiliar[k][j];
    }
    counter++;
   }
  }
  for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
  {
   delete[] coef_MOs[i];coef_MOs[i]=NULL;
   delete[] auxiliar[i];auxiliar[i]=NULL;
   delete[] DM1a[i]; DM1[i]=NULL;
  }
  delete[] coef_MOs;delete[] auxiliar;delete[] DM1a;
  coef_MOs=NULL;auxiliar=NULL;DM1a=NULL;
 }
 else
 {cout<<"Could not calculate de 'Spin SCF Density' since the fchk and/or dm1 file(s) is/are missing"<<endl;}
}
//Diagonalize DM1 with spins summed
//NO_coefficients are the coefficients that multiply MOs
//we take FCHK_for_DMN[0].nbasisf since we will sum up alpha
//and beta into a single contribution
void DMN_OPS::diagonalize()
{
 int i,j;
 if(fchk && dmn_file)
 {
  //Succeed if the fchk and dm1 files are correct (unless set or defined)
  nbasisf=FCHK_for_DMN[0].nbasisf;
  diagonalize_dm1=true;
  rho_matrix=new double*[nbasisf];
  NO_coef=new double*[nbasisf];
  for(i=0;i<nbasisf;i++)
  {
   rho_matrix[i]=new double[nbasisf];
   NO_coef[i]=new double[nbasisf];
  }
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<nbasisf;j++)
   {
    rho_matrix[i][j]=ZERO;
    NO_coef[i][j]=ZERO;
   }
  }
  Read_use_DMN();
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<nbasisf;j++)
   {
    if(i!=j)
    {rho_matrix[j][i]=rho_matrix[i][j];}
   }
  }
  jacobi(nbasisf,rho_matrix,NO_coef);
 }
 else
 {cout<<"Could not diagonalize the DM1 since the fchk and/or dm1 file(s) is/are missing"<<endl;}
}
//Diagonalize DM1 with spins separated
//rho_matrixa is assigned to alpha, rho_matrixb to beta
//NO_coefa is assigned to alpha, rho_matrixb to beta
//nbasisf=FCHK_for_DMN[0].nbasisf is used since we separate
void DMN_OPS::diagonalize_ab()
{
 int i,j;
 if(fchk && dmn_file)
 {
  //Succeed if the fchk and dm1 files are correct (unless set or defined)
  nbasisf=FCHK_for_DMN[0].nbasisf;
  diagonalize_dm1_ab=true;
  rho_matrixa=new double*[nbasisf];
  rho_matrixb=new double*[nbasisf];
  NO_coefa=new double*[nbasisf];
  NO_coefb=new double*[nbasisf];
  for(i=0;i<nbasisf;i++)
  {
   rho_matrixa[i]=new double[nbasisf];
   rho_matrixb[i]=new double[nbasisf];
   NO_coefa[i]=new double[nbasisf];
   NO_coefb[i]=new double[nbasisf];
  }
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<nbasisf;j++)
   {
    rho_matrixa[i][j]=ZERO;
    NO_coefa[i][j]=ZERO;
    rho_matrixb[i][j]=ZERO;
    NO_coefb[i][j]=ZERO;
   }
  }
  Read_use_DMN();
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<nbasisf;j++)
   {
    if(i!=j)
    {
     rho_matrixa[j][i]=rho_matrixa[i][j];
     rho_matrixb[j][i]=rho_matrixb[i][j];
    }
   }
  }
  jacobi(nbasisf,rho_matrixa,NO_coefa);
  jacobi(nbasisf,rho_matrixb,NO_coefb);
 }
 else
 {cout<<"Could not diagonalize the DM1 since the fchk and/or dm1 file(s) is/are missing"<<endl;}
}
//Build NO(i) if the fchk file and dm1 file are correctly set or defined.
void DMN_OPS::build_NO_DMN(double &NO,double point[3],bool &alpha_beta,int &numNO)
{
 int i;
 double eval=ZERO;
 double *AO,**AO_grad;
 AO_grad=new double *[3];
 for(i=0;i<3;i++)
 {AO_grad[i]=new double[FCHK_for_DMN[0].nbasisf];}
 AO=new double[FCHK_for_DMN[0].nbasisf];
 FCHK_for_DMN[0].build_AO_AOgrad2(AO,AO_grad,point);
 if(dmn_file && fchk)
 {
  if(alpha_beta)
  {
   //If the diagonalization was not done before, diagonalize the DM1
   if(!diag_ab_done)
   {diagonalize_ab();}
   if(FCHK_for_DMN[0].uhf)
   {
    //Alpha
    if(numNO%2==0)
    {
     for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
     {
      mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,2*i);
      eval=eval+NO_coefa[i][numNO/2]*mos_array[0].evaluation;
     }
    }
    //Beta
    else
    {
     for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
     {
      mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,2*i+1);
      eval=eval+NO_coefb[i][(numNO-1)/2]*mos_array[0].evaluation;
     }
    }
   }
   else
   {
    //Alpha
    if(numNO%2==0)
    {
     for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
     {
      mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,i);
      eval=eval+NO_coefa[i][numNO/2]*mos_array[0].evaluation;
     }
    }
    //Beta
    else
    {
     for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
     {
      mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,i);
      eval=eval+NO_coefb[i][(numNO-1)/2]*mos_array[0].evaluation;
     }
    }
   }
  }
  else
  {
   //If the diagonalization was not done before, diagonalize the DM1
   if(!diag_done)
   {diagonalize();}
   if(FCHK_for_DMN[0].uhf)
   {
    for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
    {
     mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,2*i);
     mos_array[1]=MO(FCHK_for_DMN[0],AO,AO_grad,2*i+1);
     eval=eval+NO_coef[i][numNO]*(mos_array[0].evaluation+mos_array[1].evaluation)/TWO;
    }
   }
   else
   {
    for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
    {
     mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,i);
     eval=eval+NO_coef[i][numNO]*mos_array[0].evaluation;
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
 NO=eval;
}
//Evaluate the DMN matrix of order n
double DMN_OPS::evaluation(double *point,double *point_prime)
{
 int i,j;
 Result=ZERO;
 if(dm1)
 {
  nbasisf=FCHK_for_DMN[0].nbasis();
 }
 else
 {cout<<"DMn of order greater than one will not be evaluated"<<endl;}
 Array_MOs=new double*[nbasisf];
 for(i=0;i<nbasisf;i++)
 {
  if(dm1){Array_MOs[i]=new double[2];}
 }
 if(dm1 && dmn_file)
 {
  density=false;
  POINT1[0]=point[0];POINT1[1]=point[1];POINT1[2]=point[2];
  POINT_PRIME1[0]=point_prime[0];POINT_PRIME1[1]=point_prime[1];POINT_PRIME1[2]=point_prime[2];
  if((POINT1[0]==POINT_PRIME1[0] && POINT1[1]==POINT_PRIME1[1]) && POINT1[2]==POINT_PRIME1[2])
  {density=true;}
  double *AO,**AO_grad;
  AO_grad=new double *[3];
  for(i=0;i<3;i++)
  {AO_grad[i]=new double[FCHK_for_DMN[0].nbasisf];}
  AO=new double[FCHK_for_DMN[0].nbasisf];
  double *AO2,**AO_grad2;
  AO_grad2=new double *[3];
  for(i=0;i<3;i++)
  {AO_grad2[i]=new double[FCHK_for_DMN[0].nbasisf];}
  AO2=new double[FCHK_for_DMN[0].nbasisf];
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++){AO[i]=ZERO;AO2[i]=ZERO;}
  for(i=0;i<3;i++)
  {for(j=0;j<FCHK_for_DMN[0].nbasisf;j++){AO_grad[i][j]=ZERO;AO_grad2[i][j]=ZERO;}}
  FCHK_for_DMN[0].build_AO_AOgrad2(AO,AO_grad,POINT1);
  FCHK_for_DMN[0].build_AO_AOgrad2(AO2,AO_grad2,POINT_PRIME1);
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<2;j++)
   {
    if(j==0)
    {
     mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,i);
     Array_MOs[i][j]=mos_array[0].evaluation;
    }
    else
    {
     mos_array[1]=MO(FCHK_for_DMN[0],AO2,AO_grad2,i);
     Array_MOs[i][j]=mos_array[1].evaluation;
    }
   }
  }
  delete[] AO;AO=NULL;
  for(i=0;i<3;i++){delete[] AO_grad[i];AO_grad[i]=NULL;}
  delete[] AO_grad;AO_grad=NULL;
  delete[] AO2;AO2=NULL;
  for(i=0;i<3;i++){delete[] AO_grad2[i];AO_grad2[i]=NULL;}
  delete[] AO_grad2;AO_grad2=NULL;
  evaluate=true;
  if(DMNinCORE)
  {
   evaluate=false;
   if(FCHK_for_DMN[0].uhf)
   {
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<=i;j++)
     {
      Result=Result+Array_MOs[i][0]*Array_MOs[j][1]*DM1[i][j];
      if(i!=j)
      {
       if(density)
       {Result=Result+Array_MOs[i][0]*Array_MOs[j][1]*DM1[i][j];}
       else
       {Result=Result+Array_MOs[j][0]*Array_MOs[i][1]*DM1[i][j];}
      }
     }
    }
   }
   else
   {
    for(i=0;i<2*nbasisf;i++)
    {
     for(j=0;j<=i;j++)
     {
      //r1'
      if(i%2==0)
      {element[0]=i/2;}
      else
      {element[0]=(i-1)/2;}
      //r1
      if(j%2==0)
      {element_prime[0]=j/2;}
      else
      {element_prime[0]=(j-1)/2;}
      Result=Result+Array_MOs[element[0]][0]*Array_MOs[element_prime[0]][1]*DM1[i][j];
      if(element[0]!=element_prime[0])
      {
       if(density)
       {Result=Result+Array_MOs[element_prime[0]][0]*Array_MOs[element[0]][0]*DM1[j][i];}
       else
       {Result=Result+Array_MOs[element_prime[0]][0]*Array_MOs[element[0]][1]*DM1[j][i];}
      }
     }
    }
   }
   density=false;
  }
  else
  {
   Read_use_DMN();
  }
 }
 else{} //DM2, DM3, DM4, etc Punctual evaluations are not needed!
 for(i=0;i<nbasisf;i++)
 {delete[] Array_MOs[i];Array_MOs[i]=NULL;}
 delete[] Array_MOs;
 Array_MOs=NULL;
 return Result;
}
//Gradient of the density and density at point
void DMN_OPS::grad_rho_r(double point[3],double grad[3],double &density)
{
 if(dm1 && dmn_file)
 {
  int i,j,k;
  for(i=0;i<3;i++){Grad[i]=ZERO;}
  Result=ZERO;
  nbasisf=FCHK_for_DMN[0].nbasis();
  Array_MOs=new double*[nbasisf];
  Grad_MOs=new double*[nbasisf];
  for(i=0;i<nbasisf;i++)
  {
   Grad_MOs[i]=new double[3];
   Array_MOs[i]=new double[1];
  }
  POINT1[0]=point[0];POINT1[1]=point[1];POINT1[2]=point[2];
  double *AO,**AO_grad;
  AO_grad=new double *[3];
  for(i=0;i<3;i++)
  {AO_grad[i]=new double[FCHK_for_DMN[0].nbasisf];}
  AO=new double[FCHK_for_DMN[0].nbasisf];
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++){AO[i]=ZERO;}
  for(i=0;i<3;i++)
  {for(j=0;j<FCHK_for_DMN[0].nbasisf;j++){AO_grad[i][j]=ZERO;}}
  FCHK_for_DMN[0].build_AO_AOgrad2(AO,AO_grad,POINT1);
  for(i=0;i<nbasisf;i++)
  {
   mos_array[0]=MO(FCHK_for_DMN[0],AO,AO_grad,i);
   Array_MOs[i][0]=mos_array[0].evaluation;
   Grad_MOs[i][0]=mos_array[0].grad[0];
   Grad_MOs[i][1]=mos_array[0].grad[1];
   Grad_MOs[i][2]=mos_array[0].grad[2];
  }
  delete[] AO;
  AO=NULL;
  for(i=0;i<3;i++)
  {delete[] AO_grad[i];AO_grad[i]=NULL;}
  delete[] AO_grad;
  AO_grad=NULL;
  grad_density=true;
  if(DMNinCORE)
  {
   if(FCHK_for_DMN[0].uhf)
   {
    for(i=0;i<nbasisf;i++)
    {
     for(j=0;j<=i;j++)
     {
      for(k=0;k<3;k++)
      {
       Grad[k]=Grad[k]+Grad_MOs[i][k]*Array_MOs[j][0]*DM1[i][j];
       Grad[k]=Grad[k]+Array_MOs[i][0]*Grad_MOs[j][k]*DM1[i][j];
       if(k==0)
       {
        Result=Result+Array_MOs[i][0]*Array_MOs[j][0]*DM1[i][j];
        if(i!=j)
        {
         Result=Result+Array_MOs[i][0]*Array_MOs[j][0]*DM1[i][j];
        }
       }
       if(i!=j)
       {
        Grad[k]=Grad[k]+Grad_MOs[j][k]*Array_MOs[i][0]*DM1[j][i];
        Grad[k]=Grad[k]+Array_MOs[j][0]*Grad_MOs[i][k]*DM1[j][i];
       }
      }
     }
    }
   }
   else
   {
    for(i=0;i<2*nbasisf;i++)
    {
     for(j=0;j<=i;j++)
     {
      //r1'
      if(i%2==0)
      {element[0]=i/2;}
      else
      {element[0]=(i-1)/2;}
      //r1
      if(j%2==0)
      {element_prime[0]=j/2;}
      else
      {element_prime[0]=(j-1)/2;}
      for(k=0;k<3;k++)
      {
       Grad[k]=Grad[k]+Grad_MOs[element[0]][k]*Array_MOs[element_prime[0]][0]*DM1[i][j];
       Grad[k]=Grad[k]+Array_MOs[element[0]][0]*Grad_MOs[element_prime[0]][k]*DM1[i][j];
       if(k==0)
       {
        Result=Result+Array_MOs[element[0]][0]*Array_MOs[element_prime[0]][0]*DM1[i][j];
        if(element[0]!=element_prime[0])
        {
         Result=Result+Array_MOs[element[0]][0]*Array_MOs[element_prime[0]][0]*DM1[i][j];
        }
       }
       if(element[0]!=element_prime[0])
       {
        Grad[k]=Grad[k]+Grad_MOs[element_prime[0]][k]*Array_MOs[element[0]][0]*DM1[j][i];
        Grad[k]=Grad[k]+Array_MOs[element_prime[0]][0]*Grad_MOs[element[0]][k]*DM1[j][i];
       }
      }
     }
    }
   }
   grad_density=false;
  }
  else
  {
   Read_use_DMN();
  }
  for(i=0;i<3;i++){grad[i]=Grad[i];}
  density=Result;
  for(i=0;i<nbasisf;i++)
  {
   delete[] Grad_MOs[i];Grad_MOs[i]=NULL;
   delete[] Array_MOs[i];Array_MOs[i]=NULL;
  }
  delete[] Grad_MOs;
  Grad_MOs=NULL;
  Array_MOs=NULL;
 }
 else
 {
  cout<<"Warning! Unable to compute the density from a DMN file of order greater than 1"<<endl;
  cout<<"or the DM1 file is missing."<<endl;
  grad[0]=ZERO;grad[1]=ZERO;grad[2]=ZERO;
 }
}
//Evaluate the DM2 at r1=r1' and r2=r2'
double DMN_OPS::dm2_coalesc(double Point[3],double Point2[3])
{
 int i,j;
 double eval=ZERO;
 if(DMNinCORE)
 {
  nbasisf=FCHK_for_DMN[0].nbasisf;
  //Array MOs will be equal to the number of AOs
  Array_MOs=new double*[nbasisf];
  for(i=0;i<nbasisf;i++){Array_MOs[i]=new double[4];}
  double *AO,**AO_grad;
  AO_grad=new double *[3];
  for(i=0;i<3;i++){AO_grad[i]=new double[FCHK_for_DMN[0].nbasisf];}
  AO=new double[FCHK_for_DMN[0].nbasisf];
  double *AO2,**AO_grad2;
  AO_grad2=new double *[3];
  for(i=0;i<3;i++){AO_grad2[i]=new double[FCHK_for_DMN[0].nbasisf];}
  AO2=new double[FCHK_for_DMN[0].nbasisf];
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++){AO[i]=ZERO;AO2[i]=ZERO;}
  for(i=0;i<3;i++)
  {for(j=0;j<FCHK_for_DMN[0].nbasisf;j++){AO_grad[i][j]=ZERO;AO_grad2[i][j]=ZERO;}}
  FCHK_for_DMN[0].build_AO_AOgrad2(AO,AO_grad,Point);
  FCHK_for_DMN[0].build_AO_AOgrad2(AO2,AO_grad2,Point2);
  for(i=0;i<nbasisf;i++)
  {
   if(FCHK_for_DMN[0].uhf)
   {
    mos_array[0]=MO(FCHK_for_DMN[0],Point,2*i);
    Array_MOs[i][0]=mos_array[0].evaluation;
    Array_MOs[i][2]=mos_array[0].evaluation;
    mos_array[1]=MO(FCHK_for_DMN[0],Point2,2*i);
    Array_MOs[i][1]=mos_array[1].evaluation;
    Array_MOs[i][3]=mos_array[1].evaluation;
   }
   else
   {
    mos_array[0]=MO(FCHK_for_DMN[0],Point,i);
    Array_MOs[i][0]=mos_array[0].evaluation;
    Array_MOs[i][2]=mos_array[0].evaluation;
    mos_array[1]=MO(FCHK_for_DMN[0],Point2,i);
    Array_MOs[i][1]=mos_array[1].evaluation;
    Array_MOs[i][3]=mos_array[1].evaluation;
   }
  }
  delete[] AO;AO=NULL;
  delete[] AO2;AO2=NULL;
  for(i=0;i<3;i++)
  {
   delete[] AO_grad2[i];AO_grad2[i]=NULL;
   delete[] AO_grad[i];AO_grad[i]=NULL;
  }
  delete[] AO_grad;AO_grad=NULL;
  delete[] AO_grad2;AO_grad2=NULL;
  for(i=0;i<nterms;i++)
  {
   element[0]=DM2[i].indexes[0];element[1]=DM2[i].indexes[1];
   element_prime[0]=DM2[i].indexes[2];element_prime[1]=DM2[i].indexes[3];
   if((element[0]%2==element_prime[0]%2) && (element[1]%2==element_prime[1]%2))
   {
    if(element[0]%2==element[1]%2)
    {
     if(element[0]%2==0)
     {
      element[0]=(element[0]/2)-1;
      element[1]=(element[1]/2)-1;
      element_prime[0]=(element_prime[0]/2)-1;
      element_prime[1]=(element_prime[1]/2)-1;
     }
     else
     {
      element[0]=(element[0]-1)/2;
      element[1]=(element[1]-1)/2;
      element_prime[0]=(element_prime[0]-1)/2;
      element_prime[1]=(element_prime[1]-1)/2;
     }
     eval=eval+DM2[i].Dijkl*Array_MOs[element[0]][0]*Array_MOs[element[1]][1]
                           *Array_MOs[element_prime[0]][2]*Array_MOs[element_prime[1]][3];
     eval=eval-DM2[i].Dijkl*Array_MOs[element[1]][0]*Array_MOs[element[0]][1]
                          *Array_MOs[element_prime[0]][2]*Array_MOs[element_prime[1]][3];
     eval=eval+DM2[i].Dijkl*Array_MOs[element[1]][0]*Array_MOs[element[0]][1]
                           *Array_MOs[element_prime[1]][2]*Array_MOs[element_prime[0]][3];
     eval=eval-DM2[i].Dijkl*Array_MOs[element[0]][0]*Array_MOs[element[1]][1]
                           *Array_MOs[element_prime[1]][2]*Array_MOs[element_prime[0]][3];
     if(element[0]!=element_prime[0]||element[1]!=element_prime[1])
     {
      eval=eval+DM2[i].Dijkl*Array_MOs[element_prime[0]][0]*Array_MOs[element_prime[1]][1]
                            *Array_MOs[element[0]][2]*Array_MOs[element[1]][3];
      eval=eval-DM2[i].Dijkl*Array_MOs[element_prime[1]][0]*Array_MOs[element_prime[0]][1]
                            *Array_MOs[element[0]][2]*Array_MOs[element[1]][3];
      eval=eval+DM2[i].Dijkl*Array_MOs[element_prime[1]][0]*Array_MOs[element_prime[0]][1]
                            *Array_MOs[element[1]][2]*Array_MOs[element[0]][3];
      eval=eval-DM2[i].Dijkl*Array_MOs[element_prime[0]][0]*Array_MOs[element_prime[1]][1]
                            *Array_MOs[element[1]][2]*Array_MOs[element[0]][3];
     }
    }
    else
    {
     if(element[0]%2==0)
     {
      element[0]=(element[0]/2)-1;
      element[1]=(element[1]-1)/2;
      element_prime[0]=(element_prime[0]/2)-1;
      element_prime[1]=(element_prime[1]-1)/2;
     }
     else
     {
      element[0]=(element[0]-1)/2;
      element[1]=(element[1]/2)-1;
      element_prime[0]=(element_prime[0]-1)/2;
      element_prime[1]=(element_prime[1]/2)-1;
     }
     eval=eval+DM2[i].Dijkl*Array_MOs[element[0]][0]*Array_MOs[element[1]][1]
                           *Array_MOs[element_prime[0]][2]*Array_MOs[element_prime[1]][3];
     eval=eval+DM2[i].Dijkl*Array_MOs[element[1]][0]*Array_MOs[element[0]][1]
                           *Array_MOs[element_prime[1]][2]*Array_MOs[element_prime[0]][3];
     if(element[0]!=element_prime[0]||element[1]!=element_prime[1])
     {
      eval=eval+DM2[i].Dijkl*Array_MOs[element_prime[0]][0]*Array_MOs[element_prime[1]][1]
                            *Array_MOs[element[0]][2]*Array_MOs[element[1]][3];
      eval=eval+DM2[i].Dijkl*Array_MOs[element_prime[1]][0]*Array_MOs[element_prime[0]][1]
                            *Array_MOs[element[1]][2]*Array_MOs[element[0]][3];
     }
    }
   }
  }
  for(i=0;i<nbasisf;i++)
  {
   delete[] Array_MOs[i];Array_MOs[i]=NULL;
  }
  delete[] Array_MOs;Array_MOs=NULL;
 }
 else
 {
  cout<<"Unable to procced since the DM2 matrix is not stored!"<<endl;
  cout<<"The evaluation has been set to zero"<<endl;
  eval=ZERO;
 }
 return eval;
}
//Read binary file .dmn
void DMN_OPS::Read_use_DMN()
{
 int i,pivot;
 int elementcpp[2],elementcpp2[2]; //we need elementcpp[2] and elementcpp2[2] because of the storing process.
 ifstream input_data(name_file.c_str(), ios::binary);
 if(dm1 && dmn_file)
 {
  cout<<"Storing/Using the DM1 matrix with threshold "<<threshold<<endl;
  element[0]=10;element_prime[0]=10;
  // open file in binary mode
  if(input_data.good())
  {
   while(element[0]!=0 || element_prime[0]!=0)
   {
    input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
    input_data.read((char*) &element[0], sizeof(element[0]));
    input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
    input_data.read((char*) &Dij, sizeof(Dij));
    input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
    if(abs(Dij)>=threshold)
    {
     //cout<<element[0]<<" "<<element_prime[0]<<" "<<Dij<<endl;
     if(trace)
     {
      if(element[0]==element_prime[0])
      {Trace=Trace+Dij;}
     }
     if(evaluate)
     {
      if(!FCHK_for_DMN[0].uhf)
      {
       //r1'
       if(element[0]%2==0)
       {element[0]=element[0]/2;}
       else
       {element[0]=(element[0]+1)/2;}
       //r1
       if(element_prime[0]%2==0)
       {element_prime[0]=element_prime[0]/2;}
       else
       {element_prime[0]=(element_prime[0]+1)/2;}
      }
      Result=Result+Array_MOs[element[0]-1][0]*Array_MOs[element_prime[0]-1][1]*Dij;
      if(element[0]!=element_prime[0])
      {
       if(density)
       {
        Result=Result+Array_MOs[element[0]-1][0]*Array_MOs[element_prime[0]-1][1]*Dij;
       }
       else
       {
        Result=Result+Array_MOs[element_prime[0]-1][0]*Array_MOs[element[0]-1][1]*Dij;
       }
      }
     }
     if(grad_density)
     {
      if(!FCHK_for_DMN[0].uhf)
      {
       //r1'
       if(element[0]%2==0)
       {element[0]=element[0]/2;}
       else
       {element[0]=(element[0]+1)/2;}
       //r1
       if(element_prime[0]%2==0)
       {element_prime[0]=element_prime[0]/2;}
       else
       {element_prime[0]=(element_prime[0]+1)/2;}
      }
      for(i=0;i<3;i++)
      {
       Grad[i]=Grad[i]+Grad_MOs[element[0]-1][i]*Array_MOs[element_prime[0]-1][0]*Dij;
       Grad[i]=Grad[i]+Array_MOs[element[0]-1][0]*Grad_MOs[element_prime[0]-1][i]*Dij;
       if(i==0)
       {Result=Result+Array_MOs[element[0]-1][0]*Array_MOs[element_prime[0]-1][0]*Dij;}
       if(element[0]!=element_prime[0])
       {
       Grad[i]=Grad[i]+Grad_MOs[element_prime[0]-1][i]*Array_MOs[element[0]-1][0]*Dij;
       Grad[i]=Grad[i]+Array_MOs[element_prime[0]-1][0]*Grad_MOs[element[0]-1][i]*Dij;
       if(i==0)
       {Result=Result+Array_MOs[element_prime[0]-1][0]*Array_MOs[element[0]-1][0]*Dij;}
       }
      }
     }
     if(diagonalize_dm1)
     {
      //I need to take as many as FCHK_for_DMN[0].nbasisf/2
      //So i devide the entire basis by two summing up alpha and beta
      //r1'
      if(element[0]%2==0)
      {element[0]=element[0]/2;}
      else
      {element[0]=(element[0]+1)/2;}
      //r1
      if(element_prime[0]%2==0)
      {element_prime[0]=element_prime[0]/2;}
      else
      {element_prime[0]=(element_prime[0]+1)/2;}
      if(rho_matrix[element[0]-1][element_prime[0]-1]==ZERO)
      {rho_matrix[element[0]-1][element_prime[0]-1]=Dij;}
      else
      {rho_matrix[element[0]-1][element_prime[0]-1]=rho_matrix[element[0]-1][element_prime[0]-1]+Dij;}
     }
     if(diagonalize_dm1_ab)
     {
      if(element[0]%2==0 && element_prime[0]%2==0)
      {
       element[0]=element[0]/2; element_prime[0]=element_prime[0]/2;
       if(rho_matrixb[element[0]-1][element_prime[0]-1]==ZERO)
       {rho_matrixb[element[0]-1][element_prime[0]-1]=Dij;}
       else
       {rho_matrixb[element[0]-1][element_prime[0]-1]=rho_matrixb[element[0]-1][element_prime[0]-1]+Dij;}
      }
      else if(element[0]%2!=0 && element_prime[0]%2!=0)
      {
       element[0]=(element[0]+1)/2; element_prime[0]=(element_prime[0]+1)/2;
       if(rho_matrixa[element[0]-1][element_prime[0]-1]==ZERO)
       {rho_matrixa[element[0]-1][element_prime[0]-1]=Dij;}
       else
       {rho_matrixa[element[0]-1][element_prime[0]-1]=rho_matrixa[element[0]-1][element_prime[0]-1]+Dij;}
      }
      else
      {}
     }
     if(DMNinCORE)
     {
      if(DM1[element[0]-1][element_prime[0]-1]==ZERO)
      {DM1[element[0]-1][element_prime[0]-1]=Dij;}
      else
      {DM1[element[0]-1][element_prime[0]-1]=DM1[element[0]-1][element_prime[0]-1]+Dij;}
      DM1[element_prime[0]-1][element[0]-1]=DM1[element[0]-1][element_prime[0]-1];
     }
    }
   }
  }
  else
  {}
 }
 else if(dm2 && dmn_file)
 {
  element[0]=10;element[1]=10;element_prime[0]=10;element_prime[1]=10;
  // open file in binary mode
  if(input_data.good())
  {
   while((element[0]!=0 || element_prime[0]!=0) || (element[1]!=0 || element_prime[1]!=0))
   {
    input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
    input_data.read((char*) &element[0], sizeof(element[0]));
    input_data.read((char*) &element[1], sizeof(element[1]));
    input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
    input_data.read((char*) &element_prime[1], sizeof(element_prime[1]));
    input_data.read((char*) &Dijkl, sizeof(Dijkl));
    input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
    if(abs(Dijkl)>=threshold)
    {
     if(trace)
     {
      nterms++;
      if(element[0]==element_prime[0] && element[1]==element_prime[1])
      {Trace=Trace+Dijkl;}
     }
     if(dm2store)
     {
      if(element[0]%2!=element_prime[0]%2 && element[0]%2==element_prime[1]%2)
      {
       pivot=element_prime[0];
       element_prime[0]=element_prime[1];
       element_prime[1]=pivot;
       Dijkl=-Dijkl;
      }
      nterms++;
      elementcpp[0]=element[0];elementcpp[1]=element[1];
      elementcpp2[0]=element_prime[0];elementcpp2[1]=element_prime[1];
      DM2.push_back(indexes_and_elements(elementcpp,elementcpp2,Dijkl));
     }
    }
   }
  }
 }
 else if (dm3 && dmn_file)
 {
  element[0]=10;element[1]=10;element[2]=10;element_prime[0]=10;element_prime[1]=10;
  element_prime[2]=10;
  // open file in binary mode
  if(input_data.good())
  {
   while(((element[0]!=0 || element_prime[0]!=0) || (element[1]!=0 || element_prime[1]!=0)) ||
        (element[2]!=0 || element_prime[2]!=0))
   {
    input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
    input_data.read((char*) &element[0], sizeof(element[0]));
    input_data.read((char*) &element[1], sizeof(element[1]));
    input_data.read((char*) &element[2], sizeof(element[2]));
    input_data.read((char*) &element_prime[0], sizeof(element_prime[0]));
    input_data.read((char*) &element_prime[1], sizeof(element_prime[1]));
    input_data.read((char*) &element_prime[2], sizeof(element_prime[2]));
    input_data.read((char*) &Dijklmn, sizeof(Dijklmn));
    input_data.seekg(RECORD_DELIMITER_LENGTH, ios::cur);
    if(abs(Dijklmn)>=threshold)
    {
     if(trace)
     {
      if((element[0]==element_prime[0] && element[1]==element_prime[1]) && (element[2]==element_prime[2]))
      {Trace=Trace+Dijklmn;}
     }
    }
   }
  }
 }
 else
 {}
 if(trace)
 {
  if(dm2) Trace=Trace*TWO;
  if(dm3) Trace=Trace*SIX;
 }
 if(diagonalize_dm1){diagonalize_dm1=false; diag_done=true;}
 if(diagonalize_dm1_ab){diagonalize_dm1_ab=false; diag_ab_done=true;}
 trace=false;density=false;grad_density=false;evaluate=false;dm2store=false;
 input_data.close();
}
