#include"DMN_ops_p_class.h"
DMN_P_OPS::DMN_P_OPS(){}
DMN_P_OPS::DMN_P_OPS(string name_dmn_file,int order)
{
 dmn_file=false;
 ifstream check;
 check.open((name_dmn_file).c_str());
 if(check.good())
 {
  dm1=false;dm2=false;dm3=false;dm4=false;fchk=false;
  evaluate=false;trace=false;grad_density=false;diagonalize_dm1=false;
  diagonalize_dm1_ab=false;diag_done=false;diag_ab_done=false;DMNinCORE=false;
  name_file=name_dmn_file;
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
DMN_P_OPS::DMN_P_OPS(const DMN_P_OPS &DMN)
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
DMN_P_OPS::~DMN_P_OPS()
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
//Make a copy of the fchk file for DMN (get the MOps coefficients)
void DMN_P_OPS::set_fchk(string name_fchk,string name_log,bool wfn_fchk,bool log,bool cas,int multiplicity)
{
 ifstream test;
 test.open((name_fchk).c_str());
 if(test.good())
 {
  test.close();
  fchk=true;
  FCHK_for_DMN=new READ_FCHK_WFN(name_fchk,name_log,wfn_fchk,log,cas,multiplicity);
 }
}
//Set a threshold
void DMN_P_OPS::set_thershold(double THRESHOLD)
{
 threshold=THRESHOLD;
 send_threshold=threshold;
}
//Size of the basis
int DMN_P_OPS::nbasis()
{
 int n=0;
 if(fchk)
 {n=FCHK_for_DMN[0].nbasisf;}
 return n;
}
//Set store or not the DMN
void DMN_P_OPS::dmnincore(bool store_or_not,double mem)
{
 int total_terms,storable_terms,i,j;
 if(dmn_file)
 {
  if(fchk)
  {
   if(store_or_not)
   {
    total_terms=2*FCHK_for_DMN[0].nbasisf*2*FCHK_for_DMN[0].nbasisf;
    storable_terms=(int)(mem*BILLION/(EIGHT));
    if(storable_terms>=total_terms)
    {
     DMNinCORE=true;
     if(dm1)
     {
      //We always have alpha and beta separated from Edu's program!!
      DM1=new double*[2*FCHK_for_DMN[0].nbasisf];
      for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
      {
       DM1[i]=new double[2*FCHK_for_DMN[0].nbasisf];
      }
      for(i=0;i<2*FCHK_for_DMN[0].nbasisf;i++)
      {
       for(j=0;j<2*FCHK_for_DMN[0].nbasisf;j++)
       {
        DM1[i][j]=ZERO;
       }
      }
      Read_use_DMN();
     }
    }
    else
    {
     DMNinCORE=false;
     cout<<"Warning! Unable to store the DM1 matrix due to the small ram memory given."<<endl;
     cout<<"Calculations will be preformed without storing the DM1 matrix (on the fly)"<<endl;
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
double DMN_P_OPS::calc_p_tr()
{
 Trace=ZERO;
 if(dmn_file)
 {
  trace=true;
  Read_use_DMN();
 }
 return Trace;
}
//Diagonalize DM1 with spins summed
void DMN_P_OPS::diagonalize()
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
void DMN_P_OPS::diagonalize_ab()
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
//Build NOp(i) if the fchk file and dm1 file are correctly set or defined.
void DMN_P_OPS::build_NOp_DMN(complex<double> &NOp,double point[3],bool &alpha_beta,int &numNO)
{
 int i;
 complex<double> eval=(ZERO,ZERO);
 complex<double> **AOp_grad;
 double **AOp;
 AOp=new double*[2];
 AOp_grad=new complex<double>*[3];
 for(i=0;i<3;i++)
 {AOp_grad[i]=new complex<double>[FCHK_for_DMN[0].nbasisf];}
 for(i=0;i<2;i++)
 {AOp[i]=new double[FCHK_for_DMN[0].nbasisf];}
 for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
 {AOp[0][i]=ZERO;AOp[1][i]=ZERO;AOp_grad[0][i]=ZERO;AOp_grad[1][i]=ZERO;AOp_grad[2][i]=ZERO;}
 FCHK_for_DMN[0].build_AOp2(AOp,point);
 FCHK_for_DMN[0].build_AOp_grad2(AOp_grad,point);
 if(fchk && dmn_file)
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
      mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,2*i);
      eval=eval+NO_coefa[i][numNO/2]*mos_array[0].evaluation;
     }
    }
    //Beta
    else
    {
     for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
     {
      mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,2*i+1);
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
      mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,i);
      eval=eval+NO_coefa[i][numNO/2]*mos_array[0].evaluation;
     }
    }
    //Beta
    else
    {
     for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
     {
      mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,i);
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
     mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,2*i);
     mos_array[1]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,2*i+1);
     eval=eval+NO_coef[i][numNO]*(mos_array[0].evaluation+mos_array[1].evaluation)/TWO;
    }
   }
   else
   {
    for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
    {
     mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,i);
     eval=eval+NO_coef[i][numNO]*mos_array[0].evaluation;
    }
   }
  }
 }
 for(i=0;i<2;i++)
 {delete[] AOp[i];AOp[i]=NULL;}
 for(i=0;i<3;i++)
 {delete[] AOp_grad[i];AOp_grad[i]=NULL;}
 delete[] AOp;
 delete[] AOp_grad;
 AOp=NULL;
 AOp_grad=NULL;
 NOp=eval;
}
//Evaluate the DMN matrix of order n
double DMN_P_OPS::evaluation_p(double *point_p,double *point_prime_p)
{
 int i,j;
 sum=(ZERO,ZERO);
 nbasisf=FCHK_for_DMN[0].nbasis();
 Array_MOs=new complex<double>*[nbasisf];
 for(i=0;i<nbasisf;i++)
 {
  if(dm1){Array_MOs[i]=new complex<double>[2];}
  else if(dm2){Array_MOs[i]=new complex<double>[4];}
  else if(dm3){Array_MOs[i]=new complex<double>[6];}
  else{Array_MOs[i]=new complex<double> [8];}
 }
 if(dm1 && dmn_file)
 {
  POINT1[0]=point_p[0];POINT1[1]=point_p[1];POINT1[2]=point_p[2];
  POINT_PRIME1[0]=point_prime_p[0];POINT_PRIME1[1]=point_prime_p[1];POINT_PRIME1[2]=point_prime_p[2];
  complex<double> **AOp_grad;
  double **AOp;
  AOp=new double*[2];
  AOp_grad=new complex<double>*[3];
  for(i=0;i<3;i++)
  {AOp_grad[i]=new complex<double>[FCHK_for_DMN[0].nbasisf];}
  for(i=0;i<2;i++)
  {AOp[i]=new double[FCHK_for_DMN[0].nbasisf];}
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
  {AOp[0][i]=ZERO;AOp[1][i]=ZERO;AOp_grad[0][i]=ZERO;AOp_grad[1][i]=ZERO;AOp_grad[2][i]=ZERO;}
  FCHK_for_DMN[0].build_AOp2(AOp,POINT1);
  FCHK_for_DMN[0].build_AOp_grad2(AOp_grad,POINT1);
  complex<double> **AOp_grad2;
  double **AOp2;
  AOp2=new double*[2];
  AOp_grad2=new complex<double>*[3];
  for(i=0;i<3;i++)
  {AOp_grad2[i]=new complex<double>[FCHK_for_DMN[0].nbasisf];}
  for(i=0;i<2;i++)
  {AOp2[i]=new double[FCHK_for_DMN[0].nbasisf];}
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
  {AOp2[0][i]=ZERO;AOp2[1][i]=ZERO;AOp_grad2[0][i]=ZERO;AOp_grad2[1][i]=ZERO;AOp_grad2[2][i]=ZERO;}
  FCHK_for_DMN[0].build_AOp2(AOp2,POINT_PRIME1);
  FCHK_for_DMN[0].build_AOp_grad2(AOp_grad2,POINT_PRIME1);
  for(i=0;i<nbasisf;i++)
  {
   for(j=0;j<2;j++)
   {
    if(j==0)
    {
     mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,i);
     Array_MOs[i][j]=mos_array[0].evaluation;
    }
    else
    {
     mos_array[1]=MOp(FCHK_for_DMN[0],AOp2,AOp_grad2,i);
     Array_MOs[i][j]=mos_array[1].evaluation;
    }
   }
  }
  for(i=0;i<2;i++)
  {delete[] AOp[i];AOp[i]=NULL;}
  for(i=0;i<3;i++)
  {delete[] AOp_grad[i];AOp_grad[i]=NULL;}
  delete[] AOp;
  delete[] AOp_grad;
  AOp=NULL;
  AOp_grad=NULL;
  for(i=0;i<2;i++)
  {delete[] AOp2[i];AOp2[i]=NULL;}
  for(i=0;i<3;i++)
  {delete[] AOp_grad2[i];AOp_grad2[i]=NULL;}
  delete[] AOp2;
  delete[] AOp_grad2;
  AOp2=NULL;
  AOp_grad2=NULL;
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
      sum=sum+conj(Array_MOs[i][0])*Array_MOs[j][1]*DM1[i][j];
      if(i!=j)
      {
       sum=sum+conj(Array_MOs[j][0])*Array_MOs[i][1]*DM1[i][j];
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
      sum=sum+conj(Array_MOs[element[0]][0])*Array_MOs[element_prime[0]][1]*DM1[i][j];
      if(element[0]!=element_prime[0])
      {
       sum=sum+conj(Array_MOs[element_prime[0]][0])*Array_MOs[element[0]][1]*DM1[j][i];
      }
     }
    }
   }
  }
  else
  {
   Read_use_DMN();
  }
 }
 else //The rest of the DMNs not needed punctual evals
 {}
 for(i=0;i<nbasisf;i++)
 {delete[] Array_MOs[i];Array_MOs[i]=NULL;}
 delete[] Array_MOs;
 Array_MOs=NULL;
 return real(sum);
}
//Gradient of the density and density at point
void DMN_P_OPS::grad_rho_p(double point_p[3],double grad[3],double &density)
{
 if(dm1 && dmn_file)
 {
  int i,j,k;
  sum=(ZERO,ZERO);
  for(i=0;i<3;i++)
  {Grad[i]=(ZERO,ZERO);}
  nbasisf=FCHK_for_DMN[0].nbasis();
  Array_MOs=new complex<double>*[nbasisf];
  Grad_MOs=new complex<double>*[nbasisf];
  POINT1[0]=point_p[0];POINT1[1]=point_p[1];POINT1[2]=point_p[2];
  complex<double> **AOp_grad;
  double **AOp;
  AOp=new double*[2];
  AOp_grad=new complex<double>*[3];
  for(i=0;i<3;i++)
  {AOp_grad[i]=new complex<double>[FCHK_for_DMN[0].nbasisf];}
  for(i=0;i<2;i++)
  {AOp[i]=new double[FCHK_for_DMN[0].nbasisf];}
  for(i=0;i<FCHK_for_DMN[0].nbasisf;i++)
  {AOp[0][i]=ZERO;AOp[1][i]=ZERO;AOp_grad[0][i]=ZERO;AOp_grad[1][i]=ZERO;AOp_grad[2][i]=ZERO;}
  FCHK_for_DMN[0].build_AOp2(AOp,POINT1);
  FCHK_for_DMN[0].build_AOp_grad2(AOp_grad,POINT1);
  for(i=0;i<nbasisf;i++)
  {
   Grad_MOs[i]=new complex<double>[3];
   Array_MOs[i]=new complex<double>[1];
  }
  for(i=0;i<nbasisf;i++)
  {
   mos_array[0]=MOp(FCHK_for_DMN[0],AOp,AOp_grad,i);
   Array_MOs[i][0]=mos_array[0].evaluation;
   Grad_MOs[i][0]=mos_array[0].grad[0];
   Grad_MOs[i][1]=mos_array[0].grad[1];
   Grad_MOs[i][2]=mos_array[0].grad[2];
  }
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
       Grad[k]=Grad[k]+conj(Grad_MOs[i][k])*Array_MOs[j][0]*DM1[i][j];
       Grad[k]=Grad[k]+conj(Array_MOs[i][0])*Grad_MOs[j][k]*DM1[i][j];
       if(k==0)
       {
        sum=sum+conj(Array_MOs[i][0])*Array_MOs[j][0]*DM1[i][j];
        if(i!=j)
        {
         sum=sum+conj(Array_MOs[i][0])*Array_MOs[j][0]*DM1[i][j];
        }
       }
       if(i!=j)
       {
        Grad[k]=Grad[k]+conj(Grad_MOs[j][k])*Array_MOs[i][0]*DM1[j][i];
        Grad[k]=Grad[k]+conj(Array_MOs[j][0])*Grad_MOs[i][k]*DM1[j][i];
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
       Grad[k]=Grad[k]+conj(Grad_MOs[element[0]][k])*Array_MOs[element_prime[0]][0]*DM1[i][j];
       Grad[k]=Grad[k]+conj(Array_MOs[element[0]][0])*Grad_MOs[element_prime[0]][k]*DM1[i][j];
       if(k==0)
       {
        sum=sum+conj(Array_MOs[element[0]][0])*Array_MOs[element_prime[0]][0]*DM1[i][j];
        if(element[0]!=element_prime[0])
        {
         sum=sum+conj(Array_MOs[element[0]][0])*Array_MOs[element_prime[0]][0]*DM1[i][j];
        }
       }
       if(element[0]!=element_prime[0])
       {
        Grad[k]=Grad[k]+conj(Grad_MOs[element_prime[0]][k])*Array_MOs[element[0]][0]*DM1[j][i];
        Grad[k]=Grad[k]+conj(Array_MOs[element_prime[0]][0])*Grad_MOs[element[0]][k]*DM1[j][i];
       }
      }
     }
    }
   }
   grad_density=false;
  }
  else
  {
   grad_density=true;
   Read_use_DMN();
  }
  for(i=0;i<3;i++)
  {grad[i]=real(Grad[i]);}
  density=real(sum);
  for(i=0;i<2;i++)
  {delete[] AOp[i];AOp[i]=NULL;}
  for(i=0;i<3;i++)
  {delete[] AOp_grad[i];AOp_grad[i]=NULL;}
  delete[] AOp;
  delete[] AOp_grad;
  AOp=NULL;
  AOp_grad=NULL;
  for(i=0;i<nbasisf;i++)
  {
   delete[] Grad_MOs[i];Grad_MOs[i]=NULL;
   delete[] Array_MOs[i];Array_MOs[i]=NULL;
  }
  delete[] Array_MOs;
  delete[] Grad_MOs;
  Array_MOs=NULL;
  Grad_MOs=NULL;
 }
 else
 {
  cout<<"Warning! Unable to compute the density from a DMN file of order greater than 1"<<endl;
  cout<<"or the DM1 file is missing."<<endl;
  grad[0]=ZERO;grad[1]=ZERO;grad[2]=ZERO;
 }
}
//Read binary file .dmn
void DMN_P_OPS::Read_use_DMN()
{
 int i;
 ifstream input_data(name_file.c_str(), ios::binary);
 if(dm1 && dmn_file)
 {
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
     if(trace)
     {
      if(element[0]==element_prime[0])
      {Trace=Trace+Dij;}
     }
     if(evaluate)
     {
      if(!FCHK_for_DMN[0].uhf)
      {
       //p1'
       if(element[0]%2==0)
       {element[0]=element[0]/2;}
       else
       {element[0]=(element[0]+1)/2;}
       //p1
       if(element_prime[0]%2==0)
       {element_prime[0]=element_prime[0]/2;}
       else
       {element_prime[0]=(element_prime[0]+1)/2;}
      }
      COMPLXconj(conjugated,Array_MOs[element[0]-1][0]);
      sum=sum+conjugated*Array_MOs[element_prime[0]-1][1]*Dij;
      if(element[0]!=element_prime[0])
      {
       COMPLXconj(conjugated,Array_MOs[element_prime[0]-1][0]);
       sum=sum+conjugated*Array_MOs[element[0]-1][1]*Dij;
      }
     }
     if(grad_density)
     {
      if(!FCHK_for_DMN[0].uhf)
      {
       //p1'
       if(element[0]%2==0)
       {element[0]=element[0]/2;}
       else
       {element[0]=(element[0]+1)/2;}
       //p1
       if(element_prime[0]%2==0)
       {element_prime[0]=element_prime[0]/2;}
       else
       {element_prime[0]=(element_prime[0]+1)/2;}
      }
      for(i=0;i<3;i++)
      {
       COMPLXconj(conjugated,Grad_MOs[element[0]-1][i]);
       Grad[i]=Grad[i]+conjugated*Array_MOs[element_prime[0]-1][0]*Dij;
       COMPLXconj(conjugated,Array_MOs[element[0]-1][0]);
       Grad[i]=Grad[i]+conjugated*Grad_MOs[element_prime[0]-1][i]*Dij;
       if(i==0)
       {
        COMPLXconj(conjugated,Array_MOs[element[0]-1][0]);
        sum=sum+conjugated*Array_MOs[element_prime[0]-1][0]*Dij;
       }
       if(element[0]!=element_prime[0])
       {
        COMPLXconj(conjugated,Grad_MOs[element_prime[0]-1][i]);
        Grad[i]=Grad[i]+conjugated*Array_MOs[element[0]-1][0]*Dij;
        COMPLXconj(conjugated,Array_MOs[element_prime[0]-1][0]);
        Grad[i]=Grad[i]+conjugated*Grad_MOs[element[0]-1][i]*Dij;
        if(i==0)
        {
         COMPLXconj(conjugated,Array_MOs[element_prime[0]-1][0]);
         sum=sum+conjugated*Array_MOs[element[0]-1][0]*Dij;
        }
       }
      }
     }
     if(diagonalize_dm1)
     {
      //Reduce indexes since we sum up alpha and beta
      //p1'
      if(element[0]%2==0)
      {element[0]=element[0]/2;}
      else
      {element[0]=(element[0]+1)/2;}
      //p1
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
      if(element[0]==element_prime[0] && element[1]==element_prime[1])
      {Trace=Trace+Dijkl;}
     }
    }
   }
  }
 }
 else if(dm3 && dmn_file)
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
 trace=false;grad_density=false;evaluate=false;
 input_data.close();
}
