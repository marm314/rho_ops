#include "utils_IO.h"

//Print INT file
void print_int(READ_FCHK_WFN &Read_fchk_wfn,string name_file,double **Sij,int nbasis,double &rho,double &rhoa,
double &rhob,string region)
{
 int i,j,k;
 ofstream file_int;
 file_int.open(name_file);
 file_int<<setprecision(14)<<fixed<<scientific;
 file_int<<" MOLECULAR SCF ENERGY (AU)  =        0.00000000000"<<endl;
 file_int<<endl;
 file_int<<" INTEGRATION IS OVER ATOM   X    "<<region<<endl;
 file_int<<" RESULTS OF THE INTEGRATION"<<endl;
 file_int<<"              N   "<<rho<<"  NET CHARGE"<<endl;
 file_int<<"              G   "<<endl;
 file_int<<"              K   0.00000000000000E+00        E(ATOM)   0.00000000000000E+00"<<endl;
 file_int<<"              L   0.00000000000000E+00"<<endl;
 file_int<<endl;
 file_int<<" The Atomic Overlap Matrix:"<<endl;
 file_int<<endl;
 if(!Read_fchk_wfn.wfn)
 {
  if(Read_fchk_wfn.uhf)
  {
   file_int<<"Unrestricted ";
   if(Read_fchk_wfn.open_shell)
   {file_int<<"Open-Shell Wavefunction"<<endl;}
   else
   {file_int<<"Closed-Shell Wavefunction"<<endl;}
  }
  else
  {
   file_int<<"Restricted Closed-Shell Wavefunction"<<endl;
  }
 }
 else
 {
  if(!Read_fchk_wfn.correlated)
  {
   if(Read_fchk_wfn.uhf)
   {
    file_int<<"Unrestricted ";
    if(Read_fchk_wfn.open_shell && Read_fchk_wfn.multiplicity!=1)
    {
     file_int<<"Open-Shell Wavefunction";
    }
    else
    {
     file_int<<"Closed-Shell Wavefunction";
    }
    file_int<<endl;
   }
   else
   {
    file_int<<"Restricted Closed-Shell Wavefunction"<<endl;
   }
  file_int<<endl;
  }
  else
  {
   file_int<<"Restricted Closed-Shell Wavefunction"<<endl; //Asummed! .wfn information
   file_int<<endl;                                         //is missing.
  }
 }
 file_int<<endl;
 k=0;
 for(i=0;i<nbasis;i++)
 {
  for(j=0;j<=i;j++)
  {
   file_int<<setw(25)<<Sij[i][j];
   k++;
   if(k==4){k=0;file_int<<endl;}
  }
 }
 file_int<<endl;
 file_int<<endl;
 file_int<<" ALPHA ELECTRONS (NA)   "<<rhoa<<endl;
 file_int<<" BETA ELECTRONS (NB)   "<<rhob<<endl;
 file_int<<endl;
 file_int<<" NORMAL TERMINATION OF PROAIMV"<<endl;
 file_int.close();
}

//Change SIJ in MOs to SIJ in NOs using DMN information
void mos_to_nos_dmn_sij(double **SIJ,READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands,string name_file,bool &wfn_fchk)
{
 int i,j,nbasis;
 double **aux,**aux1,**NO_coef;
 nbasis=Read_fchk_wfn.nbasis();
 aux=new double*[nbasis];
 aux1=new double*[nbasis];
 NO_coef=new double*[nbasis];
 for(i=0;i<nbasis;i++)
 {
  aux[i]=new double[nbasis];
  aux1[i]=new double[nbasis];
  NO_coef[i]=new double[nbasis];
 }
 for(i=0;i<nbasis;i++)
 {
  for(j=0;j<=i;j++)
  {
   aux[i][j]=ZERO;
   aux[j][i]=ZERO;
   NO_coef[i][j]=ZERO;
   NO_coef[j][i]=ZERO;
   aux1[i][j]=SIJ[i][j];
   aux1[j][i]=aux1[i][j];
  }
 }
 DMN_OPS dmn(Input_commands.name_dm1,1);
 dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
 dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
 if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
 if(Read_fchk_wfn.uhf)
 {
  for(i=0;i<nbasis;i++)
  {
   for(j=0;j<nbasis/2;j++)
   {
    if(i%2==0)
    {
     NO_coef[2*j][i]=dmn.NO_coefa[j][i/2];
    }
    else
    {
     NO_coef[2*j+1][i]=dmn.NO_coefb[j][(i-1)/2];
    }
   }
  }
  dmn.diagonalize_ab();
  matmul(nbasis,aux1,NO_coef,aux);
  mat_transpose(nbasis,NO_coef,aux1);
  mat_equal(nbasis,aux1,NO_coef);
  matmul(nbasis,NO_coef,aux,aux1);
 }
 else
 {
  dmn.diagonalize();
  matmul(nbasis,aux1,dmn.NO_coef,aux);
  mat_transpose(nbasis,dmn.NO_coef,NO_coef);
  matmul(nbasis,NO_coef,aux,aux1);
 }
 for(i=0;i<nbasis;i++)
 {
  for(j=0;j<=i;j++)
  {
   SIJ[i][j]=aux1[i][j];
  }
 }
 for(i=0;i<nbasis;i++)
 {
  delete[] NO_coef[i];NO_coef[i]=NULL;
  delete[] aux[i];aux[i]=NULL;
  delete[] aux1[i];aux1[i]=NULL;
 }
 delete[] NO_coef;NO_coef=NULL;
 delete[] aux;aux=NULL;
 delete[] aux1;aux1=NULL;
}

//Create a cube file
void cube_file(READ_FCHK_WFN &Read_fchk_wfn,string name_file,string op_cube,double cubex,double cubey, double cubez,double stepx,
double stepy,double stepz)
{
 int i,j,k,termsx,termsy,termsz;
 double eval,point[3],grad[3];
 double Density=ZERO,Density_alpha,Density_beta;
 double tauW_alpha,tauW_beta,tau_alpha,tau_beta,tcurr_alpha,tcurr_beta;
 double IND_alpha,IND_beta,ID_alpha,ID_beta;
 double coef_elf=(THREE/FIVE)*pow(SIX*PI*PI,TWO/THREE);
 bool not_elf,not_indic,tmp_false=false;
 ofstream cube;
 ifstream read_cube;
 string name_cube,aux;
 MO mo;
 NO no;
 not_elf=false;
 not_indic=false;
 termsx=(int)cubex;termsy=(int)cubey;termsz=(int)cubez;
 FILE *pFile;
 pFile=fopen((name_file+".cub_tmp").c_str(),"w");
 fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f \n",termsx,stepx,ZERO,ZERO);
 fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f \n",termsy,ZERO,stepy,ZERO);
 fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f \n",termsz,ZERO,ZERO,stepz);
 if(op_cube!="densityp" && op_cube!="shannonp" && op_cube!="fisherp")
 {
  for(i=0;i<Read_fchk_wfn.natoms;i++)
  {
   point[0]=Read_fchk_wfn.Cartesian_Coor[i][0]+cubex*stepx/TWO;
   point[1]=Read_fchk_wfn.Cartesian_Coor[i][1]+cubey*stepy/TWO;
   point[2]=Read_fchk_wfn.Cartesian_Coor[i][2]+cubez*stepz/TWO;
   fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f    %8.6f \n",(int)Read_fchk_wfn.Nu_charge[i],ZERO,
   point[0],point[1],point[2]);
  }
 }
 else
 {
  point[0]=cubex*stepx/TWO;
  point[1]=cubey*stepy/TWO;
  point[2]=cubez*stepz/TWO;
  fprintf(pFile,"%4d    %8.6f    %8.6f    %8.6f    %8.6f \n",1,ZERO,point[0],point[1],point[2]);
 }
 point[0]=-cubex*stepx/TWO;
 for(i=0;i<termsx;i++)
 {
  point[1]=-cubey*stepy/TWO;
  for(j=0;j<termsy;j++)
  {
   point[2]=-cubez*stepz/TWO;
   for(k=0;k<termsz;k++)
   {
    if(op_cube=="density")
    {Read_fchk_wfn.rho_eval(point,eval);}
    if(op_cube=="laplacian")
    {Read_fchk_wfn.rho_lapl(point,eval);}
    if(op_cube=="neglaplacian")
    {
     Read_fchk_wfn.rho_lapl(point,eval);
     eval=-eval;
    }
    if(op_cube=="shannon")
    {
     Read_fchk_wfn.rho_eval(point,eval);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=-eval/(double)Read_fchk_wfn.nelectrons*log(eval/(double)Read_fchk_wfn.nelectrons);
    }
    if(op_cube=="fisher")
    {
     Read_fchk_wfn.rho_eval(point,eval);
     Read_fchk_wfn.rho_grad(point,grad);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=pow(norm3D(grad),TWO)/(eval*(double)Read_fchk_wfn.nelectrons);
    }
    if(op_cube=="elf" || op_cube=="elfa" || op_cube=="elfb")
    {
     if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
     {
      not_elf=true;
      cout<<"Unable to compute the ELF from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
      cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
      cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
     }
     if(!not_elf)
     {
      pre_elf(Read_fchk_wfn,point,Density_alpha,Density_beta,tauW_alpha,tauW_beta,tau_alpha,tau_beta,tcurr_alpha,tcurr_beta);
      Density=Density_alpha+Density_beta;
     }
     if(op_cube=="elf")
     {
      if(pow(Density,FIVE/THREE)>pow(TEN,-TEN))
      {
       eval=TWO*(tau_alpha+tau_beta-tauW_alpha-tauW_beta-tcurr_alpha-tcurr_beta)/(coef_elf*pow(Density,FIVE/THREE));
       eval=ONE/(ONE+eval*eval);
      }
      else
      {
       eval=ZERO;
      }
     }
     else if(op_cube=="elfa")
     {
      if(pow(Density_alpha,FIVE/THREE)>pow(TEN,-TEN))
      {
       eval=TWO*(tau_alpha-tauW_alpha-tcurr_alpha)/(coef_elf*pow(Density_alpha,FIVE/THREE));
       eval=ONE/(ONE+eval*eval);
      }
      else
      {
       eval=ZERO;
      }
     }
     else
     {
      if(pow(Density_beta,FIVE/THREE)>pow(TEN,-TEN))
      {
       eval=TWO*(tau_beta-tauW_beta-tcurr_beta)/(coef_elf*pow(Density_beta,FIVE/THREE));
       eval=ONE/(ONE+eval*eval);
      }
      else
      {
       eval=ZERO;
      }
     }
    }
    if(op_cube=="id" || op_cube=="ind" || op_cube=="ida" ||  op_cube=="idb" ||  op_cube=="inda" ||  op_cube=="indb" || op_cube=="idnd")
    {
     if(Read_fchk_wfn.wfn && (Read_fchk_wfn.multiplicity!=1 && Read_fchk_wfn.correlated))
     {
      not_indic=true;
      cout<<"Unable to compute the INDICATORS from a correlated WFN file or the multiplicity was not correctly stated"<<endl;
      cout<<"remember to include the $MULTIPLICITY keyword or the $LOG file for calculations using WFN files"<<endl;
      cout<<"Nevertheless, do not forget that correlated WFN files have no SPIN components."<<endl;
     }
     if(!not_indic)
     {
      ID_IND_local(Read_fchk_wfn,point,ID_alpha,ID_beta,IND_alpha,IND_beta,tmp_false);
     }
     if(op_cube=="id")
     {
      eval=ID_alpha+ID_beta;
     }
     else if(op_cube=="ind")
     {
      eval=IND_alpha+IND_beta;
     }
     else if(op_cube=="inda")
     {
      eval=IND_alpha;
     }
     else if(op_cube=="indb")
     {
      eval=IND_beta;
     }
     else if(op_cube=="ida")
     {
      eval=ID_alpha;
     }
     else if(op_cube=="idb")
     {
      eval=ID_beta;
     }
     else
     {
      eval=ID_alpha+ID_beta+IND_alpha+IND_beta;
     }
    }
    if(op_cube=="mo")
    {
     mo=MO(Read_fchk_wfn,point,Read_fchk_wfn.Pair[0]);
     eval=mo.evaluation;
    }
    if(op_cube=="no")
    {
     no=NO(Read_fchk_wfn,point,Read_fchk_wfn.Pair[0]);
     eval=no.evaluation;
    }
    if(op_cube=="densityp")
    {Read_fchk_wfn.rho_p_eval(point,eval);}
    if(op_cube=="shannonp")
    {
     Read_fchk_wfn.rho_p_eval(point,eval);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=-eval/(double)Read_fchk_wfn.nelectrons*log(eval/(double)Read_fchk_wfn.nelectrons);
    }
    if(op_cube=="fisherp")
    {
     Read_fchk_wfn.rho_p_eval(point,eval);
     Read_fchk_wfn.rho_p_grad(point,grad);
     if(eval<pow(TEN,-EIGHT)){eval=pow(TEN,-EIGHT);}
     eval=pow(norm3D(grad),TWO)/(eval*(double)Read_fchk_wfn.nelectrons);
    }
    if((k%6==0))
    {
     if(eval>=ZERO)
     {fprintf(pFile," %11.5E ",eval);}
     else
     {fprintf(pFile,"%11.5E ",eval);}
    }
    else
    {
     if(eval>=ZERO)
     {fprintf(pFile," %11.5E ",eval);}
     else
     {fprintf(pFile,"%11.5E ",eval);}
    }
    if((k+1)%6==0)
    {fprintf(pFile,"\n");}
    point[2]=point[2]+stepz;
   }
   fprintf(pFile,"\n");
   point[1]=point[1]+stepy;
  }
  point[0]=point[0]+stepx;
 }
 fclose(pFile);
 read_cube.open((name_file+".cub_tmp").c_str());
 if(Read_fchk_wfn.wfn)
 {
  name_cube=name_file.substr(0,(name_file.length()-4));
 }
 else
 {
  name_cube=name_file.substr(0,(name_file.length()-5));
 }
 cube.open((name_cube+".cube").c_str());
 cube<<"RHO_OPS CUBE FILE."<<endl;
 cube<<"OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"<<endl;
 if(op_cube!="densityp" && op_cube!="shannonp" && op_cube!="fisherp")
 {
  //Only up to 99 atoms. I think is enough for this code =P
  cube<<"  ";
  if(Read_fchk_wfn.natoms<10){cube<<" ";}
  cube<<Read_fchk_wfn.natoms;
  cube<<"    0.000000    0.000000    0.000000"<<endl;
 }
 else
 {
  cube<<"  ";
  if(Read_fchk_wfn.natoms<10){cube<<" ";}
  cube<<1;
  cube<<"    0.000000    0.000000    0.000000"<<endl;
 }
 while(getline(read_cube,aux))
 {
  cube<<aux<<endl;
 }
 cube.close();
 read_cube.close();
 system(("rm "+name_file+".cub_tmp").c_str());
}

//Compute intermediate quantities used by ELF
void pre_elf(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &Density_alpha,double &Density_beta,
double &tauW_alpha,double &tauW_beta,double &tau_alpha,double &tau_beta,double &tcurr_alpha,double &tcurr_beta)
{
 if(Read_fchk_wfn.wfn)
 {
  int nbasis,i,j,k;
  double Grad_norm_alpha,shift_min=pow(TEN,-TWO*TEN);
  nbasis=Read_fchk_wfn.nbasis();
  Density_alpha=ZERO;Density_beta=ZERO;tauW_alpha=ZERO;tauW_beta=ZERO,tau_alpha=ZERO;tau_beta=ZERO;tcurr_alpha=ZERO;tcurr_beta=ZERO;
  if(Read_fchk_wfn.im_wfn_wfx)
  {
   complex<double>rhoa(ZERO,ZERO);
   complex<double>ztmp0(ZERO,ZERO);
   complex<double>ztmpI(ZERO,ONE);
   complex<double> **NO_orb_grad,AUX[3],gradA[3],grad_currA[3];
   NO_orb_grad=new complex<double>*[4];
   for(i=0;i<3;i++)
   {
    gradA[i]=ztmp0;grad_currA[i]=ztmp0;
   }
   for(i=0;i<4;i++)
   {
    NO_orb_grad[i]=new complex<double>[nbasis];
    for(j=0;j<nbasis;j++)
    {NO_orb_grad[i][j]=ztmp0;}
   }
   Read_fchk_wfn.orb_gradCC(Point,NO_orb_grad);
   if(Read_fchk_wfn.open_shell)
   {
    cout<<"Warning! pre_elf not prepared for open shell"<<endl;
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     rhoa=rhoa+(HALF*Read_fchk_wfn.Ocupation[i])*conj(NO_orb_grad[0][i])*NO_orb_grad[0][i];
     for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
     tau_alpha=tau_alpha+(HALF*Read_fchk_wfn.Ocupation[i])*pow(norm3DCC(AUX),TWO);
     for(j=0;j<3;j++)
     {
      gradA[j]=gradA[j]+(HALF*Read_fchk_wfn.Ocupation[i])
              *(conj(NO_orb_grad[0][i])*NO_orb_grad[j+1][i]+conj(NO_orb_grad[j+1][i])*NO_orb_grad[0][i]);
      grad_currA[j]=grad_currA[j]-ztmpI*HALF*(HALF*Read_fchk_wfn.Ocupation[i])
              *(conj(NO_orb_grad[0][i])*NO_orb_grad[j+1][i]-conj(NO_orb_grad[j+1][i])*NO_orb_grad[0][i]);
     }
    }
    Grad_norm_alpha=norm3DCC(gradA);
    Density_alpha=real(rhoa);Density_beta=Density_alpha;
    tau_alpha=HALF*tau_alpha;tau_beta=tau_alpha;
    tauW_alpha=pow(Grad_norm_alpha,TWO)/(EIGHT*Density_alpha+shift_min);tauW_beta=tauW_alpha;
    tcurr_alpha=HALF*pow(norm3DCC(grad_currA),TWO)/(Density_alpha+shift_min);tcurr_beta=tcurr_alpha;
   }
   //delete dynamic arrays
   for(i=0;i<4;i++)
   {delete[] NO_orb_grad[i];NO_orb_grad[i]=NULL;}
   delete[] NO_orb_grad;NO_orb_grad=NULL;
  }
  else
  {
   double Grad_alpha[3]={ZERO},Grad_beta[3]={ZERO};
   double Grad_norm_beta;
   double **NO_orb_grad,AUX[3]={ZERO};
   NO_orb_grad=new double*[4];
   for(i=0;i<4;i++)
   {
    NO_orb_grad[i]=new double[nbasis];
    for(j=0;j<nbasis;j++)
    {NO_orb_grad[i][j]=ZERO;}
   }
   if(Read_fchk_wfn.overlap || Read_fchk_wfn.wfn){Read_fchk_wfn.orb_grad(Point,NO_orb_grad);}//Get NOs and NO gradients
   if(!Read_fchk_wfn.correlated)
   {
    if(Read_fchk_wfn.open_shell)
    {
     for(i=0;i<nbasis;i++)
     {
      for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
      if(i%2==0)
      {
       for(k=0;k<3;k++)
       {
        Grad_alpha[k]=Grad_alpha[k]+TWO*NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       }
       Density_alpha=Density_alpha+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i];
       tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
      }
      else
      {
       for(k=0;k<3;k++)
       {
        Grad_beta[k]=Grad_beta[k]+TWO*NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       }
       Density_beta=Density_beta+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i];
       tau_beta=tau_beta+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO);
      }
     }
    }
    else
    {
     for(i=0;i<nbasis;i++)
     {
      for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
      for(k=0;k<3;k++)
      {
       Grad_alpha[k]=Grad_alpha[k]+NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       Grad_beta[k]=Grad_alpha[k];
      }
      Density_alpha=Density_alpha+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
      tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO)/TWO;
     }
     tau_beta=tau_alpha;
     Density_beta=Density_alpha;
    }
   }
   else
   {
    if(Read_fchk_wfn.open_shell)
    {
     Density_alpha=ZERO;
     Density_beta=ZERO;
     tau_alpha=ZERO;
     tau_beta=ZERO;
    }
    else
    {
     for(i=0;i<nbasis;i++)
     {
      for(j=1;j<4;j++){AUX[j-1]=NO_orb_grad[j][i];}
      for(k=0;k<3;k++)
      {
       Grad_alpha[k]=Grad_alpha[k]+NO_orb_grad[0][i]*Read_fchk_wfn.Ocupation[i]*NO_orb_grad[k+1][i];
       Grad_beta[k]=Grad_alpha[k];
      }
      Density_alpha=Density_alpha+Read_fchk_wfn.Ocupation[i]*NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
      tau_alpha=tau_alpha+Read_fchk_wfn.Ocupation[i]*pow(norm3D(AUX),TWO)/TWO;
     }
     tau_beta=tau_alpha;
     Density_beta=Density_alpha;
    }
   }
   Grad_norm_alpha=norm3D(Grad_alpha);
   Grad_norm_beta=norm3D(Grad_beta);
   tau_alpha=HALF*tau_alpha;
   tau_beta=HALF*tau_beta;
   tauW_alpha=pow(Grad_norm_alpha,TWO)/(EIGHT*Density_alpha+shift_min);
   tauW_beta=pow(Grad_norm_beta,TWO)/(EIGHT*Density_beta+shift_min);
   //delete dynamic arrays
   for(i=0;i<4;i++)
   {delete[] NO_orb_grad[i];NO_orb_grad[i]=NULL;}
   delete[] NO_orb_grad;NO_orb_grad=NULL;
  }
 }
 else
 {
  cout<<"pre_elf function not implemented for FCHK files"<<endl;
 }
}

//Change MOs to NOs for quadrature using the DMN information
void mos_to_nos_int_fchk_dm1(READ_FCHK_WFN &Read_fchk_wfn,Input &Input_commands,string name_file,double **ORBITALS,int &total_grid,int &nbasis,bool wfn_fchk)
{
 //Change MOs to NOs in ORBITALS using the DM1 file.
 int i,j,k,counter;
 double *NO_at_r;
 NO_at_r=new double[nbasis];
 DMN_OPS dmn(Input_commands.name_dm1,1);
 dmn.set_fchk(name_file,Input_commands.name_log,wfn_fchk,Input_commands.log,Input_commands.cas,Input_commands.cm,Input_commands.multiplicity);
 dmn.dmnincore(Input_commands.store_dmn,Input_commands.mem);
 if(Input_commands.dmn_thresh){dmn.set_thershold(Input_commands.dmn_threshold);}
 dmn.diagonalize_ab();
 for(i=0;i<total_grid;i++)
 {
  counter=0;
  if(Read_fchk_wfn.uhf)
  {
   for(j=0;j<nbasis;j=j+2)
   {
    NO_at_r[j]=ZERO;
    NO_at_r[j+1]=ZERO;
    for(k=0;k<Read_fchk_wfn.nbasisf;k++)
    {
     NO_at_r[j]=NO_at_r[j]+dmn.NO_coefa[k][counter]*ORBITALS[i][2*k];
     NO_at_r[j+1]=NO_at_r[j+1]+dmn.NO_coefb[k][counter]*ORBITALS[i][2*k+1];
    }
    counter++;
   }
  }
  else
  {
   for(j=0;j<nbasis;j++)
   {
    NO_at_r[j]=ZERO;
    for(k=0;k<Read_fchk_wfn.nbasisf;k++)
    {
     NO_at_r[j]=NO_at_r[j]+dmn.NO_coefa[k][j]*ORBITALS[i][k];
    }
   }
  }
  for(j=0;j<nbasis;j++)
  {
   ORBITALS[i][j]=NO_at_r[j];
  }
 }
 delete[] NO_at_r; NO_at_r=NULL;
}

