#include"Integrals_atomic.h"

int nrad_becke,nang_becke,n_bcps=0;
double *r_becke,*r_real_becke,*w_radial_becke,*w_theta_phi_becke,***wA;
double *x_becke,*y_becke,*z_becke;
//////////////////////////
//Functions description //
//////////////////////////
void Grid_atomic(READ_FCHK_WFN &Rho,string name,int &natom, int &nradial,int &nang,int &stiff,string partition)
{
 string name_bcp;
 bool becke=false,tfvc=false,becke_original=false,ssf=false;
 int i,j,k,l,m,ZB,ZC;
 double r_inf,r_sup,Point[3],Diff_Point[3],rB,rC,mu_BC,s_BC=ONE,nu_BC=ZERO,xi_BC,a_BC,u_BC,RBC,Sum_PB,a_ssf=0.64e0;
 double **Xi_XY_mat,**S_X,*P_X,Xi_rad;
 // Select the partition to use
 if(partition=="becke_original"){becke_original=true;}
 if(partition=="becke"){becke=true;}
 if(partition=="tfvc"){tfvc=true;}
 if(partition=="ssf"){ssf=true;}
 if(tfvc)
 {
  if((name[name.length()-1]=='n' || name[name.length()-1]=='N')||(name[name.length()-1]=='x' || name[name.length()-1]=='X')) // WFN or WFX
  {
   name_bcp=name.substr(0,(name.length()-4))+".bcp";
  }
  else // FCHK
  {
   name_bcp=name.substr(0,(name.length()-5))+".bcp";
  }
  ofstream init_out(name_bcp);
  init_out<<endl;
  init_out<<" Searching for BCP in "<<name<<endl;
  init_out<<endl;
  init_out.close();
 }
 //
 // Prepare the quadrature grid for each atom
 //
 //Grid for r
 nrad_becke=nradial;
 ifstream wr_read,r_read;
 r_becke=new double[nrad_becke];
 r_real_becke=new double[nrad_becke];
 w_radial_becke=new double[nrad_becke];
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(name,nrad_becke,r_inf,r_sup);
 wr_read.open((name+"_w.txt").c_str());
 r_read.open((name+"_x.txt").c_str());
 for(i=0;i<nrad_becke;i++)
 {
  r_read>>r_becke[i];
  wr_read>>w_radial_becke[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<nrad_becke;i++)
 {
  //r_real_becke is for entire space (0 to inf)
  r_real_becke[i]=r_becke[i]/(ONE-r_becke[i]);
 }
 //Grid for theta, phi
 Grid_avail_atomic(nang);
 nang_becke=nang;
 x_becke=new double[nang_becke];
 y_becke=new double[nang_becke];
 z_becke=new double[nang_becke];
 w_theta_phi_becke=new double[nang_becke];
 ld_by_order(nang_becke,x_becke,y_becke,z_becke,w_theta_phi_becke);
 //
 // Prepare the quadrature weights for each atom
 //
 P_X=new double[natom];
 S_X=new double*[natom];
 Xi_XY_mat=new double*[natom];
 wA=new double**[natom];
 for(i=0;i<natom;i++)
 {
  S_X[i]=new double[natom];
  Xi_XY_mat[i]=new double[natom];
  for(j=0;j<natom;j++)
  {
   S_X[i][j]=ONE;
   Xi_XY_mat[i][j]=ONE;
  }
  if(tfvc)
  {
   for(j=i+1;j<natom;j++)
   {
    Xi_XY_bcp(Rho,name_bcp,i,j,natom,Xi_rad);
    Xi_XY_mat[i][j]=Xi_rad;
   }
  }
 }
 if(tfvc)
 {
  ofstream write_bcp(name_bcp,std::ios_base::app);
  write_bcp<<endl;
  write_bcp<<" Number of BCP "<<setw(10)<<n_bcps<<endl;
  write_bcp<<endl;
  write_bcp.close();
  for(i=0;i<natom;i++)
  {
   for(j=0;j<=i;j++)
   {
    Xi_XY_mat[i][j]=ONE/Xi_XY_mat[j][i];
   }
  }
 }
 for(i=0;i<natom;i++)
 {
  wA[i]=new double*[nrad_becke];
  for(j=0;j<nrad_becke;j++)
  {
   wA[i][j]=new double[nang_becke];
   for(k=0;k<nang_becke;k++)
   {
    Point[0]=r_real_becke[j]*x_becke[k]+Rho.Cartesian_Coor[i][0];
    Point[1]=r_real_becke[j]*y_becke[k]+Rho.Cartesian_Coor[i][1];
    Point[2]=r_real_becke[j]*z_becke[k]+Rho.Cartesian_Coor[i][2];
    for(l=0;l<natom;l++)
    {
     // Distance r to R_B
     Diff_Point[0]=Point[0]-Rho.Cartesian_Coor[l][0];
     Diff_Point[1]=Point[1]-Rho.Cartesian_Coor[l][1];
     Diff_Point[2]=Point[2]-Rho.Cartesian_Coor[l][2];
     rB=norm3D(Diff_Point);
     ZB=(int)Rho.Nu_charge[l];
     for(m=l+1;m<natom;m++)
     {
      // Distance r to R_C
      Diff_Point[0]=Point[0]-Rho.Cartesian_Coor[m][0];
      Diff_Point[1]=Point[1]-Rho.Cartesian_Coor[m][1];
      Diff_Point[2]=Point[2]-Rho.Cartesian_Coor[m][2];
      rC=norm3D(Diff_Point);
      ZC=(int)Rho.Nu_charge[m];
      // Distance R_B to R_C
      Diff_Point[0]=Rho.Cartesian_Coor[m][0]-Rho.Cartesian_Coor[l][0];
      Diff_Point[1]=Rho.Cartesian_Coor[m][1]-Rho.Cartesian_Coor[l][1];
      Diff_Point[2]=Rho.Cartesian_Coor[m][2]-Rho.Cartesian_Coor[l][2];
      RBC=norm3D(Diff_Point);
      // Compute nu_BC
      mu_BC=(rB-rC)/RBC;
      // becke_original
      if(becke_original)
      {
       s_BC=smooth_stiff(mu_BC,stiff);
      }
      // Becke
      if(becke)
      {
       xi_BC=Xi_XY_table(ZB,ZC);
       u_BC=(xi_BC-ONE)/(xi_BC+ONE); 
       a_BC=u_BC/(u_BC*u_BC-ONE);
       if(a_BC<-HALF){a_BC=-HALF;}
       if(a_BC>HALF){a_BC=HALF;}
       nu_BC=mu_BC+a_BC*(ONE-mu_BC*mu_BC);
       s_BC=smooth_stiff(nu_BC,stiff);
      }
      // SSF Stratmann, Scuseria, Frisch, Chem. Phys. Lett. 257, 213 (1996)
      if(ssf)
      {
       if(mu_BC<-a_ssf){s_BC=-ONE;}
       else if(mu_BC>a_ssf){s_BC=ONE;}
       else
       {
        nu_BC=mu_BC/a_ssf;
        s_BC=(3.5e1*nu_BC-3.5e1*pow(nu_BC,3.0e0)+2.1e1*pow(nu_BC,5.0e0)-5.0e0*pow(nu_BC,7.0e0))/1.6e1;
       }
      }
      // TFVC
      if(tfvc)
      {
       xi_BC=Xi_XY_mat[l][m];  
       nu_BC=(ONE+mu_BC-xi_BC*(ONE-mu_BC))/(ONE+mu_BC+xi_BC*(ONE-mu_BC));
       s_BC=smooth_stiff(nu_BC,stiff);
      }
      S_X[m][l]=HALF*(ONE-s_BC);
      S_X[l][m]=HALF*(ONE+s_BC);
     } 
    }
    for(l=0;l<natom;l++)
    {
     P_X[l]=ONE;
     for(m=0;m<natom;m++)
     {
      P_X[l]=P_X[l]*S_X[m][l];
     }
    }
    Sum_PB=ZERO;
    for(l=0;l<natom;l++)
    {
     Sum_PB+=P_X[l];
    }
    wA[i][j][k]=P_X[i]/Sum_PB;
   }
  }
 }
 for(i=0;i<natom;i++)
 {
  delete[] S_X[i];S_X[i]=NULL;
  delete[] Xi_XY_mat[i];Xi_XY_mat[i]=NULL;
 }
 delete[] P_X;P_X=NULL;
 delete[] S_X;S_X=NULL;
 delete[] Xi_XY_mat;Xi_XY_mat=NULL;
}

void Integrate_atomic(READ_FCHK_WFN &Rho,double *res_integration)
{
 int i,j,k,nprops=7,natoms=Rho.natoms;
 double Point[3],normPoint,density,fact_jacob_weight,FIFTEEN=THREE*FIVE;
#ifdef HAVE_LIBXC
 nprops++;                           // Add Exc integration as an extra property.
 xc_func_type func;
 int vmajor,vminor,vmicro,func_id=1; // func_id=1 is LDAx
 double rho[1],e_xc[1];              // sigma[1]. No sigma vector needed because it is LDA example (LDAx Slater)
 xc_version(&vmajor, &vminor, &vmicro);
 printf("Libxc version: %d.%d.%d\n",vmajor,vminor,vmicro);
 if(xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0)
 {
  fprintf(stderr, "Functional '%d' not found\n", func_id);
 }
#endif
 // Molecular init.
 for(i=0;i<nprops*natoms+nprops;i++)
 {
  res_integration[i]=ZERO;
 }
 // Calc. integrals
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<nrad_becke;j++)
  { 
   for(k=0;k<nang_becke;k++)
   {
    Point[0]=r_real_becke[j]*x_becke[k]+Rho.Cartesian_Coor[i][0];
    Point[1]=r_real_becke[j]*y_becke[k]+Rho.Cartesian_Coor[i][1];
    Point[2]=r_real_becke[j]*z_becke[k]+Rho.Cartesian_Coor[i][2];
    normPoint=norm3D(Point);
    Rho.rho_eval(Point,density);
    fact_jacob_weight=wA[i][j][k]*w_theta_phi_becke[k]*w_radial_becke[j]*pow(r_real_becke[j],TWO)/pow(ONE-r_becke[j],TWO);
    res_integration[i*nprops]+=density*fact_jacob_weight;
    res_integration[i*nprops+1]+=density*fact_jacob_weight*Point[0];            // rho \times x
    res_integration[i*nprops+2]+=density*fact_jacob_weight*Point[1];            // rho \times y
    res_integration[i*nprops+3]+=density*fact_jacob_weight*Point[2];            // rho \times z
    res_integration[i*nprops+4]+=density*fact_jacob_weight*normPoint;           // rho \times r
    res_integration[i*nprops+5]+=density*fact_jacob_weight*normPoint*normPoint; // rho \times r^2
    if(density>pow(TEN,-FIFTEEN)) // Negative densities from transtion densities are omitted
    {
     res_integration[i*nprops+6]+=pow(density,FIVE/THREE)*fact_jacob_weight;    // rho^5/3 
#ifdef HAVE_LIBXC
     rho[0]=density;
     xc_lda_exc(&func,1,rho,e_xc);
     res_integration[i*nprops+7]+=density*e_xc[0]*fact_jacob_weight;             // rho \times exc^LDA[rho] 
#endif
    }
   }
  }
 }
#ifdef HAVE_LIBXC
 xc_func_end(&func);
#endif
 // Multiply by 4 Pi and sum atomic contrib. to  
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<nprops;j++)
  {
   // Times 4 Pi
   res_integration[i*nprops+j]=FOUR*PI*res_integration[i*nprops+j];   // For each atom
   // Sum atoms
   res_integration[nprops*natoms+j]+=res_integration[i*nprops+j]; // The whole molecule
  }
 }
}

void Integrate_atomic_paral(vector<READ_FCHK_WFN> Rho,double *res_integration,int &nprocs)
{
 int i,j,k,nprops=7,natoms=Rho[0].natoms;
 double Point[3],normPoint,density,fact_jacob_weight,FIFTEEN=THREE*FIVE;
 // Molecular init.
 for(i=0;i<nprops*natoms+nprops;i++)
 {
  res_integration[i]=ZERO;
 }
 // Calc. integrals
 #pragma omp parallel num_threads(nprocs) \
 private(i,j,k,Point,normPoint,density,fact_jacob_weight) \
 shared(r_becke,r_real_becke,x_becke,y_becke,z_becke,natoms,nrad_becke,nang_becke,nprops,\
 wA,w_theta_phi_becke,w_radial_becke)
 {
  int nth=omp_get_num_threads();
  int ith=omp_get_thread_num();
  double *res_integration_th;
  res_integration_th=new double[nprops*natoms+nprops];
  for(i=0;i<nprops*natoms+nprops;i++){res_integration_th[i]=ZERO;}
  for(i=ith;i<natoms;i=i+nth)
  {
   for(j=0;j<nrad_becke;j++)
   { 
    for(k=0;k<nang_becke;k++)
    {
     Point[0]=r_real_becke[j]*x_becke[k]+Rho[ith].Cartesian_Coor[i][0];
     Point[1]=r_real_becke[j]*y_becke[k]+Rho[ith].Cartesian_Coor[i][1];
     Point[2]=r_real_becke[j]*z_becke[k]+Rho[ith].Cartesian_Coor[i][2];
     normPoint=norm3D(Point);
     Rho[ith].rho_eval(Point,density);
     fact_jacob_weight=wA[i][j][k]*w_theta_phi_becke[k]*w_radial_becke[j]*pow(r_real_becke[j],TWO)/pow(ONE-r_becke[j],TWO);
     res_integration_th[i*nprops]+=density*fact_jacob_weight;
     res_integration_th[i*nprops+1]+=density*fact_jacob_weight*Point[0];            // rho \times x
     res_integration_th[i*nprops+2]+=density*fact_jacob_weight*Point[1];            // rho \times y
     res_integration_th[i*nprops+3]+=density*fact_jacob_weight*Point[2];            // rho \times z
     res_integration_th[i*nprops+4]+=density*fact_jacob_weight*normPoint;           // rho \times r
     res_integration_th[i*nprops+5]+=density*fact_jacob_weight*normPoint*normPoint; // rho \times r^2
     if(density>pow(TEN,-FIFTEEN)) // Negative densities from transtion densities are omitted
     {
      res_integration_th[i*nprops+6]+=pow(density,FIVE/THREE)*fact_jacob_weight;    // rho^5/3 
     }
    }
   }
  }
  #pragma omp barrier
  #pragma omp critical
  {
   for(i=0;i<nprops*natoms+nprops;i++){res_integration[i]+=res_integration_th[i];}
  } 
  delete[] res_integration_th;res_integration_th=NULL;
 }
 // Multiply by 4 Pi and sum atomic contrib. to  
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<nprops;j++)
  {
   // Times 4 Pi
   res_integration[i*nprops+j]=FOUR*PI*res_integration[i*nprops+j];   // For each atom
   // Sum atoms
   res_integration[nprops*natoms+j]+=res_integration[i*nprops+j]; // The whole molecule
  }
 }
}

void Integrate_atomic_sij(vector<READ_FCHK_WFN> Rho,double ***S_Aij,double *Density,int &nbasis,int &nprocs)
{
 int i,j,k,l,m,natoms=Rho[0].natoms;
 double Point[3],fact_jacob_weight,density;
 // Calc. integrals
 #pragma omp parallel num_threads(nprocs) \
 private(i,j,k,l,m,Point,density,fact_jacob_weight) \
 shared(r_becke,r_real_becke,x_becke,y_becke,z_becke,natoms,nrad_becke,nang_becke,nbasis,\
 wA,w_theta_phi_becke,w_radial_becke)
 {
  int nth=omp_get_num_threads();
  int ith=omp_get_thread_num();
  double ***S_Aij_th,*Density_th,*NOs,**NOs_grad;
  Density_th=new double[natoms];
  S_Aij_th=new double**[natoms];
  for(i=0;i<natoms;i++)
  {
   Density_th[i]=ZERO;
   S_Aij_th[i]=new double*[nbasis];
   for(j=0;j<nbasis;j++)
   {
    S_Aij_th[i][j]=new double[j+1];
    for(k=0;k<=j;k++){S_Aij_th[i][j][k]=ZERO;}
   }
  } 
  NOs=new double[nbasis];
  NOs_grad=new double*[3];
  for(i=0;i<3;i++){NOs_grad[i]=new double[nbasis];}
  // Calc. integrals
  for(i=ith;i<natoms;i=i+nth)
  {
   Density[i]=ZERO;
   for(j=0;j<nrad_becke;j++)
   { 
    for(k=0;k<nang_becke;k++)
    {
     Point[0]=r_real_becke[j]*x_becke[k]+Rho[ith].Cartesian_Coor[i][0];
     Point[1]=r_real_becke[j]*y_becke[k]+Rho[ith].Cartesian_Coor[i][1];
     Point[2]=r_real_becke[j]*z_becke[k]+Rho[ith].Cartesian_Coor[i][2];
     Rho[ith].build_NO_grad_wfn_all(NOs,NOs_grad,Point);
     fact_jacob_weight=wA[i][j][k]*w_theta_phi_becke[k]*w_radial_becke[j]*pow(r_real_becke[j],TWO)/pow(ONE-r_becke[j],TWO);
     density=ZERO;
     for(l=0;l<nbasis;l++)
     {
      if(Rho[ith].Ocupation[l]!=ZERO)
      {
       density=density+NOs[l]*NOs[l]*Rho[ith].Ocupation[l];
      }
      for(m=0;m<=l;m++)
      {
       S_Aij_th[i][l][m]+=NOs[l]*NOs[m]*fact_jacob_weight;    
      }
     }
     // Int. electronic density
     Density_th[i]+=density*fact_jacob_weight;
    }
   }
  }
  #pragma omp barrier
  #pragma omp critical
  {
   for(i=0;i<natoms;i++)
   {
    for(j=0;j<nbasis;j++)
    {
     for(k=0;k<=j;k++)
     {
      S_Aij[i][j][k]+=S_Aij_th[i][j][k];
     }
    }	    
    Density[i]+=Density_th[i];
   }
  } 
  for(i=0;i<natoms;i++)
  {
   for(j=0;j<nbasis;j++)
   {
    delete[] S_Aij_th[i][j];S_Aij_th[i][j]=NULL;
   }
   delete[] S_Aij_th[i];S_Aij_th[i]=NULL;
  } 
  for(i=0;i<3;i++)
  { 
   delete[] NOs_grad[i]; NOs_grad[i]=NULL;
  }
  delete[] S_Aij_th; S_Aij_th=NULL;
  delete[] Density_th; Density_th=NULL;
  delete[] NOs; NOs=NULL;
  delete[] NOs_grad; NOs_grad=NULL;
 }
 // Multiply by 4 Pi   
 for(i=0;i<natoms;i++)
 {
  for(l=0;l<nbasis;l++)
  {
   for(m=0;m<=l;m++)
   {
    S_Aij[i][l][m]=FOUR*PI*S_Aij[i][l][m];    
   }
  }
  Density[i]=FOUR*PI*Density[i];
 }
}

// Save serial version for backup
/*
void Integrate_atomic_sij(READ_FCHK_WFN &Rho,double ***S_Aij,double *Density,int &nbasis)
{
 int i,j,k,l,m,natoms=Rho.natoms;
 double Point[3],fact_jacob_weight,density;
 double *NOs,**NOs_grad;
 NOs=new double[nbasis];
 NOs_grad=new double*[3];
 for(i=0;i<3;i++)
 {
  NOs_grad[i]=new double[nbasis];
 }
 // Calc. integrals
 for(i=0;i<natoms;i++)
 {
  Density[i]=ZERO;
  for(j=0;j<nrad_becke;j++)
  { 
   for(k=0;k<nang_becke;k++)
   {
    Point[0]=r_real_becke[j]*x_becke[k]+Rho.Cartesian_Coor[i][0];
    Point[1]=r_real_becke[j]*y_becke[k]+Rho.Cartesian_Coor[i][1];
    Point[2]=r_real_becke[j]*z_becke[k]+Rho.Cartesian_Coor[i][2];
    Rho.build_NO_grad_wfn_all(NOs,NOs_grad,Point);
    fact_jacob_weight=wA[i][j][k]*w_theta_phi_becke[k]*w_radial_becke[j]*pow(r_real_becke[j],TWO)/pow(ONE-r_becke[j],TWO);
    density=ZERO;
    for(l=0;l<nbasis;l++)
    {
     if(Rho.Ocupation[l]!=ZERO)
     {
      density=density+NOs[l]*NOs[l]*Rho.Ocupation[l];
     }
     for(m=0;m<=l;m++)
     {
      S_Aij[i][l][m]+=NOs[l]*NOs[m]*fact_jacob_weight;    
     }
    }
    // Int. electronic density
    Density[i]+=density*fact_jacob_weight;
   }
  }
 }
 // Multiply by 4 Pi   
 for(i=0;i<natoms;i++)
 {
  for(l=0;l<nbasis;l++)
  {
   for(m=0;m<=l;m++)
   {
    S_Aij[i][l][m]=FOUR*PI*S_Aij[i][l][m];    
   }
  }
  Density[i]=FOUR*PI*Density[i];
 }
 for(i=0;i<3;i++)
 { 
  delete[] NOs_grad[i]; NOs_grad[i]=NULL;
 }
 delete[] NOs; NOs=NULL;
 delete[] NOs_grad; NOs_grad=NULL;
}
*/

void Integrate_atomic_mescal(READ_FCHK_WFN &Rho,double **F_QM,double *V_QM,double &Density,double **Coords_pdb,int &natoms_pdb)
{
 int i,j,k,l,m,natoms=Rho.natoms;
 double Point[3],diff_xyz[3],r,r3,fact_jacob_weight,density,tol10=pow(TEN,-TEN);
 // Calc. integrals
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<nrad_becke;j++)
  { 
   for(k=0;k<nang_becke;k++)
   {
    Point[0]=r_real_becke[j]*x_becke[k]+Rho.Cartesian_Coor[i][0];
    Point[1]=r_real_becke[j]*y_becke[k]+Rho.Cartesian_Coor[i][1];
    Point[2]=r_real_becke[j]*z_becke[k]+Rho.Cartesian_Coor[i][2];
    Rho.rho_eval(Point,density);
    fact_jacob_weight=wA[i][j][k]*w_theta_phi_becke[k]*w_radial_becke[j]*pow(r_real_becke[j],TWO)/pow(ONE-r_becke[j],TWO);
    for(l=0;l<natoms_pdb;l++)
    {
     r=ZERO;
     for(m=0;m<3;m++)
     {
      diff_xyz[m]=Coords_pdb[l][m]-Point[m];
      r+=diff_xyz[m]*diff_xyz[m];
     }
     r=pow(r,HALF);
     if(r<tol10){r=tol10;}
     r3=pow(r,THREE);
     // Field and potential contributions
     for(m=0;m<3;m++){F_QM[l][m]+=fact_jacob_weight*density*diff_xyz[m]/r3;} 
     V_QM[l]+=density*fact_jacob_weight/r;
    }
    // Int. electronic density
    Density+=density*fact_jacob_weight;
   }
  }
 }
 // Multiply by 4 Pi   
 for(i=0;i<natoms_pdb;i++)
 {
  for(j=0;j<3;j++)
  {
   F_QM[i][j]=FOUR*PI*F_QM[i][j];
  }
  V_QM[i]=FOUR*PI*V_QM[i];
 }
 Density=FOUR*PI*Density;
}

void Clean_quadrature_atomic(string name,int &natoms)
{
 int i,j;
 system(("rm "+name+"_w.txt").c_str());
 system(("rm "+name+"_x.txt").c_str());
 system(("rm "+name+"_r.txt").c_str());
 delete[] x_becke; x_becke=NULL;
 delete[] y_becke; y_becke=NULL;
 delete[] z_becke; z_becke=NULL;
 delete[] r_becke; r_becke=NULL;
 delete[] r_real_becke;r_real_becke=NULL;
// delete[] theta_becke; theta_becke=NULL;
// delete[] phi_becke; phi_becke=NULL;
 delete[] w_radial_becke; w_radial_becke=NULL;
 delete[] w_theta_phi_becke; w_theta_phi_becke=NULL;
 for(i=0;i<natoms;i++)
 {
  for(j=0;j<nrad_becke;j++)
  {
   delete[] wA[i][j];wA[i][j]=NULL;
  }
  delete[] wA[i];wA[i]=NULL;
 }
 delete[] wA;wA=NULL;
}

//Check for available order in Lebedev
void Grid_avail_atomic(int & Order)
{ 
 if(Order<=6)
 {
  Order =    6;
 }
 else if(Order<=14)
 {
  Order =   14;
 }
 else if(Order<=26)
 {
  Order =   26;
 }
 else if(Order<=38)
 {
  Order =   38;
 }
 else if(Order<=50)
 {
  Order =   50;
 }
 else if(Order<=74)
 {
  Order =   74;
 }
 else if(Order<=86)
 {
  Order =   86;
 }
 else if(Order<=110)
 {
  Order =  110;
 }
 else if(Order<=146)
 {
  Order =  146;
 }
 else if(Order<=170)
 {
  Order =  170;
 }
 else if(Order<=194)
 {
  Order =  194;
 }
 else if(Order<=230)
 {
  Order =  230;
 }
 else if(Order<=266)
 {
  Order =  266;
 }
 else if(Order<=302)
 {
  Order =  302;
 }
 else if(Order<=350)
 {
  Order =  350;
 }
 else if(Order<=434)
 {
  Order =  434;
 }
 else if(Order<=590)
 {
  Order =  590;
 }
 else if(Order<=770)
 {
  Order =  770;
 }
 else if(Order<=974)
 {
  Order =  974;
 }
 else if(Order<=1202)
 {
  Order = 1202;
 }
 else if(Order<=1454)
 {
  Order = 1454;
 }
 else if(Order<=1730)
 {
  Order = 1730;
 }
 else if(Order<=2030)
 {
  Order = 2030;
 }
 else if(Order<=2354)
 {
  Order = 2354;
 }
 else if(Order<=2702)
 {
  Order = 2702;
 }
 else if(Order<=3074)
 {
  Order = 3074;
 }
 else if(Order<=3470)
 {
  Order = 3470;
 }
 else if(Order<=3890)
 {
  Order = 3890;
 }
 else if(Order<=4334)
 {
  Order = 4334;
 }
 else if(Order<=4802)
 {
  Order = 4802;
 }
 else if(Order<=5294)
 {
  Order = 5294;
 }
 else
 {
  Order = 5810;
 }
}

void Xi_XY_bcp(READ_FCHK_WFN &Rho,string name,int &iatom,int &jatom,int &natom,double &Xi_rad)
{
 int icoord,katom,iter=0,info;
 double R_dist,R_dist2,rX=ONE,rY=ONE,Diff_Point[3],Diff_Point2[3],Middle[3],a,b,x1[3],x2[3],x3[3],x4[3],f1,f2,f3,f4;
 double tol8=pow(TEN,-EIGHT),tol6=pow(TEN,-SIX),norm_grad,delta;
 double **Hess,**Inv_Hess,**Eigenvect,*Grad;
 for(icoord=0;icoord<3;icoord++)
 {
  Diff_Point[icoord]=Rho.Cartesian_Coor[iatom][icoord]-Rho.Cartesian_Coor[jatom][icoord];
 }
 R_dist=norm3D(Diff_Point);
 for(icoord=0;icoord<3;icoord++)
 {
  Middle[icoord]=Rho.Cartesian_Coor[jatom][icoord]+HALF*Diff_Point[icoord]/(R_dist+tol8); 
 }
 for(katom=0;katom<natom;katom++)
 {
  if(katom!=iatom && katom!=jatom)
  {
   for(icoord=0;icoord<3;icoord++)
   {
    Diff_Point2[icoord]=Rho.Cartesian_Coor[katom][icoord]-Middle[icoord];
   }
   R_dist2=norm3D(Diff_Point2);
   if(R_dist2<HALF*R_dist)
   { 
    R_dist=TEN;
    katom=natom;
   }
  }
 }	 
 ofstream write_bcp(name,std::ios_base::app);
 if(R_dist<SIX && iatom!=jatom) // 6 au ( >3 Angstrom ) to search for BCPs
 {
  write_bcp<<endl;
  // Find CP, set radious rX and rY
  for(icoord=0;icoord<3;icoord++)
  {
   x1[icoord]=Rho.Cartesian_Coor[iatom][icoord];
   x3[icoord]=Rho.Cartesian_Coor[jatom][icoord];
  }
  Rho.rho_eval(x1,f1);
  Rho.rho_eval(x3,f3);
  do
  {
   for(icoord=0;icoord<3;icoord++)
   {
    Diff_Point[icoord]=x3[icoord]-x1[icoord];
   }
   R_dist=norm3D(Diff_Point);
   // a + b = R_dist
   // b/a = phi = 1.618033899
   // then
   // a + a phi = a ( 1 + phi ) = R_dist -> a = R_dist / ( 1 + phi )
   if(R_dist>tol6)
   { 
    a=R_dist/(golden_ratio+ONE);
    b=R_dist-a;
    for(icoord=0;icoord<3;icoord++)
    {
     x2[icoord]=x1[icoord]+a*Diff_Point[icoord]/R_dist;
     x4[icoord]=x1[icoord]+b*Diff_Point[icoord]/R_dist;
    }
    Rho.rho_eval(x2,f2);
    Rho.rho_eval(x4,f4);
    if(f4<f2)
    {
     for(icoord=0;icoord<3;icoord++)
     {
      x1[icoord]=x2[icoord];
     }
     f1=f2;
    }
    else
    {
     for(icoord=0;icoord<3;icoord++)
     {
      x3[icoord]=x4[icoord];
     }
     f3=f4;
    }
   }
   iter++;   
  }while((abs(f1-f3)>tol8 || R_dist>tol6) && iter<1000);
  if(iter==1000){write_bcp<<" Warning! Failed to find min between atoms "<<iatom+1<<","<<jatom<<"."<<endl;}
  Hess=new double*[3];
  Inv_Hess=new double*[3];
  Eigenvect=new double*[3];
  Grad=new double[3];
  for(icoord=0;icoord<3;icoord++)
  {
   Hess[icoord]=new double[3];
   Inv_Hess[icoord]=new double[3];
   Eigenvect[icoord]=new double[3];
  }
  for(icoord=0;icoord<3;icoord++)
  {
   x1[icoord]=HALF*(x1[icoord]+x3[icoord]);
   x2[icoord]=x1[icoord];
  }
  write_bcp<<setprecision(10)<<fixed<<scientific;
  write_bcp<<"Terminated the search for a minimun along the internuclear axis"<<endl;
  write_bcp<<"Starting the search for a BCP between atoms "<<iatom+1<<" and "<<jatom+1<<endl;
  iter=0;delta=ONE;
  do
  {
   info=0;
   Rho.rho_hessian(x1,Hess,Grad,f1);
   mat_inverse2(3,Hess,Inv_Hess,info);
   if(info==0)
   {
    norm_grad=pow(Grad[0]*Grad[0]+Grad[1]*Grad[1]+Grad[2]*Grad[2],HALF);
    if(iter==0)
    {
     write_bcp<<"Current values for the density, grad, and laplacian"<<endl;
     write_bcp<<"Density     : "<<setw(17)<<f1<<endl;
     write_bcp<<"|Grad|      : "<<setw(17)<<norm_grad<<endl;
     jacobi(3,Hess,Eigenvect);
     write_bcp<<"Laplacian   : "<<setw(17)<<Hess[0][0]+Hess[1][1]+Hess[2][2]<<endl;
     write_bcp<<"Eigenval_1  : "<<setw(17)<<Hess[0][0]<<endl;
     write_bcp<<"Eigenval_2  : "<<setw(17)<<Hess[1][1]<<endl;
     write_bcp<<"Eigenval_3  : "<<setw(17)<<Hess[2][2]<<endl;
    }
    x1[0]=x1[0]-delta*(Inv_Hess[0][0]*Grad[0]+Inv_Hess[0][1]*Grad[1]+Inv_Hess[0][2]*Grad[2]); 
    x1[1]=x1[1]-delta*(Inv_Hess[1][0]*Grad[0]+Inv_Hess[1][1]*Grad[1]+Inv_Hess[1][2]*Grad[2]); 
    x1[2]=x1[2]-delta*(Inv_Hess[2][0]*Grad[0]+Inv_Hess[2][1]*Grad[1]+Inv_Hess[2][2]*Grad[2]); 
    for(icoord=0;icoord<3;icoord++){x3[icoord]=x2[icoord]-x1[icoord];}
    if(norm3D(x3)>ONE_THIRD){iter=50;norm_grad=pow(TEN,TEN);}
    iter++;
   }
   else{iter=50;norm_grad=pow(TEN,TEN);}  
  }while(iter<50 && norm_grad>tol6);
  for(icoord=0;icoord<3;icoord++){x2[icoord]=x2[icoord]-x1[icoord];}
  R_dist=norm3D(x2);
  if(norm_grad<tol6 && iter<50 && R_dist<ONE_THIRD)
  { 
   for(icoord=0;icoord<3;icoord++){if(abs(x1[icoord])<tol8){x1[icoord]=ZERO;}}
   Rho.rho_hessian(x1,Hess,Grad,f1);
   write_bcp<<"Final values for the density, grad, and laplacian after "<<iter<<" iterations."<<endl;
   write_bcp<<"BCP  x_coord: "<<setw(17)<<x1[0]<<endl;
   write_bcp<<"BCP  y_coord: "<<setw(17)<<x1[1]<<endl;
   write_bcp<<"BCP  z_coord: "<<setw(17)<<x1[2]<<endl;
   write_bcp<<"Density     : "<<setw(17)<<f1<<endl;
   write_bcp<<"|Grad|      : "<<setw(17)<<norm_grad<<endl;
   jacobi(3,Hess,Eigenvect);
   write_bcp<<"Laplacian   : "<<setw(17)<<Hess[0][0]+Hess[1][1]+Hess[2][2]<<endl;
   write_bcp<<"Eigenval_1  : "<<setw(17)<<Hess[0][0]<<endl;
   write_bcp<<"Eigenval_2  : "<<setw(17)<<Hess[1][1]<<endl;
   write_bcp<<"Eigenval_3  : "<<setw(17)<<Hess[2][2]<<endl;
   for(icoord=0;icoord<3;icoord++)
   {
    x2[icoord]=x1[icoord]-Rho.Cartesian_Coor[iatom][icoord];
    x4[icoord]=x1[icoord]-Rho.Cartesian_Coor[jatom][icoord];
   }
   rX=norm3D(x2);
   rY=norm3D(x4);
   n_bcps++;
  }
  else{rX=ONE,rY=ONE;write_bcp<<"No BCP found between atoms "<<iatom+1<<" and "<<jatom+1<<endl;}
  for(icoord=0;icoord<3;icoord++)
  {
   delete[] Hess[icoord];Hess[icoord]=NULL;
   delete[] Inv_Hess[icoord];Inv_Hess[icoord]=NULL;
   delete[] Eigenvect[icoord];Eigenvect[icoord]=NULL;
  }
  delete[] Hess;Hess=NULL;
  delete[] Inv_Hess;Inv_Hess=NULL;
  delete[] Eigenvect;Eigenvect=NULL;
  delete[] Grad;Grad=NULL;
 }
 write_bcp.close();
 Xi_rad=rX/rY;
} 

double Xi_XY_table(int &Z1, int &Z2)
{
 return set_radii(Z1)/set_radii(Z2);
} 

// Radii taken from APOST-3D code
double set_radii(int &Z)
{
 double val=ZERO;
 if(Z==1){val=0.327;}
 else if(Z==2){val= 1.000;}
 else if(Z==3){val= 1.219;}
 else if(Z==4){val= 0.911;}
 else if(Z==5){val= 0.793;}
 else if(Z==6){val= 0.766;}
 else if(Z==7){val= 0.699;}
 else if(Z==8){val= 0.658;}
 else if(Z==9){val= 0.900;}
 else if(Z==10){val=1.000;}
 else if(Z==11){val=1.545;}
 else if(Z==12){val=1.333;}
 else if(Z==13){val=1.199;}
 else if(Z==14){val=1.123;}
 else if(Z==15){val=1.110;}
 else if(Z==16){val=1.071;}
 else if(Z==17){val=1.039;}
 else if(Z==18){val=1.000;}
 else if(Z==19){val=1.978;}
 else if(Z==20){val=1.745;}
 else if(Z==21){val=1.337;}
 else if(Z==22){val=1.274;}
 else if(Z==23){val=1.236;}
 else if(Z==24){val=1.128;}
 else if(Z==25){val=1.180;}
 else if(Z==26){val=1.091;}
 else if(Z==27){val=1.089;}
 else if(Z==28){val=1.077;}
 else if(Z==29){val=1.146;}
 else if(Z==30){val=1.187;}
 else if(Z==31){val=1.199;}
 else if(Z==32){val=1.179;}
 else if(Z==33){val=1.209;}
 else if(Z==34){val=1.201;}
 else if(Z==35){val=1.201;}
 else if(Z==36){val=1.000;}
 else if(Z==37){val=2.217;}
 else if(Z==38){val=1.928;}
 else if(Z==39){val=1.482;}
 else if(Z==40){val=1.377;}
 else if(Z==41){val=1.353;}
 else if(Z==42){val=1.240;}
 else if(Z==43){val=1.287;}
 else if(Z==44){val=1.212;}
 else if(Z==45){val=1.229;}
 else if(Z==46){val=1.240;}
 else if(Z==47){val=1.362;}
 else if(Z==48){val=1.429;}
 else if(Z==49){val=1.385;}
 else if(Z==50){val=1.380;}
 else if(Z==51){val=1.421;}
 else if(Z==52){val=1.400;}
 else if(Z==53){val=1.397;}
 else if(Z==54){val=1.000;}
 else if(Z==55){val=2.442;}
 else if(Z==56){val=2.149;}
 else if(Z==57){val=1.653;}
 else if(Z==58){val=1.500;}
 else if(Z==59){val=1.500;}
 else if(Z==60){val=1.500;}
 else if(Z==61){val=1.500;}
 else if(Z==62){val=1.500;}
 else if(Z==63){val=1.500;}
 else if(Z==64){val=1.500;}
 else if(Z==65){val=1.500;}
 else if(Z==66){val=1.500;}
 else if(Z==67){val=1.500;}
 else if(Z==68){val=1.500;}
 else if(Z==69){val=1.500;}
 else if(Z==70){val=1.500;}
 else if(Z==71){val=1.500;}
 else if(Z==72){val=1.364;}
 else if(Z==73){val=1.346;}
 else if(Z==74){val=1.256;}
 else if(Z==75){val=1.258;}
 else if(Z==76){val=1.222;}
 else if(Z==77){val=1.227;}
 else if(Z==78){val=1.227;}
 else if(Z==79){val=1.273;}
 else if(Z==80){val=1.465;}
 else if(Z==81){val=1.531;}
 else if(Z==82){val=1.434;}      
 else if(Z==83){val=1.496;}
 else if(Z==84){val=1.500;}
 else if(Z==85){val=1.500;}	
 else if(Z==86){val=1.500;}      
 else if(Z==87){val=1.500;}      		
 else if(Z==88){val=1.500;}
 else if(Z==89){val=1.500;}
 else if(Z==90){val=1.500;}
 else if(Z==91){val=1.500;}
 else if(Z==92){val=1.500;}
 else
 {
  cout<<"Warning! Atom with Z "<<Z<<" does not have a Bohr radius available."<<endl;
  cout<<"Bohr radius set to 1.60."<<endl;
  val=1.60;
 }
 return val;
}

double smooth_stiff(double &mu,int stiff)
{
 int i;
 double val=mu;
 for(i=0;i<stiff;i++)
 {
  val=p_mu(val);
 }
 return val;
}

double p_mu(double &mu)
{
 return HALF*mu*(THREE-pow(mu,TWO));
}

// See J. C. Slater, J. Chem. Phys., 41, 3199 (1964)
/*
double set_radii(int &Z)
{
 double val=ZERO;
 if(Z==1){val=0.25;}
 else if(Z==2){val= 0.85;} // Avg 
 else if(Z==3){val= 1.45;}    
 else if(Z==4){val= 1.05;}
 else if(Z==5){val= 0.85;} 
 else if(Z==6){val= 0.70;} 
 else if(Z==7){val= 0.65;} 
 else if(Z==8){val= 0.60;} 
 else if(Z==9){val= 0.50;} 
 else if(Z==10){val=1.15;} // Avg 
 else if(Z==11){val=1.80;}
 else if(Z==12){val=1.50;}
 else if(Z==13){val=1.25;}
 else if(Z==14){val=1.10;} 
 else if(Z==15){val=1.00;} 
 else if(Z==16){val=1.00;} 
 else if(Z==17){val=1.00;} 
 else if(Z==18){val=1.60;} // Avg 
 else if(Z==19){val=2.20;}
 else if(Z==20){val=1.80;}
 else if(Z==21){val=1.60;}
 else if(Z==22){val=1.40;}
 else if(Z==23){val=1.35;}
 else if(Z==24){val=1.40;}
 else if(Z==25){val=1.40;}
 else if(Z==26){val=1.40;} 
 else if(Z==27){val=1.35;}
 else if(Z==28){val=1.35;}
 else if(Z==29){val=1.35;}
 else if(Z==30){val=1.35;}
 else if(Z==31){val=1.30;}
 else if(Z==32){val=1.25;}
 else if(Z==33){val=1.15;}
 else if(Z==34){val=1.15;}
 else if(Z==35){val=1.15;} 	
 else if(Z==36){val=1.75;} // Avg
 else if(Z==37){val=2.35;}
 else if(Z==38){val=2.00;}
 else if(Z==39){val=1.80;}
 else if(Z==40){val=1.55;}
 else if(Z==41){val=1.45;}
 else if(Z==42){val=1.45;}
 else if(Z==43){val=1.35;}
 else if(Z==44){val=1.30;}
 else if(Z==45){val=1.35;}
 else if(Z==46){val=1.40;}
 else if(Z==47){val=1.60;}
 else if(Z==48){val=1.55;}
 else if(Z==49){val=1.55;}
 else if(Z==50){val=1.45;}
 else if(Z==51){val=1.45;}
 else if(Z==52){val=1.40;}
 else if(Z==53){val=1.40;}
 else if(Z==54){val=2.00;} // Avg
 else if(Z==55){val=2.60;}
 else if(Z==56){val=2.15;}
 else if(Z==57){val=1.95;}
 else if(Z==58){val=1.85;}
 else if(Z==59){val=1.85;}
 else if(Z==60){val=1.85;}
 else if(Z==61){val=1.85;}
 else if(Z==62){val=1.85;} 
 else if(Z==63){val=1.85;}
 else if(Z==64){val=1.80;}
 else if(Z==65){val=1.75;}
 else if(Z==66){val=1.75;}
 else if(Z==67){val=1.75;}
 else if(Z==68){val=1.75;}
 else if(Z==69){val=1.75;}
 else if(Z==70){val=1.75;}
 else if(Z==71){val=1.75;}
 else if(Z==72){val=1.55;}
 else if(Z==73){val=1.45;}
 else if(Z==74){val=1.35;} 
 else if(Z==75){val=1.35;}
 else if(Z==76){val=1.30;}
 else if(Z==77){val=1.35;}
 else if(Z==78){val=1.35;}
 else if(Z==79){val=1.35;}
 else if(Z==80){val=1.50;}
 else if(Z==81){val=1.90;}
 else if(Z==82){val=1.75;} // Avg	
 else if(Z==83){val=1.60;}
 else if(Z==84){val=1.90;}
 else if(Z==85){val=2.03;} // Avg 	
 else if(Z==86){val=2.03;} // Avg	
 else if(Z==87){val=2.03;} // Avg			
 else if(Z==88){val=2.15;}
 else if(Z==89){val=1.95;}
 else if(Z==90){val=1.80;}
 else if(Z==91){val=1.80;}
 else if(Z==92){val=1.75;}
 else if(Z==93){val=1.75;}
 else if(Z==94){val=1.75;}
 else if(Z==95){val=1.75;}
 else
 {
  cout<<"Warning! Atom with Z "<<Z<<" does not have a Bohr radius available."<<endl;
  cout<<"Bohr radius set to 1.80."<<endl;
  val=1.80;
 }
 return val;
}
*/

