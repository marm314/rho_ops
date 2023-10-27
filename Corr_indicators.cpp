#include"Corr_indicators.h"

//Calculate deviation from idempotency
double Deviation_idemp(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double I2=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {I2=ZERO;}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i]/TWO,TWO);}
     I2=ONE-I2/((double)Read_fchk_wfn.nelectrons*HALF);
    }
    else
    {
     cout<<"Warning! Unable to compute I2 open-shell correlated wfn file!"<<endl;
     I2=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i]/TWO,TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons*HALF);
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i],TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons);
   }
  }
  if(I2<pow(TEN,-FIVE)){I2=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i],TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons);
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++){I2=I2+pow(Read_fchk_wfn.Ocupation[i]/TWO,TWO);}
    I2=ONE-I2/((double)Read_fchk_wfn.nelectrons*HALF);
   }
  }
  else
  {
   cout<<"Warning! Unable to compute I2 without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(I2 will be 0.0e0)."<<endl;
   I2=ZERO;
  }
  if(I2<pow(TEN,-FIVE)){I2=ZERO;}
 }
 return I2;
}
//Calculate ID
double ID_ni(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double ID=ZERO,occ;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {ID=ZERO;}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
      ID=ID-HALF*abs(occ)*(ONE-abs(occ));
     }
     ID=TWO*ID;
    }
    else
    {
     cout<<"Warning! Unable to compute ID open-shell correlated wfn file!"<<endl;
     ID=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
    ID=TWO*ID;
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i])<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     if(abs(Read_fchk_wfn.Ocupation[i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
    ID=ID;
   }
  }
  if(ID<pow(TEN,-FIVE)){ID=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i])<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     if(abs(Read_fchk_wfn.Ocupation[i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i];}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
     ID=ID-HALF*abs(occ)*(ONE-abs(occ));
    }
    ID=ID*TWO;
   }
   if(ID<pow(TEN,-FIVE)){ID=ZERO;}
  }
  else
  {
   cout<<"Warning!  Unable to compute ID without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(ID will be 0.0e0)."<<endl;
   ID=ZERO;
  }
 }
 return ID;
}
//Calculate IND
double IND_ni(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double IND=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {IND=ZERO;}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      IND=IND+HALF*(Read_fchk_wfn.Ocupation[i]/TWO)*(ONE-Read_fchk_wfn.Ocupation[i]/TWO);
     }
     IND=IND*TWO;
    }
    else
    {
     cout<<"Warning! Unable to compute IND open-shell correlated wfn file!"<<endl;
     IND=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*(Read_fchk_wfn.Ocupation[i]/TWO)*(ONE-Read_fchk_wfn.Ocupation[i]/TWO);
    }
    IND=IND*TWO;
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*(Read_fchk_wfn.Ocupation[i])*(ONE-Read_fchk_wfn.Ocupation[i]);
    }
    IND=IND;
   }
  }
  if(IND<pow(TEN,-FIVE)){IND=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]);
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     IND=IND+HALF*(Read_fchk_wfn.Ocupation[i]/TWO)*(ONE-Read_fchk_wfn.Ocupation[i]/TWO);
    }
    IND=IND*TWO;
   }
   if(IND<pow(TEN,-FIVE)){IND=ZERO;}
  }
  else
  {
   cout<<"Warning! Unable to compute IND without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(IND will be 0.0e0)."<<endl;
   IND=ZERO;
  }
 }
 return IND;
}
//Calculate Shannon of ocupations
double Shannon_ni(READ_FCHK_WFN &Read_fchk_wfn)
{
 int i;
 double S=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {S=(double)Read_fchk_wfn.nelectrons*log((double)Read_fchk_wfn.nelectrons);}
   else
   {
    if(!Read_fchk_wfn.open_shell)
    {
     for(i=0;i<Read_fchk_wfn.nbasis();i++)
     {
      if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SIX))
      {
       S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons*TWO));
      }
      else
      {}
     }
    }
    else
    {
     cout<<"Warning! Unable to compute S(ni) open-shell correlated wfn file!"<<endl;
     S=ZERO;
    }
   }
  }
  else
  {
   if(!Read_fchk_wfn.open_shell)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SIX))
     {
      S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons*TWO));
     }
     else
     {}
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SIX))
     {
      S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons));
     }
     else
     {}
    }
   }
  }
  if(S<ZERO)
  {S=ZERO;}
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SEVEN))
     {
      S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/(double)Read_fchk_wfn.nelectrons);
     }
     else
     {}
    }
   }
   else
   {
    for(i=0;i<Read_fchk_wfn.nbasis();i++)
    {
     if(Read_fchk_wfn.Ocupation[i]>=pow(TEN,-SEVEN))
     {S=S-Read_fchk_wfn.Ocupation[i]*log(Read_fchk_wfn.Ocupation[i]/((double)Read_fchk_wfn.nelectrons*TWO));}
     else
     {}
    }
   }
  }
  else
  {
   cout<<"Warning! Unable to compute S(ni) without the S^1/2 P S^1/2 matrix"<<endl;
   cout<<"(S(ni) will be 0.0e0)."<<endl;
   S=ZERO;
  }
  if(S<ZERO)
  {S=ZERO;}
 }
 return S/(double)Read_fchk_wfn.nelectrons;
}
//Calculate deviation from idempotency from dmn
double Deviation_idemp_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double I2=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  I2=I2+pow(DMN.rho_matrixa[i][i],TWO)+pow(DMN.rho_matrixb[i][i],TWO);
 }
 I2=ONE-I2/N;
 return I2;
}
//Calculate ID
double ID_ni_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double occ,ID=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  if(abs(DMN.rho_matrixa[i][i])<pow(TEN,-SIX)){occ=ZERO;}
  else{occ=DMN.rho_matrixa[i][i];}
  if(abs(DMN.rho_matrixa[i][i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
  else{occ=DMN.rho_matrixa[i][i];}
  ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
  ID=ID-HALF*abs(occ)*(ONE-abs(occ));
  if(abs(DMN.rho_matrixb[i][i])<pow(TEN,-SIX)){occ=ZERO;}
  else{occ=DMN.rho_matrixb[i][i];}
  if(abs(DMN.rho_matrixb[i][i]-ONE)<pow(TEN,-SIX)){occ=ONE;}
  else{occ=DMN.rho_matrixb[i][i];}
  ID=ID+HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF);
  ID=ID-HALF*abs(occ)*(ONE-abs(occ));
 }
 return ID;
}
//Calculate IND
double IND_ni_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double IND=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  IND=IND+HALF*DMN.rho_matrixa[i][i]*(ONE-DMN.rho_matrixa[i][i]);
  IND=IND+HALF*DMN.rho_matrixb[i][i]*(ONE-DMN.rho_matrixb[i][i]);
 }
 return IND;
}
//Calculate Shannon of ocupations
double Shannon_ni_dmn(DMN_OPS &DMN, double &N)
{
 int i;
 double S=ZERO;
 for(i=0;i<DMN.nbasis();i++)
 {
  if(DMN.rho_matrixa[i][i]>=pow(TEN,-SEVEN))
  {
   S=S-DMN.rho_matrixa[i][i]*log(DMN.rho_matrixa[i][i]/N);
  }
  else
  {}
  if(DMN.rho_matrixb[i][i]>=pow(TEN,-SEVEN))
  {
   S=S-DMN.rho_matrixb[i][i]*log(DMN.rho_matrixb[i][i]/N);
  }
  else
  {}
 }
 return S/N;
}

//Local version for Nondynamic and Dynamic corr indcators
void ID_IND_local(READ_FCHK_WFN &Read_fchk_wfn,double Point[3],double &ID_alpha,double &ID_beta,double &IND_alpha,double &IND_beta,
bool &firstcall)
{
 int nbasis,i,j;
 double **NO_orb_grad,occ;
 nbasis=Read_fchk_wfn.nbasis();
 NO_orb_grad=new double*[4];
 for(i=0;i<4;i++)
 {NO_orb_grad[i]=new double[nbasis];}
 for(i=0;i<4;i++)
 {
  for(j=0;j<nbasis;j++)
  {NO_orb_grad[i][j]=ZERO;}
 }
////////////////////////////////////////////
 if(Read_fchk_wfn.overlap || Read_fchk_wfn.wfn){Read_fchk_wfn.orb_grad(Point,NO_orb_grad);}//Get NOs and NO gradients
//////////////////////////////////////////
 if((Read_fchk_wfn.correlated && Read_fchk_wfn.open_shell) && (Read_fchk_wfn.wfn && firstcall))
 {cout<<"Warning open-shell correlated wfn file!"<<endl;}
 ID_alpha=ZERO;
 ID_beta=ZERO;
 IND_alpha=ZERO;
 IND_beta=ZERO;
 if(Read_fchk_wfn.wfn)
 {
  if(!Read_fchk_wfn.wfx)
  {
   if(!Read_fchk_wfn.correlated)
   {
    ID_alpha=ZERO;
    ID_beta=ZERO;
    IND_alpha=ZERO;
    IND_beta=ZERO;
   }
   else
   {
    if(Read_fchk_wfn.open_shell)
    {
     if(firstcall)
     {
      cout<<"Warning! Unable to compute local indicators with open-shell"<<endl;
      cout<<"correlated wfn file!"<<endl;
      firstcall=false;
     }
     ID_alpha=ZERO;
     ID_beta=ZERO;
     IND_alpha=ZERO;
     IND_beta=ZERO;
    }
    else
    {
     for(i=0;i<nbasis;i++)
     {
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
      else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
      ID_alpha=ID_alpha+(HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF)-HALF*abs(occ)*(ONE-abs(occ)))
              *NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]/TWO)
               *NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
     }
     ID_beta=ID_alpha;
     IND_beta=IND_alpha;
    }
   }
  }
  else
  {
   if(Read_fchk_wfn.open_shell)
   {
    for(i=0;i<nbasis;i++)
    {
     if(i%2==0)
     {
      ID_alpha=ID_alpha+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
              -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
               *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
     else
     {
      ID_beta=ID_beta+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
             -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_beta=IND_beta+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
              *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
    }
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID_alpha=ID_alpha+(HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF)-HALF*abs(occ)*(ONE-abs(occ)))
             *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]/TWO)
              *NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
    }
    ID_beta=ID_alpha;
    IND_beta=IND_alpha;
   }
  }
 }
 else
 {
  if(Read_fchk_wfn.overlap)
  {
   if(Read_fchk_wfn.uhf)
   {
    for(i=0;i<nbasis;i++)
    {
     if(i%2==0)
     {
      ID_alpha=ID_alpha+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
              -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
               *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
     else
     {
      ID_beta=ID_beta+(HALF*HALF*pow(abs(Read_fchk_wfn.Ocupation[i])*abs(ONE-Read_fchk_wfn.Ocupation[i]),HALF)
             -HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]))*NO_orb_grad[0][i]*NO_orb_grad[0][i];
      IND_beta=IND_beta+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i])
              *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     }
    }
   }
   else
   {
    for(i=0;i<nbasis;i++)
    {
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO)<pow(TEN,-SIX)){occ=ZERO;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     if(abs(Read_fchk_wfn.Ocupation[i]/TWO-ONE)<pow(TEN,-SIX)){occ=ONE;}
     else{occ=Read_fchk_wfn.Ocupation[i]/TWO;}
     ID_alpha=ID_alpha+(HALF*HALF*pow(abs(occ),HALF)*pow(abs(ONE-abs(occ)),HALF)-HALF*abs(occ)*(ONE-abs(occ)))
             *NO_orb_grad[0][i]*NO_orb_grad[0][i];
     IND_alpha=IND_alpha+HALF*Read_fchk_wfn.Ocupation[i]*(ONE-Read_fchk_wfn.Ocupation[i]/TWO)
              *NO_orb_grad[0][i]*NO_orb_grad[0][i]/TWO;
    }
    ID_beta=ID_alpha;
    IND_beta=IND_alpha;
   }
  }
  else if(firstcall)
  {
   cout<<"Warning! Unable to compute local indicators without the S^1/2 P S^1/2 matrix"<<endl;
   firstcall=false;
   ID_alpha=ZERO;
   ID_beta=ZERO;
  }
  else
  {
   firstcall=false;
   ID_alpha=ZERO;
   ID_beta=ZERO;
   IND_alpha=ZERO;
   IND_beta=ZERO;
  }
 }
 //delete dynamic arrays
 for(i=0;i<4;i++)
 {delete[] NO_orb_grad[i];}
 delete[] NO_orb_grad;
}

