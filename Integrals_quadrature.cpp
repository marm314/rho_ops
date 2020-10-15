#include"Integrals_quadrature.h"
//theta from 0 to PI,
//phi from -PI  to PI
//r and r_real are not position space are RADIAL
double *r,*r_real,*w_radial,*w_theta_phi,*theta,*phi,*XX,*YY,*ZZ;
/////////////////////////////////////
//Quantities dependent of densities//
/////////////////////////////////////
//Generate entire grid for position space dependent densities
void integrate_quadrature(void *data,string name,bool dmn,int order,int order2,bool fisher,double **density_evals,double **density_grads,int mode)
{
 int i,j;
 double r_inf,r_sup,*x,*y,*z,Grad[3];
 double Point[3];
 ifstream wr_read,r_read;
 r=new double[order];
 r_real=new double[order];
 w_radial=new double[order];
 if(mode!=3)
 {
  //Grid for r
  r_inf=ZERO;
  r_sup=ONE;
  legendre_quadrature(name,order,r_inf,r_sup);
  wr_read.open((name+"_w.txt").c_str());
  r_read.open((name+"_x.txt").c_str());
  for(i=0;i<order;i++)
  {
   r_read>>r[i];
   wr_read>>w_radial[i];
  }
  wr_read.close();
  r_read.close();
 }
 else
 {
  r[0]=ZERO;
  w_radial[0]=ONE;
  ofstream nn((name+"_w.txt").c_str());
  ofstream nn1((name+"_r.txt").c_str());
  ofstream nn2((name+"_x.txt").c_str());
  nn.close();
  nn1.close();
  nn2.close();
 }
 for(i=0;i<order;i++)
 {
  //r_real is for entire space (0 to inf)
  r_real[i]=r[i]/(ONE-r[i]);
 }
 //Grid for theta, phi
 grid_avail(order2);
 x=new double[order2];
 y=new double[order2];
 z=new double[order2];
 if(mode!=0)
 {
  XX=new double[order2];
  YY=new double[order2];
  ZZ=new double[order2];
 }
 w_theta_phi=new double[order2];
 theta=new double[order2];
 phi=new double[order2];
 if(mode!=3)
 {
  ld_by_order(order2,x,y,z,w_theta_phi);
 }
 else
 {
  x[0]=ZERO;y[0]=ZERO;z[0]=ZERO;w_theta_phi[0]=ONE; 
 }
 double *temp,*temp1;
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order2;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  if(mode!=0)
  {
   XX[i]=x[i];
   YY[i]=y[i];
   ZZ[i]=z[i];
  }
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 if(mode==0)
 {
  //Loop over radial points
  for(i=0;i<order;i++)
  {//Loop over angular points
   for(j=0;j<order2;j++)
   {
    Point[0]=r_real[i]*x[j];
    Point[1]=r_real[i]*y[j];
    Point[2]=r_real[i]*z[j];
    if(!dmn)
    {
     ((READ_FCHK_WFN*)data)[0].rho_eval(Point,density_evals[i][j]);
     if(fisher)
     {
      ((READ_FCHK_WFN*)data)[0].rho_grad(Point,Grad);
      density_grads[i][j]=pow(norm3D(Grad),TWO);
     }
    }
    else
    {
     density_evals[i][j]=((DMN_OPS*)data)[0].evaluation(Point,Point);
     if(fisher)
     {
      ((DMN_OPS*)data)[0].grad_rho_r(Point,Grad,density_evals[i][j]);
      density_grads[i][j]=pow(norm3D(Grad),TWO);
     }
    }
   }
  }
 }
 delete[] x; x=NULL;
 delete[] y; y=NULL;
 delete[] z; z=NULL;
}
//Generate entire grid for momentum space dependent densities
void integrate_quadrature_p(void *data,string name,bool dmn,int order,int order2,bool fisher,double **density_evals,double **density_grads)
{
 int i,j;
 double r_inf,r_sup,*x,*y,*z,Grad[3];
 double Point[3];
 ifstream wr_read,r_read;
 r=new double[order];
 r_real=new double[order];
 w_radial=new double[order];
 //Grid for r
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(name,order,r_inf,r_sup);
 wr_read.open((name+"_w.txt").c_str());
 r_read.open((name+"_x.txt").c_str());
 for(i=0;i<order;i++)
 {
  r_read>>r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<order;i++)
 {
  //r_real is for entire space (0 to inf)
  r_real[i]=r[i]/(ONE-r[i]);
 }
 //Grid for theta, phi
 grid_avail(order2);
 x=new double[order2];
 y=new double[order2];
 z=new double[order2];
 w_theta_phi=new double[order2];
 theta=new double[order2];
 phi=new double[order2];
 ld_by_order(order2,x,y,z,w_theta_phi);
 double *temp,*temp1;
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order2;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 for(i=0;i<order;i++)
 {//Loop over angular points
  for(j=0;j<order2;j++)
  {
   Point[0]=r_real[i]*x[j];
   Point[1]=r_real[i]*y[j];
   Point[2]=r_real[i]*z[j];
   if(!dmn)
   {
    ((READ_FCHK_WFN*)data)[0].rho_p_eval(Point,density_evals[i][j]);
    if(fisher)
    {
     ((READ_FCHK_WFN*)data)[0].rho_p_grad(Point,Grad);
     density_grads[i][j]=pow(norm3D(Grad),TWO);
    }
   }
   else
   {
    density_evals[i][j]=((DMN_P_OPS*)data)[0].evaluation_p(Point,Point);
    if(fisher)
    {
     ((DMN_P_OPS*)data)[0].grad_rho_p(Point,Grad,density_evals[i][j]);
     density_grads[i][j]=pow(norm3D(Grad),TWO);
    }
   }
  }
 }
 delete[] x; x=NULL;
 delete[] y; y=NULL;
 delete[] z; z=NULL;
}
//Evaluate all integrals in position space
void integral_calc(string operation,int &order_r,int &order_theta_phi,double **density_evals,double **density_grads,
double interval[6],double &result_integration)
{
 int i,j;
 double eval,theta_inf,theta_sup,phi_inf,phi_sup,r_inf,r_sup;
 bool interval_quad=false,jump;
 #ifndef RSUP_QUAD
 #define RSUP_QUAD 1e99
 #endif // RSUP_QUAD
 if((interval[0]==ZERO && interval[1]==RSUP_QUAD))
 {
  if((interval[2]==ZERO && abs(interval[3]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
  {
   if((abs(interval[4]+TWO*NINE*TEN)<pow(TEN,-EIGHT) && abs(interval[5]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
   {interval_quad=false;}
   else{interval_quad=true;}
  }
  else{interval_quad=true;}
 }
 else{interval_quad=true;}
 if(interval_quad)
 {
  r_inf=interval[0];
  r_sup=interval[1];
  theta_inf=interval[2];
  theta_sup=interval[3];
  phi_inf=interval[4];
  phi_sup=interval[5];
 }
 result_integration=ZERO;
 //Loop over radial points
 for(i=0;i<order_r;i++)
 {//Loop over angular points
  for(j=0;j<order_theta_phi;j++)
  {
   if(interval_quad)
   {
    if((theta[j]>=theta_inf && theta[j]<=theta_sup)&&(phi[j]>phi_inf&&phi[j]<=phi_sup)&&(r_real[i]>=r_inf&&r_real[i]<=r_sup))
    {
     eval=density_evals[i][j];
     if(operation=="density")
     {
      //Density:
      eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="rho")
     {
      //Rho2:
      eval=eval*eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="shannon")
     {
      //Shannon:
      jump=false;
      if(eval<=ZERO){eval=pow(TEN,-TWO*TEN);}
      if(eval<pow(TEN,-TWO*TEN)){jump=true;}
      if(!jump)
      {
       eval=-eval*log(eval);
      }
      eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="fisher" || operation=="tw" )
     {
      //Fisher or Weiszaker:
      if(eval==ZERO){eval=pow(TEN,-TWO*TEN);}
      eval=density_grads[i][j]/eval;
      eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="ttf")
     {
      //Thomas Fermi:
      if(eval>pow(TEN,-THREE*FIVE))
      {
       eval=pow(eval,FIVE/THREE)*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
      }
      else
      {eval=ZERO;}
     }
     if(operation=="r1")
     {
      //<r>:
      eval=eval*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="r2")
     {
      //<r2>:
      eval=eval*r_real[i]*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Icmx")
     {
      eval=eval*r_real[i]*sin(theta_rad(theta[j]))*cos(phi[j])*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Icmy")
     {
      eval=eval*r_real[i]*sin(theta_rad(theta[j]))*sin(phi[j])*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Icmz")
     {
      eval=eval*r_real[i]*cos(theta_rad(theta[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Ixx")
     {
      eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO)+pow(r_real[i]*cos(theta[j]),TWO))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Ixy")
     {
      eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*sin(theta_rad(theta[j]))
          *sin(phi_rad(phi[j])))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Iyy")
     {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)
         +pow(r_real[i]*cos(theta_rad(theta[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Ixz")
     {
      eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Iyz")
     {
      eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Izz")
     {
      eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)
          +pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
    }
    else{eval=ZERO;}
   }
   else
   {
    eval=density_evals[i][j];
    if(operation=="density")
    {
     //Density:
     eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="rho")
    {
     //Rho2:
     eval=eval*eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="shannon")
    {
     //Shannon:
     jump=false;
     if(eval<=ZERO){eval=pow(TEN,-TWO*TEN);}
     if(eval<pow(TEN,-TWO*TEN)){jump=true;}
     if(!jump)
     {
      eval=-eval*log(eval);
     }
     eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="fisher" || operation=="tw" )
    {
     //Fisher or Weiszaker:
     if(eval==ZERO){eval=pow(TEN,-TWO*TEN);}
     eval=density_grads[i][j]/eval;
     eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="ttf")
    {
     //Thomas Fermi:
     if(eval>pow(TEN,-THREE*FIVE))
     {
      eval=pow(eval,FIVE/THREE)*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     else
     {eval=ZERO;}
    }
    if(operation=="r1")
    {
     //<r>:
     eval=eval*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="r2")
    {
     //<r2>:
     eval=eval*r_real[i]*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Icmx")
    {
     eval=eval*r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Icmy")
    {
     eval=eval*r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Icmz")
    {
     eval=eval*r_real[i]*cos(theta_rad(theta[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Ixx")
    {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO)
         +pow(r_real[i]*cos(theta_rad(theta[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Ixy")
    {
     eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])))
     *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Iyy")
    {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)
         +pow(r_real[i]*cos(theta_rad(theta[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Ixz")
    {
     eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
         *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Iyz")
    {
     eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
         *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Izz")
    {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)
         +pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
   }
   result_integration=result_integration+FOUR*PI*w_theta_phi[j]*w_radial[i]*eval;
  }
 }
}
//Evaluate all integrals in momentum space
void integral_calc_p(string operation,int &order_r,int &order_theta_phi,double **density_evals,double **density_grads,
double interval[6],double &result_integration)
{
 int i,j;
 double eval,theta_inf,theta_sup,phi_inf,phi_sup,r_inf,r_sup;
 bool interval_quad=false,jump;
 #ifndef RSUP_QUAD
 #define RSUP_QUAD 1e99
 #endif // RSUP_QUAD
 if((interval[0]==ZERO && interval[1]==RSUP_QUAD))
 {
  if((interval[2]==ZERO && abs(interval[3]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
  {
   if((abs(interval[4]+TWO*NINE*TEN)<pow(TEN,-EIGHT) && abs(interval[5]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
   {interval_quad=false;}
   else{interval_quad=true;}
  }
  else{interval_quad=true;}
 }
 else{interval_quad=true;}
 if(interval_quad)
 {
  r_inf=interval[0];
  r_sup=interval[1];
  theta_inf=interval[2];
  theta_sup=interval[3];
  phi_inf=interval[4];
  phi_sup=interval[5];
 }
 result_integration=ZERO;
 //Loop over radial points
 for(i=0;i<order_r;i++)
 {//Loop over angular points
  for(j=0;j<order_theta_phi;j++)
  {
   if(interval_quad)
   {
    if((theta[j]>=theta_inf && theta[j]<=theta_sup)&&(phi[j]>phi_inf&&phi[j]<=phi_sup)&&(r_real[i]>=r_inf&&r_real[i]<=r_sup))
    {
     eval=density_evals[i][j];
     if(operation=="densityp")
     {
      //Densityp:
      eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="shannonp")
     {
      //Shannonp:
      jump=false;
      if(eval<=ZERO){eval=pow(TEN,-TWO*TEN);}
      if(eval<pow(TEN,-TWO*TEN)){jump=true;}
      if(!jump)
      {
       eval=-eval*log(eval);
      }
      eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="fisherp")
     {
      //Fisherp:
      if(eval==ZERO){eval=pow(TEN,-TWO*TEN);}
      eval=density_grads[i][j]/eval;
      eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="p1")
     {
      //<p>:
      eval=eval*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="p2")
     {
      //<p2>:
      eval=eval*r_real[i]*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Icmx")
     {
      eval=eval*r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Icmy")
     {
      eval=eval*r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Icmz")
     {
      eval=eval*r_real[i]*cos(theta_rad(theta[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Ixx")
     {
      eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO)+pow(r_real[i]*cos(theta_rad(theta[j])),TWO))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Ixy")
     {
      eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Iyy")
     {
      eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)+pow(r_real[i]*cos(theta_rad(theta[j])),TWO))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Ixz")
     {
      eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Iyz")
     {
      eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
          *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
     if(operation=="Izz")
     {
      eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)
          +pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
     }
    }
    else
    {eval=ZERO;}
   }
   else
   {
    eval=density_evals[i][j];
    if(operation=="densityp")
    {
     //Densityp:
     eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="shannonp")
    {
     //Shannonp:
     jump=false;
     if(eval<=ZERO){eval=pow(TEN,-TWO*TEN);}
     if(eval<pow(TEN,-TWO*TEN)){jump=true;}
     if(!jump)
     {
      eval=-eval*log(eval);
     }
     eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="fisherp")
    {
     //Fisherp:
     if(eval==ZERO){eval=pow(TEN,-TWO*TEN);}
     eval=density_grads[i][j]/eval;
     eval=eval*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="p1")
    {
     //<p>:
     eval=eval*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="p2")
    {
     //<p2>:
     eval=eval*r_real[i]*r_real[i]*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Icmx")
    {
     eval=eval*r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Icmy")
    {
     eval=eval*r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Icmz")
    {
     eval=eval*r_real[i]*cos(theta_rad(theta[j]))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Ixx")
    {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO)
          +pow(r_real[i]*cos(theta_rad(theta[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Ixy")
    {
     eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*sin(theta_rad(theta[j]))
          *sin(phi_rad(phi[j])))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Iyy")
    {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)+pow(r_real[i]
         *cos(theta_rad(theta[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Ixz")
    {
     eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
         *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Iyz")
    {
     eval=-eval*(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])))*(r_real[i]*cos(theta_rad(theta[j])))
         *pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
    if(operation=="Izz")
    {
     eval=eval*(pow(r_real[i]*sin(theta_rad(theta[j]))*cos(phi_rad(phi[j])),TWO)
         +pow(r_real[i]*sin(theta_rad(theta[j]))*sin(phi_rad(phi[j])),TWO))*pow(r_real[i],TWO)/pow(ONE-r[i],TWO);
    }
   }
   result_integration=result_integration+FOUR*PI*w_theta_phi[j]*w_radial[i]*eval;
  }
 }
}
////////////////////////////////////
//Quantities dependent of Overlaps//
////////////////////////////////////
void build_nos_quadrature(void *data,string name,double **ORBITALS,int order,int order2)
{
 int i,j,k,l,counter=0;
 double r_inf,r_sup,*x,*y,*z;
 double Point[3];
 NO no;
 ifstream wr_read,r_read;
 r=new double[order];
 r_real=new double[order];
 w_radial=new double[order];
 //Grid for r
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(name,order,r_inf,r_sup);
 wr_read.open((name+"_w.txt").c_str());
 r_read.open((name+"_x.txt").c_str());
 for(i=0;i<order;i++)
 {
  r_read>>r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<order;i++)
 {
  //r_real is for entire space (0 to inf)
  r_real[i]=r[i]/(ONE-r[i]);
 }
 //Grid for theta, phi
 grid_avail(order2);
 x=new double[order2];
 y=new double[order2];
 z=new double[order2];
 w_theta_phi=new double[order2];
 theta=new double[order2];
 phi=new double[order2];
 ld_by_order(order2,x,y,z,w_theta_phi);
 double *temp,*temp1;
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order2;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 double *AO,**AO_grad;
 AO_grad=new double *[3];
 for(i=0;i<3;i++)
 {AO_grad[i]=new double[((READ_FCHK_WFN*)data)[0].nbasisf];}
 AO=new double[((READ_FCHK_WFN*)data)[0].nbasisf];
 for(i=0;i<order;i++)
 {//Loop over angular points
  for(j=0;j<order2;j++)
  {
   Point[0]=r_real[i]*x[j];
   Point[1]=r_real[i]*y[j];
   Point[2]=r_real[i]*z[j];
   for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasisf;k++){AO[k]=ZERO;}
   for(k=0;k<3;k++)
   {for(l=0;l<((READ_FCHK_WFN*)data)[0].nbasisf;l++){AO_grad[k][l]=ZERO;}}
   if(!((READ_FCHK_WFN*)data)[0].wfn)
   {
    ((READ_FCHK_WFN*)data)[0].build_AO_AOgrad2(AO,AO_grad,Point);
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     no=NO(((READ_FCHK_WFN*)data)[0],AO,AO_grad,k);
     ORBITALS[counter][k]=no.evaluation;
    }
   }
   else
   {
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     no=NO(((READ_FCHK_WFN*)data)[0],Point,k);
     ORBITALS[counter][k]=no.evaluation;
    }
   }
   counter++;
  }
 }
 delete[] AO;
 AO=NULL;
 for(i=0;i<3;i++)
 {delete[] AO_grad[i];AO_grad[i]=NULL;}
 delete[] AO_grad;
 AO_grad=NULL;
 delete[] x; x=NULL;
 delete[] y; y=NULL;
 delete[] z; z=NULL;
}
void build_nos_quadrature_rot(void *data,string name,double **ORBITALS,int order,int order2,double rot)
{
 int i,j,k,l,counter=0;
 double r_inf,r_sup,*x,*y,*z;
 double Point[3];
 NO no;
 ifstream wr_read,r_read;
 r=new double[order];
 r_real=new double[order];
 w_radial=new double[order];
 //Grid for r
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(name,order,r_inf,r_sup);
 wr_read.open((name+"_w.txt").c_str());
 r_read.open((name+"_x.txt").c_str());
 for(i=0;i<order;i++)
 {
  r_read>>r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<order;i++)
 {
  //r_real is for entire space (0 to inf)
  r_real[i]=r[i]/(ONE-r[i]);
 }
 //Grid for theta, phi
 grid_avail(order2);
 x=new double[order2];
 y=new double[order2];
 z=new double[order2];
 w_theta_phi=new double[order2];
 theta=new double[order2];
 phi=new double[order2];
 ld_by_order(order2,x,y,z,w_theta_phi);
 double *temp,*temp1;
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order2;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
  if((phi[i]+rot)<=TWO*NINE*TEN)
  {
   phi[i]=phi[i]+rot;
  }
  else
  {
   phi[i]=-TWO*NINE*TEN+((phi[i]+rot)-TWO*NINE*TEN);
  }
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 double *AO,**AO_grad;
 AO_grad=new double *[3];
 for(i=0;i<3;i++)
 {AO_grad[i]=new double[((READ_FCHK_WFN*)data)[0].nbasisf];}
 AO=new double[((READ_FCHK_WFN*)data)[0].nbasisf];
 for(i=0;i<order;i++)
 {//Loop over angular points
  for(j=0;j<order2;j++)
  {
   Point[0]=r_real[i]*x[j];
   Point[1]=r_real[i]*y[j];
   Point[2]=r_real[i]*z[j];
   rotate_point(Point,rot);
   for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasisf;k++){AO[k]=ZERO;}
   for(k=0;k<3;k++)
   {for(l=0;l<((READ_FCHK_WFN*)data)[0].nbasisf;l++){AO_grad[k][l]=ZERO;}}
   if(!((READ_FCHK_WFN*)data)[0].wfn)
   {
    ((READ_FCHK_WFN*)data)[0].build_AO_AOgrad2(AO,AO_grad,Point);
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     no=NO(((READ_FCHK_WFN*)data)[0],AO,AO_grad,k);
     ORBITALS[counter][k]=no.evaluation;
    }
   }
   else
   {
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     no=NO(((READ_FCHK_WFN*)data)[0],Point,k);
     ORBITALS[counter][k]=no.evaluation;
    }
   }
   counter++;
  }
 }
 delete[] AO;
 AO=NULL;
 for(i=0;i<3;i++)
 {delete[] AO_grad[i];AO_grad[i]=NULL;}
 delete[] AO_grad;
 AO_grad=NULL;
 delete[] x; x=NULL;
 delete[] y; y=NULL;
 delete[] z; z=NULL;
}
void build_mos_quadrature(void *data,string name,double **ORBITALS,int order,int order2)
{
 int i,j,k,l,counter=0;
 double r_inf,r_sup,*x,*y,*z;
 double Point[3];
 MO mo;
 ifstream wr_read,r_read;
 r=new double[order];
 r_real=new double[order];
 w_radial=new double[order];
 //Grid for r
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(name,order,r_inf,r_sup);
 wr_read.open((name+"_w.txt").c_str());
 r_read.open((name+"_x.txt").c_str());
 for(i=0;i<order;i++)
 {
  r_read>>r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<order;i++)
 {
  //r_real is for entire space (0 to inf)
  r_real[i]=r[i]/(ONE-r[i]);
 }
 //Grid for theta, phi
 grid_avail(order2);
 x=new double[order2];
 y=new double[order2];
 z=new double[order2];
 w_theta_phi=new double[order2];
 theta=new double[order2];
 phi=new double[order2];
 ld_by_order(order2,x,y,z,w_theta_phi);
 double *temp,*temp1;
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order2;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 double *AO,**AO_grad;
 AO_grad=new double *[3];
 for(i=0;i<3;i++)
 {AO_grad[i]=new double[((READ_FCHK_WFN*)data)[0].nbasisf];}
 AO=new double[((READ_FCHK_WFN*)data)[0].nbasisf];
 for(i=0;i<order;i++)
 {//Loop over angular points
  for(j=0;j<order2;j++)
  {
   Point[0]=r_real[i]*x[j];
   Point[1]=r_real[i]*y[j];
   Point[2]=r_real[i]*z[j];
   for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasisf;k++){AO[k]=ZERO;}
   for(k=0;k<3;k++)
   {for(l=0;l<((READ_FCHK_WFN*)data)[0].nbasisf;l++){AO_grad[k][l]=ZERO;}}
   if(!((READ_FCHK_WFN*)data)[0].wfn)
   {
    ((READ_FCHK_WFN*)data)[0].build_AO_AOgrad2(AO,AO_grad,Point);
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     mo=MO(((READ_FCHK_WFN*)data)[0],AO,AO_grad,k);
     ORBITALS[counter][k]=mo.evaluation;
    }
   }
   else
   {
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     mo=MO(((READ_FCHK_WFN*)data)[0],Point,k);
     ORBITALS[counter][k]=mo.evaluation;
    }
   }
   counter++;
  }
 }
 delete[] AO;
 AO=NULL;
 for(i=0;i<3;i++)
 {delete[] AO_grad[i];AO_grad[i]=NULL;}
 delete[] AO_grad;
 AO_grad=NULL;
 delete[] x; x=NULL;
 delete[] y; y=NULL;
 delete[] z; z=NULL;
}
void build_mos_quadrature_rot(void *data,string name,double **ORBITALS,int order,int order2,double rot)
{
 int i,j,k,l,counter=0;
 double r_inf,r_sup,*x,*y,*z;
 double Point[3];
 MO mo;
 ifstream wr_read,r_read;
 r=new double[order];
 r_real=new double[order];
 w_radial=new double[order];
 //Grid for r
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(name,order,r_inf,r_sup);
 wr_read.open((name+"_w.txt").c_str());
 r_read.open((name+"_x.txt").c_str());
 for(i=0;i<order;i++)
 {
  r_read>>r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<order;i++)
 {
  //r_real is for entire space (0 to inf)
  r_real[i]=r[i]/(ONE-r[i]);
 }
 //Grid for theta, phi
 grid_avail(order2);
 x=new double[order2];
 y=new double[order2];
 z=new double[order2];
 w_theta_phi=new double[order2];
 theta=new double[order2];
 phi=new double[order2];
 ld_by_order(order2,x,y,z,w_theta_phi);
 double *temp,*temp1;
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order2;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
  if((phi[i]+rot)<=TWO*NINE*TEN)
  {
   phi[i]=phi[i]+rot;
  }
  else
  {
   phi[i]=-TWO*NINE*TEN+((phi[i]+rot)-TWO*NINE*TEN);
  }
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 double *AO,**AO_grad;
 AO_grad=new double *[3];
 for(i=0;i<3;i++)
 {AO_grad[i]=new double[((READ_FCHK_WFN*)data)[0].nbasisf];}
 AO=new double[((READ_FCHK_WFN*)data)[0].nbasisf];
 for(i=0;i<order;i++)
 {//Loop over angular points
  for(j=0;j<order2;j++)
  {
   Point[0]=r_real[i]*x[j];
   Point[1]=r_real[i]*y[j];
   Point[2]=r_real[i]*z[j];
   rotate_point(Point,rot);
   for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasisf;k++){AO[k]=ZERO;}
   for(k=0;k<3;k++)
   {for(l=0;l<((READ_FCHK_WFN*)data)[0].nbasisf;l++){AO_grad[k][l]=ZERO;}}
   if(!((READ_FCHK_WFN*)data)[0].wfn)
   {
    ((READ_FCHK_WFN*)data)[0].build_AO_AOgrad2(AO,AO_grad,Point);
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     mo=MO(((READ_FCHK_WFN*)data)[0],AO,AO_grad,k);
     ORBITALS[counter][k]=mo.evaluation;
    }
   }
   else
   {
    for(k=0;k<((READ_FCHK_WFN*)data)[0].nbasis();k++)
    {
     mo=MO(((READ_FCHK_WFN*)data)[0],Point,k);
     ORBITALS[counter][k]=mo.evaluation;
    }
   }
   counter++;
  }
 }
 delete[] AO;
 AO=NULL;
 for(i=0;i<3;i++)
 {delete[] AO_grad[i];AO_grad[i]=NULL;}
 delete[] AO_grad;
 AO_grad=NULL;
 delete[] x; x=NULL;
 delete[] y; y=NULL;
 delete[] z; z=NULL;
}
//Calculate the Sij matrix
void calc_sij_mat(double **Sij,double **ORBITALS,double interval[6],int &order_r,int &order_ang,int &nbasis)
{
 int i,j,k,l,m;
 double eval,theta_inf,theta_sup,phi_inf,phi_sup,r_inf,r_sup;
 bool interval_quad=false;
 #ifndef RSUP_QUAD
 #define RSUP_QUAD 1e99
 #endif // RSUP_QUAD
 if((interval[0]==ZERO && interval[1]==RSUP_QUAD))
 {
  if((interval[2]==ZERO && abs(interval[3]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
  {
   if((abs(interval[4]+TWO*NINE*TEN)<pow(TEN,-EIGHT) && abs(interval[5]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
   {interval_quad=false;}
   else{interval_quad=true;}
  }
  else{interval_quad=true;}
 }
 else{interval_quad=true;}
 if(interval_quad)
 {
  r_inf=interval[0];
  r_sup=interval[1];
  theta_inf=interval[2];
  theta_sup=interval[3];
  phi_inf=interval[4];
  phi_sup=interval[5];
 }
 //Loop over radial points
 for(i=0;i<nbasis;i++)
 {//Loop over angular points
  for(j=0;j<=i;j++)
  {
   eval=ZERO;
   m=0;
   for(k=0;k<order_r;k++)
   {
    for(l=0;l<order_ang;l++)
    {
     if(interval_quad)
     {
      if((theta[l]>=theta_inf && theta[l]<=theta_sup)&&(phi[l]>phi_inf&&phi[l]<=phi_sup)&&(r_real[k]>=r_inf&&r_real[k]<=r_sup))
      {
       eval=eval+ORBITALS[m][i]*ORBITALS[m][j]*FOUR*PI*w_theta_phi[l]*w_radial[k]*pow(r_real[k],TWO)/pow(ONE-r[k],TWO);
      }
      else
      {
       eval=eval+ZERO;
      }
     }
     else
     {
      eval=eval+ORBITALS[m][i]*ORBITALS[m][j]*FOUR*PI*w_theta_phi[l]*w_radial[k]*pow(r_real[k],TWO)/pow(ONE-r[k],TWO);
     }
     m++;
    }
   }
   Sij[i][j]=eval;
   if(abs(Sij[i][j])<pow(TEN,-EIGHT)){Sij[i][j]=ZERO;}
  }
 }
}
//Calc Mij matrix (r1 moment of MOs)
void calc_mij_mat(double **Sij,double **ORBITALS,double interval[6],int &order_r,int &order_ang,int &nbasis)
{
 int i,j,k,l,m;
 double eval,theta_inf,theta_sup,phi_inf,phi_sup,r_inf,r_sup;
 bool interval_quad=false;
 #ifndef RSUP_QUAD
 #define RSUP_QUAD 1e99
 #endif // RSUP_QUAD
 if((interval[0]==ZERO && interval[1]==RSUP_QUAD))
 {
  if((interval[2]==ZERO && abs(interval[3]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
  {
   if((abs(interval[4]+TWO*NINE*TEN)<pow(TEN,-EIGHT) && abs(interval[5]-TWO*NINE*TEN)<pow(TEN,-EIGHT)))
   {interval_quad=false;}
   else{interval_quad=true;}
  }
  else{interval_quad=true;}
 }
 else{interval_quad=true;}
 if(interval_quad)
 {
  r_inf=interval[0];
  r_sup=interval[1];
  theta_inf=interval[2];
  theta_sup=interval[3];
  phi_inf=interval[4];
  phi_sup=interval[5];
 }
 //Loop over radial points
 for(i=0;i<nbasis;i++)
 {//Loop over angular points
  for(j=0;j<=i;j++)
  {
   eval=ZERO;
   m=0;
   for(k=0;k<order_r;k++)
   {
    for(l=0;l<order_ang;l++)
    {
     if(interval_quad)
     {
      if((theta[l]>=theta_inf && theta[l]<=theta_sup)&&(phi[l]>phi_inf&&phi[l]<=phi_sup)&&(r_real[k]>=r_inf&&r_real[k]<=r_sup))
      {
       eval=eval+ORBITALS[m][i]*ORBITALS[m][j]*FOUR*PI*w_theta_phi[l]*w_radial[k]*pow(r_real[k],THREE)/pow(ONE-r[k],TWO);
      }
      else
      {
       eval=eval+ZERO;
      }
     }
     else
     {
      eval=eval+ORBITALS[m][i]*ORBITALS[m][j]*FOUR*PI*w_theta_phi[l]*w_radial[k]*pow(r_real[k],THREE)/pow(ONE-r[k],TWO);
     }
     m++;
    }
   }
   Sij[i][j]=eval;
   if(abs(Sij[i][j])<pow(TEN,-EIGHT)){Sij[i][j]=ZERO;}
  }
 }
}
// Integrate intracule like for 1RDM
void integrate_intra_coord(double **Intrac,double Dij,double exp_i,double exp_j,double Atom1[3],double Atom2[3],int nx_exp[2],int ny_exp[2],
int nz_exp[2],int Nroot_Lmax_plus_1,double *r_gauss,double *w_gauss,int order_r, int order_ang,bool last)
{
 int i,j,k;
 double zeta_ij,zeta_ij_m1h,zeta_ij_m1h_s,Aij,alpha_ij,alpha_ij_m_12,alpha_ij_p_12,Xij,Yij,Zij,Coef_ij,Point[3],RpRimRj[3],VijX,VijY,VijZ,Iu;
 zeta_ij=exp_i+exp_j;
 zeta_ij_m1h=pow(zeta_ij,-HALF);
 Aij=pow(zeta_ij_m1h,THREE)*Dij;
 alpha_ij=HALF*(exp_i-exp_j)/zeta_ij;
 alpha_ij_m_12=alpha_ij-HALF;
 alpha_ij_p_12=alpha_ij+HALF;
 Xij=(exp_i*Atom1[0]+exp_j*Atom2[0])/zeta_ij;
 Yij=(exp_i*Atom1[1]+exp_j*Atom2[1])/zeta_ij;
 Zij=(exp_i*Atom1[2]+exp_j*Atom2[2])/zeta_ij;
 for(i=0;i<order_r;i++)
 {
  for(j=0;j<order_ang;j++)
  {
   Point[0]=r_real[i]*XX[j];
   Point[1]=r_real[i]*YY[j];
   Point[2]=r_real[i]*ZZ[j];
   RpRimRj[0]=Point[0]+Atom1[0]-Atom2[0];
   RpRimRj[1]=Point[1]+Atom1[1]-Atom2[1];
   RpRimRj[2]=Point[2]+Atom1[2]-Atom2[2];
   Coef_ij=Aij*exp(-(exp_i*exp_j/zeta_ij)*(RpRimRj[0]*RpRimRj[0]+RpRimRj[1]*RpRimRj[1]+RpRimRj[2]*RpRimRj[2]));
   VijX=ZERO;VijY=ZERO;VijZ=ZERO;
   for(k=0;k<Nroot_Lmax_plus_1;k++)
   {
    zeta_ij_m1h_s=zeta_ij_m1h*r_gauss[k];
    VijX=VijX+w_gauss[k]*pow(zeta_ij_m1h_s+alpha_ij_m_12*Point[0]+Xij-Atom1[0],double(nx_exp[0]))
                        *pow(zeta_ij_m1h_s+alpha_ij_p_12*Point[0]+Xij-Atom2[0],double(nx_exp[1]));
    VijY=VijY+w_gauss[k]*pow(zeta_ij_m1h_s+alpha_ij_m_12*Point[1]+Yij-Atom1[1],double(ny_exp[0]))
                        *pow(zeta_ij_m1h_s+alpha_ij_p_12*Point[1]+Yij-Atom2[1],double(ny_exp[1]));
    VijZ=VijZ+w_gauss[k]*pow(zeta_ij_m1h_s+alpha_ij_m_12*Point[2]+Zij-Atom1[2],double(nz_exp[0]))
                        *pow(zeta_ij_m1h_s+alpha_ij_p_12*Point[2]+Zij-Atom2[2],double(nz_exp[1]));
   }
   Intrac[i][j]=Intrac[i][j]+Coef_ij*VijX*VijY*VijZ;
  }
 }
 if(last)
 {
  for(i=0;i<order_r;i++)
  {
   Iu=ZERO;
   for(j=0;j<order_ang;j++)
   {
    // Do not use 4*PI as we want I(0) = Nelectrons. Is like, spherical-averaged quantity.
    Iu=Iu+w_theta_phi[j]*Intrac[i][j];
   }
   if(abs(Iu)<pow(TEN,-EIGHT)){Iu=ZERO;}
   Intrac[i][0]=r_real[i];
   Intrac[i][1]=Iu;
   Intrac[i][2]=Iu*pow(r_real[i],TWO);
  }
 } 
}
//Clean info
void clean_quadrature(string name,int mode)
{
 system(("rm "+name+"_w.txt").c_str());
 system(("rm "+name+"_x.txt").c_str());
 system(("rm "+name+"_r.txt").c_str());
 delete[] r; r=NULL;
 delete[] r_real;r_real=NULL;
 delete[] theta; theta=NULL;
 delete[] phi; phi=NULL;
 delete[] w_radial; w_radial=NULL;
 delete[] w_theta_phi; w_theta_phi=NULL;
 if(mode!=0)
 {
  delete[] XX;XX=NULL;
  delete[] YY;YY=NULL;
  delete[] ZZ;ZZ=NULL;
 }
}
//Check for available order in Lebedev
void grid_avail(int & Order)
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
//Adjust interval from 0 to 180
//to 0 to 3.141592654...
double theta_rad(double &theta)
{
 return (theta*PI)/(TWO*NINE*TEN);
}
//Adjust interval from -180 to 180
//to 0 to 6.283185307
double phi_rad(double &phi)
{
 return ((phi+TWO*NINE*TEN)*PI)/(TWO*NINE*TEN);
}
//Rotate in azimuthal angle
void rotate_point(double Point[3],double rot)
{
 int i;
 double AUX[3],rot_radians;
 rot_radians=rot;
 rot_radians=rot_radians*PI/(TWO*NINE*TEN);
 for(i=0;i<3;i++){AUX[i]=Point[i];}
 Point[0]=AUX[0]*cos(rot_radians)-AUX[1]*sin(rot_radians);
 Point[1]=AUX[0]*sin(rot_radians)+AUX[1]*cos(rot_radians);
 Point[2]=Point[2];
}
