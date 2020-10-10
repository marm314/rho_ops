#include"Critical_Points.h"
CP::CP(){cout<<"Not allowed default constructor"<<endl;}
CP::CP(READ_FCHK_WFN &Rho_r,int AtomA,int AtomB)
{
 int i,j;
 double h,step,steps,res,oldres,Point[3],Grad[3];
 step=pow(TEN,-ONE);
 h=step;
 atomic_num1=(int)Rho_r.Nu_charge[AtomA];
 atomic_num2=(int)Rho_r.Nu_charge[AtomB];
 for(i=0;i<3;i++)
 {
  atom1[i]=Rho_r.Cartesian_Coor[AtomA][i];
  atom2[i]=Rho_r.Cartesian_Coor[AtomB][i];
  direction[i]=atom2[i]-atom1[i];
  Point[i]=ZERO;
  Grad[i]=ZERO;
 }
 res=norm3D(direction);
 for(i=0;i<3;i++){direction[i]=direction[i]/res;}
 steps=res/step;
 res=ZERO;
 oldres=pow(TEN,TEN);
 while(step<steps && oldres>res)
 {
  if(h!=step){oldres=res;}
  for(i=0;i<3;i++)
  {Point[i]=atom1[i]+h*direction[i];}
  h=h+step;
  Rho_r.rho_eval(Point,res);
 }
 for(i=0;i<3;i++){Point[i]=atom1[i]+(h-THREE*step)*direction[i];}
 for(i=0;i<6;i++)
 {
  step=pow(TEN,-(TWO+((double)i)*ONE));
  finesearch(Rho_r,Point,step);
 }
 Rho_r.rho_grad(Point,Grad);
 grad_norm=norm3D(Grad);
 if(grad_norm>pow(TEN,-SIX))
 {
  for(j=0;j<2;j++)
  {
   for(i=0;i<3;i++){Point[i]=Point[i]+Grad[i];}
   Rho_r.rho_grad(Point,Grad);
   grad_norm=norm3D(Grad);
  }
 }
 for(i=0;i<3;i++){BCP[i]=Point[i];}
}
CP::~CP()
{}

void CP::finesearch(READ_FCHK_WFN &Rho_r,double Point[3],double step)
{
 int i;
 double h,res,oldres,AUX[3];
 h=step;
 for(i=0;i<3;i++){AUX[i]=Point[i];}
 res=ZERO;
 oldres=pow(TEN,TEN);
 while(oldres>res)
 {
  if(h!=step){oldres=res;}
  for(i=0;i<3;i++)
  {Point[i]=AUX[i]+h*direction[i];}
  h=h+step;
  Rho_r.rho_eval(Point,res);
 }
 for(i=0;i<3;i++){Point[i]=AUX[i]+(h-THREE*step)*direction[i];}
}
