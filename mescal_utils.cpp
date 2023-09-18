#include"mescal.h"

// Jacobi diag procedure
//Matrix diagonalization and calculation of eigenvectors
void MESCAL::jacobi(int n, double **m, double **v)
{
 int i,j,k,ip,iq,maxiter=10000;
 double tol = pow(10.0,-12.0),apq,t,alpha,c,s,tau,d,delta,temp1,temp2;
 //Define v(matrix eigenvectors) as identity matrix
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {v[i][j]=0.0e0;}
   v[i][i] = 1.0e0;
 }
 //Calculate Delta
 delta=0.0e0;
 for(i=0;i<n;i++)
 {
  for(j=i+1;j<n;j++)
  {delta=delta+m[i][j]*m[i][j];}
 }
 //Iterations cycle
 for(i=0;i<maxiter;i++)
 {
  apq=0.0e0;
  ip=0;
  iq=1;
  for(j=0;j<n;j++)
  {
   for(k=j+1;k<n;k++)
   {
    if(abs(m[j][k])>apq)
    {
     ip = j;
     iq = k;
     apq = abs(m[ip][iq]);
    }
   }
  }
  //Determine c, s, c, tau
  t=1.0e0;
  if(m[ip][ip]!=m[iq][iq])
  {
   alpha=(m[iq][iq]-m[ip][ip])/(2.0e0*m[ip][iq]);
   t=-alpha+alpha/abs(alpha)*sqrt(alpha*alpha+1.0e0);
  }
  c=1.0e0/sqrt(t*t+1.0e0);
  s=t*c;
  tau=s/(1.0e0+c);
  d=m[ip][iq];
  //Update matrix m en and the upper diagonal part
  m[ip][ip]=m[ip][ip]-t*m[ip][iq];
  m[iq][iq]=m[iq][iq]+t*m[ip][iq];
  for(j=0;j<ip;j++)
  {
   temp1=m[j][ip];
   temp2=m[j][iq];
   m[j][ip]=temp1-s*(temp2+tau*temp1);
   m[j][iq]=temp2+s*(temp1-tau*temp2);
  }
  for(j=ip+1;j<iq;j++)
  {
   temp1=m[ip][j];
   temp2=m[j][iq];
   m[j][iq]=temp2+s*(temp1-tau*temp2);
   m[ip][j]=temp1-s*(temp2+tau*temp1);
  }
  for(j=iq+1;j<n;j++)
  {
   temp1=m[ip][j];
   temp2=m[iq][j];
   m[iq][j]=temp2+s*(temp1-tau*temp2);
   m[ip][j]=temp1-s*(temp2+tau*temp1);
  }
  m[ip][iq]=0.0e0;
  //Update v
  for(j=0;j<n;j++)
  {
   temp1=v[j][ip];
   temp2=v[j][iq];
   v[j][ip]=c*temp1-s*temp2;
   v[j][iq]=s*temp1+c*temp2;
  }
  //Update delta
  delta=delta-d*d;
  //If it has converge, update the lower diagonal part of m and finish jacobi()
  if(abs(delta)<=tol)
  {
   for(j=0;j<n;j++)
   {
    for(k=j+1;k<n;k++)
    {m[k][j] = m[j][k];}
   }
   return;
  }
 }
}

// Bellow this point we have functions that only work as libraries to asing values. 
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
// Z -> mass  
double MESCAL::Z2mass(int &Z)
{
 double mass=0.0e0;
 if(Z==1){mass=1.00079;}
 else if(Z==2 ){mass=4.0026;}
 else if(Z==3 ){mass=6.941;}
 else if(Z==4 ){mass=9.0122;}
 else if(Z==5 ){mass=10.811;}
 else if(Z==6 ){mass=12.011;}
 else if(Z==7 ){mass=14.007;}
 else if(Z==8 ){mass=15.999;}
 else if(Z==9 ){mass=18.998;}
 else if(Z==10){mass=20.180;}
 else if(Z==11){mass=22.990;}
 else if(Z==12){mass=24.305;}
 else if(Z==13){mass=26.982;}
 else if(Z==14){mass=28.086;}
 else if(Z==15){mass=30.974;}
 else if(Z==16){mass=32.065;}
 else if(Z==17){mass=35.453;}
 else if(Z==18){mass=39.948;}
 else if(Z==19){mass=39.098;}
 else if(Z==20){mass=40.078;}
 else if(Z==21){mass=44.956;}
 else if(Z==22){mass=47.867;}
 else if(Z==23){mass=50.942;}
 else if(Z==24){mass=51.996;}
 else if(Z==25){mass=54.938;}
 else if(Z==26){mass=55.845;}
 else if(Z==27){mass=58.933;}
 else if(Z==28){mass=58.693;}
 else if(Z==29){mass=63.546;}
 else if(Z==30){mass=65.39;}
 else if(Z==31){mass=69.723;}
 else if(Z==32){mass=72.61;}
 else if(Z==33){mass=74.922;}
 else if(Z==34){mass=78.96;}
 else if(Z==35){mass=79.904;}
 else if(Z==36){mass=83.80;}
 else if(Z==37){mass=85.468;}
 else if(Z==38){mass=87.62;}
 else if(Z==39){mass=89.906;}
 else if(Z==40){mass=91.224;}
 else if(Z==41){mass=92.906;}
 else if(Z==42){mass=95.94;}
 else if(Z==43){mass=98;}
 else if(Z==44){mass=101.07;}
 else if(Z==45){mass=102.91;}
 else if(Z==46){mass=106.42;}
 else if(Z==47){mass=107.87;}
 else if(Z==48){mass=112.41;}
 else if(Z==49){mass=114.82;}
 else if(Z==50){mass=118.71;}
 else if(Z==51){mass=121.76;}
 else if(Z==52){mass=127.60;}
 else if(Z==53){mass=126.90;}
 else if(Z==54){mass=131.29;}
 else if(Z==55){mass=132.91;}
 else if(Z==56){mass=137.33;}
 else if(Z==57){mass=138.91;}
 else if(Z==58){mass=140.12;}
 else if(Z==59){mass=140.91;}
 else if(Z==60){mass=144.24;}
 else if(Z==61){mass=145;}
 else if(Z==62){mass=150.36;}
 else if(Z==63){mass=151.96;}
 else if(Z==64){mass=157.25;}
 else if(Z==65){mass=158.93;}
 else if(Z==66){mass=162.50;}
 else if(Z==67){mass=164.93;}
 else if(Z==68){mass=167.26;}
 else if(Z==69){mass=168.93;}
 else if(Z==70){mass=173.04;}
 else if(Z==71){mass=174.97;}
 else if(Z==72){mass=178.49;}
 else if(Z==73){mass=180.95;}
 else if(Z==74){mass=183.84;}
 else if(Z==75){mass=186.21;}
 else if(Z==76){mass=190.23;}
 else if(Z==77){mass=192.22;}
 else if(Z==78){mass=195.08;}
 else if(Z==79){mass=196.97;}
 else if(Z==80){mass=200.59;}
 else if(Z==81){mass=204.38;}
 else if(Z==82){mass=207.2;}
 else if(Z==83){mass=208.98;}
 else if(Z==84){mass=209;}
 else if(Z==85){mass=210;}
 else if(Z==86){mass=222;}
 else if(Z==87){mass=223;}
 else if(Z==88){mass=226;}
 else if(Z==89){mass=227;}
 else if(Z==90){mass=232.04;}
 else if(Z==91){mass=231.04;}
 else if(Z==92){mass=238.03;}
 else if(Z==93){mass=237;}
 else if(Z==94){mass=244;}
 else if(Z==95){mass=243;}
 else if(Z==96){mass=247;}
 else if(Z==97){mass=247;}
 else if(Z==98){mass=251;}
 else if(Z==99){mass=252;}
 else if(Z==100){mass=257;}
 else if(Z==101){mass=258;}
 else if(Z==102){mass=259;}
 else if(Z==103){mass=262;}
 else if(Z==104){mass=261;}
 else if(Z==105){mass=262;}
 else if(Z==106){mass=266;}
 else if(Z==107){mass=264;}
 else if(Z==108){mass=269;}
 else if(Z==109){mass=268;}
 else{cout<<"Warning! Mass not found for Atomic Number "<<Z<<endl;}
 return mass;
} 
// Z -> valence electrons
int MESCAL::Z2val_electrons(int &Z)
{
 int val_elect=0;
 if(Z==1){val_elect=1;}
 else if(Z==2 ){val_elect=2;}
 else if(Z==3 ){val_elect=1;}
 else if(Z==4 ){val_elect=2;}
 else if(Z==5 ){val_elect=3;}
 else if(Z==6 ){val_elect=4;}
 else if(Z==7 ){val_elect=5;}
 else if(Z==8 ){val_elect=6;}
 else if(Z==9 ){val_elect=7;}
 else if(Z==10){val_elect=8;}
 else if(Z==11){val_elect=1;}
 else if(Z==12){val_elect=2;}
 else if(Z==13){val_elect=3;}
 else if(Z==14){val_elect=4;}
 else if(Z==15){val_elect=5;}
 else if(Z==16){val_elect=6;}
 else if(Z==17){val_elect=7;}
 else if(Z==18){val_elect=8;}
 else if(Z==19){val_elect=1;}
 else if(Z==20){val_elect=2;}
 else if(Z==21){val_elect=3;}
 else if(Z==22){val_elect=4;}
 else if(Z==23){val_elect=5;}
 else if(Z==24){val_elect=6;}
 else if(Z==25){val_elect=7;}
 else if(Z==26){val_elect=8;}
 else if(Z==27){val_elect=9;}
 else if(Z==28){val_elect=10;}
 else if(Z==29){val_elect=11;}
 else if(Z==30){val_elect=12;}
 else if(Z==31){val_elect=13;}
 else if(Z==32){val_elect=14;}
 else if(Z==33){val_elect=15;}
 else if(Z==34){val_elect=16;}
 else if(Z==35){val_elect=17;}
 else if(Z==36){val_elect=18;}
 else if(Z==37){val_elect=1;}
 else if(Z==38){val_elect=2;}
 else if(Z==39){val_elect=3;}
 else if(Z==40){val_elect=4;}
 else if(Z==41){val_elect=5;}
 else if(Z==42){val_elect=6;}
 else if(Z==43){val_elect=7;}
 else if(Z==44){val_elect=8;}
 else if(Z==45){val_elect=9;}
 else if(Z==46){val_elect=10;}
 else if(Z==47){val_elect=11;}
 else if(Z==48){val_elect=12;}
 else if(Z==49){val_elect=13;}
 else if(Z==50){val_elect=14;}
 else if(Z==51){val_elect=15;}
 else if(Z==52){val_elect=16;}
 else if(Z==53){val_elect=17;}
 else if(Z==54){val_elect=18;}
 else if(Z==55){val_elect=1;}
 else if(Z==56){val_elect=2;}
 else if(Z==57){val_elect=3;}
 else if(Z==58){val_elect=4;}
 else if(Z==59){val_elect=5;}
 else if(Z==60){val_elect=6;}
 else if(Z==61){val_elect=7;}
 else if(Z==62){val_elect=8;}
 else if(Z==63){val_elect=9;}
 else if(Z==64){val_elect=10;}
 else if(Z==65){val_elect=11;}
 else if(Z==66){val_elect=12;}
 else if(Z==67){val_elect=13;}
 else if(Z==68){val_elect=14;}
 else if(Z==69){val_elect=15;}
 else if(Z==70){val_elect=16;}
 else if(Z==71){val_elect=17;}
 else if(Z==72){val_elect=18;}
 else if(Z==73){val_elect=19;}
 else if(Z==74){val_elect=20;}
 else if(Z==75){val_elect=21;}
 else if(Z==76){val_elect=22;}
 else if(Z==77){val_elect=23;}
 else if(Z==78){val_elect=24;}
 else if(Z==79){val_elect=25;}
 else if(Z==80){val_elect=26;}
 else if(Z==81){val_elect=27;}
 else if(Z==82){val_elect=28;}
 else if(Z==83){val_elect=29;}
 else if(Z==84){val_elect=30;}
 else if(Z==85){val_elect=31;}
 else if(Z==86){val_elect=32;}
 else if(Z==87){val_elect=1;}
 else if(Z==88){val_elect=2;}
 else if(Z==89){val_elect=3;}
 else if(Z==90){val_elect=4;}
 else if(Z==91){val_elect=5;}
 else if(Z==92){val_elect=6;}
 else if(Z==93){val_elect=7;}
 else if(Z==94){val_elect=8;}
 else if(Z==95){val_elect=9;}
 else if(Z==96){val_elect=10;}
 else if(Z==97){val_elect=11;}
 else if(Z==98){val_elect=12;}
 else if(Z==99){val_elect=13;}
 else if(Z==100){val_elect=14;}
 else if(Z==101){val_elect=15;}
 else if(Z==102){val_elect=16;}
 else if(Z==103){val_elect=17;}
 else if(Z==104){val_elect=18;}
 else if(Z==105){val_elect=19;}
 else if(Z==106){val_elect=20;}
 else if(Z==107){val_elect=21;}
 else if(Z==108){val_elect=22;}
 else if(Z==109){val_elect=23;}
 else{cout<<"Warning! Valence electrons number not found for Atomic Number "<<Z<<endl;}
 return val_elect;
} 
