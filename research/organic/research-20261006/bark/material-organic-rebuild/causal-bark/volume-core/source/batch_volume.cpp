// Own implementation of the same compressible neo-Hookean material operator.
// 18 geometric nodes/cell; all arrays double, contiguous, in declared SI order.
#include <cmath>
#include <algorithm>
static double det(const double*a){return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);}
extern "C" int volume_batch(int Cn,int Qn,const double*x,const double*G,const double*w,const double*mu,const double*la,double*gradient,double*energies,double*volumes,double*minimumJ){
 if(Cn<1||Qn<1)return -1;std::fill(gradient,gradient+Cn*54,0.);*minimumJ=1e300;
 for(int c=0;c<Cn;c++){
  double mean[3]={0,0,0};for(int i=0;i<18;i++)for(int a=0;a<3;a++)mean[a]+=x[c*54+i*3+a]/18.;double E=0,V=0;
  for(int q=0;q<Qn;q++){
   const double*g=G+(c*Qn+q)*54;double F[9]={0};for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int i=0;i<18;i++)F[a*3+b]+=(x[c*54+i*3+a]-mean[a])*g[i*3+b];double J=det(F);if(!(J>0)||!std::isfinite(J))return -2;*minimumJ=std::min(*minimumJ,J);
   double invT[9]={F[4]*F[8]-F[5]*F[7],F[5]*F[6]-F[3]*F[8],F[3]*F[7]-F[4]*F[6],F[2]*F[7]-F[1]*F[8],F[0]*F[8]-F[2]*F[6],F[1]*F[6]-F[0]*F[7],F[1]*F[5]-F[2]*F[4],F[2]*F[3]-F[0]*F[5],F[0]*F[4]-F[1]*F[3]};for(double &v:invT)v/=J;
   double A[9];double tr=0,p2=0,norm2=0;for(int a=0;a<3;a++)for(int b=0;b<3;b++){double v=0;for(int k=0;k<3;k++)v+=F[k*3+a]*F[k*3+b];A[a*3+b]=v-(a==b?1.:0.);norm2+=A[a*3+b]*A[a*3+b];}tr=A[0]+A[4]+A[8];for(int a=0;a<3;a++)for(int b=0;b<3;b++)p2+=A[a*3+b]*A[b*3+a];
   double lnJ,psi;
   if(norm2<.0625){
    // Cayley-Hamilton power sums evaluate tr(A)-log(det(I+A))
    // without subtracting two first-order quantities at small strain.
    double I2=.5*(tr*tr-p2),I3=det(A),p0=3,p1=tr,p=p2,diff=.5*p2;
    for(int k=3;k<=24;k++){double next=tr*p-I2*p1+I3*p0;diff+=(k%2?-1.:1.)*next/k;p0=p1;p1=p;p=next;}
    lnJ=.5*(tr-diff);psi=.5*mu[c]*diff+.5*la[c]*lnJ*lnJ;
   }else{lnJ=std::log(J);psi=.5*mu[c]*tr-mu[c]*lnJ+.5*la[c]*lnJ*lnJ;}
   double P[9];for(int k=0;k<9;k++)P[k]=mu[c]*(F[k]-invT[k])+la[c]*lnJ*invT[k];double weight=w[c*Qn+q];E+=weight*psi;V+=weight*J;
   for(int i=0;i<18;i++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)gradient[c*54+i*3+a]+=weight*g[i*3+b]*P[a*3+b];
  }
  // Exact derivative of centering the coordinates in F.
  for(int a=0;a<3;a++){double s=0;for(int i=0;i<18;i++)s+=gradient[c*54+i*3+a]/18.;for(int i=0;i<18;i++)gradient[c*54+i*3+a]-=s;}
  energies[c]=E;volumes[c]=V;
 }
 return 0;
}
