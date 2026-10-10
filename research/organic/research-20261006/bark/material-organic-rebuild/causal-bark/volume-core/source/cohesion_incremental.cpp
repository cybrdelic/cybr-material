// Incremental-origin variant: origin interfaces are bonded/coincident.
#include <cmath>
#include <algorithm>
static void cross(const double*a,const double*b,double*c){c[0]=a[1]*b[2]-a[2]*b[1];c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0];}
extern "C" int cohesive_incremental(int Cn,int N,const int*pairs,const double*q,const double*tref,const double*jc,const double*ac,const double*bc,const double*area,const double*kt,const double*kn,const double*damage,const double*sign,double*energy,double*gradient,double*opening){
 if(Cn<1||N<1)return -1;*energy=0;std::fill(gradient,gradient+Cn*54,0.);
 for(int p=0;p<N;p++){
  int ids[36];for(int side=0;side<2;side++){int c=pairs[2*p+side];if(c<0||c>=Cn)return -2;for(int i=0;i<18;i++)ids[side*18+i]=c*18+i;}if(!(area[p]>0&&kt[p]>0&&kn[p]>0&&damage[p]>=0&&damage[p]<=1&&(sign[p]==1||sign[p]==-1)))return -3;
  double jump[3]={0},t1[3],t2[3];for(int a=0;a<3;a++){t1[a]=tref[p*6+a];t2[a]=tref[p*6+3+a];}for(int i=0;i<36;i++)for(int a=0;a<3;a++){double x=q[ids[i]*3+a]-q[ids[0]*3+a];jump[a]+=jc[p*36+i]*x;t1[a]+=ac[p*36+i]*x;t2[a]+=bc[p*36+i]*x;}
  double v[3];cross(t1,t2,v);double length=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);if(!(length>1e-18))return -4;double n[3],dn=0,j2=0;for(int a=0;a<3;a++){n[a]=sign[p]*v[a]/length;dn+=jump[a]*n[a];j2+=jump[a]*jump[a];}double w=1-damage[p],closed=std::min(dn,0.);*energy+=.5*area[p]*(w*kt[p]*j2+w*(kn[p]-kt[p])*dn*dn+damage[p]*kn[p]*closed*closed);opening[p]=std::sqrt(std::max(0.,j2-dn*dn+kn[p]/kt[p]*std::max(dn,0.)*std::max(dn,0.)));
  double b[3];for(int a=0;a<3;a++)b[a]=sign[p]*(jump[a]-dn*n[a])/length;double v1[3],v2[3];cross(t2,b,v1);cross(b,t1,v2);double c=w*(kn[p]-kt[p])*dn+damage[p]*kn[p]*closed;
  double sum[3]={0};for(int i=0;i<36;i++)for(int a=0;a<3;a++){double g=area[p]*(w*kt[p]*jc[p*36+i]*jump[a]+c*(jc[p*36+i]*n[a]+ac[p*36+i]*v1[a]+bc[p*36+i]*v2[a]));gradient[ids[i]*3+a]+=g;sum[a]+=g;}for(int a=0;a<3;a++)gradient[ids[0]*3+a]-=sum[a];
 }
 return 0;
}
