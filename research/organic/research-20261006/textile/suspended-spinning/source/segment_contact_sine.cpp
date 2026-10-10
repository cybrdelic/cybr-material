// Same Abel-normalized material-pair potential, integrated over exact compact
// segment-pair support. No global point lattice is used to discover contact.
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <limits>
struct Vec{double x[3];};
Vec sub(Vec a,Vec b){for(int c=0;c<3;c++)a.x[c]-=b.x[c];return a;}
Vec add(Vec a,Vec b,double s){for(int c=0;c<3;c++)a.x[c]+=s*b.x[c];return a;}
double dot(Vec a,Vec b){double s=0;for(int c=0;c<3;c++)s+=a.x[c]*b.x[c];return s;}
struct Interval{double lo,hi;bool valid;};
Interval intersect(Interval a,Interval b){double l=std::max(a.lo,b.lo),h=std::min(a.hi,b.hi);return {l,h,a.valid&&b.valid&&h>l};}
Interval ball(Vec o,Vec v,double R){double A=dot(v,v),B=dot(o,v),C=dot(o,o)-R*R;if(A<1e-40)return {0,1,C<0};double D=B*B-A*C;if(D<=0)return {0,0,false};double root=std::sqrt(D);return intersect({(-B-root)/A,(-B+root)/A,true},{0,1,true});}
struct Seg{Vec a,v;double ds,s,r,mn[3],mx[3];int fiber,index;};
double segment_distance2(const Seg&a,const Seg&b){
 Vec w=sub(a.a,b.a);double A=dot(a.v,a.v),B=dot(a.v,b.v),C=dot(b.v,b.v),D=dot(a.v,w),F=dot(b.v,w),best=1e300;
 auto put=[&](double u,double v){Vec q=sub(add(a.a,a.v,u),add(b.a,b.v,v));best=std::min(best,dot(q,q));};
 put(0,std::clamp(F/C,0.,1.));put(1,std::clamp((F+B)/C,0.,1.));put(std::clamp(-D/A,0.,1.),0);put(std::clamp((B-D)/A,0.,1.),1);
 double den=A*C-B*B;if(den>1e-14*A*C){double u=(B*F-C*D)/den,v=(A*F-B*D)/den;if(u>=0&&u<=1&&v>=0&&v<=1)put(u,v);}return best;
}
void gauss(int n,std::vector<double>&x,std::vector<double>&w){x.resize(n);w.resize(n);for(int i=0;i<(n+1)/2;i++){double z=std::cos(M_PI*(i+.75)/(n+.5)),pp=0;for(int it=0;it<30;it++){double p1=1,p2=0;for(int j=1;j<=n;j++){double p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j;}pp=n*(z*p1-p2)/(z*z-1);double next=z-p1/pp;if(std::abs(next-z)<2e-15){z=next;break;}z=next;}x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=2/((1-z*z)*pp*pp);}}
extern "C" double segment_pair_contact_sine(int N,int K,const double*p,const double*r,const double*ref,double kc,int order,int self_contact,double exclusion,double*g,double*stats){
 std::fill(g,g+N*K*3,0.);std::fill(stats,stats+7,0.);std::vector<Seg> segs;double rmax=*std::max_element(r,r+N);
 for(int i=0;i<N;i++){double s=0;for(int a=0;a<K-1;a++){Seg z;z.fiber=i;z.index=a;z.r=r[i];z.ds=ref[i*(K-1)+a];z.s=s;for(int c=0;c<3;c++){z.a.x[c]=p[(i*K+a)*3+c];z.v.x[c]=p[(i*K+a+1)*3+c]-z.a.x[c];z.mn[c]=std::min(z.a.x[c],z.a.x[c]+z.v.x[c]);z.mx[c]=std::max(z.a.x[c],z.a.x[c]+z.v.x[c]);}segs.push_back(z);s+=z.ds;}}
 std::sort(segs.begin(),segs.end(),[](const Seg&a,const Seg&b){return a.mn[0]<b.mn[0];});std::vector<double>gx,gw;gauss(order,gx,gw);double E=0;uint64_t count=0,points=0;
 for(size_t ia=0;ia<segs.size();ia++)for(size_t ib=ia+1;ib<segs.size();ib++){Seg a=segs[ia],b=segs[ib];if(b.mn[0]>a.mx[0]+2*rmax)break;if(a.fiber>b.fiber)std::swap(a,b);if(a.fiber==b.fiber){if(!self_contact)continue;if(a.s>b.s)std::swap(a,b);if(b.s+b.ds-a.s<=exclusion*a.r)continue;}
 double R=a.r+b.r;bool skip=false;for(int c=0;c<3;c++)if(a.mn[c]>b.mx[c]+R||b.mn[c]>a.mx[c]+R)skip=true;if(skip)continue;
 Vec rel=sub(a.a,b.a);double b2=dot(b.v,b.v);Interval s0=ball(rel,a.v,R),s1=ball(sub(rel,b.v),a.v,R);double beta=dot(a.v,b.v)/b2,alpha=dot(rel,b.v)/b2;Vec ap=add(a.v,b.v,-beta),rp=add(rel,b.v,-alpha);Interval cyl=ball(rp,ap,R),proj{0,1,true};if(std::abs(beta)<1e-14){proj.valid=alpha>=0&&alpha<=1;}else{double u0=-alpha/beta,u1=(1-alpha)/beta;proj=intersect({std::min(u0,u1),std::max(u0,u1),true},proj);}cyl=intersect(cyl,proj);
 std::vector<double> cuts;for(Interval v:{s0,s1,cyl})if(v.valid){cuts.push_back(v.lo);cuts.push_back(v.hi);}if(cuts.empty())continue;
 if((a.fiber!=b.fiber||b.s-a.s-a.ds>=exclusion*a.r)&&segment_distance2(a,b)<=R*R*1e-24){stats[6]++;return NAN;}
 double lo=*std::min_element(cuts.begin(),cuts.end()),hi=*std::max_element(cuts.begin(),cuts.end());cuts.push_back(lo);cuts.push_back(hi);
 // Split also where a fixed material exclusion meets target endpoints.
 if(a.fiber==b.fiber){for(double z:{(b.s-a.s-exclusion*a.r)/a.ds,(b.s+b.ds-a.s-exclusion*a.r)/a.ds})if(z>lo&&z<hi)cuts.push_back(z);}
 std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<1e-14;}),cuts.end());count++;if(count>2000000)return NAN;
 for(size_t cut=0;cut+1<cuts.size();cut++){double ul=cuts[cut],uh=cuts[cut+1];if(uh<=ul)continue;for(int qi=0;qi<order;qi++){double u=.5*((uh-ul)*gx[qi]+uh+ul),wu=.5*(uh-ul)*gw[qi]*a.ds;Vec x=add(a.a,a.v,u);Interval v=ball(sub(b.a,x),b.v,R);if(!v.valid)continue;if(a.fiber==b.fiber)v=intersect(v,{(a.s+u*a.ds+exclusion*a.r-b.s)/b.ds,1,true});if(!v.valid)continue;
 double center=-dot(sub(b.a,x),b.v)/b2;Vec perpendicular=add(sub(b.a,x),b.v,center);double support_radius=std::sqrt(std::max(0.,(R*R-dot(perpendicular,perpendicular))/b2));if(support_radius<=0)continue;double theta_lo=std::asin(std::clamp((v.lo-center)/support_radius,-1.,1.)),theta_hi=std::asin(std::clamp((v.hi-center)/support_radius,-1.,1.));
 for(int qj=0;qj<order;qj++){double theta=.5*((theta_hi-theta_lo)*gx[qj]+theta_hi+theta_lo),vv=center+support_radius*std::sin(theta),weight=wu*.5*(theta_hi-theta_lo)*gw[qj]*b.ds*support_radius*std::cos(theta);Vec delta=sub(x,add(b.a,b.v,vv));double dd=dot(delta,delta);if(dd>=R*R)continue;if(dd<=R*R*1e-24){stats[6]++;return NAN;}double dist=std::sqrt(dd),root=std::sqrt(R*R-dd),t=root/R,term;if(t<.02){double t2=t*t;term=R*t*t2*(1./3+t2*(1./5+t2*(1./7+t2/9)));}else term=R*std::acosh(R/dist)-root;E+=kc/M_PI*term*weight;double coeff=-kc/M_PI*root/dd*weight;stats[0]=std::max(stats[0],R-dist);stats[4]+=-coeff*dist;points++;
 for(int c=0;c<3;c++){double vgrad=coeff*delta.x[c];g[(a.fiber*K+a.index)*3+c]+=(1-u)*vgrad;g[(a.fiber*K+a.index+1)*3+c]+=u*vgrad;g[(b.fiber*K+b.index)*3+c]-=(1-vv)*vgrad;g[(b.fiber*K+b.index+1)*3+c]-=vv*vgrad;}
 }}}
 }
 stats[1]=count;stats[2]=points;stats[3]=segs.size();return E;
}
