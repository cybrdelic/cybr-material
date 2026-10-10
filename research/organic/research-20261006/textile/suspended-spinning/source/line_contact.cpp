// Symmetric closest-centreline contact integral. Input/output units are SI.
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cstdint>
struct Cell {long long x,y,z; bool operator==(const Cell&o)const{return x==o.x&&y==o.y&&z==o.z;}};
struct Hash {size_t operator()(const Cell&a)const{return uint64_t(a.x)*73856093ULL ^ uint64_t(a.y)*19349663ULL ^ uint64_t(a.z)*83492791ULL;}};
static const double gx[4][4]={{.5,0,0,0},{.21132486540518713,.7886751345948129,0,0},{.1127016653792583,.5,.8872983346207417,0},{.06943184420297371,.33000947820757187,.6699905217924281,.9305681557970262}};
static const double gw[4][4]={{1,0,0,0},{.5,.5,0,0},{.2777777777777778,.4444444444444444,.2777777777777778,0},{.17392742256872693,.32607257743127305,.32607257743127305,.17392742256872693}};
extern "C" {
double line_contact_impl(int N,int K,const double*p,const double*r,const double*ref,double kc,int order,int current_measure,int self_contact,double exclusion_radii,int grid_mode,double*g,double*stats,int M,const double*phi,double*H){
 if(N<1||K<2||order<1||order>4||kc<0)return NAN;
 const int S=K-1;const int ndof=N*M*2;if(H)std::fill(H,H+ndof*ndof,0.);std::fill(g,g+N*K*3,0.);std::fill(stats,stats+7,0.);double energy=0,rmax=0;
 std::vector<double> arc(N*K,0),lens(N*S);for(int i=0;i<N;i++){rmax=std::max(rmax,r[i]);for(int k=0;k<S;k++){double l2=0;for(int d=0;d<3;d++){double q=p[3*(i*K+k+1)+d]-p[3*(i*K+k)+d];l2+=q*q;}lens[i*S+k]=std::sqrt(l2);arc[i*K+k+1]=arc[i*K+k]+ref[i*S+k];}}
 if(rmax<=0)return NAN;double cellsize=2*rmax;
 auto cell=[&](const double*x){return Cell{(long long)std::floor(x[0]/cellsize),(long long)std::floor(x[1]/cellsize),(long long)std::floor(x[2]/cellsize)};};
 std::unordered_map<Cell,std::vector<int>,Hash> grid;bool usegrid=grid_mode==1||(grid_mode==2&&N>8);size_t inserted=0;
 if(usegrid){grid.reserve(N*S*3);for(int i=0;i<N&&usegrid;i++)for(int k=0;k<S&&usegrid;k++){
  const double*A=p+3*(i*K+k),*B=A+3;double lo[3],hi[3];for(int d=0;d<3;d++){lo[d]=std::min(A[d],B[d])-r[i]-rmax;hi[d]=std::max(A[d],B[d])+r[i]+rmax;}Cell a=cell(lo),b=cell(hi);
  double needed=double(b.x-a.x+1)*double(b.y-a.y+1)*double(b.z-a.z+1);if(needed+inserted>5000000){usegrid=false;grid.clear();break;}
  for(long long x=a.x;x<=b.x;x++)for(long long y=a.y;y<=b.y;y++)for(long long z=a.z;z<=b.z;z++){grid[Cell{x,y,z}].push_back(i*S+k);inserted++;}
 }}
 std::vector<int> all(N*S);for(int i=0;i<N*S;i++)all[i]=i;
 std::vector<int> stamp(N,-1),bestSeg(N),touched;std::vector<double> bestD2(N),bestBeta(N);int query=0;
 for(int i=0;i<N;i++)for(int k=0;k<S;k++)for(int iq=0;iq<order;iq++,query++){
  const double*A=p+3*(i*K+k),*B=A+3;double a=gx[order-1][iq],sp[3];for(int d=0;d<3;d++)sp[d]=(1-a)*A[d]+a*B[d];double sarc=arc[i*K+k]+a*ref[i*S+k];
  const std::vector<int>*cand=&all;
  if(usegrid){auto found=grid.find(cell(sp));if(found==grid.end())continue;cand=&found->second;}
  touched.clear();
  for(int sid:*cand){int j=sid/S,l=sid%S;if(i==j&&!self_contact)continue;double R=r[i]+r[j];const double*C=p+3*(j*K+l),*D=C+3;
   bool far=false;for(int d=0;d<3;d++)if(sp[d]<std::min(C[d],D[d])-R||sp[d]>std::max(C[d],D[d])+R)far=true;if(far)continue;
   double v[3],len2=0,dot=0;for(int d=0;d<3;d++){v[d]=D[d]-C[d];len2+=v[d]*v[d];dot+=(sp[d]-C[d])*v[d];}if(len2<=1e-30)continue;
   auto candidate=[&](double lo,double hi){if(hi<lo)return;double b=std::clamp(dot/len2,lo,hi),d2=0;for(int d=0;d<3;d++){double q=sp[d]-C[d]-b*v[d];d2+=q*q;}if(d2>=R*R)return;
    if(stamp[j]!=query){stamp[j]=query;bestD2[j]=d2;bestSeg[j]=l;bestBeta[j]=b;touched.push_back(j);}else if(d2<bestD2[j]){bestD2[j]=d2;bestSeg[j]=l;bestBeta[j]=b;}};
   if(i!=j)candidate(0.,1.);else{
    double e=exclusion_radii*r[i],r0=arc[j*K+l],dl=ref[j*S+l];
    if(sarc-e>=r0)candidate(0.,std::min(1.,(sarc-e-r0)/dl));
    if(sarc+e<=r0+dl)candidate(std::max(0.,(sarc+e-r0)/dl),1.);
   }
  }
  for(int j:touched){int l=bestSeg[j];double b=bestBeta[j],dist=std::sqrt(bestD2[j]),pen=r[i]+r[j]-dist;const double*C=p+3*(j*K+l),*D=C+3;
   double measure=current_measure?lens[i*S+k]:ref[i*S+k],weight=.5*gw[order-1][iq]*measure,E=.5*kc*weight*pen*pen;energy+=E;stats[6]+=kc*weight*pen;stats[1]=std::max(stats[1],pen);stats[2]++;
   if(dist<=1e-15){stats[3]++;continue;}
   if(H){int ids[32],nn=0;double J[32];for(int m=0;m<M;m++)for(int c=0;c<2;c++){
    int d=c==0?0:2;double n=(sp[d]-(1-b)*C[d]-b*D[d])/dist,ps=(1-a)*phi[k*M+m]+a*phi[(k+1)*M+m],pt=(1-b)*phi[l*M+m]+b*phi[(l+1)*M+m];
    ids[nn]=(i*M+m)*2+c;J[nn]=-n*ps+(current_measure?.5*pen/lens[i*S+k]*(B[d]-A[d])/lens[i*S+k]*(phi[(k+1)*M+m]-phi[k*M+m]):0.);nn++;
    ids[nn]=(j*M+m)*2+c;J[nn]=n*pt;nn++;
   }for(int a=0;a<nn;a++)for(int b=0;b<nn;b++)H[ids[a]*ndof+ids[b]]+=kc*weight*J[a]*J[b];}
   for(int d=0;d<3;d++){
    double normal=(sp[d]-(1-b)*C[d]-b*D[d])/dist,q=-kc*weight*pen*normal;
    g[3*(i*K+k)+d]+=(1-a)*q;g[3*(i*K+k+1)+d]+=a*q;g[3*(j*K+l)+d]-=(1-b)*q;g[3*(j*K+l+1)+d]-=b*q;
    if(current_measure){double v=E/lens[i*S+k]*(B[d]-A[d])/lens[i*S+k];g[3*(i*K+k)+d]-=v;g[3*(i*K+k+1)+d]+=v;}
   }
  }
 }
 stats[0]=energy;stats[4]=usegrid;stats[5]=inserted;return energy;
}
// Exact integral on each linear segment of quadratic one-sided plane penalty.
double plane_contact_impl(int N,int K,const double*p,const double*r,const double*ref,double kc,double top,double bottom,int current_measure,double*g,double*stats,double*top_node_forces,int M,const double*phi,double*H){
 const int ndof=N*M*2;if(H)std::fill(H,H+ndof*ndof,0.);std::fill(g,g+N*K*3,0.);std::fill(stats,stats+5,0.);std::fill(top_node_forces,top_node_forces+N*K,0.);double Etotal=0;
 for(int i=0;i<N;i++)for(int k=0;k<K-1;k++){
  const double*A=p+3*(i*K+k),*B=A+3;double len2=0;for(int d=0;d<3;d++){double v=B[d]-A[d];len2+=v*v;}double len=std::sqrt(len2),measure=current_measure?len:ref[i*(K-1)+k];
  for(int side=0;side<2;side++){
   double sign=side==0?1.:-1.,a=side==0?A[2]+r[i]-top:bottom+r[i]-A[2],b=sign*(B[2]-A[2]);
   stats[side+2]=std::max(stats[side+2],std::max(a,a+b));double lo=0,hi=1;
   if(std::abs(b)<1e-30){if(a<=0)continue;}else if(b>0)lo=std::max(0.,-a/b);else hi=std::min(1.,-a/b);
   lo=std::clamp(lo,0.,1.);hi=std::clamp(hi,0.,1.);if(hi<=lo)continue;
   double width=hi-lo,pl=std::max(0.,a+b*lo),ph=std::max(0.,a+b*hi);
   double I0=width*(pl+ph)/2,I1=lo*I0+width*width*(pl+2*ph)/6,I2=width*(pl*pl+pl*ph+ph*ph)/3;
   double E=.5*kc*measure*I2,fa=kc*measure*(I0-I1),fb=kc*measure*I1;Etotal+=E;stats[side]+=fa+fb;
   if(H){for(int iq=0;iq<3;iq++){double s=lo+(hi-lo)*gx[2][iq],pen=std::max(0.,a+b*s),qw=(hi-lo)*gw[2][iq]*measure;int ids[16],nn=0;double J[16];
    for(int m=0;m<M;m++)for(int c=0;c<2;c++){int d=c==0?0:2;ids[nn]=(i*M+m)*2+c;J[nn]=(c==1?sign*((1-s)*phi[k*M+m]+s*phi[(k+1)*M+m]):0.)+(current_measure?.5*pen/len*(B[d]-A[d])/len*(phi[(k+1)*M+m]-phi[k*M+m]):0.);nn++;}
    for(int a=0;a<nn;a++)for(int b=0;b<nn;b++)H[ids[a]*ndof+ids[b]]+=kc*qw*J[a]*J[b];
   }}
   g[3*(i*K+k)+2]+=sign*fa;g[3*(i*K+k+1)+2]+=sign*fb;
   if(side==0){top_node_forces[i*K+k]+=fa;top_node_forces[i*K+k+1]+=fb;}
   if(current_measure&&len>1e-15)for(int d=0;d<3;d++){double q=E/len*(B[d]-A[d])/len;g[3*(i*K+k)+d]-=q;g[3*(i*K+k+1)+d]+=q;}
  }
 }
 stats[4]=Etotal;return Etotal;
}
}
extern "C" double line_contact(int N,int K,const double*p,const double*r,const double*ref,double kc,int order,int current_measure,int self_contact,double exclusion_radii,int grid_mode,double*g,double*stats){return line_contact_impl(N,K,p,r,ref,kc,order,current_measure,self_contact,exclusion_radii,grid_mode,g,stats,0,nullptr,nullptr);}
extern "C" double line_contact_modal(int N,int K,const double*p,const double*r,const double*ref,double kc,int order,int current_measure,int self_contact,double exclusion_radii,int grid_mode,double*g,double*stats,int M,const double*phi,double*H){if(M<1||M>4||N*M*2>2048)return NAN;return line_contact_impl(N,K,p,r,ref,kc,order,current_measure,self_contact,exclusion_radii,grid_mode,g,stats,M,phi,H);}
extern "C" double plane_contact(int N,int K,const double*p,const double*r,const double*ref,double kc,double top,double bottom,int current_measure,double*g,double*stats,double*node){return plane_contact_impl(N,K,p,r,ref,kc,top,bottom,current_measure,g,stats,node,0,nullptr,nullptr);}
extern "C" double plane_contact_modal(int N,int K,const double*p,const double*r,const double*ref,double kc,double top,double bottom,int current_measure,double*g,double*stats,double*node,int M,const double*phi,double*H){if(M<1||M>4||N*M*2>2048)return NAN;return plane_contact_impl(N,K,p,r,ref,kc,top,bottom,current_measure,g,stats,node,M,phi,H);}
extern "C" void capsule_audit(int N,int K,const double*p,const double*r,const double*ref,int self_contact,double exclusion_radii,double*out){
 const int S=K-1;std::vector<double> arc(N*K,0),dy(N,0);for(int i=0;i<N;i++){for(int k=0;k<S;k++)arc[i*K+k+1]=arc[i*K+k]+ref[i*S+k];double d=p[3*(i*K+1)+1]-p[3*i*K+1];bool uniform=d>0;for(int k=1;k<K;k++)if(std::abs(p[3*(i*K+k)+1]-p[3*i*K+1]-k*d)>1e-12)uniform=false;dy[i]=uniform?d:0;}
 double maxpen=0,mind=1e100;int count=0;
 for(int i=0;i<N;i++)for(int j=i;j<N;j++){if(i==j&&!self_contact)continue;double R=r[i]+r[j];
  for(int ka=0;ka<S;ka++){
   const double*A=p+3*(i*K+ka);int lo=i==j?ka+1:0,hi=S-1;
   if(dy[j]>0){double y0=p[3*j*K+1];lo=std::max(lo,int(std::floor((std::min(A[1],A[4])-R-y0)/dy[j]))-1);hi=std::min(hi,int(std::floor((std::max(A[1],A[4])+R-y0)/dy[j]))+1);}
   for(int kb=lo;kb<=hi;kb++){
    const double*B=p+3*(j*K+kb);bool far=false;for(int d=0;d<3;d++)if(std::min(A[d],A[3+d])-std::max(B[d],B[3+d])>=R||std::min(B[d],B[3+d])-std::max(A[d],A[3+d])>=R)far=true;if(far)continue;
    double u[3],v[3],w[3],aa=0,bb=0,cc=0,dd=0,ee=0;for(int d=0;d<3;d++){u[d]=A[3+d]-A[d];v[d]=B[3+d]-B[d];w[d]=A[d]-B[d];aa+=u[d]*u[d];bb+=u[d]*v[d];cc+=v[d]*v[d];dd+=u[d]*w[d];ee+=v[d]*w[d];}
    if(aa<=1e-30||cc<=1e-30)continue;double den=aa*cc-bb*bb,s=den>1e-30?std::clamp((bb*ee-cc*dd)/den,0.,1.):0.,t=(bb*s+ee)/cc;
    if(t<0){t=0;s=std::clamp(-dd/aa,0.,1.);}else if(t>1){t=1;s=std::clamp((bb-dd)/aa,0.,1.);}
    if(i==j){double e=exclusion_radii*r[i],a0=arc[i*K+ka],b0=arc[j*K+kb],la=ref[i*S+ka],lb=ref[j*S+kb];if(b0+lb-a0<e)continue;
     if(b0+t*lb-a0-s*la<e){double low=std::max(0.,(b0-a0-e)/la),high=std::min(1.,(b0+lb-a0-e)/la);if(high<low)continue;double alpha=la/lb,beta=(e+a0-b0)/lb,dp=0,qq=0;
      for(int d=0;d<3;d++){double vv=u[d]-alpha*v[d],ww=w[d]-beta*v[d];dp+=ww*vv;qq+=vv*vv;}s=qq>1e-30?std::clamp(-dp/qq,low,high):low;t=alpha*s+beta;
     }
    }
    double dist=0;for(int d=0;d<3;d++){double q=w[d]+s*u[d]-t*v[d];dist+=q*q;}dist=std::sqrt(dist);mind=std::min(mind,dist);if(dist<R){count++;maxpen=std::max(maxpen,R-dist);}
   }
  }
 }
 out[0]=maxpen;out[1]=mind;out[2]=count;
}
