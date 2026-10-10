// Isolated double-precision port of frozen incremental_envelope.py.
// No fast math, constitutive fitting, history variables or production callers.
// Row-major inputs; material axes are columns. Energy is per reference volume.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {
constexpr int PARAMS = 33;
constexpr int DIAG = 6; // J, lambda_r, lambda_t, lambda_a, axial W, angular W
bool finite(const double* x, std::int64_t n) {
    for (std::int64_t i=0; i<n; ++i) if (!std::isfinite(x[i])) return false;
    return true;
}
double det(const double* f) {
    return f[0]*(f[4]*f[8]-f[5]*f[7])-f[1]*(f[3]*f[8]-f[5]*f[6])+f[2]*(f[3]*f[7]-f[4]*f[6]);
}
int parameters(const double* p) {
    if (!finite(p, PARAMS)) return -3;
    for (int b=0; b<6; ++b) {
        const double* v=p+5*b;
        if (!(v[0]>0 && v[1]>=0 && v[2]>0 && v[3]>=0 && v[4]>0 && v[4]<1)) return -3;
        if (v[3] && (b<3 || v[2]>=v[4])) return -3;
    }
    return p[30]>0 && p[31]>0 && p[31]<1 && p[32]>1 ? 0 : -3;
}
int frame(const double* a) {
    if (!finite(a,9)) return -4;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        double s=0; for(int k=0;k<3;++k) s+=a[3*k+i]*a[3*k+j];
        if (std::abs(s-(i==j?1.:0.))>1e-12) return -4;
    }
    return det(a)>0 ? 0 : -4;
}
void branch(double s, const double* p, double& energy, double& stress) {
    const double first=std::min(s,p[2]), post=std::max(s-p[2],0.);
    energy=.5*p[0]*first*first+p[0]*p[2]*post+.5*p[1]*post*post;
    stress=p[0]*first+p[1]*post;
    if (p[3]) {
        const double t=std::max(s-p[4],0.), a=1-p[4], u=t/a;
        double remainder;
        if (u<.001) {
            double power=u*u*u; remainder=0;
            for(int k=3;k<=10;++k) { remainder+=power/k; power*=u; }
        } else remainder=-std::log1p(-u)-u-.5*u*u;
        energy+=p[3]*a*a*remainder;
        stress+=p[3]*t*t/(a-t);
    }
}
// Partial-pivoted LU solve follows the frozen NumPy Gram-inverse operation.
// Using 1+det_delta as the inverse denominator loses stress accuracy in valid,
// nearly coplanar frames even before the declared Gram-domain gate rejects them.
bool inverse3(const double* off, double* inverse) {
    double lu[9];int permutation[3]={0,1,2};
    for(int i=0;i<9;++i) lu[i]=off[i]+(i%4==0?1.:0.);
    for(int k=0;k<3;++k) {
        int pivot=k;for(int i=k+1;i<3;++i)if(std::abs(lu[3*i+k])>std::abs(lu[3*pivot+k]))pivot=i;
        if(!(std::abs(lu[3*pivot+k])>0))return false;
        if(pivot!=k){for(int j=0;j<3;++j)std::swap(lu[3*k+j],lu[3*pivot+j]);std::swap(permutation[k],permutation[pivot]);}
        for(int i=k+1;i<3;++i){lu[3*i+k]/=lu[3*k+k];for(int j=k+1;j<3;++j)lu[3*i+j]-=lu[3*i+k]*lu[3*k+j];}
    }
    for(int column=0;column<3;++column){
        double rhs[3];for(int i=0;i<3;++i)rhs[i]=permutation[i]==column?1.:0.;
        for(int i=0;i<3;++i)for(int k=0;k<i;++k)rhs[i]-=lu[3*i+k]*rhs[k];
        for(int i=2;i>=0;--i){for(int k=i+1;k<3;++k)rhs[i]-=lu[3*i+k]*rhs[k];rhs[i]/=lu[3*i+i];}
        for(int i=0;i<3;++i)inverse[3*i+column]=rhs[i];
    }
    return finite(inverse,9);
}
int point(const double* h,const double* a,const double* p,double& w,double* stress,double* d) {
    if (!finite(h,9)) return -2;
    int rc=frame(a); if(rc) return rc;
    double f[9],ha[9]={},b[9],l[3],e[3];
    for(int i=0;i<9;++i) f[i]=h[i]+(i%4==0?1.:0.);
    const double j=det(f); if (!(j>0) || !std::isfinite(j)) return -5;
    for(int i=0;i<3;++i) for(int k=0;k<3;++k) for(int v=0;v<3;++v) ha[3*i+k]+=h[3*i+v]*a[3*v+k];
    for(int i=0;i<9;++i) b[i]=a[i]+ha[i];
    double axial=0, signed_p[3];
    for(int i=0;i<3;++i) {
        double linear=0,squared=0;
        for(int k=0;k<3;++k) { linear+=a[3*k+i]*ha[3*k+i]; squared+=ha[3*k+i]*ha[3*k+i]; }
        const double delta=2*linear+squared;
        if (!(delta>-1) || !std::isfinite(delta)) return -6;
        l[i]=std::sqrt(1+delta); e[i]=delta/(l[i]+1);
        if (l[i]<p[31] || l[i]>p[32]) return -7;
        double ut,pt,uc,pc;
        branch(std::max(e[i],0.),p+5*i,ut,pt);
        branch(std::max(-e[i],0.),p+15+5*i,uc,pc);
        axial+=ut+uc; signed_p[i]=pt-pc;
    }
    double off[9]={};
    // Reference A^T A is the identity by contract, as in the frozen correction.
    for(int i=0;i<3;++i) for(int j2=0;j2<3;++j2) if(i!=j2) {
        double v1=0,v2=0,v3=0;
        for(int k=0;k<3;++k) {
            v1+=a[3*k+i]*ha[3*k+j2];
            v2+=ha[3*k+i]*a[3*k+j2];
            v3+=ha[3*k+i]*ha[3*k+j2];
        }
        off[3*i+j2]=(v1+v2+v3)/(l[i]*l[j2]);
    }
    const double x=off[1],y=off[2],z=off[5];
    const double delta_det=2*x*y*z-x*x-y*y-z*z;
    if (!(delta_det>-1) || !std::isfinite(delta_det)) return -8;
    const double angular=-.5*p[30]*std::log1p(delta_det);
    double inv[9];if(!inverse3(off,inv))return -8;
    // Form off*inverse, never I-inverse (which cancels near identity).
    double product[9]={},material_p[9]={};
    for(int i=0;i<3;++i) for(int j2=0;j2<3;++j2) for(int k=0;k<3;++k)
        product[3*i+j2]+=off[3*i+k]*inv[3*k+j2];
    for(int i=0;i<3;++i) for(int j2=0;j2<3;++j2) {
        double v=0;for(int k=0;k<3;++k) v+=(b[3*i+k]/l[k])*product[3*k+j2];
        material_p[3*i+j2]=b[3*i+j2]*(signed_p[j2]/l[j2])+p[30]*v/l[j2];
    }
    for(int i=0;i<3;++i) for(int j2=0;j2<3;++j2) {
        stress[3*i+j2]=0;for(int k=0;k<3;++k) stress[3*i+j2]+=material_p[3*i+k]*a[3*j2+k];
    }
    w=axial+angular; d[0]=j;for(int i=0;i<3;++i)d[1+i]=l[i];d[4]=axial;d[5]=angular;
    return std::isfinite(w) && finite(stress,9) && finite(d,DIAG) ? 0 : -9;
}
}

extern "C" int directional_points(std::int64_t n,const double* h,const double* axes,
    const double* params,double* w,double* p,double* diagnostics,std::int64_t* bad) {
    *bad=-1;if(n<1)return -1;int rc=parameters(params);if(rc)return rc;
    for(std::int64_t i=0;i<n;++i) {
        rc=point(h+9*i,axes+9*i,params,w[i],p+9*i,diagnostics+DIAG*i);
        if(rc){*bad=i;return rc;}
    }
    return 0;
}

extern "C" int directional_prisms(std::int64_t cells,std::int64_t quadrature,
    const double* du,const double* h0,const double* axes,const double* grad,
    const double* weights,const double* params,double* nodal_gradient,
    double* energies,double* volumes,double* diagnostics,std::int64_t* bad) {
    *bad=-1;if(cells<1||quadrature<1)return -1;int rc=parameters(params);if(rc)return rc;
    for(std::int64_t c=0;c<cells;++c) {
        const double* u=du+54*c;
        if(!finite(u,54)){*bad=c*quadrature;return -2;}
        double mean[3]={};for(int i=0;i<18;++i)for(int a=0;a<3;++a)mean[a]+=u[3*i+a]/18.;
        double* ng=nodal_gradient+54*c;std::fill(ng,ng+54,0.);double energy=0,volume=0;
        for(std::int64_t q=0;q<quadrature;++q) {
            const std::int64_t index=c*quadrature+q;const double* g=grad+54*index;
            const double weight=weights[index];double h[9],p[9],w;
            if(!finite(g,54)||!(weight>0)||!std::isfinite(weight)){*bad=index;return -10;}
            std::copy(h0+9*index,h0+9*index+9,h);
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int i=0;i<18;++i)
                h[3*a+b]+=(u[3*i+a]-mean[a])*g[3*i+b];
            double* d=diagnostics+DIAG*index;
            rc=point(h,axes+9*index,params,w,p,d);if(rc){*bad=index;return rc;}
            energy+=weight*w;volume+=weight*d[0];
            for(int i=0;i<18;++i)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                ng[3*i+a]+=weight*g[3*i+b]*p[3*a+b];
        }
        // Exact derivative of the centered local displacement input.
        for(int a=0;a<3;++a){double m=0;for(int i=0;i<18;++i)m+=ng[3*i+a]/18.;for(int i=0;i<18;++i)ng[3*i+a]-=m;}
        energies[c]=energy;volumes[c]=volume;
        if(!std::isfinite(energy)||!std::isfinite(volume)||!finite(ng,54)){*bad=c*quadrature;return -9;}
    }
    return 0;
}
