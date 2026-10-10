"""Two closed planar cells coupled by shared nodes and junction rotations.

Axes: x = transverse/circumferential section coordinate, z = radial coordinate.
Positions are stored in material-node order, never as a single-valued x(z).
The SI energy is a finite-dimensional extensible Euler--Bernoulli strip energy.
End half-cell bending quadrature couples each wall to a common junction rotation.
Contact is checked after equilibrium, not solved. This is a mechanism prototype.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import eigh, solve
from corrugated_wall_frozen import Parameters, wrap


@dataclass(frozen=True)
class NetworkParameters(Parameters):
    cell_width_m: float = 30e-6
    junction_contact_radius_m: float = 2e-6


class TwoCellNetwork:
    def __init__(self, radial_segments=24, transverse_segments=8,
                 parameters=NetworkParameters()):
        self.p = parameters
        self.L = parameters.height_m
        self.A = parameters.wall_thickness_m * parameters.strip_width_m
        self.I = parameters.strip_width_m * parameters.wall_thickness_m**3 / 12
        self.EA = parameters.wall_modulus_Pa * self.A
        self.EI = parameters.wall_modulus_Pa * self.I
        self.unit = self.EI / self.L
        width = parameters.cell_width_m / self.L
        vertices = np.array([[0, 0], [width, 0], [2*width, 0],
                             [0, 1], [width, 1], [2*width, 1]], dtype=float)
        positions = list(vertices)
        self.walls = []
        for name, a, b, radial in [('left',0,3,True), ('shared',1,4,True),
                ('right',2,5,True), ('bottom_left',0,1,False),
                ('bottom_right',1,2,False), ('top_left',3,4,False),
                ('top_right',4,5,False)]:
            n = radial_segments if radial else transverse_segments
            if n < 4: raise ValueError('Each wall needs at least four segments')
            u = np.linspace(0,1,n+1)
            X = (1-u[:,None])*vertices[a] + u[:,None]*vertices[b]
            if radial:
                X[:,0] += parameters.corrugation_amplitude_m/self.L * np.sin(
                    2*np.pi*parameters.corrugation_cycles*u)
            ids = [a]
            for point in X[1:-1]:
                ids.append(len(positions)); positions.append(point)
            ids.append(b)
            self.walls.append(dict(name=name, junctions=[a,b],
                                   nodes=np.array(ids), radial=radial))
        self.X = np.array(positions)
        self.nn = len(self.X); self.ndof = 2*self.nn+6
        self.q0 = np.r_[self.X.ravel(), np.zeros(6)]
        segment_nodes = []
        for wall in self.walls:
            wall['segments'] = np.arange(len(segment_nodes),
                                         len(segment_nodes)+len(wall['nodes'])-1)
            segment_nodes.extend(zip(wall['nodes'][:-1],wall['nodes'][1:]))
        self.segment_nodes = np.array(segment_nodes)
        self.ns = len(segment_nodes)
        self.D = np.zeros((2*self.ns,2*self.nn))
        for i,(a,b) in enumerate(segment_nodes):
            self.D[2*i:2*i+2,2*a:2*a+2] = -np.eye(2)
            self.D[2*i:2*i+2,2*b:2*b+2] = np.eye(2)
        d = (self.D@self.X.ravel()).reshape(-1,2)
        self.l0 = np.linalg.norm(d,axis=1)
        self.theta0 = np.arctan2(d[:,1],d[:,0])
        self.axial = self.EA*self.L**2/self.EI/self.l0
        theta_rows=[]; phi_rows=[]; stiffness=[]; locations=[]
        for wi,wall in enumerate(self.walls):
            ss=wall['segments']; vv=wall['junctions']
            for i,j in zip(ss[:-1],ss[1:]):
                row=np.zeros(self.ns); row[i]=-1; row[j]=1
                theta_rows.append(row); phi_rows.append(np.zeros(6))
                stiffness.append(2/(self.l0[i]+self.l0[j]))
                locations.append(('interior',wi,int(i)))
            for i,v in zip(ss[[0,-1]],vv):
                row=np.zeros(self.ns); row[i]=1
                prow=np.zeros(6); prow[v]=-1
                theta_rows.append(row); phi_rows.append(prow)
                stiffness.append(2/self.l0[i])
                locations.append(('junction',wi,int(v)))
        self.B=np.array(theta_rows); self.P=np.array(phi_rows)
        self.k=np.array(stiffness); self.locations=locations
        self.b0=self.B@self.theta0
        self.C=self.B.T@(self.k[:,None]*self.B)
        self.Cphi=self.B.T@(self.k[:,None]*self.P)
        self.Hphi=self.P.T@(self.k[:,None]*self.P)
        # Six junction z positions receive prescribed radial displacement.
        # These are point supports, not continuous contacting flat platens.
        # A single x anchor removes rigid translation; x and all rotations
        # otherwise remain free. Wall interiors have no displacement control.
        self.fixed=np.array([0]+[2*i+1 for i in range(6)])
        self.free=np.setdiff1d(np.arange(self.ndof),self.fixed)
        self.top_y=np.array([7,9,11]); self.bottom_y=np.array([1,3,5])
        self.cells=[[(3,1),(1,1),(5,-1),(0,-1)],
                    [(4,1),(2,1),(6,-1),(1,-1)]]

    def evaluate(self,q,hessian=False):
        q=np.asarray(q); x=q[:2*self.nn].reshape(-1,2)
        phi=q[2*self.nn:]
        d=(self.D@x.ravel()).reshape(-1,2)
        ell=np.linalg.norm(d,axis=1)
        if np.min(ell)<1e-8: raise ValueError('Degenerate wall segment')
        t=d/ell[:,None]; theta=np.arctan2(d[:,1],d[:,0])
        r=wrap(self.B@theta+self.P@phi-self.b0)
        stretch=ell-self.l0
        ea=.5*np.dot(self.axial,stretch**2)
        eb=.5*np.dot(self.k,r**2)
        moment=self.B.T@(self.k*r)
        jt=np.c_[-d[:,1],d[:,0]]/ell[:,None]**2
        gd=self.axial[:,None]*stretch[:,None]*t+moment[:,None]*jt
        gx=self.D.T@gd.ravel(); gp=self.P.T@(self.k*r)
        g=np.r_[gx,gp]
        bend_strain=np.abs(r*self.k)*self.p.wall_thickness_m/(2*self.L)
        parts=dict(axial_J=float(ea*self.unit),bending_J=float(eb*self.unit),
                   max_axial_strain=float(np.max(abs(stretch/self.l0))),
                   max_outer_fibre_bending_increment=float(bend_strain.max()),
                   max_junction_rotation_rad=float(np.max(abs(phi))))
        if not hessian: return float(ea+eb),g,parts
        Hd=np.zeros((2*self.ns,2*self.ns)); Jt=np.zeros((self.ns,2*self.ns))
        for i in range(self.ns):
            a,b=d[i]; l=ell[i]
            ht=np.array([[2*a*b,b*b-a*a],[b*b-a*a,-2*a*b]])/l**4
            ha=self.axial[i]*((1-self.l0[i]/l)*np.eye(2)
                 +self.l0[i]/l*np.outer(t[i],t[i]))
            Hd[2*i:2*i+2,2*i:2*i+2]=ha+moment[i]*ht
            Jt[i,2*i:2*i+2]=jt[i]
        H=np.zeros((self.ndof,self.ndof))
        angle_gradient=Jt@self.D
        H[:2*self.nn,:2*self.nn]=self.D.T@Hd@self.D+angle_gradient.T@self.C@angle_gradient
        cross=angle_gradient.T@self.Cphi
        H[:2*self.nn,2*self.nn:]=cross
        H[2*self.nn:,:2*self.nn]=cross.T
        H[2*self.nn:,2*self.nn:]=self.Hphi
        return float(ea+eb),g,H,parts

    def control(self,q,strain):
        out=q.copy(); x=out[:2*self.nn].reshape(-1,2)
        old=np.mean(x[3:6,1])-np.mean(x[:3,1])
        x[:,1] *= (1+strain)/old
        x[:3,1]=0; x[3:6,1]=1+strain; x[0,0]=0
        return out

    def reduced_gradient(self,g):
        return g[self.free]

    def reduced_hessian(self,H):
        return H[np.ix_(self.free,self.free)]

    def increment(self,q,dq,scale=1.):
        out=q.copy();out[self.free]+=scale*dq
        return out

    def support_gradient(self,g):
        support=np.zeros_like(g);support[self.fixed]=g[self.fixed]
        return support

    def solve(self,strain,initial=None,maxiter=90):
        q=self.control(self.q0 if initial is None else initial,strain)
        trace=[]; f=self.free; passed=False
        for it in range(maxiter):
            e,g,H,parts=self.evaluate(q,hessian=True)
            gf=self.reduced_gradient(g); hf=self.reduced_hessian(H)
            eig=eigh(hf,eigvals_only=True,subset_by_index=[0,0])[0]
            res=float(np.max(abs(gf)))
            row=dict(iteration=it,residual=res,minimum_eigenvalue=float(eig));trace.append(row)
            if res<2e-8 and eig>0: passed=True;break
            if eig<0 and res<1e-5:
                _,V=eigh(hf,subset_by_index=[0,0]); v=V[:,0]
                step=.003/max(np.max(abs(v)),1e-30)
                escaped=False
                for ls in range(24):
                    candidates=[]
                    for sign in (-1.,1.):
                        trial=self.increment(q,v,sign*step)
                        try: et=self.evaluate(trial)[0]
                        except ValueError: et=np.inf
                        candidates.append((et,trial,sign))
                    et,trial,sign=min(candidates,key=lambda a:a[0])
                    if et<e-1e-13:
                        q=trial;escaped=True
                        row.update(negative_curvature_step=True,sign=sign)
                        break
                    step*=.5
                if escaped: continue
            scale=max(float(np.max(abs(np.diag(hf)))),1.)
            shift=-eig+1e-8*scale if eig<=0 else 0.
            dq=solve(hf+shift*np.eye(len(f)),-gf,assume_a='pos')
            descent=gf@dq; alpha=min(1.,.02/max(np.max(abs(dq)),1e-30))
            accepted=False
            for ls in range(36):
                trial=self.increment(q,dq,alpha)
                try: et=self.evaluate(trial)[0]
                except ValueError: et=np.inf
                if et<=e+1e-4*alpha*descent+1e-15: accepted=True;break
                alpha*=.5
            row.update(metric_shift=float(shift),backtracks=ls,step_fraction=float(alpha))
            if not accepted:break
            q=trial
        e,g,H,parts=self.evaluate(q,hessian=True)
        eig=eigh(self.reduced_hessian(H),eigvals_only=True,subset_by_index=[0,0])[0]
        x=q[:2*self.nn].reshape(-1,2)
        support=self.support_gradient(g)
        sf=support[:2*self.nn].reshape(-1,2)*self.unit/self.L
        moment=np.sum((x[:,0]*sf[:,1]-x[:,1]*sf[:,0])*self.L)
        report=dict(engineering_strain=float(strain),passed=bool(passed),iterations=len(trace),
                    energy_J=float(e*self.unit),radial_reaction_N=float(np.sum(g[self.top_y])*self.unit/self.L),
                    dimensionless_residual=float(np.max(abs(self.reduced_gradient(g)))),
                    minimum_free_Hessian_eigenvalue=float(eig),
                    support_force_balance_N=np.sum(sf,axis=0).tolist(),
                    support_moment_balance_Nm=float(moment),
                    max_free_junction_moment_Nm=float(np.max(abs(g[2*self.nn:]))*self.unit),
                    max_free_position_residual_N=float(np.max(abs((g-support)[:2*self.nn]))*self.unit/self.L),
                    junction_rotations_rad=q[2*self.nn:].tolist(),**parts,trace=trace)
        return q,report

    def cell_polygons(self,q):
        x=q[:2*self.nn].reshape(-1,2)
        return [np.vstack([x[self.walls[wi]['nodes']][::direction][:-1]
                  for wi,direction in boundary]) for boundary in self.cells]

    def topology(self):
        return dict(parameters=asdict(self.p),axes=['transverse_x','radial_z'],
                    unique_wall_count=len(self.walls),junction_count=6,
                    position_node_count=self.nn,degrees_of_freedom=self.ndof,
                    free_degrees_of_freedom=len(self.free),
                    walls=[dict(name=w['name'],junctions=w['junctions'],
                                nodes=w['nodes'].tolist(),segments=w['segments'].tolist()) for w in self.walls],
                    cells=self.cells,shared_wall_index=1,
                    boundary_conditions='Six prescribed radial junction coordinates; one x translation anchor; all other positions and six junction rotations free; no continuous platen contact')
