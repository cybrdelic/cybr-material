"""SI Q4 plane-stress reference, analytic elastic energy/forces, no damage.

Engineering strain = [u_x,x, u_y,y, u_x,y+u_y,x].
Plane-stress D is diagonal for the deliberately zero-Poisson fixtures.
The 2x2 rule integrates the affine rectangular Q4 stiffness exactly.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def shape(xi, eta):
    signs = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    N = (1 + signs[:, 0]*xi)*(1 + signs[:, 1]*eta)/4
    dN = np.column_stack((signs[:, 0]*(1 + signs[:, 1]*eta),
                          signs[:, 1]*(1 + signs[:, 0]*xi)))/4
    return N, dN


def element_rule(hx, hy, D, thickness):
    """Return exact integrated stiffness and physical Gauss gradients."""
    K = np.zeros((8, 8))
    rule = []
    for xi in [-1/np.sqrt(3), 1/np.sqrt(3)]:
        for eta in [-1/np.sqrt(3), 1/np.sqrt(3)]:
            N, dN = shape(xi, eta)
            grad = dN @ np.diag([2/hx, 2/hy])
            B = np.zeros((3, 8))
            B[0, 0::2] = grad[:, 0]
            B[1, 1::2] = grad[:, 1]
            B[2, 0::2] = grad[:, 1]
            B[2, 1::2] = grad[:, 0]
            area_weight = hx*hy/4
            K += B.T @ D @ B * area_weight * thickness
            rule.append((N, grad, B, area_weight))
    return K, rule


@dataclass
class Model:
    n: int
    width: float
    half_height: float
    thickness: float
    D: np.ndarray

    def __post_init__(self):
        self.hx = self.width/self.n
        self.hy = self.half_height/self.n
        xx, yy = np.meshgrid(np.linspace(0, self.width, self.n+1),
                             np.linspace(0, self.half_height, self.n+1))
        self.xy = np.column_stack([xx.ravel(), yy.ravel()])
        ids = np.arange((self.n+1)**2).reshape(self.n+1, self.n+1)
        self.conn = np.column_stack([ids[:-1, :-1].ravel(), ids[:-1, 1:].ravel(),
                                     ids[1:, 1:].ravel(), ids[1:, :-1].ravel()])
        self.dofs = np.stack([2*self.conn, 2*self.conn+1], axis=-1).reshape(-1, 8)
        self.ke, self.rule = element_rule(self.hx, self.hy, self.D, self.thickness)
        rows = np.repeat(self.dofs, 8, axis=1).ravel()
        cols = np.tile(self.dofs, (1, 8)).ravel()
        values = np.broadcast_to(self.ke.ravel(), (len(self.conn), 64)).ravel()
        ndof = 2*len(self.xy)
        self.K = coo_matrix((values, (rows, cols)), shape=(ndof, ndof)).tocsr()
        self.top = ids[-1, :]
        self.bottom = ids[0, :]

    def solve(self, crack_index, delta):
        """Half symmetry model; tip is the first intact bottom node."""
        fixed_y = 2*np.concatenate([self.top, self.bottom[crack_index:]])+1
        # One harmless translation gauge: ux at top right. No side clamping.
        anchor = 2*self.top[-1]
        fixed = np.concatenate([fixed_y, [anchor]])
        free_mask = np.ones(self.K.shape[0], dtype=bool)
        free_mask[fixed] = False
        free = np.flatnonzero(free_mask)
        u = np.zeros(self.K.shape[0])
        u[2*self.top+1] = delta
        u[free] = spsolve(self.K[free][:, free], -(self.K @ u)[free], permc_spec='COLAMD')
        forces = self.K @ u
        ue = u[self.dofs]
        energy_half = float(0.5*np.einsum('ei,ij,ej->', ue, self.ke, ue))
        P = float(forces[2*self.top+1].sum())
        bottom_R = float(forces[2*self.bottom[crack_index:]+1].sum())
        # Full reflected coupon has two moving grips and separation Delta=2delta.
        return {"u": u, "forces": forces, "crack_index": crack_index,
                "a_m": crack_index*self.hx, "delta_m": delta,
                "P_N": P, "Delta_m": 2*delta,
                "compliance_m_N": 2*delta/P,
                "energy_half_J": energy_half, "energy_full_J": 2*energy_half,
                "ramp_work_full_J": P*delta,
                "bottom_half_reaction_N": bottom_R,
                "horizontal_anchor_reaction_N": float(forces[anchor]),
                "relative_free_residual": float(np.linalg.norm(forces[free])/abs(P)),
                "relative_force_balance": abs(P+bottom_R)/abs(P),
                "relative_work_error": abs(2*energy_half-P*delta)/(P*delta),
                "matrix_energy_full_J": float(u @ (self.K @ u))}

    def domain_J(self, state, outer_radii, inner_fraction=0.5):
        """J=2 int_upper (sigma_ij u_i,1-W delta_1j) q_,j dA.

        q is one near the tip and zero outside the chosen domain; no thickness
        factor occurs: stress [Pa] times area [m2] times grad(q) [1/m] gives
        J [N/m]=[J/m2]. The 2 restores the lower reflected half.
        """
        ue = state['u'][self.dofs].reshape(-1, 4, 2)
        distance = np.linalg.norm(self.xy-[state['a_m'], 0.0], axis=1)
        Js = []
        for radius in outer_radii:
            assert radius < min(state['a_m'], self.width-state['a_m'], self.half_height)
            inner = inner_fraction*radius
            qn = np.clip((radius-distance)/(radius-inner), 0.0, 1.0)
            qe = qn[self.conn]
            Jhalf = 0.0
            for N, grad, B, area in self.rule:
                # du[e,i,j] = partial u_i / partial x_j.
                du = np.einsum('eai,aj->eij', ue, grad)
                eps = np.stack([du[:, 0, 0], du[:, 1, 1], du[:, 0, 1]+du[:, 1, 0]], axis=1)
                sig = eps @ self.D.T
                W = 0.5*np.einsum('ei,ei->e', eps, sig)
                # FULL contraction over displacement component i; no trace substitute.
                A1 = sig[:, 0]*du[:, 0, 0] + sig[:, 2]*du[:, 1, 0] - W
                A2 = sig[:, 2]*du[:, 0, 0] + sig[:, 1]*du[:, 1, 0]
                dq = qe @ grad
                Jhalf += float(np.sum(A1*dq[:, 0]+A2*dq[:, 1])*area)
            Js.append(2*Jhalf)
        return np.array(Js)


def serial_state(state):
    return {k: v for k, v in state.items() if k not in ('u', 'forces')}


def relative(a, b):
    return abs(a-b)/abs(b)
