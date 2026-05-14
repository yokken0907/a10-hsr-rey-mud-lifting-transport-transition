from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any
import numpy as np
from scipy.optimize import root_scalar


@dataclass
class TruthParams:
    """
    Axisymmetric 2D truth-model (radial BVP form) for a fully developed,
    steady vertical pipe flow with shear-induced migration of coarse BCP
    particles in a fine Herschel-Bulkley matrix.

    This v1 implementation is the *base closure generator* for a_perp = 0.
    The unknowns are the radial axial velocity profile u(r) and the radial
    coarse-particle profile phi_c(r). The solver alternates between:
      1) momentum solve for u(r) at fixed phi_c(r)
      2) pseudo-time relaxation of phi_c(r) at fixed u(r)

    The outputs are the 1.5D closure quantities:
      - phi_bar   : mean coarse-particle fraction
      - beta_eq   : segregation amplitude (projection on Psi = 1 - 2 r^2/R^2)
      - tau_w     : wall shear stress
      - r_p       : plug radius

    Notes
    -----
    - The constitutive model is regularized Herschel-Bulkley.
    - Gravity is not solved as an absolute-pressure problem here; instead,
      a dynamic pressure gradient G = -dp_dyn/dz is adjusted so that the
      imposed bulk velocity is satisfied.
    - The migration model follows a Phillips-type diffusive-flux form.
    - This code is intended as a robust base generator for closure tables,
      not as a final certification-grade solver.
    """

    # Geometry / mesh
    R: float = 0.125
    Nr: int = 301

    # Mixture composition inputs
    phi_tot: float = 0.28          # total solids fraction (fine + coarse)
    phi_c_bulk_target: float = 0.12  # target bulk coarse fraction

    # Densities
    rho_l: float = 1030.0
    rho_c: float = 2700.0
    rho_fstar: float = 1600.0

    # Environment
    T: float = 2.0                 # degC
    P: float = 55e6                # Pa
    d_p: float = 80e-6             # representative coarse BCP size [m]

    # HB base matrix parameters (at reference composition)
    tau_y_m0: float = 8.0          # Pa
    K_m0: float = 0.35             # Pa.s^n
    n_m0: float = 0.60
    T_ref: float = 20.0

    # Composition sensitivity for the fine matrix
    a_tau_f: float = 7.0
    a_K_f: float = 4.0
    a_n_f: float = 0.14
    b_tau_T: float = 0.05
    b_K_T: float = 0.03

    # Coarse-particle amplification of matrix rheology
    y1: float = 2.5
    y2: float = 10.0
    k1: float = 1.5
    k2: float = 4.0
    n1: float = 0.12

    # Migration coefficients
    k_c: float = 0.41
    k_mu: float = 0.62
    k_phi: float = 0.20

    # Regularization and numerics
    eps_reg: float = 1e-6
    pseudo_dt: float = 1e-4
    pseudo_diff: float = 2e-9      # extra numerical diffusion for phi_c
    max_iter: int = 4000
    tol: float = 1e-8
    under_relax_phi: float = 0.65

    # Search range for dynamic pressure gradient G [Pa/m]
    G_min: float = 1e1
    G_max: float = 1e7

    def __post_init__(self):
        if self.Nr < 11:
            raise ValueError("Nr must be >= 11")
        if not (0.0 < self.phi_c_bulk_target < self.phi_tot):
            raise ValueError("phi_c_bulk_target must lie in (0, phi_tot)")
        self.r = np.linspace(0.0, self.R, self.Nr)
        self.dr = self.r[1] - self.r[0]


def chi_y(phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    return 1.0 + prm.y1 * phi_c + prm.y2 * phi_c**2


def chi_K(phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    return 1.0 + prm.k1 * phi_c + prm.k2 * phi_c**2


def chi_n(phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    return prm.n1 * phi_c


def phi_f(phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    return np.clip(prm.phi_tot - phi_c, 0.0, prm.phi_tot)


def hb_params(phi_c: np.ndarray, prm: TruthParams):
    pf = phi_f(phi_c, prm)
    tau_y_m = prm.tau_y_m0 * np.exp(prm.a_tau_f * pf + prm.b_tau_T * (prm.T_ref - prm.T))
    K_m = prm.K_m0 * np.exp(prm.a_K_f * pf + prm.b_K_T * (prm.T_ref - prm.T))
    n_m = np.clip(prm.n_m0 - prm.a_n_f * pf, 0.20, 0.95)

    tau_y = tau_y_m * chi_y(phi_c, prm)
    K = K_m * chi_K(phi_c, prm)
    n = np.clip(n_m - chi_n(phi_c, prm), 0.20, 0.95)
    return tau_y, K, n


def mixture_density(phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    pf = phi_f(phi_c, prm)
    return ((1.0 - prm.phi_tot) * prm.rho_l
            + pf * prm.rho_fstar
            + phi_c * prm.rho_c)


def regularized_shear_rate(du_dr: np.ndarray, prm: TruthParams) -> np.ndarray:
    return np.sqrt(du_dr * du_dr + prm.eps_reg**2)


def mu_eff_from_profile(du_dr: np.ndarray, phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    gdot = regularized_shear_rate(du_dr, prm)
    tau_y, K, n = hb_params(phi_c, prm)
    mu = tau_y / gdot + K * np.power(gdot, n - 1.0)
    return np.clip(mu, 1e-6, 1e9)


def integrate_tau_from_G(G: float, r: np.ndarray) -> np.ndarray:
    """
    For fully developed axisymmetric flow with dynamic gradient G = -dp_dyn/dz,
    the radial shear stress satisfies d(r tau)/dr = r G, so tau = G r / 2.
    This exact form is used for the base generator.
    """
    return 0.5 * G * r


def reconstruct_velocity_from_tau(tau: np.ndarray, phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    r = prm.r
    dr = prm.dr
    tau_y, K, n = hb_params(phi_c, prm)

    gdot = np.zeros_like(r)
    yielded = np.abs(tau) > tau_y
    gdot[yielded] = np.power(
        (np.abs(tau[yielded]) - tau_y[yielded]) / np.maximum(K[yielded], 1e-14),
        1.0 / np.maximum(n[yielded], 1e-14),
    )

    # du/dr is negative in upward flow when tau>0
    u = np.zeros_like(r)
    u[-1] = 0.0  # no-slip at the wall
    for j in range(len(r) - 2, -1, -1):
        u[j] = u[j + 1] + 0.5 * (gdot[j + 1] + gdot[j]) * dr
    return u


def bulk_velocity(u: np.ndarray, prm: TruthParams) -> float:
    return 2.0 / prm.R**2 * np.trapezoid(u * prm.r, prm.r)


def solve_velocity_given_phi(phi_c: np.ndarray, U_bulk_target: float, prm: TruthParams):
    def residual(G):
        tau = integrate_tau_from_G(G, prm.r)
        u = reconstruct_velocity_from_tau(tau, phi_c, prm)
        return bulk_velocity(u, prm) - U_bulk_target

    fmin = residual(prm.G_min)
    fmax = residual(prm.G_max)
    if fmin * fmax > 0.0:
        # bracket expansion
        Gmin, Gmax = prm.G_min, prm.G_max
        for _ in range(20):
            Gmin *= 0.5
            Gmax *= 2.0
            fmin = residual(Gmin)
            fmax = residual(Gmax)
            if fmin * fmax <= 0.0:
                break
        else:
            raise RuntimeError("Could not bracket pressure gradient G for target bulk velocity")
    else:
        Gmin, Gmax = prm.G_min, prm.G_max

    sol = root_scalar(residual, bracket=[Gmin, Gmax], method="bisect", xtol=1e-8, rtol=1e-8)
    if not sol.converged:
        raise RuntimeError("root_scalar did not converge for G")

    G_star = sol.root
    tau = integrate_tau_from_G(G_star, prm.r)
    u = reconstruct_velocity_from_tau(tau, phi_c, prm)
    Ub = bulk_velocity(u, prm)
    return G_star, u, tau, Ub


def migration_flux(phi_c: np.ndarray, u: np.ndarray, prm: TruthParams) -> np.ndarray:
    r, dr = prm.r, prm.dr
    du_dr = np.gradient(u, dr, edge_order=2)
    gdot = regularized_shear_rate(du_dr, prm)
    mu = mu_eff_from_profile(du_dr, phi_c, prm)

    dgdot_dr = np.gradient(gdot, dr, edge_order=2)
    dlnmu_dr = np.gradient(np.log(np.maximum(mu, 1e-20)), dr, edge_order=2)
    dphi_dr = np.gradient(phi_c, dr, edge_order=2)

    J_coll = -prm.k_c * prm.d_p**2 * phi_c**2 * dgdot_dr
    J_mu = -prm.k_mu * prm.d_p**2 * phi_c**2 * gdot * dlnmu_dr
    J_phi = -prm.k_phi * prm.d_p**2 * phi_c * gdot * dphi_dr

    J = J_coll + J_mu + J_phi

    # symmetry at center and no-flux at wall
    J[0] = 0.0
    J[-1] = 0.0
    return J


def enforce_bulk_coarse_fraction(phi_c: np.ndarray, prm: TruthParams) -> np.ndarray:
    target = prm.phi_c_bulk_target
    current = 2.0 / prm.R**2 * np.trapezoid(phi_c * prm.r, prm.r)
    corrected = phi_c + (target - current)
    corrected = np.clip(corrected, 0.0, prm.phi_tot)

    # one more conservative correction after clipping
    current2 = 2.0 / prm.R**2 * np.trapezoid(corrected * prm.r, prm.r)
    corrected = np.clip(corrected + (target - current2), 0.0, prm.phi_tot)
    return corrected


def update_phi_pseudotime(phi_c: np.ndarray, u: np.ndarray, prm: TruthParams) -> np.ndarray:
    r, dr = prm.r, prm.dr
    J = migration_flux(phi_c, u, prm)

    rJ = r * J
    drJ_dr = np.gradient(rJ, dr, edge_order=2)
    divJ = np.zeros_like(phi_c)
    divJ[1:] = drJ_dr[1:] / np.maximum(r[1:], 1e-20)
    # regular center limit
    divJ[0] = 2.0 * (J[1] - J[0]) / dr

    # light radial numerical diffusion for robustness
    dphi_dr = np.gradient(phi_c, dr, edge_order=2)
    rdphi = r * dphi_dr
    drdphi_dr = np.gradient(rdphi, dr, edge_order=2)
    lap_axi = np.zeros_like(phi_c)
    lap_axi[1:] = drdphi_dr[1:] / np.maximum(r[1:], 1e-20)
    lap_axi[0] = 2.0 * (dphi_dr[1] - dphi_dr[0]) / dr

    phi_trial = phi_c - prm.pseudo_dt * divJ + prm.pseudo_dt * prm.pseudo_diff * lap_axi
    phi_trial = np.clip(phi_trial, 0.0, prm.phi_tot)
    phi_trial = enforce_bulk_coarse_fraction(phi_trial, prm)

    phi_new = prm.under_relax_phi * phi_trial + (1.0 - prm.under_relax_phi) * phi_c
    return np.clip(phi_new, 0.0, prm.phi_tot)


def extract_closures(phi_c: np.ndarray, u: np.ndarray, tau: np.ndarray, prm: TruthParams) -> Dict[str, Any]:
    r = prm.r
    phi_bar = 2.0 / prm.R**2 * np.trapezoid(phi_c * r, r)

    Psi = 1.0 - 2.0 * (r / prm.R) ** 2
    beta_eq = (6.0 / prm.R**2) * np.trapezoid((phi_c - phi_bar) * Psi * r, r)

    du_dr = np.gradient(u, prm.dr, edge_order=2)
    mu = mu_eff_from_profile(du_dr, phi_c, prm)
    tau_num = mu * du_dr
    tau_w = abs(tau_num[-1])

    tau_y, _, _ = hb_params(phi_c, prm)
    yielded = np.abs(tau) > tau_y
    if np.all(yielded):
        r_p = 0.0
    else:
        r_p = r[np.where(~yielded)[0][-1]]

    return {
        "phi_bar": float(phi_bar),
        "beta_eq": float(beta_eq),
        "tau_w": float(tau_w),
        "r_p": float(r_p),
        "phi_center": float(phi_c[0]),
        "phi_wall": float(phi_c[-1]),
        "u_center": float(u[0]),
        "u_wall": float(u[-1]),
        "yielded_fraction": float(np.mean(yielded.astype(float))),
    }


def solve_truth_case(U_bulk_target: float, prm: TruthParams) -> Dict[str, Any]:
    """
    Main base-truth solve.

    Parameters
    ----------
    U_bulk_target : float
        Target bulk axial velocity [m/s].
    prm : TruthParams
        Solver / material parameters.

    Returns
    -------
    dict
        Contains closures and full radial profiles.
    """
    # initialize around the target bulk coarse loading with a weak center bias
    rr = prm.r / prm.R
    phi_c = prm.phi_c_bulk_target * (1.0 + 0.05 * (1.0 - 2.0 * rr**2))
    phi_c = np.clip(phi_c, 0.0, prm.phi_tot)
    phi_c = enforce_bulk_coarse_fraction(phi_c, prm)

    last_err = np.inf
    history = []

    for it in range(prm.max_iter):
        G, u, tau, Ub = solve_velocity_given_phi(phi_c, U_bulk_target, prm)
        phi_new = update_phi_pseudotime(phi_c, u, prm)
        err = np.linalg.norm(phi_new - phi_c, ord=np.inf)
        history.append(err)
        phi_c = phi_new

        if err < prm.tol:
            break
        if err > 10.0 * last_err and it > 10:
            # defensive dt reduction if pseudo-time becomes too aggressive
            prm.pseudo_dt *= 0.5
        last_err = err
    else:
        raise RuntimeError("Truth-case iteration did not converge within max_iter")

    closures = extract_closures(phi_c, u, tau, prm)
    closures.update({
        "G": float(G),
        "U_bulk": float(Ub),
        "r": prm.r.copy(),
        "u_profile": u.copy(),
        "phi_profile": phi_c.copy(),
        "tau_profile": tau.copy(),
        "iter_history": np.array(history),
        "params": prm,
    })
    return closures


def demo_case():
    prm = TruthParams(
        max_iter=10000,         # 最大ループ回数を4000から10000に増やす
        tol=1e-5,               # 許容誤差を少し現実的なラインに緩める
        under_relax_phi=0.2,    # 1回の計算での変化量を抑え、振動を防ぐ
        pseudo_dt=5e-5          # 擬似時間の進み方を小さくして慎重に計算する
    )
    out = solve_truth_case(U_bulk_target=1.20, prm=prm)

    print("Axisymmetric truth-model v1 complete")
    print(f"bulk coarse fraction phi_bar: {out['phi_bar']:.6f}")
    print(f"segregation amplitude beta_eq: {out['beta_eq']:.6f}")
    print(f"wall shear stress tau_w: {out['tau_w']:.3f} Pa")
    print(f"plug radius r_p: {out['r_p']:.6f} m")
    print(f"center coarse fraction: {out['phi_center']:.6f}")
    print(f"wall coarse fraction: {out['phi_wall']:.6f}")
    print(f"required dynamic gradient G: {out['G']:.3f} Pa/m")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 3, figsize=(13, 4))
        axs[0].plot(out['r'], out['u_profile'])
        axs[0].set_xlabel('r [m]')
        axs[0].set_ylabel('u(r) [m/s]')
        axs[0].set_title('Velocity profile')

        axs[1].plot(out['r'], out['phi_profile'])
        axs[1].set_xlabel('r [m]')
        axs[1].set_ylabel('phi_c(r) [-]')
        axs[1].set_title('Coarse-particle profile')

        axs[2].semilogy(out['iter_history'])
        axs[2].set_xlabel('iteration')
        axs[2].set_ylabel('||phi^(n+1)-phi^n||_inf')
        axs[2].set_title('Convergence history')
        plt.tight_layout()
        plt.savefig('axisym_truth_generator_v1_demo.png', dpi=160)
        print('saved demo plot: axisym_truth_generator_v1_demo.png')
    except Exception:
        pass


if __name__ == "__main__":
    demo_case()
