from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class ModelParams15D:
    # ------------------------
    # Geometry / mesh
    # ------------------------
    L: float = 5500.0
    N: int = 550
    D: float = 0.25

    # ------------------------
    # Fine/coarse phase densities
    # ------------------------
    rho_l: float = 1030.0          # interstitial liquid
    rho_fstar: float = 1600.0      # effective fine-matrix density
    rho_cstar: float = 2700.0      # coarse BCP density
    d_rep: float = 80e-6           # representative coarse BCP size [m]
    T: float = 2.0                 # degC

    # ------------------------
    # Background fine matrix loading
    # ------------------------
    phi_f_bg: float = 0.10         # constant fine matrix for v3 skeleton
    phi_c_inlet: float = 0.18      # transported coarse-phase inlet fraction
    phi_max: float = 0.60

    # ------------------------
    # Acoustic / transient
    # ------------------------
    a_wave: float = 700.0
    unsteady_damping_beta: float = 2.0

    # ------------------------
    # Steady operating point
    # ------------------------
    u_ss: float = 1.60
    t_blackout: float = 5.0
    tau_coastdown: float = 2.0
    p_out_dyn: float = 0.0
    check_valve: bool = True

    # ------------------------
    # Herschel-Bulkley baseline (matrix + coarse correction)
    # ------------------------
    phi_ref: float = 0.28
    tau_y_ref: float = 12.0
    K_ref: float = 0.45
    n_ref: float = 0.55
    a_tau: float = 7.0
    a_K: float = 3.0
    a_n: float = 0.10
    b_tau_T: float = 0.06
    b_K_T: float = 0.03
    T_ref: float = 20.0

    # coarse-phase correction multipliers
    chi_tau_c: float = 2.0
    chi_K_c: float = 1.3
    chi_n_c: float = 0.08

    # ------------------------
    # Segregation closures
    # ------------------------
    tau_m: float = 3.0             # migration relaxation time [s]
    D_beta: float = 0.03           # axial diffusion of segregation mode
    kappa_a: float = 8.0e-4        # coupling from transverse acceleration to beta
    beta_max: float = 0.12
    beta_relax_gain: float = 0.08
    beta_Bn_gain: float = 0.02
    beta_acc_gain: float = 0.01

    # ------------------------
    # Critical-velocity closures
    # ------------------------
    Re_transition: float = 3000.0
    Y_crit: float = 0.20
    hinder_exp: float = 4.65
    slip_sf: float = 1.15
    tau_mob_fac: float = 1.8
    ucrit_margin: float = 1.15

    # ------------------------
    # Deposition / erosion
    # ------------------------
    phi_c_eq_wall: float = 0.06
    U_ero: float = 1.8
    k_dep: float = 4.5e-4
    k_ero: float = 1.0e-2
    c_dep: float = 0.15
    c_ero: float = 0.03
    delta_max_frac: float = 0.45

    # ------------------------
    # Structure / VIV modal surrogate
    # ------------------------
    M: int = 2
    mode_freq_hz: tuple[float, ...] = (0.18, 0.52)
    mode_zeta: tuple[float, ...] = (0.03, 0.04)
    T0: float = 2.5e5             # baseline effective top tension proxy [N]
    m_struct: float = 180.0       # lineic structural mass [kg/m]
    m_added: float = 120.0        # lineic added mass [kg/m]
    c_rho_modal: float = 5.0e-4   # placeholder coupling gains
    c_p_modal: float = 1.5e-7
    c_u_modal: float = 5.0e-3
    Uo_ext: float = 0.55          # external cross-flow current [m/s]
    St: float = 0.20
    C_L0: float = 0.8
    wake_eps: float = 0.25
    wake_A: float = 0.8
    viv_force_gain: float = 0.6

    # ------------------------
    # Numerics
    # ------------------------
    t_end: float = 40.0
    cfl: float = 0.40
    dt_max: float = 0.01
    c_art: float = 0.05
    out_absorb: float = 1.0
    in_impedance_fac: float = 1.2
    outlet_relax: float = 0.18
    sponge_n: int = 60

    # ------------------------
    # Monitoring section
    # ------------------------
    i_mon0: int = 0
    i_mon1: int = 80

    def __post_init__(self):
        self.dx = self.L / self.N
        self.x = np.linspace(0.5*self.dx, self.L - 0.5*self.dx, self.N)
        self.A0 = 0.25*np.pi*self.D**2
        self.delta_max = self.delta_max_frac * 0.5 * self.D
        self.R0 = 0.5*self.D
        self.mode_omega = 2.0*np.pi*np.asarray(self.mode_freq_hz, dtype=float)
        self.mode_zeta = np.asarray(self.mode_zeta, dtype=float)
        assert len(self.mode_omega) == self.M
        assert len(self.mode_zeta) == self.M

        # simple sine modes on [0, L]
        self.phi_modes = np.zeros((self.M, self.N))
        for k in range(self.M):
            n = k + 1
            self.phi_modes[k, :] = np.sin(n*np.pi*self.x / self.L)
        # L2 normalize on the discrete grid
        for k in range(self.M):
            norm = np.sqrt(np.sum(self.phi_modes[k]**2) * self.dx)
            self.phi_modes[k] /= max(norm, 1.0e-12)

        self.p_in0: float | None = None


@dataclass
class State15D:
    p: np.ndarray
    U: np.ndarray
    phi_c_bar: np.ndarray
    beta: np.ndarray
    delta: np.ndarray
    q: np.ndarray
    qdot: np.ndarray
    w: np.ndarray
    wdot: np.ndarray


# ----------------------------
# Basic field utilities
# ----------------------------

def clip_phi(phi: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    return np.clip(phi, 0.0, prm.phi_max)


def wall_center_concentrations(phi_bar: np.ndarray, beta: np.ndarray, prm: ModelParams15D):
    phi0 = clip_phi(phi_bar + beta, prm)
    phiw = clip_phi(phi_bar - beta, prm)
    return phi0, phiw


def effective_diameter(delta: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    return np.maximum(prm.D - 2.0*delta, 1.0e-4)


def effective_area_ratio(delta: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    return np.clip((effective_diameter(delta, prm) / prm.D)**2, 1.0e-4, 1.0)


def mixture_density(phi_c_bar: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    phi_tot = np.clip(prm.phi_f_bg + phi_c_bar, 0.0, prm.phi_max)
    return ((1.0 - phi_tot)*prm.rho_l + prm.phi_f_bg*prm.rho_fstar + phi_c_bar*prm.rho_cstar)


# ----------------------------
# HB closures with segregation
# ----------------------------

def hb_matrix_params(phi_tot: np.ndarray, prm: ModelParams15D):
    tau_y = prm.tau_y_ref * np.exp(prm.a_tau*(phi_tot - prm.phi_ref) + prm.b_tau_T*(prm.T_ref - prm.T))
    K = prm.K_ref * np.exp(prm.a_K*(phi_tot - prm.phi_ref) + prm.b_K_T*(prm.T_ref - prm.T))
    n = np.clip(prm.n_ref - prm.a_n*(phi_tot - prm.phi_ref), 0.20, 0.95)
    return tau_y, K, n


def hb_wall_params(phiw: np.ndarray, prm: ModelParams15D):
    phi_tot_w = np.clip(prm.phi_f_bg + phiw, 0.0, prm.phi_max)
    tau_m, K_m, n_m = hb_matrix_params(phi_tot_w, prm)

    # coarse-phase multipliers at wall
    chi_tau = 1.0 + prm.chi_tau_c * phiw
    chi_K = 1.0 + prm.chi_K_c * phiw
    chi_n = prm.chi_n_c * phiw

    tau_yw = tau_m * chi_tau
    K_w = K_m * chi_K
    n_w = np.clip(n_m - chi_n, 0.15, 0.95)
    return tau_yw, K_w, n_w


def apparent_viscosity_hb_local(U: np.ndarray, D_eff: np.ndarray, phiw: np.ndarray, prm: ModelParams15D):
    tau_yw, K_w, n_w = hb_wall_params(phiw, prm)
    gamma_w = 8.0*np.abs(U) / np.maximum(D_eff, 1.0e-6)
    mu = tau_yw/(gamma_w + 1.0e-9) + K_w*np.power(gamma_w + 1.0e-9, n_w - 1.0)
    return np.clip(mu, 1.0e-4, 1.0e6)


def wall_shear_hb(U: np.ndarray, D_eff: np.ndarray, phiw: np.ndarray, prm: ModelParams15D):
    tau_yw, K_w, n_w = hb_wall_params(phiw, prm)
    gamma_w = 8.0*np.abs(U) / np.maximum(D_eff, 1.0e-6)
    tauw = tau_yw + K_w*np.power(gamma_w + 1.0e-9, n_w)
    return tauw, gamma_w, tau_yw


def plug_radius(phi0: np.ndarray, tauw: np.ndarray, prm: ModelParams15D):
    phi_tot0 = np.clip(prm.phi_f_bg + phi0, 0.0, prm.phi_max)
    tau_y0, _, _ = hb_matrix_params(phi_tot0, prm)
    rp = prm.R0 * np.clip(tau_y0 / np.maximum(tauw, 1.0e-9), 0.0, 1.0)
    return rp, tau_y0


# ----------------------------
# Segregation and Ucrit closures
# ----------------------------

def beta_equilibrium(phi_bar: np.ndarray, U: np.ndarray, gamma_w: np.ndarray, tau_yw: np.ndarray,
                     a_perp: np.ndarray, prm: ModelParams15D):
    Bn = tau_yw / np.maximum(gamma_w, 1.0e-9)
    raw = prm.beta_relax_gain*(U - 1.0) + prm.beta_Bn_gain*Bn + prm.beta_acc_gain*a_perp
    return np.clip(prm.beta_max * np.tanh(raw), -prm.beta_max, prm.beta_max)


def hindered_settling_velocity(phi_bar: np.ndarray, phi0: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    # coarse-core based surrogate; replace with 2D/column calibration later
    phi_tot0 = np.clip(prm.phi_f_bg + phi0, 0.0, prm.phi_max)
    tau_y0, K0, n0 = hb_matrix_params(phi_tot0, prm)
    gamma_settle = 1.0
    mu_settle = tau_y0/gamma_settle + K0*(gamma_settle**(n0 - 1.0))
    Y = tau_y0 / np.maximum((prm.rho_cstar - prm.rho_l)*9.81*prm.d_rep, 1.0e-9)
    trap = np.clip(1.0 - Y/prm.Y_crit, 0.0, 1.0)
    w0 = (prm.rho_cstar - prm.rho_l)*9.81*(prm.d_rep**2) / np.maximum(18.0*mu_settle, 1.0e-9)
    wh = w0 * np.power(np.clip(1.0 - phi_bar, 1.0e-6, 1.0), prm.hinder_exp) * trap
    return np.clip(wh, 0.0, 10.0)


def u_transition_segment(D_eff: np.ndarray, rho_m: np.ndarray, mu_w: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    Utr = prm.Re_transition * mu_w / np.maximum(rho_m * D_eff, 1.0e-9)
    return np.clip(Utr, 1.0e-3, 20.0)


def u_mobilization_segment(D_eff: np.ndarray, tauw: np.ndarray, tau_y0: np.ndarray, rp: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    # higher center yield / larger plug radius => harder mobilization
    tau_mob = prm.tau_mob_fac * tau_y0 * (1.0 + 0.8*rp / np.maximum(0.5*D_eff, 1.0e-9))
    gamma_req = np.power(np.maximum(tau_mob - tau_y0, 0.0) / np.maximum(0.45, 1.0e-12), 1.0/0.55)
    U = 0.125 * D_eff * gamma_req
    return np.clip(U, 0.0, 20.0)


def critical_transport_velocity_segment(phi_bar: np.ndarray, beta: np.ndarray, delta: np.ndarray,
                                        U: np.ndarray, rho_m: np.ndarray, prm: ModelParams15D):
    D_eff = effective_diameter(delta, prm)
    phi0, phiw = wall_center_concentrations(phi_bar, beta, prm)
    tauw, gamma_w, tau_yw = wall_shear_hb(U, D_eff, phiw, prm)
    mu_w = apparent_viscosity_hb_local(U, D_eff, phiw, prm)
    rp, tau_y0 = plug_radius(phi0, tauw, prm)

    Utr = u_transition_segment(D_eff, rho_m, mu_w, prm)
    Usl = prm.slip_sf * hindered_settling_velocity(phi_bar, phi0, prm)
    Umb = u_mobilization_segment(D_eff, tauw, tau_y0, rp, prm)
    Ucrit = prm.ucrit_margin * np.maximum.reduce([Utr, Usl, Umb])
    return Ucrit, tauw, gamma_w, tau_yw, mu_w, rp, tau_y0


# ----------------------------
# Boundary conditions
# ----------------------------

def pump_pressure_target(t: float, prm: ModelParams15D) -> float:
    if t < prm.t_blackout:
        return prm.p_in0
    return prm.p_in0 * np.exp(-(t - prm.t_blackout)/prm.tau_coastdown)


def quasi_steady_velocity_target(t: float, prm: ModelParams15D) -> float:
    return prm.u_ss * pump_pressure_target(t, prm) / max(prm.p_in0, 1.0e-9)


def inlet_bc(pL: float, uL: float, rhoL: float, t: float, prm: ModelParams15D):
    Z = max(rhoL*prm.a_wave, 1.0e-6)
    w_minus = uL - pL/Z
    p_tgt = pump_pressure_target(t, prm)
    u_ref = quasi_steady_velocity_target(t, prm)
    R = prm.in_impedance_fac * Z
    u_b = (w_minus + p_tgt/Z + (R/Z)*u_ref) / (1.0 + R/Z)
    p_b = p_tgt - R*(u_b - u_ref)
    if prm.check_valve and (t >= prm.t_blackout) and (u_b < 0.0):
        u_b = 0.0
        p_b = Z*(u_b - w_minus)
    return p_b, u_b


def outlet_bc(pR: float, uR: float, rhoR: float, t: float, prm: ModelParams15D):
    Z = max(rhoR*prm.a_wave, 1.0e-6)
    w_plus = uR + pR/Z
    u_ref = quasi_steady_velocity_target(t, prm)
    w_minus_bc = u_ref - prm.p_out_dyn/Z
    u_b = 0.5*(w_plus + w_minus_bc)
    p_b = 0.5*Z*(w_plus - w_minus_bc)
    p_b *= (1.0 - 0.35*prm.out_absorb)
    u_b = (1.0 - prm.outlet_relax)*u_b + prm.outlet_relax*u_ref
    return p_b, u_b


def extend_with_ghosts(p: np.ndarray, U: np.ndarray, phi_bar: np.ndarray, rho: np.ndarray, t: float, prm: ModelParams15D):
    N = prm.N
    p_ext = np.empty(N + 2)
    u_ext = np.empty(N + 2)
    phi_ext = np.empty(N + 2)
    rho_ext = np.empty(N + 2)

    p_ext[1:-1] = p
    u_ext[1:-1] = U
    phi_ext[1:-1] = phi_bar
    rho_ext[1:-1] = rho

    p_ext[0], u_ext[0] = inlet_bc(p[0], U[0], rho[0], t, prm)
    phi_ext[0] = prm.phi_c_inlet
    rho_ext[0] = rho[0]

    p_ext[-1], u_ext[-1] = outlet_bc(p[-1], U[-1], rho[-1], t, prm)
    phi_ext[-1] = phi_bar[-1]
    rho_ext[-1] = rho[-1]
    return p_ext, u_ext, phi_ext, rho_ext


def apply_artificial_viscosity(q: np.ndarray, eps: float) -> np.ndarray:
    qn = q.copy()
    qn[1:-1] += eps * (q[:-2] - 2.0*q[1:-1] + q[2:])
    return qn


# ----------------------------
# Fluid step
# ----------------------------

def fluid_step(p: np.ndarray, U: np.ndarray, phi_bar: np.ndarray, beta: np.ndarray, delta: np.ndarray,
               dt: float, t: float, prm: ModelParams15D):
    rho = mixture_density(phi_bar, prm)
    D_eff = effective_diameter(delta, prm)

    Ucrit, tauw, gamma_w, tau_yw, mu_w, rp, tau_y0 = critical_transport_velocity_segment(phi_bar, beta, delta, U, rho, prm)
    p_ext, u_ext, _, rho_ext = extend_with_ghosts(p, U, phi_bar, rho, t, prm)

    qLp = p_ext[:-1]
    qRp = p_ext[1:]
    qLu = u_ext[:-1]
    qRu = u_ext[1:]
    rho_face = 0.5*(rho_ext[:-1] + rho_ext[1:])
    a = prm.a_wave

    FL0 = rho_face * a*a * qLu
    FR0 = rho_face * a*a * qRu
    FL1 = qLp / np.maximum(rho_face, 1.0e-9)
    FR1 = qRp / np.maximum(rho_face, 1.0e-9)

    smax = 1.15 * (a + np.maximum(np.abs(qLu), np.abs(qRu)))
    flux_p = 0.5*(FL0 + FR0) - 0.5*smax*(qRp - qLp)
    flux_u = 0.5*(FL1 + FR1) - 0.5*smax*(qRu - qLu)

    p_new = p - (dt/prm.dx)*(flux_p[1:] - flux_p[:-1])
    U_new = U - (dt/prm.dx)*(flux_u[1:] - flux_u[:-1])

    # friction
    S_fric = -(4.0*tauw / np.maximum(rho*D_eff, 1.0e-9)) * np.sign(U)
    U_new += dt * S_fric

    # unresolved transient-friction surrogate: damp toward the quasi-steady backbone
    u_ref = quasi_steady_velocity_target(t, prm)
    U_new += dt * (-prm.unsteady_damping_beta * (U_new - u_ref))

    # artificial viscosity
    nu_num = prm.c_art * prm.a_wave * prm.dx
    eps = min(nu_num * dt / (prm.dx*prm.dx), 0.20)
    p_new = apply_artificial_viscosity(p_new, eps)
    U_new = apply_artificial_viscosity(U_new, eps)

    # outlet sponge
    sigma = np.zeros_like(p_new)
    sigma[-prm.sponge_n:] = np.linspace(0.0, 1.0, prm.sponge_n)
    p_new = (1.0 - 0.12*sigma)*p_new + (0.12*sigma)*prm.p_out_dyn
    U_new = (1.0 - 0.18*sigma)*U_new + (0.18*sigma)*u_ref

    return p_new, U_new, Ucrit, tauw, gamma_w, tau_yw, mu_w, rp, tau_y0


# ----------------------------
# Coarse-phase and segregation transport
# ----------------------------

def scalar_upwind_step(phi: np.ndarray, U: np.ndarray, dt: float, t: float, prm: ModelParams15D):
    rho = mixture_density(phi, prm)
    _, u_ext, phi_ext, _ = extend_with_ghosts(np.zeros_like(phi), U, phi, rho, t, prm)
    u_face = 0.5*(u_ext[:-1] + u_ext[1:])
    phi_up = np.where(u_face >= 0.0, phi_ext[:-1], phi_ext[1:])
    F = u_face * phi_up
    phi_new = phi - (dt/prm.dx)*(F[1:] - F[:-1])
    return clip_phi(phi_new, prm)


def laplacian_1d(q: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(q)
    out[1:-1] = (q[:-2] - 2.0*q[1:-1] + q[2:]) / (dx*dx)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def segregation_step(beta: np.ndarray, U: np.ndarray, beta_eq_loc: np.ndarray,
                     a_perp: np.ndarray, dt: float, prm: ModelParams15D):
    beta_adv = scalar_upwind_step(beta, U, dt, 0.0, prm)  # same scheme, no special boundary forcing
    relax = -(beta_adv - beta_eq_loc) / max(prm.tau_m, 1.0e-6)
    diff = prm.D_beta * laplacian_1d(beta_adv, prm.dx)
    src = prm.kappa_a * a_perp
    beta_new = beta_adv + dt*(relax + diff + src)
    return np.clip(beta_new, -prm.beta_max, prm.beta_max)


def coarse_fraction_step(phi_bar: np.ndarray, U: np.ndarray, beta: np.ndarray, delta: np.ndarray,
                         Ucrit: np.ndarray, dt: float, t: float, prm: ModelParams15D):
    phi_adv = scalar_upwind_step(phi_bar, U, dt, t, prm)
    _, phiw = wall_center_concentrations(phi_adv, beta, prm)
    dep_drive = np.maximum(Ucrit - np.abs(U), 0.0) / np.maximum(Ucrit, 1.0e-9)
    ero_drive = np.maximum(np.abs(U) - prm.U_ero, 0.0)
    phi_new = (
        phi_adv
        - dt * prm.c_dep * dep_drive * np.maximum(phiw - prm.phi_c_eq_wall, 0.0)
        + dt * prm.c_ero * ero_drive * delta
    )
    return clip_phi(phi_new, prm)


def deposit_step(delta: np.ndarray, U: np.ndarray, phi_bar: np.ndarray, beta: np.ndarray,
                 Ucrit: np.ndarray, dt: float, prm: ModelParams15D):
    _, phiw = wall_center_concentrations(phi_bar, beta, prm)
    dep_drive = np.maximum(Ucrit - np.abs(U), 0.0) / np.maximum(Ucrit, 1.0e-9)
    ero_drive = np.maximum(np.abs(U) - prm.U_ero, 0.0)
    delta_dot = prm.k_dep * dep_drive * np.maximum(phiw - prm.phi_c_eq_wall, 0.0) - prm.k_ero * ero_drive * delta
    delta_new = np.clip(delta + dt*delta_dot, 0.0, prm.delta_max)
    return delta_new, delta_dot


# ----------------------------
# Structural modal surrogate
# ----------------------------

def project_modal(field: np.ndarray, modes: np.ndarray, dx: float) -> np.ndarray:
    return np.array([np.sum(field * modes[k]) * dx for k in range(modes.shape[0])])


def transverse_acceleration(qddot: np.ndarray, phi_modes: np.ndarray) -> np.ndarray:
    return np.sum(qddot[:, None] * phi_modes, axis=0)


def structure_step(q: np.ndarray, qdot: np.ndarray, w: np.ndarray, wdot: np.ndarray,
                   p: np.ndarray, U: np.ndarray, rho_m: np.ndarray, prm: ModelParams15D, dt: float):
    # reduced-order placeholder consistent with the 1.5D coupling discussion
    p_ref = np.mean(p)
    rho_ref = np.mean(rho_m)
    U_ref = np.mean(U)

    p_proj = project_modal(p - p_ref, prm.phi_modes, prm.dx)
    rho_proj = project_modal(rho_m - rho_ref, prm.phi_modes, prm.dx)
    U_proj = project_modal(U - U_ref, prm.phi_modes, prm.dx)

    # current-induced shedding frequency per mode (simplified)
    omega_s = 2.0*np.pi*prm.St*np.maximum(prm.Uo_ext, 1.0e-6) / prm.D * np.ones(prm.M)

    # internal forcing: reduced surrogate for density / pressure / velocity-wave modulation
    F_int = prm.c_rho_modal * rho_proj + prm.c_p_modal * p_proj + prm.c_u_modal * U_proj

    # wake oscillator
    wddot = -prm.wake_eps*omega_s*(w*w - 1.0)*wdot - omega_s**2 * w + prm.wake_A * 0.5 * qdot
    wdot_new = wdot + dt*wddot
    w_new = w + dt*wdot_new

    # VIV forcing from wake variable
    F_viv = prm.viv_force_gain * w_new

    qddot = F_viv + F_int - 2.0*prm.mode_zeta*prm.mode_omega*qdot - (prm.mode_omega**2)*q
    qdot_new = qdot + dt*qddot
    q_new = q + dt*qdot_new

    return q_new, qdot_new, w_new, wdot_new, qddot


# ----------------------------
# Initialization
# ----------------------------

def initialize_steady_state(prm: ModelParams15D):
    U = np.full(prm.N, prm.u_ss)
    phi_c_bar = np.full(prm.N, prm.phi_c_inlet)
    beta = np.zeros(prm.N)
    delta = np.zeros(prm.N)

    # clogging nucleus as coarse-particle enrichment
    seed = slice(40, 46)
    phi_c_bar[seed] += 0.05
    phi_c_bar = clip_phi(phi_c_bar, prm)

    rho = mixture_density(phi_c_bar, prm)
    D_eff = effective_diameter(delta, prm)
    _, phiw = wall_center_concentrations(phi_c_bar, beta, prm)
    tauw, _, _ = wall_shear_hb(U, D_eff, phiw, prm)
    dpdx = 4.0*tauw / np.maximum(D_eff, 1.0e-9)

    p_dyn = np.zeros(prm.N)
    p_dyn[-1] = prm.p_out_dyn
    for i in range(prm.N - 2, -1, -1):
        p_dyn[i] = p_dyn[i + 1] + dpdx[i + 1]*prm.dx

    prm.p_in0 = float(p_dyn[0])
    q = np.zeros(prm.M)
    qdot = np.zeros(prm.M)
    w = np.zeros(prm.M)
    wdot = np.zeros(prm.M)
    return State15D(p=p_dyn, U=U, phi_c_bar=phi_c_bar, beta=beta, delta=delta,
                    q=q, qdot=qdot, w=w, wdot=wdot)


# ----------------------------
# Diagnostics
# ----------------------------

def section_residual_pressure(p: np.ndarray, U: np.ndarray, phi_bar: np.ndarray, beta: np.ndarray,
                              delta: np.ndarray, prm: ModelParams15D):
    i0, i1 = prm.i_mon0, prm.i_mon1
    rho = mixture_density(phi_bar, prm)
    Ucrit, tauw, gamma_w, tau_yw, mu_w, rp, tau_y0 = critical_transport_velocity_segment(
        phi_bar, beta, delta, U, rho, prm
    )
    p_drop_meas = p[i0] - p[i1]
    D_eff = effective_diameter(delta, prm)
    dp_fric_model = np.sum((4.0*tauw[i0:i1+1] / np.maximum(D_eff[i0:i1+1], 1.0e-9)) * prm.dx)
    r_dp = p_drop_meas - dp_fric_model
    return r_dp, p_drop_meas, dp_fric_model, Ucrit


# ----------------------------
# Stable time step
# ----------------------------

def stable_dt(state: State15D, prm: ModelParams15D) -> float:
    wave = prm.a_wave + np.max(np.abs(state.U))
    dt_wave = prm.cfl * prm.dx / max(wave, 1.0e-9)
    dt_adv = prm.cfl * prm.dx / max(np.max(np.abs(state.U)), 1.0e-9)
    dt_beta = 0.5 * prm.tau_m
    dt_struct = 0.2 / max(np.max(prm.mode_omega), 1.0e-9)

    # deposit-source limiter
    rho = mixture_density(state.phi_c_bar, prm)
    Ucrit, *_ = critical_transport_velocity_segment(state.phi_c_bar, state.beta, state.delta, state.U, rho, prm)
    _, phiw = wall_center_concentrations(state.phi_c_bar, state.beta, prm)
    dep_drive = np.maximum(Ucrit - np.abs(state.U), 0.0) / np.maximum(Ucrit, 1.0e-9)
    delta_dot_est = prm.k_dep * dep_drive * np.maximum(phiw - prm.phi_c_eq_wall, 0.0)
    if np.any(delta_dot_est > 1.0e-12):
        dt_dep = 0.2*np.min((prm.delta_max - state.delta[delta_dot_est > 1.0e-12]) / delta_dot_est[delta_dot_est > 1.0e-12])
        dt_dep = max(dt_dep, 1.0e-5)
    else:
        dt_dep = prm.dt_max

    return min(prm.dt_max, dt_wave, dt_adv, dt_beta, dt_struct, dt_dep)


# ----------------------------
# Main simulator
# ----------------------------

def run(prm: ModelParams15D):
    state = initialize_steady_state(prm)

    t = 0.0
    r_prev, _, _, _ = section_residual_pressure(state.p, state.U, state.phi_c_bar, state.beta, state.delta, prm)
    dres_filt = 0.0
    alpha_filter = 0.95

    hist = {k: [] for k in [
        "t", "U_in", "U_mid", "U_out",
        "resid_dp", "dresid_dt",
        "seed_phi", "seed_beta", "seed_delta",
        "max_q", "max_aperp", "min_ucrit_ratio"
    ]}

    while t < prm.t_end:
        dt = stable_dt(state, prm)
        if t + dt > prm.t_end:
            dt = prm.t_end - t

        # fluid update using old beta / delta
        p, U, Ucrit_old, tauw, gamma_w, tau_yw, mu_w, rp, tau_y0 = fluid_step(
            state.p, state.U, state.phi_c_bar, state.beta, state.delta, dt, t, prm
        )

        # structural update driven by current fluid state
        rho_m = mixture_density(state.phi_c_bar, prm)
        q, qdot, w, wdot, qddot = structure_step(state.q, state.qdot, state.w, state.wdot,
                                                 state.p, state.U, rho_m, prm, dt)
        a_perp = transverse_acceleration(qddot, prm.phi_modes)

        # segregation update
        beta_eq_loc = beta_equilibrium(state.phi_c_bar, state.U, gamma_w, tau_yw, a_perp, prm)
        beta = segregation_step(state.beta, state.U, beta_eq_loc, a_perp, dt, prm)

        # coarse-phase transport with updated beta
        rho_m_new = mixture_density(state.phi_c_bar, prm)
        Ucrit_seg, *_ = critical_transport_velocity_segment(state.phi_c_bar, beta, state.delta, U, rho_m_new, prm)
        phi_c_bar = coarse_fraction_step(state.phi_c_bar, U, beta, state.delta, Ucrit_seg, dt, t, prm)

        # deposit update from updated transport state
        delta, delta_dot = deposit_step(state.delta, U, phi_c_bar, beta, Ucrit_seg, dt, prm)

        # positivity / clipping
        p = np.nan_to_num(p, nan=0.0, posinf=1.0e8, neginf=-1.0e8)
        U = np.clip(np.nan_to_num(U), -5.0, 5.0)
        phi_c_bar = clip_phi(np.nan_to_num(phi_c_bar), prm)
        beta = np.clip(np.nan_to_num(beta), -prm.beta_max, prm.beta_max)
        delta = np.clip(np.nan_to_num(delta), 0.0, prm.delta_max)

        state = State15D(p=p, U=U, phi_c_bar=phi_c_bar, beta=beta, delta=delta,
                         q=q, qdot=qdot, w=w, wdot=wdot)

        t += dt

        r_dp, _, _, Ucrit_mon = section_residual_pressure(state.p, state.U, state.phi_c_bar, state.beta, state.delta, prm)
        dres = (r_dp - r_prev) / max(dt, 1.0e-9)
        dres_filt = alpha_filter*dres_filt + (1.0 - alpha_filter)*dres
        r_prev = r_dp

        hist["t"].append(t)
        hist["U_in"].append(state.U[0])
        hist["U_mid"].append(state.U[prm.N // 2])
        hist["U_out"].append(state.U[-1])
        hist["resid_dp"].append(r_dp)
        hist["dresid_dt"].append(dres_filt)
        hist["seed_phi"].append(state.phi_c_bar[42])
        hist["seed_beta"].append(state.beta[42])
        hist["seed_delta"].append(state.delta[42])
        hist["max_q"].append(np.max(np.abs(state.q)))
        hist["max_aperp"].append(np.max(np.abs(a_perp)))
        hist["min_ucrit_ratio"].append(np.min(np.abs(state.U) / np.maximum(Ucrit_mon, 1.0e-9)))

    hist = {k: np.asarray(v) for k, v in hist.items()}
    hist["x"] = prm.x.copy()
    hist["final_p"] = state.p.copy()
    hist["final_U"] = state.U.copy()
    hist["final_phi_c_bar"] = state.phi_c_bar.copy()
    hist["final_beta"] = state.beta.copy()
    hist["final_delta"] = state.delta.copy()
    hist["final_q"] = state.q.copy()
    hist["final_qdot"] = state.qdot.copy()
    return hist


if __name__ == "__main__":
    prm = ModelParams15D()
    out = run(prm)

    print(f"steady inlet dynamic pressure: {prm.p_in0:.3e} Pa")
    print(f"max |d(residual Δp)/dt|: {np.max(np.abs(out['dresid_dt'])):.3e} Pa/s")
    print(f"final max deposit thickness: {np.max(out['final_delta']):.4e} m")
    print(f"final min velocity: {np.min(out['final_U']):.4f} m/s")
    print(f"final max velocity: {np.max(out['final_U']):.4f} m/s")
    print(f"final max |beta|: {np.max(np.abs(out['final_beta'])):.4e}")
    print(f"final max |q|: {np.max(np.abs(out['final_q'])):.4e} m")

    try:
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(3, 2, figsize=(12, 10))

        axs[0, 0].plot(out['t'], out['U_in'], label='U @ inlet')
        axs[0, 0].plot(out['t'], out['U_mid'], label='U @ mid')
        axs[0, 0].plot(out['t'], out['U_out'], label='U @ outlet')
        axs[0, 0].set_title('Velocity response')
        axs[0, 0].legend()

        axs[0, 1].plot(out['t'], out['resid_dp'], label='residual Δp')
        axs[0, 1].plot(out['t'], out['dresid_dt'], label='d(residual Δp)/dt')
        axs[0, 1].set_title('Residual pressure diagnostics')
        axs[0, 1].legend()

        axs[1, 0].plot(out['t'], out['seed_phi'], label='coarse φ̄ @ seed')
        axs[1, 0].plot(out['t'], out['seed_beta'], label='β @ seed')
        axs[1, 0].plot(out['t'], out['seed_delta'], label='δ @ seed')
        axs[1, 0].set_title('Segregation / deposit growth')
        axs[1, 0].legend()

        axs[1, 1].plot(out['t'], out['max_q'], label='max |q|')
        axs[1, 1].plot(out['t'], out['max_aperp'], label='max |a_perp|')
        axs[1, 1].set_title('Structural modal response')
        axs[1, 1].legend()

        axs[2, 0].plot(out['x'], out['final_phi_c_bar'], label='final coarse φ̄')
        axs[2, 0].plot(out['x'], out['final_beta'], label='final β')
        axs[2, 0].plot(out['x'], out['final_delta'], label='final δ')
        axs[2, 0].set_title('Final 1.5D fields')
        axs[2, 0].legend()

        axs[2, 1].plot(out['t'], out['min_ucrit_ratio'], label='min |U|/Ucrit_seg')
        axs[2, 1].axhline(1.0, color='k', ls='--', lw=1.0)
        axs[2, 1].set_title('Transport margin')
        axs[2, 1].legend()

        for ax in axs.flat:
            ax.set_xlabel('t [s]' if ax in axs[:, 0].tolist() + axs[:, 1].tolist() else '')
            ax.grid(alpha=0.25)

        plt.tight_layout()
        out_png = Path(__file__).with_suffix('.png')
        plt.savefig(out_png, dpi=160, bbox_inches='tight')
        print(f"saved plot: {out_png}")
    except Exception as exc:
        print(f"Plotting skipped: {exc}")
