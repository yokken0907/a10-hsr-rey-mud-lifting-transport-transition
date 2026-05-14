
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class ModelParams:
    # ------------------------
    # Geometry / mesh
    # ------------------------
    L: float = 5500.0
    N: int = 550
    D: float = 0.25  # nominal inner diameter [m]

    # ------------------------
    # Slurry properties
    # ------------------------
    rho_f: float = 1030.0         # carrier fluid density [kg/m3]
    rho_s: float = 2700.0         # representative solid density [kg/m3]
    d_rep: float = 80e-6          # representative coarse particle size [m]
    T: float = 2.0                # deep-sea envelope [degC]

    # ------------------------
    # Acoustic / transient
    # ------------------------
    a_wave: float = 700.0         # effective wave speed [m/s]
    unsteady_damping_beta: float = 2.0  # Rayleigh-like transient damping [1/s]

    # ------------------------
    # Steady operating point
    # ------------------------
    C_inlet: float = 0.28
    u_ss: float = 1.60            # chosen steady velocity [m/s]

    # blackout schedule
    t_blackout: float = 5.0       # [s]
    tau_coastdown: float = 2.0    # [s], exponential decay of pump head
    p_out_dyn: float = 0.0        # outlet dynamic pressure [Pa]

    # ------------------------
    # Herschel-Bulkley envelope
    # ------------------------
    C_ref: float = 0.22
    tau_y_ref: float = 12.0       # Pa
    K_ref: float = 0.45           # Pa s^n
    n_ref: float = 0.55

    a_tau: float = 7.0            # concentration sensitivity
    a_K: float = 3.0
    a_n: float = 0.12
    b_tau_T: float = 0.06         # lower T => higher yield
    b_K_T: float = 0.03
    T_ref: float = 20.0

    # ------------------------
    # Critical-velocity closures
    # ------------------------
    Re_transition: float = 3000.0
    Y_crit: float = 0.20
    hinder_exp: float = 4.65
    slip_sf: float = 1.15
    tau_mob_fac: float = 1.8

    # ------------------------
    # Deposition / erosion
    # ------------------------
    C_eq: float = 0.05
    U_ero: float = 1.8
    k_dep: float = 5.0e-4         # tuned to make clog growth visible in 1D toy model
    k_ero: float = 1.0e-2
    c_dep: float = 0.20
    c_ero: float = 0.03
    delta_max_frac: float = 0.45

    # ------------------------
    # Numerics
    # ------------------------
    t_end: float = 40.0
    cfl: float = 0.40
    dt_max: float = 0.01
    c_art: float = 0.05           # artificial viscosity coefficient
    out_absorb: float = 1.0       # outlet sponge strength
    in_impedance_fac: float = 1.2 # inlet impedance matching
    outlet_relax: float = 0.18
    sponge_n: int = 60
    check_valve: bool = True

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

        # determined from steady initialization
        self.p_in0: float | None = None


# ----------------------------
# Material closures
# ----------------------------

def mixture_density(C: np.ndarray, prm: ModelParams) -> np.ndarray:
    return prm.rho_f*(1.0 - C) + prm.rho_s*C


def hb_params(C: np.ndarray, prm: ModelParams):
    tau_y = prm.tau_y_ref * np.exp(
        prm.a_tau*(C - prm.C_ref) + prm.b_tau_T*(prm.T_ref - prm.T)
    )
    K = prm.K_ref * np.exp(
        prm.a_K*(C - prm.C_ref) + prm.b_K_T*(prm.T_ref - prm.T)
    )
    n = np.clip(prm.n_ref - prm.a_n*(C - prm.C_ref), 0.20, 0.95)
    return tau_y, K, n


def effective_diameter(delta: np.ndarray, prm: ModelParams) -> np.ndarray:
    return np.maximum(prm.D - 2.0*delta, 1.0e-4)


def effective_area_ratio(delta: np.ndarray, prm: ModelParams) -> np.ndarray:
    return np.clip((effective_diameter(delta, prm) / prm.D)**2, 1.0e-4, 1.0)


def wall_shear_hb(u_loc: np.ndarray, D_eff: np.ndarray, C: np.ndarray, prm: ModelParams) -> np.ndarray:
    tau_y, K, n = hb_params(C, prm)
    gamma_w = 8.0*np.abs(u_loc) / np.maximum(D_eff, 1.0e-6)
    return tau_y + K*np.power(gamma_w + 1.0e-9, n)


def apparent_viscosity_hb(u_loc: np.ndarray, D_eff: np.ndarray, C: np.ndarray, prm: ModelParams) -> np.ndarray:
    tau_y, K, n = hb_params(C, prm)
    gamma_w = 8.0*np.abs(u_loc) / np.maximum(D_eff, 1.0e-6)
    mu = tau_y/(gamma_w + 1.0e-9) + K*np.power(gamma_w + 1.0e-9, n - 1.0)
    return np.clip(mu, 1.0e-4, 1.0e6)


def u_transition(C: np.ndarray, delta: np.ndarray, prm: ModelParams) -> np.ndarray:
    """
    Crude HB transition-based velocity:
    solve Re_HB(U) ~ Re_transition by fixed-point iteration.
    """
    D_eff = effective_diameter(delta, prm)
    rho = mixture_density(C, prm)
    U = np.full_like(C, 0.30)
    for _ in range(6):
        mu = apparent_viscosity_hb(U, D_eff, C, prm)
        U = prm.Re_transition * mu / np.maximum(rho*D_eff, 1.0e-9)
    return np.clip(U, 1.0e-3, 20.0)


def hindered_settling_velocity(C: np.ndarray, prm: ModelParams) -> np.ndarray:
    """
    Placeholder HB-matrix hindered settling closure.
    Replace with pressurized settling-column calibration when data exist.
    """
    tau_y, K, n = hb_params(C, prm)
    gamma_settle = 1.0
    mu_settle = tau_y/gamma_settle + K*(gamma_settle**(n - 1.0))

    Y = tau_y / np.maximum((prm.rho_s - prm.rho_f)*9.81*prm.d_rep, 1.0e-9)
    trap = np.clip(1.0 - Y/prm.Y_crit, 0.0, 1.0)

    w0 = (prm.rho_s - prm.rho_f)*9.81*(prm.d_rep**2) / np.maximum(18.0*mu_settle, 1.0e-9)
    wh = w0 * np.power(np.clip(1.0 - C, 1.0e-6, 1.0), prm.hinder_exp) * trap
    return np.clip(wh, 0.0, 10.0)


def u_mobilization(C: np.ndarray, delta: np.ndarray, prm: ModelParams) -> np.ndarray:
    D_eff = effective_diameter(delta, prm)
    tau_y, K, n = hb_params(C, prm)
    tau_mob = prm.tau_mob_fac * tau_y
    surplus = np.maximum(tau_mob - tau_y, 0.0)
    gamma_req = np.power(surplus / np.maximum(K, 1.0e-12), 1.0/np.maximum(n, 1.0e-6))
    U = 0.125 * D_eff * gamma_req
    return np.clip(U, 0.0, 20.0)


def critical_transport_velocity(C: np.ndarray, delta: np.ndarray, prm: ModelParams) -> np.ndarray:
    U_tr = u_transition(C, delta, prm)
    U_sl = prm.slip_sf * hindered_settling_velocity(C, prm)
    U_mb = u_mobilization(C, delta, prm)
    return np.maximum.reduce([U_tr, U_sl, U_mb])


# ----------------------------
# Stable steady initialization
# ----------------------------

def initialize_steady_state(prm: ModelParams):
    """
    Build a pressure profile that is exactly consistent with:
    - chosen steady velocity u_ss
    - local HB friction
    - outlet dynamic pressure reference

    This removes the startup shock that v1 produced.
    """
    u = np.full(prm.N, prm.u_ss)
    C = np.full(prm.N, prm.C_inlet)
    delta = np.zeros(prm.N)

    # local concentration seed: "clogging nucleus"
    seed = slice(40, 46)
    C[seed] += 0.08
    C = np.clip(C, 0.0, 0.65)

    D_eff = effective_diameter(delta, prm)
    alpha = effective_area_ratio(delta, prm)
    u_loc = u / alpha
    tau_w = wall_shear_hb(u_loc, D_eff, C, prm)
    dpdx = 4.0*tau_w / np.maximum(D_eff, 1.0e-9)   # dynamic friction gradient [Pa/m]

    p_dyn = np.zeros(prm.N)
    p_dyn[-1] = prm.p_out_dyn
    for i in range(prm.N - 2, -1, -1):
        p_dyn[i] = p_dyn[i + 1] + dpdx[i + 1]*prm.dx

    prm.p_in0 = float(p_dyn[0])
    return p_dyn, u, C, delta


def pump_pressure_target(t: float, prm: ModelParams) -> float:
    if t < prm.t_blackout:
        return prm.p_in0
    return prm.p_in0 * np.exp(-(t - prm.t_blackout)/prm.tau_coastdown)


def quasi_steady_velocity_target(t: float, prm: ModelParams) -> float:
    return prm.u_ss * pump_pressure_target(t, prm) / max(prm.p_in0, 1.0e-9)


# ----------------------------
# Characteristic boundary conditions
# ----------------------------

def inlet_bc(pL: float, uL: float, rhoL: float, t: float, prm: ModelParams):
    """
    Impedance-matched inlet boundary with a pump-pressure target.
    A check-valve option suppresses nonphysical backflow after blackout.
    """
    Z = max(rhoL*prm.a_wave, 1.0e-6)
    w_minus = uL - pL/Z

    p_tgt = pump_pressure_target(t, prm)
    u_ref = quasi_steady_velocity_target(t, prm)
    R = prm.in_impedance_fac * Z

    # p_b = p_tgt - R*(u_b - u_ref)
    # w_minus = u_b - p_b/Z
    u_b = (w_minus + p_tgt/Z + (R/Z)*u_ref) / (1.0 + R/Z)
    p_b = p_tgt - R*(u_b - u_ref)

    if prm.check_valve and (t >= prm.t_blackout) and (u_b < 0.0):
        u_b = 0.0
        p_b = Z*(u_b - w_minus)

    return p_b, u_b


def outlet_bc(pR: float, uR: float, rhoR: float, t: float, prm: ModelParams):
    """
    Approximate non-reflecting outlet:
    outgoing characteristic leaves freely,
    incoming characteristic is tied to a decaying quasi-steady target.
    """
    Z = max(rhoR*prm.a_wave, 1.0e-6)
    w_plus = uR + pR/Z
    u_ref = quasi_steady_velocity_target(t, prm)
    w_minus_bc = u_ref - prm.p_out_dyn/Z

    u_b = 0.5*(w_plus + w_minus_bc)
    p_b = 0.5*Z*(w_plus - w_minus_bc)

    # light impedance damping at outlet
    p_b *= (1.0 - 0.35*prm.out_absorb)
    u_b = (1.0 - prm.outlet_relax)*u_b + prm.outlet_relax*u_ref
    return p_b, u_b


def extend_with_ghosts(p: np.ndarray, u: np.ndarray, C: np.ndarray, rho: np.ndarray, t: float, prm: ModelParams):
    N = prm.N
    p_ext = np.empty(N + 2)
    u_ext = np.empty(N + 2)
    C_ext = np.empty(N + 2)
    rho_ext = np.empty(N + 2)

    p_ext[1:-1] = p
    u_ext[1:-1] = u
    C_ext[1:-1] = C
    rho_ext[1:-1] = rho

    p_ext[0], u_ext[0] = inlet_bc(p[0], u[0], rho[0], t, prm)
    C_ext[0] = prm.C_inlet
    rho_ext[0] = rho[0]

    p_ext[-1], u_ext[-1] = outlet_bc(p[-1], u[-1], rho[-1], t, prm)
    C_ext[-1] = C[-1]
    rho_ext[-1] = rho[-1]
    return p_ext, u_ext, C_ext, rho_ext


# ----------------------------
# Numerical helpers
# ----------------------------

def apply_artificial_viscosity(q: np.ndarray, eps: float) -> np.ndarray:
    qn = q.copy()
    qn[1:-1] += eps * (q[:-2] - 2.0*q[1:-1] + q[2:])
    return qn


# ----------------------------
# Time-marching kernels
# ----------------------------

def acoustic_momentum_step(p: np.ndarray, u: np.ndarray, C: np.ndarray, delta: np.ndarray, dt: float, t: float, prm: ModelParams):
    """
    Acoustic block:
      p_t + rho a^2 u_x = 0
      u_t + (1/rho) p_x = friction + unsteady damping

    Discretization:
    - finite-volume, first-order Rusanov for the acoustic part
    - explicit HB friction
    - explicit unsteady damping surrogate
    - KO/artificial viscosity + outlet sponge
    """
    rho = mixture_density(C, prm)
    D_eff = effective_diameter(delta, prm)
    alpha = effective_area_ratio(delta, prm)
    u_loc = u / alpha

    p_ext, u_ext, _, rho_ext = extend_with_ghosts(p, u, C, rho, t, prm)

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

    # slightly stronger than plain Rusanov to damp standing-wave growth
    smax = 1.15 * (a + np.maximum(np.abs(qLu), np.abs(qRu)))
    flux_p = 0.5*(FL0 + FR0) - 0.5*smax*(qRp - qLp)
    flux_u = 0.5*(FL1 + FR1) - 0.5*smax*(qRu - qLu)

    p_new = p - (dt/prm.dx)*(flux_p[1:] - flux_p[:-1])
    u_new = u - (dt/prm.dx)*(flux_u[1:] - flux_u[:-1])

    # HB friction
    tau_w = wall_shear_hb(u_loc, D_eff, C, prm)
    S_fric = -(4.0*tau_w / np.maximum(rho*D_eff, 1.0e-9)) * np.sign(u_loc)
    u_new += dt * S_fric

    # unresolved transient-friction surrogate: damp toward the quasi-steady backbone
    u_ref = quasi_steady_velocity_target(t, prm)
    u_new += dt * (-prm.unsteady_damping_beta * (u_new - u_ref))

    # artificial viscosity for p and u
    nu_num = prm.c_art * prm.a_wave * prm.dx
    eps = min(nu_num * dt / (prm.dx*prm.dx), 0.20)
    p_new = apply_artificial_viscosity(p_new, eps)
    u_new = apply_artificial_viscosity(u_new, eps)

    # outlet sponge to suppress repeated reflections
    sigma = np.zeros_like(p_new)
    sigma[-prm.sponge_n:] = np.linspace(0.0, 1.0, prm.sponge_n)
    p_new = (1.0 - 0.12*sigma)*p_new + (0.12*sigma)*prm.p_out_dyn
    u_new = (1.0 - 0.18*sigma)*u_new + (0.18*sigma)*u_ref

    return p_new, u_new


def scalar_deposit_step(C: np.ndarray, u: np.ndarray, delta: np.ndarray, dt: float, t: float, prm: ModelParams):
    """
    Suspended concentration transport + deposition / erosion.
    """
    rho = mixture_density(C, prm)
    _, u_ext, C_ext, _ = extend_with_ghosts(np.zeros_like(C), u, C, rho, t, prm)

    # upwind scalar advection
    u_face = 0.5*(u_ext[:-1] + u_ext[1:])
    C_up = np.where(u_face >= 0.0, C_ext[:-1], C_ext[1:])
    F_C = u_face * C_up
    C_adv = C - (dt/prm.dx)*(F_C[1:] - F_C[:-1])
    C_adv = np.clip(C_adv, 0.0, 0.65)

    # deposit / erosion source
    alpha = effective_area_ratio(delta, prm)
    u_loc = np.abs(u) / alpha
    Ucrit = critical_transport_velocity(C_adv, delta, prm)

    dep_drive = np.maximum(Ucrit - u_loc, 0.0) / np.maximum(Ucrit, 1.0e-9)
    ero_drive = np.maximum(u_loc - prm.U_ero, 0.0)

    delta_dot = (
        prm.k_dep * dep_drive * np.maximum(C_adv - prm.C_eq, 0.0)
        - prm.k_ero * ero_drive * delta
    )

    delta_new = np.clip(delta + dt*delta_dot, 0.0, prm.delta_max)

    C_new = (
        C_adv
        - dt * prm.c_dep * dep_drive * np.maximum(C_adv - prm.C_eq, 0.0)
        + dt * prm.c_ero * ero_drive * delta
    )
    C_new = np.clip(C_new, 0.0, 0.65)
    return C_new, delta_new, Ucrit, delta_dot


def stable_dt(u: np.ndarray, delta: np.ndarray, C: np.ndarray, prm: ModelParams) -> float:
    alpha = effective_area_ratio(delta, prm)
    u_loc = np.abs(u) / alpha

    # acoustic CFL
    dt_wave = prm.cfl * prm.dx / max(prm.a_wave + np.max(u_loc), 1.0e-9)

    # scalar advection CFL
    dt_adv = prm.cfl * prm.dx / max(np.max(u_loc), 1.0e-9)

    # deposit-source timescale limiter
    Ucrit = critical_transport_velocity(C, delta, prm)
    dep_drive = np.maximum(Ucrit - u_loc, 0.0) / np.maximum(Ucrit, 1.0e-9)
    delta_dot_est = prm.k_dep * dep_drive * np.maximum(C - prm.C_eq, 0.0)

    if np.any(delta_dot_est > 1.0e-12):
        dt_dep = 0.20*np.min((prm.delta_max - delta[delta_dot_est > 1.0e-12]) / delta_dot_est[delta_dot_est > 1.0e-12])
        dt_dep = max(dt_dep, 1.0e-5)
    else:
        dt_dep = prm.dt_max

    return min(prm.dt_max, dt_wave, dt_adv, dt_dep)


# ----------------------------
# Diagnostics
# ----------------------------

def section_residual_pressure(p: np.ndarray, u: np.ndarray, C: np.ndarray, delta: np.ndarray, prm: ModelParams):
    """
    Dynamic residual section pressure:
      r_dp = measured dynamic section drop - modeled HB friction drop
    """
    i0, i1 = prm.i_mon0, prm.i_mon1
    D_eff = effective_diameter(delta, prm)
    alpha = effective_area_ratio(delta, prm)
    u_loc = np.abs(u[i0:i1+1]) / alpha[i0:i1+1]

    tau_w = wall_shear_hb(u_loc, D_eff[i0:i1+1], C[i0:i1+1], prm)
    p_drop_meas = p[i0] - p[i1]
    dp_fric_model = np.sum((4.0*tau_w / np.maximum(D_eff[i0:i1+1], 1.0e-9)) * prm.dx)

    r_dp = p_drop_meas - dp_fric_model
    return r_dp, p_drop_meas, dp_fric_model


# ----------------------------
# Main driver
# ----------------------------

def run(prm: ModelParams):
    p, u, C, delta = initialize_steady_state(prm)

    t = 0.0
    r_prev, _, _ = section_residual_pressure(p, u, C, delta, prm)
    dres_filt = 0.0
    alpha_filter = 0.95

    hist = {k: [] for k in [
        "t", "u_in", "u_mid", "u_out",
        "resid_dp", "dresid_dt",
        "seed_C", "seed_delta",
        "p_in", "p_out"
    ]}

    while t < prm.t_end:
        dt = stable_dt(u, delta, C, prm)
        if t + dt > prm.t_end:
            dt = prm.t_end - t

        p, u = acoustic_momentum_step(p, u, C, delta, dt, t, prm)
        C, delta, Ucrit, delta_dot = scalar_deposit_step(C, u, delta, dt, t, prm)

        # positivity / clipping
        p = np.nan_to_num(p, nan=0.0, posinf=1.0e8, neginf=-1.0e8)
        u = np.clip(np.nan_to_num(u), -5.0, 5.0)
        C = np.clip(np.nan_to_num(C), 0.0, 0.65)
        delta = np.clip(np.nan_to_num(delta), 0.0, prm.delta_max)

        t += dt

        r_dp, _, _ = section_residual_pressure(p, u, C, delta, prm)
        dres = (r_dp - r_prev) / max(dt, 1.0e-9)
        dres_filt = alpha_filter*dres_filt + (1.0 - alpha_filter)*dres
        r_prev = r_dp

        hist["t"].append(t)
        hist["u_in"].append(u[0])
        hist["u_mid"].append(u[prm.N // 2])
        hist["u_out"].append(u[-1])
        hist["resid_dp"].append(r_dp)
        hist["dresid_dt"].append(dres_filt)
        hist["seed_C"].append(C[42])
        hist["seed_delta"].append(delta[42])
        hist["p_in"].append(p[0])
        hist["p_out"].append(p[-1])

    hist = {k: np.array(v) for k, v in hist.items()}
    state = {
        "x": prm.x.copy(),
        "p": p.copy(),
        "u": u.copy(),
        "C": C.copy(),
        "delta": delta.copy(),
    }
    return hist, state


if __name__ == "__main__":
    prm = ModelParams()
    hist, state = run(prm)

    print("Simulation complete")
    print(f"steady inlet dynamic pressure: {prm.p_in0:.3e} Pa")
    print(f"max |d(residual Δp)/dt|: {np.max(np.abs(hist['dresid_dt'])):.3e} Pa/s")
    print(f"final max deposit thickness: {np.max(state['delta']):.4e} m")
    print(f"final min velocity: {np.min(state['u']):.4f} m/s")
    print(f"final max velocity: {np.max(state['u']):.4f} m/s")

    try:
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 2, figsize=(12, 8))

        axs[0, 0].plot(hist["t"], hist["u_in"], label="u @ inlet")
        axs[0, 0].plot(hist["t"], hist["u_mid"], label="u @ mid")
        axs[0, 0].plot(hist["t"], hist["u_out"], label="u @ outlet")
        axs[0, 0].axvline(prm.t_blackout, color="k", linestyle="--", linewidth=1)
        axs[0, 0].set_title("Velocity response")
        axs[0, 0].set_xlabel("t [s]")
        axs[0, 0].set_ylabel("u [m/s]")
        axs[0, 0].legend()

        axs[0, 1].plot(hist["t"], hist["resid_dp"], label="residual Δp")
        axs[0, 1].plot(hist["t"], hist["dresid_dt"], label="d(residual Δp)/dt")
        axs[0, 1].axvline(prm.t_blackout, color="k", linestyle="--", linewidth=1)
        axs[0, 1].set_title("Residual pressure diagnostics")
        axs[0, 1].set_xlabel("t [s]")
        axs[0, 1].legend()

        axs[1, 0].plot(hist["t"], hist["seed_C"], label="C @ seed")
        axs[1, 0].plot(hist["t"], hist["seed_delta"], label="δ @ seed")
        axs[1, 0].axvline(prm.t_blackout, color="k", linestyle="--", linewidth=1)
        axs[1, 0].set_title("Local clog growth")
        axs[1, 0].set_xlabel("t [s]")
        axs[1, 0].legend()

        axs[1, 1].plot(state["x"], state["u"], label="final u")
        axs[1, 1].plot(state["x"], state["C"], label="final C")
        axs[1, 1].plot(state["x"], state["delta"], label="final δ")
        axs[1, 1].set_title("Final profiles")
        axs[1, 1].set_xlabel("x [m]")
        axs[1, 1].legend()

        plt.tight_layout()
        plt.savefig("hsr_simulation_v2.png")
        print("結果を 'hsr_simulation_v2.png' として保存しました。")

    except ImportError:
        pass
