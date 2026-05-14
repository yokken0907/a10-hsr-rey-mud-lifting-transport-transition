from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass
class ModelParams:
    L: float = 5500.0
    N: int = 550
    D: float = 0.25
    rho_f: float = 1030.0
    rho_s: float = 2700.0
    a_wave: float = 1100.0
    T: float = 2.0
    d_rep: float = 80e-6

    C_ref: float = 0.22
    tau_y_ref: float = 12.0
    K_ref: float = 0.45
    n_ref: float = 0.55

    a_tau: float = 7.0
    a_K: float = 3.0
    a_n: float = 0.12
    b_tau_T: float = 0.06
    b_K_T: float = 0.03
    T_ref: float = 20.0

    Re_transition: float = 3000.0
    Y_crit: float = 0.20
    hinder_exp: float = 4.65
    slip_sf: float = 1.15
    tau_mob_fac: float = 1.8

    C_eq: float = 0.10
    U_ero: float = 1.8
    k_dep: float = 2.5e-4
    k_ero: float = 6.0e-2
    c_dep: float = 0.15
    c_ero: float = 0.08
    delta_max_frac: float = 0.45

    H_pump0: float = 6200.0
    t_blackout: float = 5.0
    tau_coastdown: float = 0.15
    p_out_dyn: float = 0.0
    C_inlet: float = 0.28

    t_end: float = 40.0
    cfl: float = 0.45
    dt_max: float = 0.01

    i_mon0: int = 0
    i_mon1: int = 80

    def __post_init__(self):
        self.dx = self.L / self.N
        self.x = np.linspace(0.5*self.dx, self.L - 0.5*self.dx, self.N)
        self.A0 = 0.25*np.pi*self.D**2
        self.delta_max = self.delta_max_frac * 0.5 * self.D

def mixture_density(C, prm): return prm.rho_f*(1.0 - C) + prm.rho_s*C

def hb_params(C, prm):
    tau_y = prm.tau_y_ref * np.exp(prm.a_tau*(C - prm.C_ref) + prm.b_tau_T*(prm.T_ref - prm.T))
    K = prm.K_ref * np.exp(prm.a_K*(C - prm.C_ref) + prm.b_K_T*(prm.T_ref - prm.T))
    n = np.clip(prm.n_ref - prm.a_n*(C - prm.C_ref), 0.20, 0.95)
    return tau_y, K, n

def effective_diameter(delta, prm): return np.maximum(prm.D - 2.0*delta, 1.0e-4)
def effective_area_ratio(delta, prm): return np.clip((effective_diameter(delta, prm) / prm.D)**2, 1.0e-4, 1.0)

def wall_shear_hb(u_loc, D_eff, C, prm):
    tau_y, K, n = hb_params(C, prm)
    gamma_w = 8.0*np.abs(u_loc) / np.maximum(D_eff, 1.0e-6)
    return tau_y + K*np.power(gamma_w + 1.0e-9, n)

def apparent_viscosity_hb(u_loc, D_eff, C, prm):
    tau_y, K, n = hb_params(C, prm)
    gamma_w = 8.0*np.abs(u_loc) / np.maximum(D_eff, 1.0e-6)
    mu = tau_y/(gamma_w + 1.0e-9) + K*np.power(gamma_w + 1.0e-9, n - 1.0)
    return np.clip(mu, 1.0e-4, 1.0e6)

def u_transition(C, delta, prm):
    D_eff = effective_diameter(delta, prm)
    rho = mixture_density(C, prm)
    U = np.full_like(C, 0.5)
    for _ in range(6):
        mu = apparent_viscosity_hb(U, D_eff, C, prm)
        U = prm.Re_transition * mu / np.maximum(rho*D_eff, 1.0e-9)
    return np.clip(U, 1.0e-3, 20.0)

def hindered_settling_velocity(C, prm):
    tau_y, K, n = hb_params(C, prm)
    gamma_settle = 1.0
    mu_settle = tau_y/gamma_settle + K*(gamma_settle**(n - 1.0))
    Y = tau_y / np.maximum((prm.rho_s - prm.rho_f)*9.81*prm.d_rep, 1.0e-9)
    trap = np.clip(1.0 - Y/prm.Y_crit, 0.0, 1.0)
    w0 = (prm.rho_s - prm.rho_f)*9.81*(prm.d_rep**2) / np.maximum(18.0*mu_settle, 1.0e-9)
    wh = w0 * np.power(np.clip(1.0 - C, 1.0e-6, 1.0), prm.hinder_exp) * trap
    return np.clip(wh, 0.0, 10.0)

def u_mobilization(C, delta, prm):
    D_eff = effective_diameter(delta, prm)
    tau_y, K, n = hb_params(C, prm)
    tau_mob = prm.tau_mob_fac * tau_y
    surplus = np.maximum(tau_mob - tau_y, 0.0)
    gamma_req = np.power(surplus / np.maximum(K, 1.0e-12), 1.0/np.maximum(n, 1.0e-6))
    return np.clip(0.125 * D_eff * gamma_req, 0.0, 20.0)

def critical_transport_velocity(C, delta, prm):
    U_tr = u_transition(C, delta, prm)
    U_sl = prm.slip_sf * hindered_settling_velocity(C, prm)
    U_mb = u_mobilization(C, delta, prm)
    return np.maximum.reduce([U_tr, U_sl, U_mb])

def pump_head(t, prm):
    if t < prm.t_blackout: return prm.H_pump0
    return prm.H_pump0 * np.exp(-(t - prm.t_blackout)/prm.tau_coastdown)

def inlet_dynamic_pressure(t, rho_ref, prm):
    return rho_ref * 9.81 * pump_head(t, prm)

def extend_with_ghosts(p, u, C, rho, prm, t):
    N = prm.N
    p_ext, u_ext, C_ext, rho_ext = np.empty(N + 2), np.empty(N + 2), np.empty(N + 2), np.empty(N + 2)
    p_ext[1:-1], u_ext[1:-1], C_ext[1:-1], rho_ext[1:-1] = p, u, C, rho
    p_ext[0] = inlet_dynamic_pressure(t, rho[0], prm)
    u_ext[0], C_ext[0], rho_ext[0] = u[0], prm.C_inlet, rho[0]
    p_ext[-1], u_ext[-1], C_ext[-1], rho_ext[-1] = prm.p_out_dyn, u[-1], C[-1], rho[-1]
    return p_ext, u_ext, C_ext, rho_ext

def acoustic_fvm_step(p, u, C, delta, dt, prm, t):
    rho = mixture_density(C, prm)
    D_eff = effective_diameter(delta, prm)
    alpha = effective_area_ratio(delta, prm)
    u_loc = u / alpha
    p_ext, u_ext, _, rho_ext = extend_with_ghosts(p, u, C, rho, prm, t)
    
    qL_p, qR_p = p_ext[:-1], p_ext[1:]
    qL_u, qR_u = u_ext[:-1], u_ext[1:]
    rho_f = 0.5*(rho_ext[:-1] + rho_ext[1:])
    a = prm.a_wave

    FL0, FR0 = rho_f * a*a * qL_u, rho_f * a*a * qR_u
    FL1, FR1 = qL_p / np.maximum(rho_f, 1.0e-9), qR_p / np.maximum(rho_f, 1.0e-9)
    
    flux_p = 0.5*(FL0 + FR0) - 0.5*a*(qR_p - qL_p)
    flux_u = 0.5*(FL1 + FR1) - 0.5*a*(qR_u - qL_u)

    p_new = p - (dt/prm.dx)*(flux_p[1:] - flux_p[:-1])
    u_new = u - (dt/prm.dx)*(flux_u[1:] - flux_u[:-1])

    u_ext2 = np.empty(prm.N + 2)
    u_ext2[1:-1] = u_new
    u_ext2[0], u_ext2[-1] = u_new[0], u_new[-1]
    du_dx = (u_ext2[2:] - u_ext2[:-2]) / (2.0*prm.dx)

    alpha_new = effective_area_ratio(delta, prm)
    u_loc_new = u_new / alpha_new
    tau_w = wall_shear_hb(u_loc_new, D_eff, C, prm)
    S_fric = -(4.0*tau_w / np.maximum(rho*D_eff, 1.0e-9)) * np.sign(u_loc_new)

    u_new += dt * (-u_new*du_dx + S_fric)
    return p_new, u_new

def scalar_transport_step(C, u, delta, dt, prm, t):
    rho = mixture_density(C, prm)
    _, u_ext, C_ext, _ = extend_with_ghosts(np.zeros_like(C), u, C, rho, prm, t)
    u_face = 0.5*(u_ext[:-1] + u_ext[1:])
    C_up = np.where(u_face >= 0.0, C_ext[:-1], C_ext[1:])
    F_C = u_face * C_up
    C_adv = C - (dt/prm.dx)*(F_C[1:] - F_C[:-1])

    alpha = effective_area_ratio(delta, prm)
    u_loc = np.abs(u) / alpha
    Ucrit = critical_transport_velocity(C_adv, delta, prm)

    dep_drive = np.maximum(Ucrit - u_loc, 0.0) / np.maximum(Ucrit, 1.0e-9)
    ero_drive = np.maximum(u_loc - prm.U_ero, 0.0)

    delta_dot = prm.k_dep * dep_drive * np.maximum(C_adv - prm.C_eq, 0.0) - prm.k_ero * ero_drive * delta
    delta_new = np.clip(delta + dt*delta_dot, 0.0, prm.delta_max)

    C_new = C_adv - dt * prm.c_dep * dep_drive * np.maximum(C_adv - prm.C_eq, 0.0) + dt * prm.c_ero * ero_drive * delta
    return np.clip(C_new, 0.0, 0.65), delta_new, Ucrit, delta_dot

def stable_dt(u, delta, C, prm):
    alpha = effective_area_ratio(delta, prm)
    u_loc = np.abs(u) / alpha
    dt_wave = prm.cfl * prm.dx / max(prm.a_wave + np.max(u_loc), 1.0e-9)
    dt_adv  = prm.cfl * prm.dx / max(np.max(u_loc),  1.0e-9)

    Ucrit = critical_transport_velocity(C, delta, prm)
    dep_drive = np.maximum(Ucrit - u_loc, 0.0) / np.maximum(Ucrit, 1.0e-9)
    delta_dot_est = prm.k_dep * dep_drive * np.maximum(C - prm.C_eq, 0.0)
    pos = delta_dot_est > 1.0e-12
    dt_dep = max(0.25*np.min((prm.delta_max - delta[pos]) / delta_dot_est[pos]), 1.0e-5) if np.any(pos) else prm.dt_max
    return min(prm.dt_max, dt_wave, dt_adv, dt_dep)

def section_residual_pressure(p, u, C, delta, prm):
    i0, i1 = prm.i_mon0, prm.i_mon1
    D_eff = effective_diameter(delta, prm)
    alpha = effective_area_ratio(delta, prm)
    p_drop_meas = p[i0] - p[i1]
    u_loc = np.abs(u[i0:i1+1]) / alpha[i0:i1+1]
    tau_w = wall_shear_hb(u_loc, D_eff[i0:i1+1], C[i0:i1+1], prm)
    dp_fric_model = np.sum((4.0*tau_w / np.maximum(D_eff[i0:i1+1], 1.0e-9)) * prm.dx)
    return p_drop_meas - dp_fric_model, p_drop_meas, dp_fric_model

def initialize_state(prm):
    u0 = np.full(prm.N, 1.8)
    C0 = np.full(prm.N, prm.C_inlet)
    delta0 = np.zeros(prm.N)
    
    seed_i = 40
    C0[seed_i:seed_i+5] += 0.07
    C0 = np.clip(C0, 0.0, 0.65)

    D_eff0 = effective_diameter(delta0, prm)
    tau_w0 = wall_shear_hb(u0, D_eff0, C0, prm)
    dpdx0 = 4.0*tau_w0 / np.maximum(D_eff0, 1.0e-9)

    p0 = np.zeros(prm.N)
    p0[-1] = prm.p_out_dyn
    for i in range(prm.N-2, -1, -1):
        p0[i] = p0[i+1] + dpdx0[i+1]*prm.dx
    return p0, u0, C0, delta0

def run(prm):
    p, u, C, delta = initialize_state(prm)
    t = 0.0
    hist = {"t": [], "resid_dp": [], "dresid_dt": [], "u_in": [], "u_mid": [], "u_out": [], "seed_C": [], "seed_delta": []}
    
    r_prev, _, _ = section_residual_pressure(p, u, C, delta, prm)
    dresid_filt = 0.0

    while t < prm.t_end:
        dt = stable_dt(u, delta, C, prm)
        if t + dt > prm.t_end: dt = prm.t_end - t

        p, u = acoustic_fvm_step(p, u, C, delta, dt, prm, t)
        C, delta, _, _ = scalar_transport_step(C, u, delta, dt, prm, t)

        p = np.nan_to_num(p, nan=0.0, posinf=1.0e9, neginf=-1.0e9)
        u = np.nan_to_num(u, nan=0.0, posinf=50.0, neginf=-50.0)
        
        t += dt
        r_dp, _, _ = section_residual_pressure(p, u, C, delta, prm)
        dresid = (r_dp - r_prev) / max(dt, 1.0e-9)
        dresid_filt = 0.85*dresid_filt + 0.15*dresid
        r_prev = r_dp

        hist["t"].append(t); hist["resid_dp"].append(r_dp); hist["dresid_dt"].append(dresid_filt)
        hist["u_in"].append(u[0]); hist["u_mid"].append(u[prm.N//2]); hist["u_out"].append(u[-1])
        hist["seed_C"].append(C[40]); hist["seed_delta"].append(delta[40])

    hist.update({"final_p": p, "final_u": u, "final_C": C, "final_delta": delta, "x": prm.x})
    return {k: np.array(v) for k, v in hist.items()}

if __name__ == "__main__":
    print("シミュレーションを開始します...（数十秒かかる場合があります）")
    prm = ModelParams()
    out = run(prm)

    print("計算完了！")
    print(f"最大 残差圧力上昇率 = {np.max(np.abs(out['dresid_dt'])):.3e} Pa/s")
    print(f"最終的な最大堆積厚 = {np.max(out['final_delta']):.4f} m")

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    axs[0, 0].plot(out["t"], out["u_in"], label="Inlet Velocity")
    axs[0, 0].plot(out["t"], out["u_mid"], label="Mid Velocity")
    axs[0, 0].plot(out["t"], out["u_out"], label="Outlet Velocity")
    axs[0, 0].set_title("Velocity after Blackout (t=5s)")
    axs[0, 0].legend()

    axs[0, 1].plot(out["t"], out["resid_dp"], label="Residual $\\Delta p$")
    axs[0, 1].plot(out["t"], out["dresid_dt"], label="d(Residual $\\Delta p$)/dt")
    axs[0, 1].set_title("Clogging Indicator (Residual Pressure)")
    axs[0, 1].legend()

    axs[1, 0].plot(out["t"], out["seed_C"], label="Concentration at Clog Seed")
    axs[1, 0].plot(out["t"], out["seed_delta"], label="Deposit Thickness at Seed")
    axs[1, 0].set_title("Local Clog Growth")
    axs[1, 0].legend()

    axs[1, 1].plot(out["x"], out["final_u"], label="Final Velocity")
    axs[1, 1].plot(out["x"], out["final_delta"], label="Final Deposit Thickness")
    axs[1, 1].set_title("Final Pipe Profile")
    axs[1, 1].legend()

    plt.tight_layout()
    # 罠対策：画面に表示するのではなく、画像ファイルとして保存する
    plt.savefig("hsr_simulation.png")
    print("結果を 'hsr_simulation.png' として保存しました。")
