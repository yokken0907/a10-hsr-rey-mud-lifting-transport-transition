#!/usr/bin/env python3
from __future__ import annotations

# -----------------------------------------------------------------------------
# Best-effort single-thread enforcement for BLAS / OpenMP backends.
# This is useful on CPU-constrained hosts and matches the closure-builder policy.
# -----------------------------------------------------------------------------
import os
for _env_name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ.setdefault(_env_name, "1")

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, RegularGridInterpolator


# =============================================================================
# Closure-table loader / interpolator
# =============================================================================

@dataclass
class ClosureTableModel:
    csv_path: Path
    clip_to_domain: bool = True

    input_cols: tuple[str, ...] = (
        "phi_tot", "phi_c_bulk_target", "T_degC", "P_Pa", "d_p_m", "U_bulk_target_mps",
    )
    output_cols: tuple[str, ...] = (
        "beta_eq", "tau_w_Pa", "r_p_m", "phi_center", "phi_wall", "G_Pa_per_m", "U_bulk_solved_mps",
    )

    axes: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    grid_shape: tuple[int, ...] = field(init=False, default=())
    is_tensor_grid: bool = field(init=False, default=False)
    linear_rgis: Dict[str, RegularGridInterpolator] = field(init=False, default_factory=dict)
    nearest_rgis: Dict[str, RegularGridInterpolator] = field(init=False, default_factory=dict)
    point_cloud: np.ndarray = field(init=False, default_factory=lambda: np.empty((0, 6), dtype=float))
    u_axis: np.ndarray = field(init=False, default_factory=lambda: np.empty(0, dtype=float))

    def __post_init__(self) -> None:
        self.csv_path = Path(self.csv_path)
        rows = self._read_rows(self.csv_path)
        if not rows:
            raise ValueError(f"No valid closure rows found in {self.csv_path}")
        self._build_interpolators(rows)

    @staticmethod
    def _safe_float(row: Mapping[str, str], key: str) -> float:
        val = row.get(key, "")
        if val is None or str(val).strip() == "":
            raise ValueError(f"Empty value for column {key}")
        return float(str(val).strip())

    def _read_rows(self, csv_path: Path) -> List[Dict[str, float]]:
        rows: List[Dict[str, float]] = []
        with csv_path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if str(row.get("status", "ok")).strip().lower() != "ok":
                    continue
                try:
                    item = {k: self._safe_float(row, k) for k in (*self.input_cols, *self.output_cols)}
                    rows.append(item)
                except Exception:
                    continue
        return rows

    def _build_interpolators(self, rows: List[Dict[str, float]]) -> None:
        points = np.asarray([[row[k] for k in self.input_cols] for row in rows], dtype=float)
        self.point_cloud = points
        
        # 各次元のユニークな値を抽出してソート
        self.axes = {k: np.sort(np.unique(points[:, i])) for i, k in enumerate(self.input_cols)}
        self.u_axis = self.axes["U_bulk_target_mps"].copy()
        self.grid_shape = tuple(len(self.axes[k]) for k in self.input_cols)
        
        # 強制的にテンソルグリッド（超高速モード）として扱う
        self.is_tensor_grid = True
        axis_index = {k: {float(v): i for i, v in enumerate(self.axes[k])} for k in self.input_cols}
        
        field_arrays: Dict[str, np.ndarray] = {
            out_key: np.full(self.grid_shape, np.nan, dtype=float) for out_key in self.output_cols
        }
        
        for row in rows:
            idx = tuple(axis_index[k][float(row[k])] for k in self.input_cols)
            for out_key in self.output_cols:
                field_arrays[out_key][idx] = float(row[out_key])
                
        # 穴（NaN）を平均値で埋める
        for out_key, arr in field_arrays.items():
            if np.isnan(arr).any():
                mean_val = np.nanmean(arr)
                arr[np.isnan(arr)] = mean_val
                
        # 超高速な RegularGridInterpolator のみを生成
        for out_key, arr in field_arrays.items():
            self.linear_rgis[out_key] = RegularGridInterpolator(
                tuple(self.axes[k] for k in self.input_cols),
                arr, method="linear", bounds_error=False, fill_value=np.nan
            )
            self.nearest_rgis[out_key] = RegularGridInterpolator(
                tuple(self.axes[k] for k in self.input_cols),
                arr, method="nearest", bounds_error=False, fill_value=np.nan
            )

    def _clip_points(self, pts: np.ndarray) -> np.ndarray:
        if not self.clip_to_domain:
            return pts
        pts = np.asarray(pts, dtype=float).copy()
        for i, key in enumerate(self.input_cols):
            ax = self.axes[key]
            pts[:, i] = np.clip(pts[:, i], ax[0], ax[-1])
        return pts

    def evaluate_fields(self, pts: np.ndarray, fields: Sequence[str]) -> Dict[str, np.ndarray]:
        pts_eval = self._clip_points(np.atleast_2d(pts))
        out: Dict[str, np.ndarray] = {}
        for key in fields:
            vals = np.asarray(self.linear_rgis[key](pts_eval), dtype=float)
            bad = ~np.isfinite(vals)
            if np.any(bad):
                vals[bad] = np.asarray(self.nearest_rgis[key](pts_eval[bad]), dtype=float)
            out[key] = vals
        return out

    def backbone_closures(
        self, phi_tot: np.ndarray, phi_c_bulk: np.ndarray, T_degC: np.ndarray,
        P_Pa: np.ndarray, d_p_m: np.ndarray, U_bulk: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        U_mag = np.abs(np.asarray(U_bulk, dtype=float))
        pts = np.column_stack([
            np.asarray(phi_tot, dtype=float), np.asarray(phi_c_bulk, dtype=float),
            np.asarray(T_degC, dtype=float), np.asarray(P_Pa, dtype=float),
            np.asarray(d_p_m, dtype=float), U_mag,
        ])
        return self.evaluate_fields(pts, self.output_cols)

# =============================================================================
# 1.5D model parameters / state
# =============================================================================

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
    rho_l: float = 1030.0
    rho_fstar: float = 1600.0
    rho_cstar: float = 2700.0
    d_rep: float = 80e-6
    T: float = 2.0
    P: float = 55e6

    # ------------------------
    # Background fine matrix loading
    # ------------------------
    phi_f_bg: float = 0.10
    phi_c_inlet: float = 0.18
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
    t_blackout: float = 1.0e9
    tau_coastdown: float = 2.0
    p_out_dyn: float = 0.0
    check_valve: bool = True

    # ------------------------
    # Herschel-Bulkley baseline (used as off-manifold correction model)
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
    tau_m: float = 3.0
    D_beta: float = 0.03
    kappa_a: float = 8.0e-4         # structural acceleration -> beta_eq shift
    beta_max: float = 0.12
    beta_eq_relax_clip: float = 0.08

    # ------------------------
    # Critical-velocity closures
    # ------------------------
    Re_transition: float = 3000.0
    Y_crit: float = 0.20
    hinder_exp: float = 4.65
    slip_sf: float = 1.15
    tau_mob_fac: float = 1.8
    ucrit_margin: float = 1.25

    # ------------------------
    # Truth-table / off-equilibrium feedback
    # ------------------------
    closure_table_path: str = "closure_table.csv"
    tauw_beta_feedback: float = 0.35
    rp_beta_feedback: float = 0.30

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
    T0: float = 2.5e5
    m_struct: float = 180.0
    m_added: float = 120.0
    c_rho_modal: float = 5.0e-4
    c_p_modal: float = 1.5e-7
    c_u_modal: float = 5.0e-3
    Uo_ext: float = 0.55
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

    # Runtime-populated objects
    closure_model: ClosureTableModel | None = field(init=False, default=None, repr=False)

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

        self.phi_modes = np.zeros((self.M, self.N))
        for k in range(self.M):
            n = k + 1
            self.phi_modes[k, :] = np.sin(n*np.pi*self.x / self.L)
        for k in range(self.M):
            norm = np.sqrt(np.sum(self.phi_modes[k]**2) * self.dx)
            self.phi_modes[k] /= max(norm, 1.0e-12)

        self.i_mon0 = int(np.clip(self.i_mon0, 0, self.N - 1))
        self.i_mon1 = int(np.clip(self.i_mon1, self.i_mon0 + 1, self.N - 1))

        self.p_in0: float | None = None
        closure_path = Path(self.closure_table_path)
        if not closure_path.is_absolute():
            closure_path = Path(__file__).resolve().parent / closure_path
        self.closure_model = ClosureTableModel(closure_path)


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


# =============================================================================
# Basic field utilities
# =============================================================================


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


def hb_matrix_params(phi_tot: np.ndarray, prm: ModelParams15D):
    tau_y = prm.tau_y_ref * np.exp(prm.a_tau*(phi_tot - prm.phi_ref) + prm.b_tau_T*(prm.T_ref - prm.T))
    K = prm.K_ref * np.exp(prm.a_K*(phi_tot - prm.phi_ref) + prm.b_K_T*(prm.T_ref - prm.T))
    n = np.clip(prm.n_ref - prm.a_n*(phi_tot - prm.phi_ref), 0.20, 0.95)
    return tau_y, K, n


def hb_wall_params(phiw: np.ndarray, prm: ModelParams15D):
    phi_tot_w = np.clip(prm.phi_f_bg + phiw, 0.0, prm.phi_max)
    tau_m, K_m, n_m = hb_matrix_params(phi_tot_w, prm)
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


# =============================================================================
# Truth-backed closures
# =============================================================================


def _truth_input_pack(
    phi_bar: np.ndarray,
    U: np.ndarray,
    prm: ModelParams15D,
) -> Dict[str, np.ndarray]:
    phi_bar = np.asarray(phi_bar, dtype=float)
    return {
        "phi_tot": np.clip(prm.phi_f_bg + phi_bar, 0.0, prm.phi_max),
        "phi_c_bulk_target": clip_phi(phi_bar, prm),
        "T_degC": np.full_like(phi_bar, prm.T, dtype=float),
        "P_Pa": np.full_like(phi_bar, prm.P, dtype=float),
        "d_p_m": np.full_like(phi_bar, prm.d_rep, dtype=float),
        "U_bulk": np.maximum(np.abs(np.asarray(U, dtype=float)), 0.0),
    }


@dataclass
class SegmentClosures:
    Ucrit: np.ndarray
    tauw: np.ndarray
    gamma_w: np.ndarray
    tau_yw: np.ndarray
    mu_w: np.ndarray
    rp: np.ndarray
    tau_y0: np.ndarray
    beta_eq_base: np.ndarray
    beta_eq_fsi: np.ndarray
    phi0_eff: np.ndarray
    phiw_eff: np.ndarray



def beta_equilibrium_truth(phi_bar: np.ndarray, U: np.ndarray, a_perp: np.ndarray, prm: ModelParams15D) -> np.ndarray:
    assert prm.closure_model is not None
    inp = _truth_input_pack(phi_bar, U, prm)
    base = prm.closure_model.backbone_closures(
        inp["phi_tot"], inp["phi_c_bulk_target"], inp["T_degC"], inp["P_Pa"], inp["d_p_m"], inp["U_bulk"]
    )["beta_eq"]
    beta_eq = base + prm.kappa_a * np.asarray(a_perp, dtype=float)
    return np.clip(beta_eq, -prm.beta_max, prm.beta_max)



def _extract_ucrit_from_table(
    phi_bar: np.ndarray,
    beta_curr: np.ndarray,
    delta: np.ndarray,
    rho_m: np.ndarray,
    prm: ModelParams15D,
) -> np.ndarray:
    """
    Extract a transport threshold from the truth-table backbone by scanning the U-axis.

    Important limitation:
    the current closure table does not store an explicit U_crit column. This routine therefore
    derives a provisional Ucrit from truth-backed tau_w/r_p/phi_center/phi_wall along the
    tabulated U-axis, then applies the existing transport criterion on top of that backbone.
    """
    assert prm.closure_model is not None
    Ugrid = prm.closure_model.u_axis
    nseg = len(phi_bar)
    if Ugrid.size == 0:
        raise RuntimeError("Closure table has an empty U-axis")

    phi_bar = clip_phi(np.asarray(phi_bar, dtype=float), prm)
    beta_curr = np.asarray(beta_curr, dtype=float)
    D_eff = effective_diameter(delta, prm)
    rho_m = np.asarray(rho_m, dtype=float)

    # Build all query points in one vectorized batch.
    Umesh = np.tile(Ugrid[None, :], (nseg, 1))
    phi_bar_mesh = np.tile(phi_bar[:, None], (1, Ugrid.size))
    inp = _truth_input_pack(phi_bar_mesh.ravel(), Umesh.ravel(), prm)
    base = prm.closure_model.backbone_closures(
        inp["phi_tot"], inp["phi_c_bulk_target"], inp["T_degC"], inp["P_Pa"], inp["d_p_m"], inp["U_bulk"]
    )

    tauw_base = base["tau_w_Pa"].reshape(nseg, Ugrid.size)
    rp_base = base["r_p_m"].reshape(nseg, Ugrid.size)
    beta_eq_base = base["beta_eq"].reshape(nseg, Ugrid.size)
    phi_center_base = base["phi_center"].reshape(nseg, Ugrid.size)
    phi_wall_base = base["phi_wall"].reshape(nseg, Ugrid.size)

    dbeta = beta_curr[:, None] - beta_eq_base
    phi0_eff = clip_phi(phi_center_base + dbeta, prm)
    phiw_eff = clip_phi(phi_wall_base - dbeta, prm)

    tauw_eff = tauw_base * np.clip(
        1.0 + prm.tauw_beta_feedback * dbeta / max(prm.beta_max, 1.0e-9),
        0.50,
        2.00,
    )
    rp_eff = rp_base * np.clip(
        1.0 + prm.rp_beta_feedback * np.maximum(dbeta, 0.0) / max(prm.beta_max, 1.0e-9),
        0.50,
        2.00,
    )

    Ucand = Umesh
    Dm = D_eff[:, None]
    rhom = rho_m[:, None]

    mu_w = apparent_viscosity_hb_local(Ucand, Dm, phiw_eff, prm)
    tau_y0, _, _ = hb_matrix_params(np.clip(prm.phi_f_bg + phi0_eff, 0.0, prm.phi_max), prm)

    Utr = prm.Re_transition * mu_w / np.maximum(rhom * Dm, 1.0e-9)
    Utr = np.clip(Utr, 1.0e-3, 20.0)

    gamma_settle = 1.0
    mu_settle = tau_y0/gamma_settle + 0.45*(gamma_settle**(0.55 - 1.0))
    Y = tau_y0 / np.maximum((prm.rho_cstar - prm.rho_l)*9.81*prm.d_rep, 1.0e-9)
    trap = np.clip(1.0 - Y/prm.Y_crit, 0.0, 1.0)
    w0 = (prm.rho_cstar - prm.rho_l)*9.81*(prm.d_rep**2) / np.maximum(18.0*mu_settle, 1.0e-9)
    Usl = prm.slip_sf * np.clip(
        w0 * np.power(np.clip(1.0 - phi_bar_mesh, 1.0e-6, 1.0), prm.hinder_exp) * trap,
        0.0,
        10.0,
    )

    tau_mob = prm.tau_mob_fac * tau_y0 * (1.0 + 0.8*rp_eff / np.maximum(0.5*Dm, 1.0e-9))
    gamma_req = np.power(np.maximum(tau_mob - tau_y0, 0.0) / np.maximum(0.45, 1.0e-12), 1.0/0.55)
    Umb = np.clip(0.125 * Dm * gamma_req, 0.0, 20.0)

    rhs = prm.ucrit_margin * np.maximum.reduce([Utr, Usl, Umb])
    margin_ratio = Ucand / np.maximum(rhs, 1.0e-9)

    Ucrit = np.empty(nseg, dtype=float)
    for i in range(nseg):
        idx = np.flatnonzero(margin_ratio[i] >= 1.0)
        if idx.size == 0:
            Ucrit[i] = float(max(Ugrid[-1], rhs[i, -1]))
            continue
        j = int(idx[0])
        if j == 0:
            Ucrit[i] = float(max(Ugrid[0], rhs[i, 0]))
            continue

        r0 = margin_ratio[i, j - 1]
        r1 = margin_ratio[i, j]
        u0 = Ugrid[j - 1]
        u1 = Ugrid[j]
        if not np.isfinite(r0) or not np.isfinite(r1) or abs(r1 - r0) < 1.0e-12:
            Ucrit[i] = float(u1)
        else:
            alpha = np.clip((1.0 - r0) / (r1 - r0), 0.0, 1.0)
            Ucrit[i] = float(u0 + alpha*(u1 - u0))

    return Ucrit



def truth_segment_closures(
    phi_bar: np.ndarray,
    beta: np.ndarray,
    delta: np.ndarray,
    U: np.ndarray,
    rho_m: np.ndarray,
    prm: ModelParams15D,
    a_perp: np.ndarray | None = None,
) -> SegmentClosures:
    assert prm.closure_model is not None
    phi_bar = clip_phi(np.asarray(phi_bar, dtype=float), prm)
    beta = np.asarray(beta, dtype=float)
    U = np.asarray(U, dtype=float)
    rho_m = np.asarray(rho_m, dtype=float)
    D_eff = effective_diameter(delta, prm)

    inp = _truth_input_pack(phi_bar, U, prm)
    base = prm.closure_model.backbone_closures(
        inp["phi_tot"], inp["phi_c_bulk_target"], inp["T_degC"], inp["P_Pa"], inp["d_p_m"], inp["U_bulk"]
    )

    beta_eq_base = base["beta_eq"]
    if a_perp is None:
        beta_eq_fsi = beta_eq_base
    else:
        beta_eq_fsi = np.clip(beta_eq_base + prm.kappa_a*np.asarray(a_perp, dtype=float), -prm.beta_max, prm.beta_max)

    dbeta = beta - beta_eq_base
    phi0_eff = clip_phi(base["phi_center"] + dbeta, prm)
    phiw_eff = clip_phi(base["phi_wall"] - dbeta, prm)

    tauw = base["tau_w_Pa"] * np.clip(
        1.0 + prm.tauw_beta_feedback * dbeta / max(prm.beta_max, 1.0e-9),
        0.50,
        2.00,
    )
    rp = base["r_p_m"] * np.clip(
        1.0 + prm.rp_beta_feedback * np.maximum(dbeta, 0.0) / max(prm.beta_max, 1.0e-9),
        0.50,
        2.00,
    )

    tau_yw, K_w, n_w = hb_wall_params(phiw_eff, prm)
    gamma_w = 8.0*np.abs(U) / np.maximum(D_eff, 1.0e-6)
    mu_w = tau_yw/(gamma_w + 1.0e-9) + K_w*np.power(gamma_w + 1.0e-9, n_w - 1.0)
    tau_y0, _, _ = hb_matrix_params(np.clip(prm.phi_f_bg + phi0_eff, 0.0, prm.phi_max), prm)
    mu_w = np.clip(mu_w, 1.0e-4, 1.0e6)

    Ucrit = _extract_ucrit_from_table(phi_bar, beta, delta, rho_m, prm)
    return SegmentClosures(
        Ucrit=Ucrit,
        tauw=tauw,
        gamma_w=gamma_w,
        tau_yw=tau_yw,
        mu_w=mu_w,
        rp=rp,
        tau_y0=tau_y0,
        beta_eq_base=beta_eq_base,
        beta_eq_fsi=beta_eq_fsi,
        phi0_eff=phi0_eff,
        phiw_eff=phiw_eff,
    )


# =============================================================================
# Boundary conditions / basic numerics
# =============================================================================


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


def laplacian_1d(q: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(q)
    out[1:-1] = (q[:-2] - 2.0*q[1:-1] + q[2:]) / (dx*dx)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def scalar_upwind_step(phi: np.ndarray, U: np.ndarray, dt: float, t: float, prm: ModelParams15D):
    rho = mixture_density(phi, prm)
    _, u_ext, phi_ext, _ = extend_with_ghosts(np.zeros_like(phi), U, phi, rho, t, prm)
    u_face = 0.5*(u_ext[:-1] + u_ext[1:])
    phi_up = np.where(u_face >= 0.0, phi_ext[:-1], phi_ext[1:])
    F = u_face * phi_up
    phi_new = phi - (dt/prm.dx)*(F[1:] - F[:-1])
    return clip_phi(phi_new, prm)


# =============================================================================
# Fluid, segregation, deposition, structure
# =============================================================================


def fluid_step(p: np.ndarray, U: np.ndarray, phi_bar: np.ndarray, beta: np.ndarray, delta: np.ndarray,
               dt: float, t: float, prm: ModelParams15D):
    rho = mixture_density(phi_bar, prm)
    D_eff = effective_diameter(delta, prm)

    cls = truth_segment_closures(phi_bar, beta, delta, U, rho, prm)
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

    S_fric = -(4.0*cls.tauw / np.maximum(rho*D_eff, 1.0e-9)) * np.sign(U)
    U_new += dt * S_fric

    u_ref = quasi_steady_velocity_target(t, prm)
    U_new += dt * (-prm.unsteady_damping_beta * (U_new - u_ref))

    nu_num = prm.c_art * prm.a_wave * prm.dx
    eps = min(nu_num * dt / (prm.dx*prm.dx), 0.20)
    p_new = apply_artificial_viscosity(p_new, eps)
    U_new = apply_artificial_viscosity(U_new, eps)

    sigma = np.zeros_like(p_new)
    sigma[-prm.sponge_n:] = np.linspace(0.0, 1.0, prm.sponge_n)
    p_new = (1.0 - 0.12*sigma)*p_new + (0.12*sigma)*prm.p_out_dyn
    U_new = (1.0 - 0.18*sigma)*U_new + (0.18*sigma)*u_ref

    return p_new, U_new, cls



def segregation_step(beta: np.ndarray, U: np.ndarray, beta_eq_loc: np.ndarray,
                     dt: float, prm: ModelParams15D):
    beta_adv = scalar_upwind_step(beta, U, dt, 0.0, prm)
    relax = -(beta_adv - beta_eq_loc) / max(prm.tau_m, 1.0e-6)
    diff = prm.D_beta * laplacian_1d(beta_adv, prm.dx)
    beta_new = beta_adv + dt*(relax + diff)
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



def project_modal(field: np.ndarray, modes: np.ndarray, dx: float) -> np.ndarray:
    return np.array([np.sum(field * modes[k]) * dx for k in range(modes.shape[0])])


def transverse_acceleration(qddot: np.ndarray, phi_modes: np.ndarray) -> np.ndarray:
    return np.sum(qddot[:, None] * phi_modes, axis=0)



def structure_step(q: np.ndarray, qdot: np.ndarray, w: np.ndarray, wdot: np.ndarray,
                   p: np.ndarray, U: np.ndarray, rho_m: np.ndarray, prm: ModelParams15D, dt: float):
    p_ref = np.mean(p)
    rho_ref = np.mean(rho_m)
    U_ref = np.mean(U)

    p_proj = project_modal(p - p_ref, prm.phi_modes, prm.dx)
    rho_proj = project_modal(rho_m - rho_ref, prm.phi_modes, prm.dx)
    U_proj = project_modal(U - U_ref, prm.phi_modes, prm.dx)

    omega_s = 2.0*np.pi*prm.St*np.maximum(prm.Uo_ext, 1.0e-6) / prm.D * np.ones(prm.M)
    F_int = prm.c_rho_modal * rho_proj + prm.c_p_modal * p_proj + prm.c_u_modal * U_proj

    wddot = -prm.wake_eps*omega_s*(w*w - 1.0)*wdot - omega_s**2 * w + prm.wake_A * 0.5 * qdot
    wdot_new = wdot + dt*wddot
    w_new = w + dt*wdot_new

    F_viv = prm.viv_force_gain * w_new

    qddot = F_viv + F_int - 2.0*prm.mode_zeta*prm.mode_omega*qdot - (prm.mode_omega**2)*q
    qdot_new = qdot + dt*qddot
    q_new = q + dt*qdot_new
    return q_new, qdot_new, w_new, wdot_new, qddot


# =============================================================================
# Initialization, diagnostics, time step, main loop
# =============================================================================


def initialize_steady_state(prm: ModelParams15D):
    U = np.full(prm.N, prm.u_ss)
    phi_c_bar = np.full(prm.N, prm.phi_c_inlet)
    beta = np.zeros(prm.N)
    delta = np.zeros(prm.N)

    i_seed0 = min(40, max(prm.N - 2, 0))
    i_seed1 = min(46, prm.N)
    seed = slice(i_seed0, i_seed1)
    phi_c_bar[seed] += 0.05
    phi_c_bar = clip_phi(phi_c_bar, prm)

    rho = mixture_density(phi_c_bar, prm)
    cls = truth_segment_closures(phi_c_bar, beta, delta, U, rho, prm)
    D_eff = effective_diameter(delta, prm)
    dpdx = 4.0*cls.tauw / np.maximum(D_eff, 1.0e-9)

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



def section_residual_pressure(p: np.ndarray, U: np.ndarray, phi_bar: np.ndarray, beta: np.ndarray,
                              delta: np.ndarray, prm: ModelParams15D):
    i0, i1 = prm.i_mon0, prm.i_mon1
    rho = mixture_density(phi_bar, prm)
    cls = truth_segment_closures(phi_bar, beta, delta, U, rho, prm)
    p_drop_meas = p[i0] - p[i1]
    D_eff = effective_diameter(delta, prm)
    dp_fric_model = np.sum((4.0*cls.tauw[i0:i1+1] / np.maximum(D_eff[i0:i1+1], 1.0e-9)) * prm.dx)
    r_dp = p_drop_meas - dp_fric_model
    return r_dp, p_drop_meas, dp_fric_model, cls.Ucrit



def stable_dt(state: State15D, prm: ModelParams15D) -> float:
    wave = prm.a_wave + np.max(np.abs(state.U))
    dt_wave = prm.cfl * prm.dx / max(wave, 1.0e-9)
    dt_adv = prm.cfl * prm.dx / max(np.max(np.abs(state.U)), 1.0e-9)
    dt_beta = 0.5 * prm.tau_m
    dt_struct = 0.2 / max(np.max(prm.mode_omega), 1.0e-9)

    rho = mixture_density(state.phi_c_bar, prm)
    cls = truth_segment_closures(state.phi_c_bar, state.beta, state.delta, state.U, rho, prm)
    dep_drive = np.maximum(cls.Ucrit - np.abs(state.U), 0.0) / np.maximum(cls.Ucrit, 1.0e-9)
    delta_dot_est = prm.k_dep * dep_drive * np.maximum(cls.phiw_eff - prm.phi_c_eq_wall, 0.0)
    if np.any(delta_dot_est > 1.0e-12):
        dt_dep = 0.2*np.min((prm.delta_max - state.delta[delta_dot_est > 1.0e-12]) / delta_dot_est[delta_dot_est > 1.0e-12])
        dt_dep = max(dt_dep, 1.0e-5)
    else:
        dt_dep = prm.dt_max

    return max(1.0e-4, min(prm.dt_max, dt_wave, dt_adv, dt_beta, dt_struct, dt_dep))



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

        p, U, cls_old = fluid_step(state.p, state.U, state.phi_c_bar, state.beta, state.delta, dt, t, prm)

        rho_m = mixture_density(state.phi_c_bar, prm)
        q, qdot, w, wdot, qddot = structure_step(
            state.q, state.qdot, state.w, state.wdot, state.p, state.U, rho_m, prm, dt
        )
        a_perp = transverse_acceleration(qddot, prm.phi_modes)

        beta_eq_loc = beta_equilibrium_truth(state.phi_c_bar, state.U, a_perp, prm)
        beta = segregation_step(state.beta, state.U, beta_eq_loc, dt, prm)

        phi_c_bar = coarse_fraction_step(state.phi_c_bar, state.U, beta, state.delta, cls_old.Ucrit, dt, t, prm)
        delta, delta_dot = deposit_step(state.delta, state.U, state.phi_c_bar, beta, cls_old.Ucrit, dt, prm)

        p = np.nan_to_num(p, nan=0.0, posinf=1.0e8, neginf=-1.0e8)
        U = np.clip(np.nan_to_num(U), -5.0, 5.0)
        phi_c_bar = np.clip(np.nan_to_num(phi_c_bar), 0.0, prm.phi_max)
        beta = np.clip(np.nan_to_num(beta), -prm.beta_max, prm.beta_max)
        delta = np.clip(np.nan_to_num(delta), 0.0, prm.delta_max)

        state = State15D(
            p=p, U=U, phi_c_bar=phi_c_bar, beta=beta, delta=delta,
            q=np.nan_to_num(q), qdot=np.nan_to_num(qdot), w=np.nan_to_num(w), wdot=np.nan_to_num(wdot)
        )
        t += dt

        r_dp, _, _, Ucrit_now = section_residual_pressure(state.p, state.U, state.phi_c_bar, state.beta, state.delta, prm)
        dres = (r_dp - r_prev) / max(dt, 1.0e-9)
        dres_filt = alpha_filter*dres_filt + (1.0 - alpha_filter)*dres
        r_prev = r_dp

        hist["t"].append(t)
        hist["U_in"].append(state.U[0])
        hist["U_mid"].append(state.U[prm.N // 2])
        hist["U_out"].append(state.U[-1])
        hist["resid_dp"].append(r_dp)
        hist["dresid_dt"].append(dres_filt)
        i_seed_hist = min(42, prm.N - 1)
        hist["seed_phi"].append(state.phi_c_bar[i_seed_hist])
        hist["seed_beta"].append(state.beta[i_seed_hist])
        hist["seed_delta"].append(state.delta[i_seed_hist])
        hist["max_q"].append(float(np.max(np.abs(state.q))))
        hist["max_aperp"].append(float(np.max(np.abs(a_perp))))
        hist["min_ucrit_ratio"].append(float(np.min(np.abs(state.U) / np.maximum(Ucrit_now, 1.0e-9))))

        # 進捗を0.5秒ごとに画面に出力する
        if len(hist["t"]) % 500 == 0:
            print(f"Time: {t:.2f} / {prm.t_end:.2f} s (dt={dt:.2e} s)")

    out = {k: np.asarray(v) for k, v in hist.items()}
    out.update({
        "x": prm.x.copy(),
        "final_p": state.p.copy(),
        "final_U": state.U.copy(),
        "final_phi_c_bar": state.phi_c_bar.copy(),
        "final_beta": state.beta.copy(),
        "final_delta": state.delta.copy(),
        "final_q": state.q.copy(),
        "final_w": state.w.copy(),
    })
    return out


# =============================================================================
# CLI / main
# =============================================================================


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="HSR 1.5D transient solver v4 with truth-table-backed closures")
    p.add_argument("--closure-table", type=str, default="closure_table.csv", help="Path to closure_table.csv")
    p.add_argument("--t-end", type=float, default=40.0)
    p.add_argument("--n", type=int, default=550)
    p.add_argument("--dt-max", type=float, default=0.01)
    p.add_argument("--u-ss", type=float, default=1.60)
    p.add_argument("--no-plot", action="store_true")
    return p.parse_args(argv)



def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    prm = ModelParams15D(
        closure_table_path=args.closure_table,
        t_end=args.t_end,
        N=args.n,
        dt_max=args.dt_max,
        u_ss=args.u_ss,
    )
    out = run(prm)

    print("Simulation complete")
    print(f"closure table: {Path(prm.closure_table_path)}")
    print(f"closure grid type: {'tensor' if prm.closure_model and prm.closure_model.is_tensor_grid else 'scattered'}")
    print(f"steady inlet dynamic pressure: {prm.p_in0:.3e} Pa")
    print(f"max |d(residual Δp)/dt|: {np.max(np.abs(out['dresid_dt'])):.3e} Pa/s")
    print(f"final max deposit thickness: {np.max(out['final_delta']):.4e} m")
    print(f"final min velocity: {np.min(out['final_U']):.4f} m/s")
    print(f"final max velocity: {np.max(out['final_U']):.4f} m/s")
    print(f"final max |beta|: {np.max(np.abs(out['final_beta'])):.4e}")
    print(f"final max |q|: {np.max(np.abs(out['final_q'])):.4e} m")

    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
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

            fig.tight_layout()
            png_path = Path(__file__).with_name("slurry_15d_transient_v4_demo.png")
            fig.savefig(png_path, dpi=160)
            print(f"saved demo plot: {png_path}")
        except Exception as exc:
            print(f"plot skipped: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
