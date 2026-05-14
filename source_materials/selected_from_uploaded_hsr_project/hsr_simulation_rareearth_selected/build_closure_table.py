#!/usr/bin/env python3
from __future__ import annotations

# -----------------------------------------------------------------------------
# Best-effort single-thread enforcement for BLAS / OpenMP backends.
# Must happen before importing numpy/scipy or the truth-generator module.
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
import hashlib
import itertools
import json
import math
import signal
import sys
import time
import traceback
from dataclasses import asdict, is_dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Sequence, Set, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent

# Import after single-thread environment variables have been set.
try:
    if str(SCRIPT_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPT_DIR))
    import axisym_truth_generator_v1 as truthmod
except Exception as exc:  # pragma: no cover - import path failure is fatal at runtime
    raise SystemExit(
        "Failed to import axisym_truth_generator_v1.py from the same directory as this script.\n"
        f"Expected path: {SCRIPT_DIR / 'axisym_truth_generator_v1.py'}\n"
        f"Import error: {exc}"
    )


TruthParams = truthmod.TruthParams
solve_truth_case = truthmod.solve_truth_case


CSV_COLUMNS: List[str] = [
    "case_key",
    "case_hash",
    "status",
    "grid_name",
    "case_index",
    "completed_at_utc",
    "elapsed_s",
    "phi_tot",
    "phi_c_bulk_target",
    "T_degC",
    "P_Pa",
    "d_p_m",
    "U_bulk_target_mps",
    "phi_bar",
    "beta_eq",
    "tau_w_Pa",
    "r_p_m",
    "phi_center",
    "phi_wall",
    "u_center_mps",
    "u_wall_mps",
    "yielded_fraction",
    "G_Pa_per_m",
    "U_bulk_solved_mps",
    "iter_count",
    "iter_last_error",
    "iter_max_error",
    "Nr",
    "R_m",
    "tol",
    "max_iter",
    "pseudo_dt_final",
    "pseudo_diff",
    "under_relax_phi",
    "error_type",
    "error_message",
]


def _linspace(a: float, b: float, n: int) -> List[float]:
    if n <= 1:
        return [float(a)]
    step = (b - a) / (n - 1)
    return [float(a + i * step) for i in range(n)]


DEFAULT_FULL_GRID: Dict[str, Sequence[float]] = {
    # Example production-oriented grid; edit freely for your campaign.
    "phi_tot": (0.24, 0.26, 0.28, 0.30),
    "phi_c_bulk_target": tuple(round(x, 6) for x in _linspace(0.06, 0.18, 10)),
    "T": (1.0, 2.0, 4.0),
    "P": (45e6, 55e6, 65e6),
    "d_p": (50e-6, 80e-6, 120e-6),
    "U_bulk_target": tuple(round(x, 6) for x in _linspace(0.60, 1.80, 10)),
}

DEFAULT_SMALL_GRID: Dict[str, Sequence[float]] = {
    # Smoke test / constrained-resource grid: 3 x 3 = 9 cases total.
    "phi_tot": (0.28,),
    "phi_c_bulk_target": (0.08, 0.12, 0.16),
    "T": (2.0,),
    "P": (55e6,),
    "d_p": (80e-6,),
    "U_bulk_target": (0.80, 1.20, 1.60),
}

_STOP_REQUESTED = False


class GracefulStop(RuntimeError):
    pass


def install_signal_handlers() -> None:
    def _handler(signum, _frame):
        global _STOP_REQUESTED
        _STOP_REQUESTED = True
        name = signal.Signals(signum).name
        print(f"\\n[{utc_now_iso()}] Received {name}. Will stop after the current case finishes.", flush=True)

    for sig in (signal.SIGINT, signal.SIGTERM):
        try:
            signal.signal(sig, _handler)
        except Exception:
            pass


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def float_token(x: float) -> str:
    return format(float(x), ".12g")


def make_case_key(case: Dict[str, float]) -> str:
    return (
        f"phi_tot={float_token(case['phi_tot'])}|"
        f"phi_c_bulk_target={float_token(case['phi_c_bulk_target'])}|"
        f"T={float_token(case['T'])}|"
        f"P={float_token(case['P'])}|"
        f"d_p={float_token(case['d_p'])}|"
        f"U_bulk_target={float_token(case['U_bulk_target'])}"
    )


def make_case_hash(case_key: str) -> str:
    return hashlib.sha1(case_key.encode("utf-8")).hexdigest()[:16]


def normalize_grid(raw_grid: Dict[str, Sequence[float]]) -> Dict[str, List[float]]:
    expected = ["phi_tot", "phi_c_bulk_target", "T", "P", "d_p", "U_bulk_target"]
    missing = [k for k in expected if k not in raw_grid]
    if missing:
        raise ValueError(f"Grid is missing required keys: {missing}")

    grid: Dict[str, List[float]] = {}
    for key in expected:
        vals = [float(v) for v in raw_grid[key]]
        if not vals:
            raise ValueError(f"Grid entry '{key}' is empty")
        grid[key] = vals
    return grid


def generate_cases(grid: Dict[str, Sequence[float]]) -> Iterator[Dict[str, float]]:
    grid_n = normalize_grid(grid)
    keys = ["phi_tot", "phi_c_bulk_target", "T", "P", "d_p", "U_bulk_target"]
    for values in itertools.product(*(grid_n[k] for k in keys)):
        case = dict(zip(keys, values))
        # Physical consistency guard.
        if case["phi_c_bulk_target"] >= case["phi_tot"]:
            continue
        yield case


def grid_size(grid: Dict[str, Sequence[float]]) -> int:
    return sum(1 for _ in generate_cases(grid))


def load_grid_from_json(json_path: Path) -> Dict[str, List[float]]:
    with json_path.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError("Grid JSON must be an object/dict")
    return normalize_grid(payload)


def load_base_param_overrides(json_path: Path | None) -> Dict[str, Any]:
    if json_path is None:
        return {}
    with json_path.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError("Base parameter override JSON must be an object/dict")
    return payload


def make_truth_params(case: Dict[str, float], base_overrides: Dict[str, Any]) -> Any:
    base = TruthParams(**base_overrides)
    replacements = {
        "phi_tot": float(case["phi_tot"]),
        "phi_c_bulk_target": float(case["phi_c_bulk_target"]),
        "T": float(case["T"]),
        "P": float(case["P"]),
        "d_p": float(case["d_p"]),
    }
    return replace(base, **replacements)


def safe_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        if math.isfinite(float(value)):
            return float(value)
        return None
    try:
        value_f = float(value)
    except Exception:
        return None
    return value_f if math.isfinite(value_f) else None


def result_row_from_output(
    *,
    case: Dict[str, float],
    case_key: str,
    case_hash: str,
    case_index: int,
    grid_name: str,
    prm: Any,
    out: Dict[str, Any] | None,
    elapsed_s: float,
    status: str,
    error_type: str = "",
    error_message: str = "",
) -> Dict[str, Any]:
    iter_history = None if out is None else out.get("iter_history")
    iter_count = None
    iter_last_error = None
    iter_max_error = None
    if iter_history is not None:
        try:
            iter_count = int(len(iter_history))
            if iter_count > 0:
                iter_last_error = safe_float(iter_history[-1])
                iter_max_error = safe_float(max(iter_history))
        except Exception:
            pass

    row = {
        "case_key": case_key,
        "case_hash": case_hash,
        "status": status,
        "grid_name": grid_name,
        "case_index": case_index,
        "completed_at_utc": utc_now_iso(),
        "elapsed_s": round(float(elapsed_s), 6),
        "phi_tot": safe_float(case["phi_tot"]),
        "phi_c_bulk_target": safe_float(case["phi_c_bulk_target"]),
        "T_degC": safe_float(case["T"]),
        "P_Pa": safe_float(case["P"]),
        "d_p_m": safe_float(case["d_p"]),
        "U_bulk_target_mps": safe_float(case["U_bulk_target"]),
        "phi_bar": None if out is None else safe_float(out.get("phi_bar")),
        "beta_eq": None if out is None else safe_float(out.get("beta_eq")),
        "tau_w_Pa": None if out is None else safe_float(out.get("tau_w")),
        "r_p_m": None if out is None else safe_float(out.get("r_p")),
        "phi_center": None if out is None else safe_float(out.get("phi_center")),
        "phi_wall": None if out is None else safe_float(out.get("phi_wall")),
        "u_center_mps": None if out is None else safe_float(out.get("u_center")),
        "u_wall_mps": None if out is None else safe_float(out.get("u_wall")),
        "yielded_fraction": None if out is None else safe_float(out.get("yielded_fraction")),
        "G_Pa_per_m": None if out is None else safe_float(out.get("G")),
        "U_bulk_solved_mps": None if out is None else safe_float(out.get("U_bulk")),
        "iter_count": iter_count,
        "iter_last_error": iter_last_error,
        "iter_max_error": iter_max_error,
        "Nr": getattr(prm, "Nr", None),
        "R_m": safe_float(getattr(prm, "R", None)),
        "tol": safe_float(getattr(prm, "tol", None)),
        "max_iter": getattr(prm, "max_iter", None),
        "pseudo_dt_final": safe_float(getattr(prm, "pseudo_dt", None)),
        "pseudo_diff": safe_float(getattr(prm, "pseudo_diff", None)),
        "under_relax_phi": safe_float(getattr(prm, "under_relax_phi", None)),
        "error_type": error_type,
        "error_message": error_message,
    }
    return row


def ensure_csv_has_header(csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    if csv_path.exists() and csv_path.stat().st_size > 0:
        return
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        f.flush()
        os.fsync(f.fileno())


def append_csv_row(csv_path: Path, row: Dict[str, Any]) -> None:
    with csv_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writerow(row)
        f.flush()
        os.fsync(f.fileno())


def load_seen_cases(csv_path: Path, retry_errors: bool) -> Tuple[Set[str], Set[str]]:
    completed_ok: Set[str] = set()
    completed_err: Set[str] = set()
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return completed_ok, completed_err

    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            case_key = (row.get("case_key") or "").strip()
            status = (row.get("status") or "").strip().lower()
            if not case_key:
                continue
            if status == "ok":
                completed_ok.add(case_key)
            elif status == "error":
                completed_err.add(case_key)

    if retry_errors:
        completed_err.clear()
    return completed_ok, completed_err




def to_jsonable(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    return obj

def save_manifest(
    manifest_path: Path,
    *,
    csv_path: Path,
    grid_name: str,
    grid: Dict[str, Sequence[float]],
    base_overrides: Dict[str, Any],
    args: argparse.Namespace,
) -> None:
    payload = {
        "created_at_utc": utc_now_iso(),
        "csv_path": str(csv_path),
        "grid_name": grid_name,
        "grid": to_jsonable({k: list(v) for k, v in grid.items()}),
        "base_overrides": to_jsonable(base_overrides),
        "single_thread_env": {
            key: os.environ.get(key)
            for key in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "BLIS_NUM_THREADS",
            )
        },
        "truth_module": str(Path(truthmod.__file__).resolve()),
        "truth_module_name": truthmod.__name__,
        "script": str(Path(__file__).resolve()),
        "resume_policy": {
            "skip_ok": True,
            "skip_error_when_retry_errors_false": True,
        },
        "args": to_jsonable(vars(args)),
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())


def select_grid(args: argparse.Namespace) -> Tuple[str, Dict[str, Sequence[float]]]:
    if args.grid_json is not None:
        return "custom_json", load_grid_from_json(args.grid_json)
    if args.small_grid:
        return "small", normalize_grid(DEFAULT_SMALL_GRID)
    return "full", normalize_grid(DEFAULT_FULL_GRID)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sequential closure-table builder for axisym_truth_generator_v1.py. "
            "Runs one case at a time, appends one CSV row per finished case, "
            "and resumes automatically on the next launch."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIR / "closure_table.csv",
        help="Output CSV path. Existing rows are used for resume unless --overwrite is given.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Optional manifest JSON path. Default: <output_stem>_manifest.json",
    )
    parser.add_argument(
        "--small-grid",
        action="store_true",
        help="Run the built-in 3x3 smoke-test grid (9 cases total).",
    )
    parser.add_argument(
        "--grid-json",
        type=Path,
        default=None,
        help="Load a custom grid from JSON. Overrides --small-grid / built-in full grid.",
    )
    parser.add_argument(
        "--base-param-json",
        type=Path,
        default=None,
        help="JSON file with TruthParams constructor overrides (e.g. Nr, tol, max_iter).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete existing output CSV and start from scratch.",
    )
    parser.add_argument(
        "--retry-errors",
        action="store_true",
        help="Re-run rows previously recorded with status=error.",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Abort immediately after the first failing case instead of logging the error row and continuing.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of NEW cases to run in this invocation.",
    )
    parser.add_argument(
        "--sleep-sec",
        type=float,
        default=0.0,
        help="Optional sleep inserted after each finished case to reduce system pressure.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the selected grid and resume counts without solving any case.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    install_signal_handlers()

    if args.overwrite and args.output.exists():
        args.output.unlink()

    grid_name, grid = select_grid(args)
    base_overrides = load_base_param_overrides(args.base_param_json)
    manifest_path = args.manifest or args.output.with_name(args.output.stem + "_manifest.json")

    ensure_csv_has_header(args.output)
    save_manifest(
        manifest_path,
        csv_path=args.output,
        grid_name=grid_name,
        grid=grid,
        base_overrides=base_overrides,
        args=args,
    )

    completed_ok, completed_err = load_seen_cases(args.output, retry_errors=args.retry_errors)
    all_cases = list(generate_cases(grid))
    total_cases = len(all_cases)

    pending_cases: List[Tuple[int, Dict[str, float], str, str]] = []
    for idx, case in enumerate(all_cases, start=1):
        case_key = make_case_key(case)
        case_hash = make_case_hash(case_key)
        if case_key in completed_ok:
            continue
        if case_key in completed_err:
            continue
        pending_cases.append((idx, case, case_key, case_hash))

    print(f"[{utc_now_iso()}] Grid name       : {grid_name}")
    print(f"[{utc_now_iso()}] Output CSV      : {args.output}")
    print(f"[{utc_now_iso()}] Manifest JSON   : {manifest_path}")
    print(f"[{utc_now_iso()}] Total cases     : {total_cases}")
    print(f"[{utc_now_iso()}] Completed OK    : {len(completed_ok)}")
    print(f"[{utc_now_iso()}] Completed error : {0 if args.retry_errors else len(completed_err)}")
    print(f"[{utc_now_iso()}] Pending         : {len(pending_cases)}")

    if args.dry_run:
        print(f"[{utc_now_iso()}] Dry-run only. No solves executed.")
        return 0

    if not pending_cases:
        print(f"[{utc_now_iso()}] Nothing to do. Closure table is already complete for this grid.")
        return 0

    new_cases_done = 0
    t_batch0 = time.perf_counter()

    for case_index, case, case_key, case_hash in pending_cases:
        if _STOP_REQUESTED:
            break
        if args.limit is not None and new_cases_done >= args.limit:
            break

        print(
            f"[{utc_now_iso()}] START case {case_index}/{total_cases} | {case_key}",
            flush=True,
        )
        prm = make_truth_params(case, base_overrides)
        t0 = time.perf_counter()

        try:
            out = solve_truth_case(U_bulk_target=float(case["U_bulk_target"]), prm=prm)
            elapsed_s = time.perf_counter() - t0
            row = result_row_from_output(
                case=case,
                case_key=case_key,
                case_hash=case_hash,
                case_index=case_index,
                grid_name=grid_name,
                prm=prm,
                out=out,
                elapsed_s=elapsed_s,
                status="ok",
            )
            append_csv_row(args.output, row)
            new_cases_done += 1
            print(
                f"[{utc_now_iso()}] DONE  case {case_index}/{total_cases} | "
                f"elapsed={elapsed_s:.2f}s | phi_bar={row['phi_bar']:.6g} | "
                f"beta_eq={row['beta_eq']:.6g} | tau_w={row['tau_w_Pa']:.6g} | r_p={row['r_p_m']:.6g}",
                flush=True,
            )
        except Exception as exc:
            elapsed_s = time.perf_counter() - t0
            err_msg = f"{type(exc).__name__}: {exc}"
            row = result_row_from_output(
                case=case,
                case_key=case_key,
                case_hash=case_hash,
                case_index=case_index,
                grid_name=grid_name,
                prm=prm,
                out=None,
                elapsed_s=elapsed_s,
                status="error",
                error_type=type(exc).__name__,
                error_message=err_msg,
            )
            append_csv_row(args.output, row)
            print(
                f"[{utc_now_iso()}] ERROR case {case_index}/{total_cases} | elapsed={elapsed_s:.2f}s | {err_msg}",
                flush=True,
            )
            traceback.print_exc()
            if args.stop_on_error:
                print(f"[{utc_now_iso()}] Stopping because --stop-on-error is set.", flush=True)
                return 2

        if args.sleep_sec > 0.0 and not _STOP_REQUESTED:
            time.sleep(args.sleep_sec)

    batch_elapsed_s = time.perf_counter() - t_batch0
    print(
        f"[{utc_now_iso()}] Finished this invocation. New rows written: {new_cases_done}. "
        f"Wall time: {batch_elapsed_s:.2f}s.",
        flush=True,
    )
    if _STOP_REQUESTED:
        print(f"[{utc_now_iso()}] Stopped cleanly after the current case. Re-run to resume.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
