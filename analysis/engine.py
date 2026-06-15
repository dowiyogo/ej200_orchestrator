#!/usr/bin/env python3.12
"""
Shared analysis engine for EXEC_14.
Motor híbrido: uproot lee hits, PyROOT construye sidecars.

Importa common.py desde ej200_endonly (single source of truth).
Reimplementa leading_edge_time / earliest desde endonly_sum4.py sin depender
de scipy (usamos TH1::Fit vía PyROOT en lugar de curve_fit).
"""
from __future__ import annotations

# Must be set BEFORE importing numpy to avoid OpenBLAS/ROOT signal-handler conflict
import os as _os
_os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
_os.environ.setdefault('OMP_NUM_THREADS', '1')

import hashlib
import json
import math
import pathlib
import sys
import datetime
from typing import NamedTuple

import numpy as np
import uproot
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")

# ── Single-source channel map ─────────────────────────────────────────────────
sys.path.insert(0, '/home/reriosto/SHiP/ej200_endonly/analysis/exec07')
import common  # noqa: E402 — must be after sys.path modification

# ── SPR / SUM4 constants (mirroring endonly_sum4.py / congruent_sum4_timing.C) ─
SPR_RISE_NS = common.SPR_RISE_NS        # 0.5 ns
SPR_FALL_NS = common.SPR_FALL_NS        # 5.0 ns
LEADING_EDGE_THRESHOLD_PE = common.LEADING_EDGE_THRESHOLD_PE  # 4.0

# Simulation jitter per hit (SiPMSD.cc:79 / SiPMSD.hh:49)
JITTER_PER_HIT_NS = 20.0e-3   # 20 ps default, NOT overridden in scan.mac

# Python analysis post-processing readout jitter (endonly_sum4.py)
READOUT_JITTER_QUADRATURE_PS = 20.0

# Double-counting note: both the per-hit sim jitter AND the quadrature jitter
# are 20 ps and both model "readout electronics timing resolution". This flag
# documents the risk; it does NOT fix anything (per René/Gerardo decision).
JITTER_DOUBLE_COUNT_RISK = (
    "Per-hit jitter in time_ns (SiPMSD.cc fJitterSigma=20ps default, scan.mac "
    "line commented → uses default) PLUS readout jitter in quadrature "
    "(endonly_sum4.py DEFAULT_READOUT_JITTER_PS=20). Same physical value at two "
    "stages. Likely double-counting. Pending René/Gerardo decision."
)

SQRT2 = math.sqrt(2.0)

# ── Time convention (verified from SiPMSD.cc:80 and data) ────────────────────
TIME_CONVENTION = (
    "time_ns = (GetGlobalTime() + G4RandGauss(0, 20ps)) / ns. "
    "t=0 is primary muon generation (Geant4 global clock). "
    "Includes: muon transit (~ps, negligible) + scintillation emission delay "
    "+ optical propagation to SiPM + per-hit Gaussian jitter (20ps sigma). "
    "Verified: t_min per event at x=0 ~ 2.8 ns, consistent with 700mm / 27.7 cm/ns."
)


# ── SHA-256 ───────────────────────────────────────────────────────────────────
def sha256_file(path: pathlib.Path | str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify_sha(path: pathlib.Path | str, expected: str, label: str = "") -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise RuntimeError(
            f"SHA-256 mismatch for {label or path}\n"
            f"  expected: {expected}\n"
            f"  actual:   {actual}"
        )


# ── ROOT file health check ────────────────────────────────────────────────────
def verify_root_input(path: pathlib.Path | str) -> None:
    path = str(path)
    with uproot.open(path) as f:
        if 'sipm_hits' not in f:
            raise RuntimeError(f"{path}: TTree 'sipm_hits' not found")
    # Check file size > 1 KB
    if pathlib.Path(path).stat().st_size < 1024:
        raise RuntimeError(f"{path}: file too small (<1 KB), likely corrupt")


# ── Load hits (uproot, memory-efficient) ─────────────────────────────────────
def load_end_hits(path: pathlib.Path | str,
                  branches: list[str] | None = None
                  ) -> dict[str, np.ndarray]:
    """
    Load end SiPM hits only (global_id < 16).
    Returns dict with n_events and n_hits added.
    Caller must del the returned dict to free memory.
    """
    if branches is None:
        branches = ['event_id', 'global_id', 'time_ns']
    with uproot.open(str(path)) as f:
        tree = f['sipm_hits']
        n_hits = int(tree.num_entries)
        arrays = tree.arrays(branches, library='np')

    event_id = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)

    n_events = int(np.max(event_id)) + 1
    if not (1500 <= n_events <= 3000):
        raise RuntimeError(
            f"{path}: n_events={n_events} outside [1500,3000]. "
            f"Hint: n_hits={n_hits} — do NOT use num_entries as n_events."
        )

    end_mask = (global_id >= 0) & (global_id < 16)
    result = {k: arrays[k][end_mask] for k in branches}
    result['_n_events'] = np.int64(n_events)
    result['_n_hits_total'] = np.int64(n_hits)
    result['_n_hits_end'] = np.int64(int(end_mask.sum()))
    return result


def load_all_hits(path: pathlib.Path | str,
                  branches: list[str] | None = None
                  ) -> dict[str, np.ndarray]:
    """Load ALL hits (End + Top). Returns with n_events, n_hits metadata."""
    if branches is None:
        branches = ['event_id', 'global_id', 'time_ns']
    with uproot.open(str(path)) as f:
        tree = f['sipm_hits']
        n_hits = int(tree.num_entries)
        arrays = tree.arrays(branches, library='np')

    event_id = arrays['event_id'].astype(np.int64)
    n_events = int(np.max(event_id)) + 1
    if not (1500 <= n_events <= 3000):
        raise RuntimeError(
            f"{path}: n_events={n_events} outside [1500,3000] (n_hits={n_hits})"
        )

    result = dict(arrays)
    result['_n_events'] = np.int64(n_events)
    result['_n_hits_total'] = np.int64(n_hits)
    return result


# ── SUM4 leading-edge time (mirror of endonly_sum4.py / congruent_sum4_timing.C) ─

def _pulse(slow: float, fast: float, dt: float) -> float:
    peak = SPR_RISE_NS * SPR_FALL_NS / (SPR_FALL_NS - SPR_RISE_NS) * math.log(SPR_FALL_NS / SPR_RISE_NS)
    norm = 1.0 / (math.exp(-peak / SPR_FALL_NS) - math.exp(-peak / SPR_RISE_NS))
    return norm * (slow * math.exp(-dt / SPR_FALL_NS) - fast * math.exp(-dt / SPR_RISE_NS))


def leading_edge_time(arrivals: np.ndarray) -> float:
    """Exact port of congruent_sum4_timing.C LeadingEdgeTime."""
    if arrivals.size == 0:
        return math.nan
    arrivals = np.sort(arrivals)
    slow, fast, idx = 0.0, 0.0, 0
    while idx < arrivals.size:
        cur = float(arrivals[idx])
        nxt = idx
        while nxt < arrivals.size and arrivals[nxt] == cur:
            slow += 1.0; fast += 1.0; nxt += 1
        interval = float(arrivals[nxt] - cur) if nxt < arrivals.size else math.inf
        deriv0 = fast / SPR_RISE_NS - slow / SPR_FALL_NS
        if deriv0 > 0.0:
            pk = math.log(fast * SPR_FALL_NS / (slow * SPR_RISE_NS)) / (1.0 / SPR_RISE_NS - 1.0 / SPR_FALL_NS)
            reach = min(pk, interval)
            if reach >= 0.0 and _pulse(slow, fast, reach) >= LEADING_EDGE_THRESHOLD_PE:
                lo, hi = 0.0, reach
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if _pulse(slow, fast, mid) >= LEADING_EDGE_THRESHOLD_PE:
                        hi = mid
                    else:
                        lo = mid
                return cur + hi
        if nxt >= arrivals.size:
            break
        slow *= math.exp(-interval / SPR_FALL_NS)
        fast *= math.exp(-interval / SPR_RISE_NS)
        idx = nxt
    return math.nan


def earliest(a: float, b: float) -> float:
    if not math.isfinite(a): return b
    if not math.isfinite(b): return a
    return min(a, b)


# ── SUM4 per-event trigger times ─────────────────────────────────────────────
class SUM4Result(NamedTuple):
    delta_lr: np.ndarray   # ΔT_LR = t_L - t_R for events with both triggers
    t_left:   np.ndarray   # t_L per event (NaN if no trigger)
    t_right:  np.ndarray   # t_R per event (NaN if no trigger)
    n_events: int


def compute_sum4(arrays: dict[str, np.ndarray]) -> SUM4Result:
    """Compute SUM4 leading-edge trigger times from loaded End hits."""
    n_events = int(arrays['_n_events'])
    event_id = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    time_ns   = arrays['time_ns'].astype(np.float64)

    # Group by (event, cluster): cluster = global_id // 4, groups 0-3
    combined = event_id * 4 + global_id // 4
    order = np.argsort(combined, kind='stable')
    combined_s = combined[order]
    time_s = time_ns[order]

    triggers = np.full((n_events, 4), np.nan)
    unique_keys, starts = np.unique(combined_s, return_index=True)
    stops = np.r_[starts[1:], combined_s.size]

    for key, start, stop in zip(unique_keys, starts, stops):
        ev = int(key // 4)
        grp = int(key % 4)
        if 0 <= ev < n_events and 0 <= grp < 4:
            triggers[ev, grp] = leading_edge_time(time_s[start:stop])

    t_left  = np.array([earliest(triggers[i, 0], triggers[i, 1]) for i in range(n_events)])
    t_right = np.array([earliest(triggers[i, 2], triggers[i, 3]) for i in range(n_events)])
    both = np.isfinite(t_left) & np.isfinite(t_right)
    delta_lr = t_left[both] - t_right[both]

    return SUM4Result(delta_lr=delta_lr, t_left=t_left, t_right=t_right, n_events=n_events)


# ── TH1::Fit Gaussian (replaces curve_fit of endonly_sum4.py) ────────────────
class FitCoreResult(NamedTuple):
    mean_ns: float
    sigma_ns: float
    sigma_err_ns: float
    chi2_ndf: float
    rms_ns: float
    n: int
    fit_used: bool


def fit_core_root(values: np.ndarray, name: str) -> FitCoreResult:
    """
    Gaussian core fit using TH1::Fit (PyROOT).
    Mirrors congruent_sum4_timing.C FitCore: 100 bins, 4 iterative ±2σ fits.
    Uses TH1::Fit with Poisson weights (unlike endonly_sum4.py curve_fit).
    """
    n = len(values)
    rms = float(np.std(values, ddof=1)) if n > 1 else math.nan

    if n < 20:
        mean = float(np.mean(values)) if n else math.nan
        return FitCoreResult(mean, rms, rms / math.sqrt(max(2 * (n - 1), 1)), -1.0, rms, n, False)

    center = float(np.median(values))
    sigma  = max(1.4826 * float(np.median(np.abs(values - center))), 1e-4)
    lo, hi = center - 8.0 * sigma, center + 8.0 * sigma

    h = ROOT.TH1D(f"_fcr_{name}", "", 100, lo, hi)
    h.SetDirectory(ROOT.nullptr)  # Python owns; safe to Delete() after use
    for v in values:
        h.Fill(float(v))

    gaus = ROOT.TF1(f"_gaus_{name}", "gaus", lo, hi)
    fit_used = False

    for _ in range(4):
        fl, fh = max(lo, center - 2.0 * sigma), min(hi, center + 2.0 * sigma)
        gaus.SetRange(fl, fh)
        gaus.SetParameters(h.GetMaximum(), center, sigma)
        gaus.SetParLimits(2, 1e-6, hi - lo)
        result = h.Fit(gaus, "QRNS")
        if result.Status() != 0 or gaus.GetParameter(2) <= 0.0:
            break
        center = gaus.GetParameter(1)
        sigma  = abs(gaus.GetParameter(2))
        fit_used = True

    # Read all values from gaus BEFORE releasing objects.
    # CRITICAL: use `del` (Python GC), NEVER `.Delete()` (double-free → segfault).
    if fit_used:
        final_mean    = float(gaus.GetParameter(1))
        final_sigma   = abs(float(gaus.GetParameter(2)))
        final_sig_err = float(gaus.GetParError(2))
        chi2          = float(gaus.GetChisquare())
        ndf           = int(gaus.GetNDF())
        chi2_ndf      = chi2 / ndf if ndf > 0 else -1.0

    # Let Python ref-counting call destructors exactly once
    del h
    del gaus

    if not fit_used:
        mean = float(np.mean(values))
        return FitCoreResult(mean, rms, rms / math.sqrt(max(2 * (n - 1), 1)), -1.0, rms, n, False)

    return FitCoreResult(
        mean_ns=final_mean,
        sigma_ns=final_sigma,
        sigma_err_ns=final_sig_err if final_sig_err >= 0 else math.nan,
        chi2_ndf=chi2_ndf,
        rms_ns=rms,
        n=n,
        fit_used=True,
    )


# ── Meta.json writer ──────────────────────────────────────────────────────────
def write_meta(path: pathlib.Path, data: dict) -> None:
    """Write meta.json, ensuring all numpy scalars are JSON-serializable."""
    def _convert(obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray):    return obj.tolist()
        raise TypeError(f"Not JSON serializable: {type(obj)}")
    path.write_text(json.dumps(data, indent=2, default=_convert))


def write_csv_arrays(path: pathlib.Path, columns: dict[str, list | np.ndarray]) -> None:
    """Write CSV with one column per key."""
    import csv
    names = list(columns.keys())
    rows = list(zip(*[columns[k] for k in names]))
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(names)
        w.writerows(rows)
