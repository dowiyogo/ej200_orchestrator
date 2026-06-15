#!/usr/bin/env python3.12
"""
F1 sidecar generator — EXEC_14.

Histograma de tiempo relativo al primer fotón del evento (t_rel = t - t_min^ev).
0–2 ns, bin 2.0 ps (variante 2.5 ps), eje Y log.

Modelo físico de emisión de centelleo (QA-1 corrección):
  I(t) = N * (1-exp(-t/tau_rise)) * (A_fast*exp(-t/tau_fast) + A_slow*exp(-t/tau_slow))
Ajuste vía TH1::Fit en rango [0.004, 2.0] ns (salta spike de t_min a t=0).
Test de significancia: ajuste de un solo modo vs doble modo (F-test vía delta-chi2).
"""
from __future__ import annotations

import os as _os
_os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
_os.environ.setdefault('OMP_NUM_THREADS', '1')

import math
import pathlib
import sys
import datetime

import numpy as np
import uproot
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from engine import (
    sha256_file, verify_sha, verify_root_input,
    write_meta, write_csv_arrays,
    TIME_CONVENTION, JITTER_PER_HIT_NS, READOUT_JITTER_QUADRATURE_PS,
    JITTER_DOUBLE_COUNT_RISK,
)

BRANCHES_F1 = ['event_id', 'global_id', 'time_ns']

# Physical scintillation model (rise + two decay components)
# [0]=N, [1]=tau_rise, [2]=A_fast, [3]=tau_fast, [4]=A_slow, [5]=tau_slow
SCINT_2_FORMULA = "[0]*(1-exp(-x/[1]))*([2]*exp(-x/[3])+[4]*exp(-x/[5]))"
SCINT_1_FORMULA = "[0]*(1-exp(-x/[1]))*[2]*exp(-x/[3])"

# Fit range: skip first 2 bins (avoid t=0 spike from per-event t_min subtraction)
FIT_LO = 0.004   # ns (2 × 2ps)
FIT_HI = 2.0     # ns


def compute_t_rel_end(path: pathlib.Path) -> tuple[np.ndarray, int, int, np.ndarray]:
    """
    Load End SiPM hits, compute t_rel = time_ns - min(time_ns per event).
    Returns: (t_rel_0_to_2ns, n_hits_total, n_events, t_min_per_event)
    """
    with uproot.open(str(path)) as f:
        tree = f['sipm_hits']
        n_hits_total = int(tree.num_entries)
        arrays = tree.arrays(BRANCHES_F1, library='np')

    event_id = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    time_ns   = arrays['time_ns'].astype(np.float64)

    n_events = int(np.max(event_id)) + 1
    if not (1500 <= n_events <= 3000):
        raise RuntimeError(
            f"n_events={n_events} outside [1500,3000] for {path} (n_hits_total={n_hits_total})"
        )

    end_mask = (global_id >= 0) & (global_id < 16)
    ev_end  = event_id[end_mask]
    t_end   = time_ns[end_mask]
    del arrays, event_id, global_id, time_ns

    t_min = np.full(n_events, np.inf)
    np.minimum.at(t_min, ev_end, t_end)
    invalid = ~np.isfinite(t_min)
    if invalid.any():
        t_min[invalid] = np.nan

    t_min_hit = t_min[ev_end]
    valid = np.isfinite(t_min_hit)
    t_rel = t_end[valid] - t_min_hit[valid]
    t_rel = t_rel[(t_rel >= 0.0) & (t_rel < 2.0)]

    del ev_end, t_end, t_min_hit, valid
    return t_rel, n_hits_total, n_events, t_min


def _build_th1(t_rel: np.ndarray, bin_ps: float, name: str, title: str) -> ROOT.TH1D:
    bin_ns = bin_ps * 1e-3
    n_bins = int(round(2.0 / bin_ns))
    h = ROOT.TH1D(name, title, n_bins, 0.0, 2.0)
    h.SetDirectory(ROOT.nullptr)
    h.GetXaxis().SetTitle("t_{foton} - t^{ev}_{min} [ns]")
    h.GetYaxis().SetTitle(f"Hits / ({bin_ps:.1f} ps)")
    for t in t_rel:
        h.Fill(float(t))
    return h


def _read_fit_params_2(f) -> dict:
    chi2 = float(f.GetChisquare())
    ndf  = int(f.GetNDF())
    return {
        "N":          float(f.GetParameter(0)), "N_err":           float(f.GetParError(0)),
        "tau_rise_ns": float(f.GetParameter(1)), "tau_rise_err_ns": float(f.GetParError(1)),
        "A_fast":      float(f.GetParameter(2)), "A_fast_err":      float(f.GetParError(2)),
        "tau_fast_ns": float(f.GetParameter(3)), "tau_fast_err_ns": float(f.GetParError(3)),
        "A_slow":      float(f.GetParameter(4)), "A_slow_err":      float(f.GetParError(4)),
        "tau_slow_ns": float(f.GetParameter(5)), "tau_slow_err_ns": float(f.GetParError(5)),
        "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
    }


def _read_fit_params_1(f) -> dict:
    chi2 = float(f.GetChisquare())
    ndf  = int(f.GetNDF())
    return {
        "N":          float(f.GetParameter(0)), "N_err":           float(f.GetParError(0)),
        "tau_rise_ns": float(f.GetParameter(1)), "tau_rise_err_ns": float(f.GetParError(1)),
        "A_fast":      float(f.GetParameter(2)), "A_fast_err":      float(f.GetParError(2)),
        "tau_fast_ns": float(f.GetParameter(3)), "tau_fast_err_ns": float(f.GetParError(3)),
        "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
    }


def _fit_scint_models(h: ROOT.TH1D, name: str) -> tuple[ROOT.TF1, ROOT.TF1, dict, dict, dict]:
    """
    Fit double and single exponential decay models via TH1::Fit (Poisson chi2).
    Returns (f2, f1, params_2, params_1, significance_test_dict).
    Caller must use `del` (not .Delete()) on returned TF1 objects.
    """
    h_max = float(h.GetMaximum())

    # ── Double-component model ────────────────────────────────────────────────
    f2 = ROOT.TF1(f"f2scint_{name}", SCINT_2_FORMULA, FIT_LO, FIT_HI)
    f2.SetParNames("N", "tau_rise", "A_fast", "tau_fast", "A_slow", "tau_slow")
    f2.SetParameters(h_max, 0.7, 1.0, 1.8, 0.3, 7.0)
    f2.SetParLimits(0, 0.0, 1e12)
    f2.SetParLimits(1, 0.2, 1.5)
    f2.SetParLimits(2, 0.0, 1e10)
    f2.SetParLimits(3, 1.0, 3.0)
    f2.SetParLimits(4, 0.0, 1e10)
    f2.SetParLimits(5, 3.0, 20.0)
    r2 = h.Fit(f2, "QRNS", "", FIT_LO, FIT_HI)
    status2 = r2.Status() if r2 else 99
    params_2 = {**_read_fit_params_2(f2), "status": status2,
                "model": "N*(1-exp(-t/tau_rise))*(A_fast*exp(-t/tau_fast)+A_slow*exp(-t/tau_slow))"}

    # ── Single-component model ────────────────────────────────────────────────
    f1 = ROOT.TF1(f"f1scint_{name}", SCINT_1_FORMULA, FIT_LO, FIT_HI)
    f1.SetParNames("N", "tau_rise", "A_fast", "tau_fast")
    f1.SetParameters(h_max, 0.7, 1.0, 1.8)
    f1.SetParLimits(0, 0.0, 1e12)
    f1.SetParLimits(1, 0.2, 1.5)
    f1.SetParLimits(2, 0.0, 1e10)
    f1.SetParLimits(3, 1.0, 3.0)
    r1 = h.Fit(f1, "QRNS+", "", FIT_LO, FIT_HI)
    status1 = r1.Status() if r1 else 99
    params_1 = {**_read_fit_params_1(f1), "status": status1,
                "model": "N*(1-exp(-t/tau_rise))*A_fast*exp(-t/tau_fast)"}

    # ── Significance of slow component (delta-chi2 test) ─────────────────────
    chi2_2 = params_2.get("chi2", float("inf"))
    ndf_2  = params_2.get("ndf", 0)
    chi2_1 = params_1.get("chi2", float("inf"))
    ndf_1  = params_1.get("ndf", 0)
    delta_chi2 = chi2_1 - chi2_2
    delta_ndf  = ndf_1 - ndf_2
    A_slow       = params_2.get("A_slow", 0.0)
    A_slow_err   = params_2.get("A_slow_err", float("inf"))
    A_slow_signif = A_slow / A_slow_err if A_slow_err > 0 else 0.0

    second_resolved = (
        status2 == 0
        and delta_chi2 > max(4.0 * delta_ndf, 10.0)
        and A_slow_signif > 2.0
    )
    significance = {
        "delta_chi2": delta_chi2,
        "delta_ndf": delta_ndf,
        "A_slow_significance_sigma": A_slow_signif,
        "second_component_resolved": second_resolved,
        "note": (
            "Segundo modo resuelto" if second_resolved else
            "Segunda componente NO resuelta a esta posicion "
            "(propagacion o escasez de fotones enmascara el modo lento)"
        ),
    }

    return f2, f1, params_2, params_1, significance


def generate_f1_sidecar(
    *,
    fig_id: str,
    material: str,
    dataset: str,
    root_input: pathlib.Path,
    sha256_expected: str,
    x_mm: int,
    output_dir: pathlib.Path,
    bin_ps: float = 2.0,
    optica: str = "sslg4",
) -> None:
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    root_input = pathlib.Path(root_input)
    stem = output_dir / fig_id

    for ext in ('.root', '.csv', '.meta.json'):
        p = stem.with_suffix(ext)
        if p.exists(): p.unlink()

    verify_root_input(root_input)
    print(f"[{fig_id}] SHA verify: {root_input.name} ...", end=' ', flush=True)
    verify_sha(root_input, sha256_expected, fig_id)
    print("OK")

    print(f"[{fig_id}] Loading End hits (x={x_mm}mm) ...", flush=True)
    t_rel, n_hits, n_events, t_min_ev = compute_t_rel_end(root_input)
    print(f"[{fig_id}] n_events={n_events}, n_hits={n_hits}, n_trel_0-2ns={len(t_rel)}")

    alt_ps = 2.5 if bin_ps == 2.0 else 2.0
    h_main = _build_th1(t_rel, bin_ps, f"h_trel_{fig_id}",
                         f"F1 {material} x={x_mm}mm")
    h_alt  = _build_th1(t_rel, alt_ps, f"h_trel_{fig_id}_alt",
                         f"F1 {material} x={x_mm}mm {alt_ps:.1f}ps")

    print(f"[{fig_id}] Fitting scintillation models ...", flush=True)
    f2, f1, params_2, params_1, sig = _fit_scint_models(h_main, fig_id)

    print(f"[{fig_id}] chi2/ndf: 2-mode={params_2['chi2_ndf']:.2f}(st={params_2['status']}) "
          f"1-mode={params_1['chi2_ndf']:.2f}(st={params_1['status']}) "
          f"slow_resolved={sig['second_component_resolved']}")

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")
    h_main.Write(); h_alt.Write(); f2.Write(); f1.Write()
    ROOT.SetOwnership(f2, False); ROOT.SetOwnership(f1, False)
    tf.Close()

    # ── CSV ──────────────────────────────────────────────────────────────────
    n_bins = h_main.GetNbinsX()
    write_csv_arrays(stem.with_suffix('.csv'), {
        'bin_center_ns': [h_main.GetBinCenter(i+1) for i in range(n_bins)],
        'count':         [h_main.GetBinContent(i+1) for i in range(n_bins)],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    now = datetime.datetime.now(datetime.UTC).isoformat()
    meta = {
        "fig_id": fig_id, "figura": "F1", "material": material,
        "dataset": dataset, "optica": optica, "x_mm": x_mm,
        "root_input": str(root_input), "sha256_input": sha256_expected,
        "binning": f"{bin_ps} ps primary, {alt_ps} ps variant; fit range [{FIT_LO},{FIT_HI}] ns",
        "escala": "Y logaritmico (figura principal); Y lineal (h_trel_alt tambien en .root)",
        "n_hits_total": int(n_hits), "n_events": int(n_events),
        "n_trel_hits_0_2ns": int(len(t_rel)),
        "t_min_per_event_median_ns": float(np.nanmedian(t_min_ev)),
        "t_min_per_event_p5_ns":     float(np.nanpercentile(t_min_ev, 5)),
        "t_min_per_event_p95_ns":    float(np.nanpercentile(t_min_ev, 95)),
        "fit_double_decay": params_2,
        "fit_single_decay": params_1,
        "slow_component_significance": sig,
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "readout_jitter_quadrature_ps": READOUT_JITTER_QUADRATURE_PS,
        "jitter_note": JITTER_DOUBLE_COUNT_RISK,
        "caption_label": (
            f"intrinseco (Etapa 1): sin time-walk/ToT/SPTR. "
            f"t_rel = t_foton - t_min^evento (incl. jitter 20ps por hit). "
            f"Bin {bin_ps} ps; escala Y log. "
            + (f"Dos modos resueltos: tau_fast={params_2.get('tau_fast_ns',0):.2f}ns, "
               f"tau_slow={params_2.get('tau_slow_ns',0):.2f}ns."
               if sig['second_component_resolved']
               else "Segunda componente NO resuelta (propagacion domina en esta posicion).")
        ),
        "comando": f"python3.12 f1_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    print(f"[{fig_id}] Done → {root_out.name}")

    del h_main, h_alt, t_rel, t_min_ev


if __name__ == "__main__":
    generate_f1_sidecar(
        fig_id="fig01a", material="EJ-204", dataset="endonly_mylar_20260614",
        root_input=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614/photon_hits_x0mm.root"),
        sha256_expected="63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
        x_mm=0, output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
