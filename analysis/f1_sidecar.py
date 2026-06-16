#!/usr/bin/env python3.12
"""
F1 sidecar generator — EXEC_14 (QA-1b update).

Histograma t_rel = t_fotón − min(t del evento) para SiPMs End.
Dos ventanas:
  - 0–2 ns, bin 2.0 ps (principal, detail del rise/pico)
  - 0–10 ns, bin 10 ps (extendida, modo lento τ~7 ns visible)

Modelo de dos modos de emisión (rise común + τ_fast + τ_slow) ajustado sobre
[0.01, 10] ns (extendida); el mismo modelo en [0.004, 2] ns (corta, comparación).
Test de significancia via delta-chi2 / A_slow/err.
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

FIT_LO_SHORT = 0.004   # ns — skip t=0 spike (2 bins @ 2ps)
FIT_HI_SHORT = 2.0
FIT_LO_EXT   = 0.01    # ns — skip first bin (1 bin @ 10ps)
FIT_HI_EXT   = 10.0


def compute_t_rel_end(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, int, int, np.ndarray]:
    """
    Load End SiPM hits, compute t_rel = time_ns - min(time_ns per event).
    Returns: (t_rel_2ns, t_rel_10ns, n_hits_total, n_events, t_min_per_event)
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
    valid     = np.isfinite(t_min_hit)
    t_rel     = t_end[valid] - t_min_hit[valid]
    t_rel     = t_rel[t_rel >= 0.0]
    del ev_end, t_end, t_min_hit, valid

    return (
        t_rel[t_rel < 2.0],
        t_rel[t_rel < 10.0],
        n_hits_total, n_events, t_min,
    )


def _build_th1(t_rel: np.ndarray, bin_ps: float, hi_ns: float,
               name: str, title: str) -> ROOT.TH1D:
    bin_ns = bin_ps * 1e-3
    n_bins = int(round(hi_ns / bin_ns))
    h = ROOT.TH1D(name, title, n_bins, 0.0, hi_ns)
    h.SetDirectory(ROOT.nullptr)
    h.GetXaxis().SetTitle("t_{foton} - t^{ev}_{min} [ns]")
    h.GetYaxis().SetTitle(f"Hits / ({bin_ps:.1f} ps)")
    for t in t_rel:
        h.Fill(float(t))
    return h


def _read_fit_params_2(f, status: int) -> dict:
    chi2 = float(f.GetChisquare())
    ndf  = int(f.GetNDF())
    return {
        "N": float(f.GetParameter(0)), "N_err": float(f.GetParError(0)),
        "tau_rise_ns": float(f.GetParameter(1)), "tau_rise_err_ns": float(f.GetParError(1)),
        "A_fast": float(f.GetParameter(2)), "A_fast_err": float(f.GetParError(2)),
        "tau_fast_ns": float(f.GetParameter(3)), "tau_fast_err_ns": float(f.GetParError(3)),
        "A_slow": float(f.GetParameter(4)), "A_slow_err": float(f.GetParError(4)),
        "tau_slow_ns": float(f.GetParameter(5)), "tau_slow_err_ns": float(f.GetParError(5)),
        "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
        "status": status,
        "model": "N*(1-exp(-t/tau_rise))*(A_fast*exp(-t/tau_fast)+A_slow*exp(-t/tau_slow))",
    }


def _read_fit_params_1(f, status: int) -> dict:
    chi2 = float(f.GetChisquare())
    ndf  = int(f.GetNDF())
    return {
        "N": float(f.GetParameter(0)), "N_err": float(f.GetParError(0)),
        "tau_rise_ns": float(f.GetParameter(1)), "tau_rise_err_ns": float(f.GetParError(1)),
        "A_fast": float(f.GetParameter(2)), "A_fast_err": float(f.GetParError(2)),
        "tau_fast_ns": float(f.GetParameter(3)), "tau_fast_err_ns": float(f.GetParError(3)),
        "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
        "status": status,
        "model": "N*(1-exp(-t/tau_rise))*A_fast*exp(-t/tau_fast)",
    }


def _significance(params_2: dict, params_1: dict) -> dict:
    chi2_2 = params_2.get("chi2", float("inf"))
    chi2_1 = params_1.get("chi2", float("inf"))
    delta_chi2 = chi2_1 - chi2_2
    delta_ndf  = params_1.get("ndf", 0) - params_2.get("ndf", 0)
    A_slow      = params_2.get("A_slow", 0.0)
    A_slow_err  = params_2.get("A_slow_err", float("inf"))
    A_slow_sig  = A_slow / A_slow_err if A_slow_err > 0 else 0.0
    second_resolved = (
        params_2.get("status", 99) == 0
        and delta_chi2 > max(4.0 * max(delta_ndf, 1), 10.0)
        and A_slow_sig > 2.0
    )
    return {
        "delta_chi2": delta_chi2, "delta_ndf": delta_ndf,
        "A_slow_significance_sigma": A_slow_sig,
        "second_component_resolved": second_resolved,
        "note": (
            "Segundo modo resuelto" if second_resolved else
            "Segunda componente NO resuelta (propagacion/estadistica enmascarada)"
        ),
    }


def _fit_scint_models(h: ROOT.TH1D, name: str, fit_lo: float, fit_hi: float
                      ) -> tuple[ROOT.TF1, ROOT.TF1, dict, dict, dict]:
    """
    Fit double + single exponential decay via TH1::Fit (Poisson chi2).
    Caller must `del` (not .Delete()) the returned TF1 objects.
    """
    h_max = float(h.GetMaximum())

    f2 = ROOT.TF1(f"f2_{name}", SCINT_2_FORMULA, fit_lo, fit_hi)
    f2.SetParNames("N", "tau_rise", "A_fast", "tau_fast", "A_slow", "tau_slow")
    f2.SetParameters(h_max, 0.7, 1.0, 1.8, 0.3, 7.0)
    f2.SetParLimits(0, 0.0, 1e12)
    f2.SetParLimits(1, 0.2, 1.5)
    f2.SetParLimits(2, 0.0, 1e10)
    f2.SetParLimits(3, 1.0, 3.0)
    f2.SetParLimits(4, 0.0, 1e10)
    f2.SetParLimits(5, 3.0, 20.0)
    r2 = h.Fit(f2, "QRNS", "", fit_lo, fit_hi)
    params_2 = _read_fit_params_2(f2, r2.Status() if r2 else 99)

    f1 = ROOT.TF1(f"f1_{name}", SCINT_1_FORMULA, fit_lo, fit_hi)
    f1.SetParNames("N", "tau_rise", "A_fast", "tau_fast")
    f1.SetParameters(h_max, 0.7, 1.0, 1.8)
    f1.SetParLimits(0, 0.0, 1e12)
    f1.SetParLimits(1, 0.2, 1.5)
    f1.SetParLimits(2, 0.0, 1e10)
    f1.SetParLimits(3, 1.0, 3.0)
    r1 = h.Fit(f1, "QRNS+", "", fit_lo, fit_hi)
    params_1 = _read_fit_params_1(f1, r1.Status() if r1 else 99)

    sig = _significance(params_2, params_1)
    return f2, f1, params_2, params_1, sig


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
    t_rel_2, t_rel_10, n_hits, n_events, t_min_ev = compute_t_rel_end(root_input)
    print(f"[{fig_id}] n_events={n_events}, n_hits={n_hits}, "
          f"n(0-2ns)={len(t_rel_2)}, n(0-10ns)={len(t_rel_10)}")

    alt_ps = 2.5 if bin_ps == 2.0 else 2.0
    h_short = _build_th1(t_rel_2,  bin_ps,  2.0,  f"h_trel_{fig_id}",
                          f"F1 {material} x={x_mm}mm 0-2ns {bin_ps}ps/bin")
    h_alt   = _build_th1(t_rel_2,  alt_ps,  2.0,  f"h_trel_{fig_id}_alt",
                          f"F1 {material} x={x_mm}mm 0-2ns {alt_ps}ps/bin")
    h_ext   = _build_th1(t_rel_10, 10.0,   10.0,  f"h_trel_{fig_id}_ext",
                          f"F1 {material} x={x_mm}mm 0-10ns 10ps/bin")

    print(f"[{fig_id}] Fitting 0-2ns ...", flush=True)
    f2s, f1s, p2s, p1s, sig_s = _fit_scint_models(
        h_short, f"{fig_id}_short", FIT_LO_SHORT, FIT_HI_SHORT)
    print(f"[{fig_id}] 0-2ns: 2-mode chi2/ndf={p2s['chi2_ndf']:.2f}(st={p2s['status']}) "
          f"1-mode={p1s['chi2_ndf']:.2f}(st={p1s['status']}) "
          f"slow={sig_s['second_component_resolved']}")

    print(f"[{fig_id}] Fitting 0-10ns ...", flush=True)
    f2e, f1e, p2e, p1e, sig_e = _fit_scint_models(
        h_ext, f"{fig_id}_ext", FIT_LO_EXT, FIT_HI_EXT)
    print(f"[{fig_id}] 0-10ns: 2-mode chi2/ndf={p2e['chi2_ndf']:.2f}(st={p2e['status']}) "
          f"1-mode={p1e['chi2_ndf']:.2f}(st={p1e['status']}) "
          f"tau_slow={p2e['tau_slow_ns']:.2f}±{p2e['tau_slow_err_ns']:.2f}ns "
          f"A_slow_sig={sig_e['A_slow_significance_sigma']:.1f}σ "
          f"slow={sig_e['second_component_resolved']}")

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")
    h_short.Write(); h_alt.Write(); h_ext.Write()
    f2s.Write(); f1s.Write(); f2e.Write(); f1e.Write()
    ROOT.SetOwnership(f2s, False); ROOT.SetOwnership(f1s, False)
    ROOT.SetOwnership(f2e, False); ROOT.SetOwnership(f1e, False)
    tf.Close()

    # ── CSV (short window) ────────────────────────────────────────────────────
    n_bins = h_short.GetNbinsX()
    write_csv_arrays(stem.with_suffix('.csv'), {
        'bin_center_ns': [h_short.GetBinCenter(i+1) for i in range(n_bins)],
        'count':         [h_short.GetBinContent(i+1) for i in range(n_bins)],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    now = datetime.datetime.now(datetime.UTC).isoformat()
    slow_note = (
        f"Dos modos resueltos en 0-10ns: "
        f"tau_fast={p2e.get('tau_fast_ns',0):.2f}ns, "
        f"tau_slow={p2e.get('tau_slow_ns',0):.2f}±{p2e.get('tau_slow_err_ns',0):.2f}ns."
        if sig_e['second_component_resolved']
        else "Segunda componente NO resuelta incluso en 0-10ns."
    )
    meta = {
        "fig_id": fig_id, "figura": "F1", "material": material,
        "dataset": dataset, "optica": optica, "x_mm": x_mm,
        "root_input": str(root_input), "sha256_input": sha256_expected,
        "n_hits_total": int(n_hits), "n_events": int(n_events),
        "n_trel_0_2ns":  int(len(t_rel_2)),
        "n_trel_0_10ns": int(len(t_rel_10)),
        "t_min_per_event_median_ns": float(np.nanmedian(t_min_ev)),
        "t_min_per_event_p5_ns":     float(np.nanpercentile(t_min_ev, 5)),
        "t_min_per_event_p95_ns":    float(np.nanpercentile(t_min_ev, 95)),
        "fit_double_decay_2ns":  p2s,
        "fit_single_decay_2ns":  p1s,
        "slow_significance_2ns": sig_s,
        "fit_double_decay_10ns": p2e,
        "fit_single_decay_10ns": p1e,
        "slow_significance_10ns": sig_e,
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "readout_jitter_quadrature_ps": READOUT_JITTER_QUADRATURE_PS,
        "jitter_note": JITTER_DOUBLE_COUNT_RISK,
        "caption_label": (
            f"Intrinseco (Etapa 1): t_rel = t_foton - t_min^ev. "
            f"Panel izq: 0-2ns bin {bin_ps}ps (rise/pico). "
            f"Panel der: 0-10ns bin 10ps. {slow_note}"
        ),
        "comando": f"python3.12 f1_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    print(f"[{fig_id}] Done → {root_out.name}")

    del h_short, h_alt, h_ext, t_rel_2, t_rel_10, t_min_ev


if __name__ == "__main__":
    generate_f1_sidecar(
        fig_id="fig01a", material="EJ-204", dataset="endonly_mylar_20260614",
        root_input=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614/photon_hits_x0mm.root"),
        sha256_expected="63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
        x_mm=0, output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
