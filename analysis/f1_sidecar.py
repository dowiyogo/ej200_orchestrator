#!/usr/bin/env python3.12
"""
F1 sidecar generator — EXEC_14.

Histograma de tiempo relativo al primer fotón del evento.
0–2 ns, bin 2.0 ps (variante 2.5 ps), eje Y log.
Ajuste de dos gaussianas (componente rápida + lenta) via TH1::Fit.
"""
from __future__ import annotations

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


def compute_t_rel_end(path: pathlib.Path) -> tuple[np.ndarray, int, int, np.ndarray]:
    """
    Load End SiPM hits, compute t_rel = time_ns - min(time_ns per event).
    Returns: (t_rel_array, n_hits_total, n_events, t_min_per_event)
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
            f"F1: n_events={n_events} outside [1500,3000] for {path}. "
            f"n_hits_total={n_hits_total} — check hits vs events distinction."
        )

    # Filter to End SiPMs (global_id 0-15)
    end_mask = (global_id >= 0) & (global_id < 16)
    ev_end  = event_id[end_mask]
    t_end   = time_ns[end_mask]

    # Per-event minimum time across End SiPMs
    t_min = np.full(n_events, np.inf)
    np.minimum.at(t_min, ev_end, t_end)
    invalid = ~np.isfinite(t_min)
    if invalid.any():
        t_min[invalid] = np.nan

    # Relative time for each End hit
    t_min_hit = t_min[ev_end]
    valid_hit = np.isfinite(t_min_hit)
    t_rel = t_end[valid_hit] - t_min_hit[valid_hit]

    # Keep only 0–2 ns range (guard against negative values from jitter)
    t_rel = t_rel[(t_rel >= 0.0) & (t_rel < 2.0)]

    del arrays, event_id, global_id, time_ns, end_mask, ev_end, t_end, t_min_hit, valid_hit
    return t_rel, n_hits_total, n_events, t_min


def _build_th1(t_rel: np.ndarray, bin_ps: float, name: str, title: str) -> ROOT.TH1D:
    bin_ns = bin_ps * 1e-3
    n_bins = int(round(2.0 / bin_ns))
    h = ROOT.TH1D(name, title, n_bins, 0.0, 2.0)
    h.SetDirectory(ROOT.nullptr)
    h.GetXaxis().SetTitle("t_{foton} - t^{evento}_{min} [ns]")
    h.GetYaxis().SetTitle(f"Hits / ({bin_ps:.1f} ps)")
    for t in t_rel:
        h.Fill(float(t))
    return h


def _fit_two_gaussians(h: ROOT.TH1D, name: str) -> tuple[ROOT.TF1, dict]:
    """
    Fit two Gaussians to time histogram (fast + slow scintillation modes).
    Seeds the fast component from the peak, slow from a later region.
    Returns (TF1 object, dict of fit parameters).
    """
    # Seed from histogram
    peak_bin = h.GetMaximumBin()
    peak_x   = h.GetBinCenter(peak_bin)
    peak_amp = h.GetBinContent(peak_bin)

    # Estimate slow region: look for secondary structure or use 1.0 ns
    # Fast component is near the peak; slow is typically at larger t_rel
    f2g = ROOT.TF1(f"f2g_{name}", "gaus(0)+gaus(3)", 0.0, 2.0)
    f2g.SetParNames("A_fast", "mu_fast_ns", "sig_fast_ns",
                    "A_slow", "mu_slow_ns", "sig_slow_ns")

    # Initial seeds
    f2g.SetParameters(peak_amp, peak_x, 0.05,
                      peak_amp * 0.2, min(peak_x + 0.5, 1.5), 0.3)
    # Bounds to keep fast near peak and slow later
    f2g.SetParLimits(0, 0.0, 1e12)
    f2g.SetParLimits(1, 0.0, 0.8)        # fast mean 0-0.8 ns
    f2g.SetParLimits(2, 0.005, 0.5)      # fast sigma 5-500 ps
    f2g.SetParLimits(3, 0.0, 1e12)
    f2g.SetParLimits(4, 0.1, 2.0)        # slow mean
    f2g.SetParLimits(5, 0.05, 1.0)       # slow sigma

    result2g = h.Fit(f2g, "QRNS", "", 0.0, 2.0)
    status2g = result2g.Status() if result2g else 99

    params2g = {
        "A_fast":     float(f2g.GetParameter(0)),
        "mu_fast_ns": float(f2g.GetParameter(1)),
        "sig_fast_ns":float(f2g.GetParameter(2)),
        "A_slow":     float(f2g.GetParameter(3)),
        "mu_slow_ns": float(f2g.GetParameter(4)),
        "sig_slow_ns":float(f2g.GetParameter(5)),
        "chi2":       float(f2g.GetChisquare()),
        "ndf":        int(f2g.GetNDF()),
        "chi2_ndf":   float(f2g.GetChisquare() / f2g.GetNDF()) if f2g.GetNDF() > 0 else -1.0,
        "status":     status2g,
    }
    return f2g, params2g


def _fit_one_gaussian(h: ROOT.TH1D, name: str) -> tuple[ROOT.TF1, dict]:
    """Single Gaussian fit for comparison (stored in sidecar but not the primary fit)."""
    peak_bin = h.GetMaximumBin()
    peak_x   = h.GetBinCenter(peak_bin)
    peak_amp = h.GetBinContent(peak_bin)

    f1g = ROOT.TF1(f"f1g_{name}", "gaus", 0.0, 2.0)
    f1g.SetParameters(peak_amp, peak_x, 0.3)
    result1g = h.Fit(f1g, "QRNS+", "", 0.0, 2.0)  # "+" = add to existing
    status1g = result1g.Status() if result1g else 99

    params1g = {
        "amplitude":  float(f1g.GetParameter(0)),
        "mean_ns":    float(f1g.GetParameter(1)),
        "sigma_ns":   float(f1g.GetParameter(2)),
        "chi2":       float(f1g.GetChisquare()),
        "ndf":        int(f1g.GetNDF()),
        "chi2_ndf":   float(f1g.GetChisquare() / f1g.GetNDF()) if f1g.GetNDF() > 0 else -1.0,
        "status":     status1g,
    }
    return f1g, params1g


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
    """Generate fig01X sidecar: .root + .csv + .meta.json"""
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    root_input = pathlib.Path(root_input)
    stem = output_dir / fig_id

    # Clean stale outputs
    for ext in ('.root', '.csv', '.meta.json'):
        p = stem.with_suffix(ext)
        if p.exists():
            p.unlink()

    # Verify input
    verify_root_input(root_input)
    print(f"[{fig_id}] SHA-256 verify: {root_input.name} ...", end=' ', flush=True)
    verify_sha(root_input, sha256_expected, fig_id)
    print("OK")

    # Load data
    print(f"[{fig_id}] Loading End SiPM hits (x={x_mm}mm) ...", flush=True)
    t_rel, n_hits_total, n_events, t_min_per_event = compute_t_rel_end(root_input)
    print(f"[{fig_id}] n_events={n_events}, n_hits_total={n_hits_total}, "
          f"n_rel_hits (0-2ns)={len(t_rel)}")

    t_min_median = float(np.nanmedian(t_min_per_event))
    t_min_p5     = float(np.nanpercentile(t_min_per_event, 5))
    t_min_p95    = float(np.nanpercentile(t_min_per_event, 95))

    # Build histograms (primary bin_ps and variant)
    h_main = _build_th1(t_rel, bin_ps,
                         f"h_trel_{fig_id}",
                         f"F1 {material} x={x_mm}mm;t_{{foton}}-t_{{min}}^{{ev}} [ns];Hits/({bin_ps:.1f}ps)")

    alt_ps  = 2.5 if bin_ps == 2.0 else 2.0
    h_alt   = _build_th1(t_rel, alt_ps,
                          f"h_trel_{fig_id}_alt",
                          f"F1 {material} x={x_mm}mm variant {alt_ps:.1f}ps")

    # Fits on primary histogram
    f2g, params_2g = _fit_two_gaussians(h_main, fig_id)
    f1g, params_1g = _fit_one_gaussian(h_main, fig_id)

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")
    h_main.Write(h_main.GetName())
    h_alt.Write(h_alt.GetName())
    f2g.Write(f2g.GetName())
    f1g.Write(f1g.GetName())
    tf.Close()

    # ── CSV sidecar ──────────────────────────────────────────────────────────
    n_bins_main = h_main.GetNbinsX()
    bin_centers = [h_main.GetBinCenter(i + 1) for i in range(n_bins_main)]
    bin_counts  = [h_main.GetBinContent(i + 1) for i in range(n_bins_main)]
    write_csv_arrays(stem.with_suffix('.csv'), {
        'bin_center_ns': bin_centers,
        'count': bin_counts,
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    meta = {
        "fig_id":     fig_id,
        "figura":     "F1",
        "material":   material,
        "dataset":    dataset,
        "optica":     optica,
        "x_mm":       x_mm,
        "root_input": str(root_input),
        "sha256_input": sha256_expected,
        "binning":    f"{bin_ps} ps primary, {alt_ps} ps variant",
        "escala":     "Y logaritmico; X lineal 0-2 ns",
        "n_hits_total": int(n_hits_total),
        "n_events":   int(n_events),
        "t_min_per_event_median_ns": t_min_median,
        "t_min_per_event_p5_ns":    t_min_p5,
        "t_min_per_event_p95_ns":   t_min_p95,
        "fit_two_gaussians": params_2g,
        "fit_single_gaussian_comparison": params_1g,
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "readout_jitter_quadrature_ps": READOUT_JITTER_QUADRATURE_PS,
        "jitter_note": JITTER_DOUBLE_COUNT_RISK,
        "caption_label": (
            "intrinseco (Etapa 1): sin time-walk/ToT/SPTR. "
            f"t_foton = GetGlobalTime() + jitter(20ps); "
            f"t_rel = t_foton - t_min^evento. Bin {bin_ps} ps."
        ),
        "comando": f"python3.12 f1_sidecar.py fig_id={fig_id} x_mm={x_mm}",
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)

    print(f"[{fig_id}] Sidecar written: {root_out.name}, {stem.name}.csv, {stem.name}.meta.json")

    # Let Python ref-counting call destructors exactly once (never .Delete() on Python-owned ROOT objects)
    del h_main, h_alt, f2g, f1g


if __name__ == "__main__":
    # Proof run: fig01a — EJ-204, x=0
    generate_f1_sidecar(
        fig_id="fig01a",
        material="EJ-204",
        dataset="endonly_mylar_20260614",
        root_input=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614/photon_hits_x0mm.root"
        ),
        sha256_expected="63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
        x_mm=0,
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
        bin_ps=2.0,
    )
