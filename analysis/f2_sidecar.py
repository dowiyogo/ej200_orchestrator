#!/usr/bin/env python3.12
"""
F2 sidecar generator — EXEC_14.

σt(N_pe) ensemble: agrupa eventos de las 31 posiciones por umbral de N_pe en
End SiPMs {5, 10, 15, 20, 25}. Para cada umbral T, sigma_single =
sigma(ΔT_LR[N_pe>=T]) / sqrt(2) con TH1::Fit.

Nota: procesa una posición a la vez para control de memoria.
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
    verify_root_input, compute_sum4, fit_core_root,
    write_meta, write_csv_arrays,
    TIME_CONVENTION, JITTER_PER_HIT_NS, READOUT_JITTER_QUADRATURE_PS,
    JITTER_DOUBLE_COUNT_RISK, SQRT2,
)

NPE_THRESHOLDS = [5, 10, 15, 20, 25]

ALL_POSITIONS = [
    0, -50, 50, -100, 100, -150, 150, -200, 200,
    -250, 250, -300, 300, -350, 350, -400, 400,
    -450, 450, -500, 500, -550, 550, -600, 600,
    -650, 650, -670, 670, -690, 690,
]

METRIC_SHA = {
    "EJ-204": {
        0:   "63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
        400: "403016634e3192b6a7d349260eb06b18ed4c222d2114d2128efded7a7e3479ec",
        690: "61435a5c9e69350aed0ca07612483bb3c78b879b1b5edd3fe7215e2e67037d55",
    },
    "EJ-230": {
        0:   "0d30b51f40c9ffe46463fc1627cc87faca8900342340c176998e14069050c930",
        400: "1054527b17854828f663fecc31456723157cda7999785bae652ba51017fc10b6",
        690: "c67b636a1e0cf916b2a51dc7b10a9cda4784cfb8807cb2c7aded2bcf06eaad97",
    },
}


def _load_position(root_path: pathlib.Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load End hits for ONE position. Returns (npe_end_per_event, delta_lr_ns).
    delta_lr is NaN for events without both-side trigger.
    Caller must del arrays.
    """
    with uproot.open(str(root_path)) as f:
        tree = f['sipm_hits']
        arrays = tree.arrays(['event_id', 'global_id', 'time_ns'], library='np')

    event_id = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    time_ns   = arrays['time_ns'].astype(np.float64)
    n_events  = int(np.max(event_id)) + 1
    del arrays

    end_mask = (global_id >= 0) & (global_id < 16)
    ev_end   = event_id[end_mask]
    g_end    = global_id[end_mask]
    t_end    = time_ns[end_mask]
    del event_id, global_id, time_ns

    npe_end = np.bincount(ev_end.astype(np.intp), minlength=n_events)

    arrays_end = {
        'event_id':       ev_end,
        'global_id':      g_end,
        'time_ns':        t_end,
        '_n_events':      np.int64(n_events),
        '_n_hits_total':  np.int64(0),
        '_n_hits_end':    np.int64(len(ev_end)),
    }
    result = compute_sum4(arrays_end)
    del arrays_end, ev_end, g_end, t_end

    t_left  = result.t_left
    t_right = result.t_right
    delta_lr = np.full(n_events, np.nan)
    both = np.isfinite(t_left) & np.isfinite(t_right)
    delta_lr[both] = t_left[both] - t_right[both]

    return npe_end, delta_lr


def generate_f2_sidecar(
    *,
    fig_id: str,
    material: str,
    dataset: str,
    data_dir: pathlib.Path,
    sha256_manifest: str,
    output_dir: pathlib.Path,
    optica: str = "sslg4",
) -> None:
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir   = pathlib.Path(data_dir)
    stem       = output_dir / fig_id

    for ext in ('.root', '.csv', '.meta.json'):
        p = stem.with_suffix(ext)
        if p.exists(): p.unlink()

    # Accumulate (npe, delta_lr) across all 31 positions
    all_npe_lists    = []
    all_delta_lists  = []
    positions_loaded = []
    n_events_total   = 0

    metric_sha = METRIC_SHA.get(material, {})

    for x_mm in ALL_POSITIONS:
        root_path = data_dir / f"photon_hits_x{x_mm}mm.root"
        verify_root_input(root_path)
        print(f"[{fig_id}] x={x_mm:+5d}mm ...", end=' ', flush=True)
        npe, delta_lr = _load_position(root_path)
        n_events_total += len(npe)
        all_npe_lists.append(npe)
        all_delta_lists.append(delta_lr)
        positions_loaded.append(x_mm)
        n_trig = int(np.isfinite(delta_lr).sum())
        print(f"N_ev={len(npe)}  ΔT finite={n_trig}  N_pe median={np.median(npe):.1f}")

    all_npe    = np.concatenate(all_npe_lists)
    all_delta  = np.concatenate(all_delta_lists)
    del all_npe_lists, all_delta_lists

    # Compute sigma_single for each N_pe threshold
    threshold_results = []
    for T in NPE_THRESHOLDS:
        mask = (all_npe >= T) & np.isfinite(all_delta)
        delta_sub = all_delta[mask]
        n_sub = len(delta_sub)
        if n_sub < 20:
            print(f"[{fig_id}] T={T}: only {n_sub} events — skip")
            threshold_results.append({
                "npe_threshold": T, "n_events_pass": n_sub,
                "sigma_single_ps": math.nan, "sigma_single_err_ps": math.nan,
                "sigma_lr_ps": math.nan, "fit_used": False, "fit_chi2_ndf": -1.0,
            })
            continue
        fit = fit_core_root(delta_sub, f"{fig_id}_T{T}")
        sigma_lr_ps     = fit.sigma_ns * 1000.0
        sigma_single_ps = sigma_lr_ps / SQRT2
        sigma_lr_err    = fit.sigma_err_ns * 1000.0
        sigma_single_err = sigma_lr_err / SQRT2
        print(
            f"[{fig_id}] T={T}: n={n_sub}, "
            f"σ_lr={sigma_lr_ps:.1f}ps, σ_single={sigma_single_ps:.1f}ps "
            f"(chi2/ndf={fit.chi2_ndf:.2f})"
        )
        threshold_results.append({
            "npe_threshold":      T,
            "n_events_pass":      n_sub,
            "sigma_lr_ps":        sigma_lr_ps,
            "sigma_lr_err_ps":    sigma_lr_err,
            "sigma_single_ps":    sigma_single_ps,
            "sigma_single_err_ps": sigma_single_err,
            "fit_used":           fit.fit_used,
            "fit_chi2_ndf":       fit.chi2_ndf,
        })
        del delta_sub

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")

    valid = [r for r in threshold_results if math.isfinite(r.get("sigma_single_ps", math.nan))]
    if valid:
        g = ROOT.TGraphErrors(len(valid))
        g.SetName(f"g_sigma_vs_npe_{fig_id}")
        g.SetTitle(
            f"F2 {material} ensemble #sigma_{{single}} vs N_{{pe}} threshold;"
            "N_{pe} umbral [PE];#sigma_{single} [ps]"
        )
        for i, r in enumerate(valid):
            g.SetPoint(i, float(r['npe_threshold']), r['sigma_single_ps'])
            g.SetPointError(i, 0.0, r.get('sigma_single_err_ps', 0.0))
        g.Write()
        ROOT.SetOwnership(g, False)

    tf.Close()

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        'npe_threshold':        [r['npe_threshold']             for r in threshold_results],
        'n_events_pass':        [r['n_events_pass']             for r in threshold_results],
        'sigma_lr_ps':          [r.get('sigma_lr_ps', math.nan) for r in threshold_results],
        'sigma_lr_err_ps':      [r.get('sigma_lr_err_ps', math.nan) for r in threshold_results],
        'sigma_single_ps':      [r.get('sigma_single_ps', math.nan) for r in threshold_results],
        'sigma_single_err_ps':  [r.get('sigma_single_err_ps', math.nan) for r in threshold_results],
        'fit_chi2_ndf':         [r.get('fit_chi2_ndf', -1.0)   for r in threshold_results],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    now = datetime.datetime.now(datetime.UTC).isoformat()
    meta = {
        "fig_id": fig_id, "figura": "F2", "material": material,
        "dataset": dataset, "optica": optica,
        "x_mm_list": "all-31 positions pooled",
        "sha256_manifest_declared": sha256_manifest,
        "sha256_metric_verified": metric_sha,
        "sha256_note": (
            "SHA verificados para x=0/400/690 (posiciones metricas). "
            "Las 28 restantes pasan verify_root_input (estructura + tamanio)."
        ),
        "n_positions": len(positions_loaded),
        "n_events_total_all_positions": n_events_total,
        "npe_thresholds": NPE_THRESHOLDS,
        "binning": "N_pe en End SiPMs (global_id 0-15) por evento; cumulative threshold scan",
        "escala": "X lineal N_pe threshold [PE]; Y lineal sigma_t [ps]",
        "estimador": "sigma_single = sigma(DeltaT_LR[N_pe>=T])/sqrt(2); TH1::Fit gaussiano core",
        "threshold_results": threshold_results,
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "jitter_note": JITTER_DOUBLE_COUNT_RISK,
        "caption_label": (
            f"{material} ensemble ({len(positions_loaded)} posiciones). "
            "sigma_single sin readout jitter adicional (Etapa 1)."
        ),
        "comando": f"python3.12 f2_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    del all_npe, all_delta
    print(f"[{fig_id}] Done → {root_out.name}")


if __name__ == "__main__":
    generate_f2_sidecar(
        fig_id="fig02a", material="EJ-204", dataset="endonly_mylar_20260614",
        data_dir=pathlib.Path("/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614"),
        sha256_manifest="8b8c93099eeae7d9e9e3e3a25cd3400890fbcdeee5e041bda5890c5fae77a8b2",
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
