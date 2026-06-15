#!/usr/bin/env python3.12
"""
F5 sidecar generator — EXEC_14.

σt SUM4 izquierda/derecha en x=0/400/690 mm, para EJ-204 y EJ-230.
Usa TH1::Fit (PyROOT) en lugar de curve_fit de endonly_sum4.py.
Reporta sigma nuevo Y sigma del cache existente para comparación de método.
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
    load_end_hits, compute_sum4, fit_core_root,
    write_meta, write_csv_arrays,
    TIME_CONVENTION, JITTER_PER_HIT_NS, READOUT_JITTER_QUADRATURE_PS,
    JITTER_DOUBLE_COUNT_RISK, SQRT2,
)

# Known sigma_single values from endonly_sum4.py cache (curve_fit method)
# Used only for comparison in meta.json; do NOT use as hardcoded values in the deck.
CACHE_SIGMA_PS = {
    "EJ-204": {0: 141.5, 400: 209.3, 690: 372.1},
    "EJ-230": {0: 140.2, 400: 243.9, 690: 434.0},
}

METRIC_POSITIONS = [0, 400, 690]


def _process_position(root_path: pathlib.Path, x_mm: int, fig_id: str
                       ) -> tuple[dict, ROOT.TH1D]:
    """
    Load hits for one position, compute SUM4 ΔT_LR, fit, return stats + TH1.
    Caller is responsible for deleting the returned TH1.
    """
    verify_root_input(root_path)

    arrays = load_end_hits(root_path)
    n_events = int(arrays['_n_events'])
    n_hits   = int(arrays['_n_hits_total'])

    result = compute_sum4(arrays)
    del arrays  # release memory

    delta_lr = result.delta_lr
    n_triggered = len(delta_lr)
    trigger_eff = n_triggered / n_events

    fit = fit_core_root(delta_lr, f"{fig_id}_x{x_mm}")

    sigma_lr_ps    = fit.sigma_ns * 1000.0
    sigma_lr_err   = fit.sigma_err_ns * 1000.0
    sigma_single_ps= sigma_lr_ps / SQRT2
    sigma_single_err= sigma_lr_err / SQRT2
    rms_lr_ps      = fit.rms_ns * 1000.0

    # Build TH1 of ΔT_LR for sidecar
    median = float(np.median(delta_lr)) if len(delta_lr) else 0.0
    mad    = max(1.4826 * float(np.median(np.abs(delta_lr - median))), 1e-4) if len(delta_lr) else 0.001
    lo, hi = median - 8.0 * mad, median + 8.0 * mad
    h = ROOT.TH1D(
        f"h_delta_{fig_id}_x{x_mm}",
        f"#DeltaT_{{LR}} {fig_id} x={x_mm}mm;#DeltaT_{{LR}} [ns];Events/bin",
        100, lo, hi,
    )
    h.SetDirectory(ROOT.nullptr)  # Python owns; safe to Delete() after write
    for v in delta_lr:
        h.Fill(float(v))

    stats = {
        "x_mm":            x_mm,
        "n_events":        n_events,
        "n_hits_total":    n_hits,
        "n_triggered":     n_triggered,
        "trigger_efficiency": trigger_eff,
        "mean_delta_lr_ns":fit.mean_ns,
        "sigma_lr_ps":     sigma_lr_ps,
        "sigma_lr_err_ps": sigma_lr_err,
        "sigma_lr_rms_ps": rms_lr_ps,
        "sigma_single_ps": sigma_single_ps,
        "sigma_single_err_ps": sigma_single_err,
        "fit_used":        fit.fit_used,
        "fit_chi2_ndf":    fit.chi2_ndf,
        "fit_n_values":    fit.n,
    }
    return stats, h


def generate_f5_sidecar(
    *,
    fig_id: str,
    material: str,
    dataset: str,
    data_dir: pathlib.Path,
    sha256_per_position: dict[int, str],
    output_dir: pathlib.Path,
    optica: str = "sslg4",
) -> None:
    """Generate fig05X sidecar: .root + .csv + .meta.json"""
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir   = pathlib.Path(data_dir)

    stem = output_dir / fig_id

    # Clean stale outputs
    for ext in ('.root', '.csv', '.meta.json'):
        p = stem.with_suffix(ext)
        if p.exists():
            p.unlink()

    all_stats = []
    th1_handles = []

    for x_mm in METRIC_POSITIONS:
        root_path = data_dir / f"photon_hits_x{x_mm}mm.root"
        expected_sha = sha256_per_position[x_mm]

        print(f"[{fig_id}] x={x_mm}mm: SHA verify ...", end=' ', flush=True)
        verify_sha(root_path, expected_sha, f"{fig_id}_x{x_mm}")
        print("OK")

        print(f"[{fig_id}] x={x_mm}mm: computing SUM4 ...", flush=True)
        stats, h = _process_position(root_path, x_mm, fig_id)
        all_stats.append(stats)
        th1_handles.append(h)

        cache_sigma = CACHE_SIGMA_PS.get(material, {}).get(x_mm)
        diff = (stats['sigma_single_ps'] - cache_sigma) if cache_sigma else None
        print(
            f"[{fig_id}] x={x_mm}mm: "
            f"σ_single(TH1::Fit)={stats['sigma_single_ps']:.2f} ps  "
            f"cache(curve_fit)={cache_sigma} ps  "
            f"Δ={f'{diff:+.2f}' if diff is not None else 'N/A'} ps"
        )

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")
    for h in th1_handles:
        h.Write(h.GetName())

    # TGraphErrors: sigma_single vs x
    g = ROOT.TGraphErrors(len(METRIC_POSITIONS))
    g.SetName(f"g_sigma_single_{fig_id}")
    g.SetTitle(
        f"F5 {material} SUM4 #sigma_{{single}} = #sigma(#DeltaT_{{LR}})/#sqrt{{2}};"
        "Position x [mm];#sigma_{{single}} [ps]"
    )
    for i, s in enumerate(all_stats):
        g.SetPoint(i, float(s['x_mm']), s['sigma_single_ps'])
        g.SetPointError(i, 0.0, s['sigma_single_err_ps'])
    g.Write()
    ROOT.SetOwnership(g, False)  # ROOT file owns g; Python must not GC it
    tf.Close()

    # ── CSV sidecar ──────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        'x_mm':              [s['x_mm']              for s in all_stats],
        'n_events':          [s['n_events']           for s in all_stats],
        'n_triggered':       [s['n_triggered']        for s in all_stats],
        'trigger_efficiency':[s['trigger_efficiency'] for s in all_stats],
        'sigma_lr_ps':       [s['sigma_lr_ps']        for s in all_stats],
        'sigma_lr_err_ps':   [s['sigma_lr_err_ps']    for s in all_stats],
        'sigma_single_ps':   [s['sigma_single_ps']    for s in all_stats],
        'sigma_single_err_ps':[s['sigma_single_err_ps'] for s in all_stats],
        'sigma_lr_rms_ps':   [s['sigma_lr_rms_ps']    for s in all_stats],
        'cache_sigma_single_ps': [
            CACHE_SIGMA_PS.get(material, {}).get(s['x_mm'], float('nan'))
            for s in all_stats
        ],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    cache_comparison = {
        str(x): {
            "sigma_single_ps_TH1Fit": s['sigma_single_ps'],
            "sigma_single_ps_cache_curvefit": CACHE_SIGMA_PS.get(material, {}).get(x),
            "delta_ps": (
                s['sigma_single_ps'] - CACHE_SIGMA_PS.get(material, {}).get(x)
                if CACHE_SIGMA_PS.get(material, {}).get(x) is not None
                else None
            ),
        }
        for x, s in zip(METRIC_POSITIONS, all_stats)
    }

    meta = {
        "fig_id":      fig_id,
        "figura":      "F5",
        "material":    material,
        "dataset":     dataset,
        "optica":      optica,
        "x_mm_list":   METRIC_POSITIONS,
        "root_inputs": {
            str(x): str(data_dir / f"photon_hits_x{x}mm.root")
            for x in METRIC_POSITIONS
        },
        "sha256_inputs": sha256_per_position,
        "binning":     "100 bins over median±8·MAD; 4-iter Gaussian core fit (TH1::Fit)",
        "escala":      "Y lineal; tabla de 3 puntos",
        "estimador":   "sigma_single = sigma(DeltaT_LR)/sqrt(2); congruente con endonly_sum4.py",
        "positions": all_stats,
        "cache_comparison_curvefit_vs_TH1Fit": cache_comparison,
        "fit_method_note": (
            "TH1::Fit usa pesos Poisson internos (a diferencia de curve_fit sin pesos "
            "en endonly_sum4.py). Diferencias esperadas <pocos ps para distribuciones "
            "bien gaussianas."
        ),
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "readout_jitter_quadrature_ps": READOUT_JITTER_QUADRATURE_PS,
        "jitter_note": JITTER_DOUBLE_COUNT_RISK,
        "time_convention": TIME_CONVENTION,
        "caption_label": (
            "intrinseco (Etapa 1): sin time-walk/ToT/SPTR. "
            "sigma_single = sigma(DeltaT_LR)/sqrt(2). "
            "Topologia: min() gana-el-primero entre dos clusters SUM4 por extremo (hipotesis)."
        ),
        "caveats": [
            "Topologia SUM4 = hipotesis gana-el-primero; pendiente confirmacion de Gerardo.",
            "sigma_single intrinsecos; readout jitter 20ps en cuadratura pendiente revision doble conteo.",
        ],
        "comando": f"python3.12 f5_sidecar.py fig_id={fig_id} material={material}",
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)

    print(f"[{fig_id}] Sidecar written: {root_out.name}, {stem.name}.csv, {stem.name}.meta.json")

    # h: Python-owned (SetDirectory nullptr), use del (never .Delete())
    for h in th1_handles:
        del h
    th1_handles.clear()
    # g: ROOT-owned after SetOwnership(g, False); just release Python ref


if __name__ == "__main__":
    # Proof run: fig05a — EJ-204, x=0/400/690
    generate_f5_sidecar(
        fig_id="fig05a",
        material="EJ-204",
        dataset="endonly_mylar_20260614",
        data_dir=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614"
        ),
        sha256_per_position={
            0:   "63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
            400: "403016634e3192b6a7d349260eb06b18ed4c222d2114d2128efded7a7e3479ec",
            690: "61435a5c9e69350aed0ca07612483bb3c78b879b1b5edd3fe7215e2e67037d55",
        },
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
