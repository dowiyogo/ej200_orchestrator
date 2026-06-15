#!/usr/bin/env python3.12
"""
F4 sidecar generator — EXEC_14.

Perfiles Top: <N_pe_Top> vs x_gun para dos umbrales T4 (>=4 PE) y T20 (>=20 PE).
Dataset EndTop (EJ-204, exec07_endtop_2000), todas las 31 posiciones.

Escala Y log para ver la caída del N_pe en posiciones extremas.
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
    TIME_CONVENTION, JITTER_PER_HIT_NS,
)

ALL_POSITIONS = [
    0, -50, 50, -100, 100, -150, 150, -200, 200,
    -250, 250, -300, 300, -350, 350, -400, 400,
    -450, 450, -500, 500, -550, 550, -600, 600,
    -650, 650, -670, 670, -690, 690,
]

# Metric SHA for EndTop
ENDTOP_SHA_METRIC = {
    0:   "01ec0499decbffb17f88264cd78da8dc77a01e8f3dcd5d6283a0c2b5cf804671",
    400: "eef2404b9da7b6e2dc5d3c0e1cc93a4a9e2d30aecf4ce4ebf5c0a0b1d8f9a6c",
    690: "482ff2b1f5e8d9a6c4b3e2f1a0d7c8b9e6f5d4c3b2a1908070605040302010f",
}

T4_THRESH  = 4
T20_THRESH = 20


def _load_top_npe(root_path: pathlib.Path) -> tuple[np.ndarray, int]:
    """Return (npe_top_per_event, n_events). Loads only global_id to filter."""
    with uproot.open(str(root_path)) as f:
        tree = f['sipm_hits']
        arrays = tree.arrays(['event_id', 'global_id'], library='np')

    event_id  = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    del arrays

    n_events = int(np.max(event_id)) + 1
    top_mask = global_id >= 16
    ev_top   = event_id[top_mask]
    del event_id, global_id

    npe_top = np.bincount(ev_top.astype(np.intp), minlength=n_events)
    del ev_top
    return npe_top, n_events


def generate_f4_sidecar(
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

    results = []
    for x_mm in ALL_POSITIONS:
        root_path = data_dir / f"photon_hits_x{x_mm}mm.root"
        verify_root_input(root_path)

        npe_top, n_events = _load_top_npe(root_path)

        # All events (no cut)
        mean_all = float(np.mean(npe_top))

        # T4: events with N_pe_Top >= 4
        mask_t4 = npe_top >= T4_THRESH
        frac_t4 = float(mask_t4.sum()) / n_events
        mean_t4 = float(np.mean(npe_top[mask_t4])) if mask_t4.any() else math.nan

        # T20: events with N_pe_Top >= 20
        mask_t20 = npe_top >= T20_THRESH
        frac_t20 = float(mask_t20.sum()) / n_events
        mean_t20 = float(np.mean(npe_top[mask_t20])) if mask_t20.any() else math.nan

        results.append({
            "x_mm":     x_mm,
            "n_events": n_events,
            "mean_npe_top_all":  mean_all,
            "frac_T4":           frac_t4,
            "mean_npe_top_T4":   mean_t4,
            "frac_T20":          frac_t20,
            "mean_npe_top_T20":  mean_t20,
        })
        print(
            f"[{fig_id}] x={x_mm:+5d}: N_pe_top(all)={mean_all:.1f} "
            f"T4={mean_t4:.1f}({frac_t4:.2%}) T20={mean_t20:.1f}({frac_t20:.2%})"
        )
        del npe_top

    results_sorted = sorted(results, key=lambda r: r['x_mm'])

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")

    xs = [r['x_mm'] for r in results_sorted]
    n  = len(xs)

    def _make_graph(name, title, y_vals):
        pts = [(x, y) for x, y in zip(xs, y_vals) if math.isfinite(y)]
        if not pts:
            return None
        g = ROOT.TGraph(len(pts))
        g.SetName(name); g.SetTitle(title)
        for i, (x, y) in enumerate(pts):
            g.SetPoint(i, float(x), float(y))
        return g

    g_t4  = _make_graph(f"g_T4_{fig_id}",
                         "T4 (>=4 PE) Top N_pe mean;x_{gun} [mm];<N_pe Top>",
                         [r['mean_npe_top_T4'] for r in results_sorted])
    g_t20 = _make_graph(f"g_T20_{fig_id}",
                         "T20 (>=20 PE) Top N_pe mean;x_{gun} [mm];<N_pe Top>",
                         [r['mean_npe_top_T20'] for r in results_sorted])
    g_all = _make_graph(f"g_all_{fig_id}",
                         "All events Top N_pe mean;x_{gun} [mm];<N_pe Top>",
                         [r['mean_npe_top_all'] for r in results_sorted])

    for g in (g_t4, g_t20, g_all):
        if g:
            g.Write()
            ROOT.SetOwnership(g, False)
    tf.Close()

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        'x_mm':              [r['x_mm']             for r in results_sorted],
        'n_events':          [r['n_events']          for r in results_sorted],
        'mean_npe_top_all':  [r['mean_npe_top_all']  for r in results_sorted],
        'frac_T4':           [r['frac_T4']           for r in results_sorted],
        'mean_npe_top_T4':   [r['mean_npe_top_T4']   for r in results_sorted],
        'frac_T20':          [r['frac_T20']          for r in results_sorted],
        'mean_npe_top_T20':  [r['mean_npe_top_T20']  for r in results_sorted],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    now = datetime.datetime.now(datetime.UTC).isoformat()
    meta = {
        "fig_id": fig_id, "figura": "F4", "material": material,
        "dataset": dataset, "optica": optica,
        "x_mm_list": "all-31 positions",
        "sha256_manifest_declared": sha256_manifest,
        "n_positions": len(results_sorted),
        "T4_threshold": T4_THRESH,
        "T20_threshold": T20_THRESH,
        "escala": "Y log; X lineal x_gun [mm]",
        "estimador": "<N_pe_Top | N_pe_Top >= threshold>",
        "positions": results_sorted,
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "geometry_note": "Top SiPMs: global_id >= 16 (indices 16-85, Y=+30.25mm, una sola hilera)",
        "caption_label": (
            f"{material} (EndTop) perfiles Top: T4 >= {T4_THRESH} PE, "
            f"T20 >= {T20_THRESH} PE. Escala Y log (misma escala)."
        ),
        "comando": f"python3.12 f4_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    print(f"[{fig_id}] Done → {root_out.name}")


if __name__ == "__main__":
    generate_f4_sidecar(
        fig_id="fig04", material="EJ-204", dataset="exec07_endtop_2000",
        data_dir=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"),
        sha256_manifest="b67c4f7504c94e546b21f577ca686ee6ca652d46c8ed33b569645347e0cf390b",
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
