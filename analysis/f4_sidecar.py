#!/usr/bin/env python3.12
"""
F4 sidecar generator — EXEC_14 (QA-1b rewrite).

Perfiles de N_pe Top vs x_gun. Cuatro curvas en escala log-Y:
  - "total": <N_pe_Top total> (suma de todos los SiPMs Top, media por evento)
  - "nearest": <N_pe del SiPM Top más cercano al muón> via nearest_top_ids(x,1)
  - "T4": <N_pe_Top total | N_pe_Top >= 4> — idem pero solo eventos que pasan corte
  - "T20": <N_pe_Top total | N_pe_Top >= 20> — ídem con corte más alto

Si T4 y T20 coinciden con "total" al 100%, se documenta que el Top nunca baja de
20 PE en el rango simulado (señal robusta; corte no discrimina).

PENDIENTE DE CONFIRMACIÓN CON GERARDO: ¿T4/T20 era corte en N_pe_total (implementado)
o eficiencia de trigger SUM4 Top (que satura al 100% según datos)?
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
    verify_root_input, write_meta, write_csv_arrays,
    TIME_CONVENTION, JITTER_PER_HIT_NS,
)

sys.path.insert(0, '/home/reriosto/SHiP/ej200_endonly/analysis/exec07')
import common  # nearest_top_ids, TOP_LEFT_IDS, TOP_RIGHT_IDS

ALL_POSITIONS = [
    0, -50, 50, -100, 100, -150, 150, -200, 200,
    -250, 250, -300, 300, -350, 350, -400, 400,
    -450, 450, -500, 500, -550, 550, -600, 600,
    -650, 650, -670, 670, -690, 690,
]

T4_THRESH  = 4
T20_THRESH = 20


def _load_top_per_event(root_path: pathlib.Path, x_mm_gun: int
                         ) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Returns (npe_total_top, npe_nearest_top, n_events).
      npe_total_top[ev]   = sum of N_pe across all Top SiPMs for event ev
      npe_nearest_top[ev] = N_pe of the single nearest Top SiPM for event ev
    """
    with uproot.open(str(root_path)) as f:
        tree = f['sipm_hits']
        arrays = tree.arrays(['event_id', 'global_id'], library='np')

    event_id  = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    del arrays

    n_events = int(np.max(event_id)) + 1
    top_mask = global_id >= 16
    ev_top   = event_id[top_mask]
    g_top    = global_id[top_mask]
    del event_id, global_id

    npe_total = np.bincount(ev_top.astype(np.intp), minlength=n_events)

    nearest_id = int(common.nearest_top_ids(x_mm_gun, 1)[0])
    near_mask  = g_top == nearest_id
    ev_near    = ev_top[near_mask]
    npe_nearest = np.bincount(ev_near.astype(np.intp), minlength=n_events)

    del ev_top, g_top, ev_near
    return npe_total, npe_nearest, n_events


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

        npe_total, npe_nearest, n_events = _load_top_per_event(root_path, x_mm)

        mean_total   = float(np.mean(npe_total))
        mean_nearest = float(np.mean(npe_nearest))

        mask_t4  = npe_total >= T4_THRESH
        frac_t4  = float(mask_t4.sum()) / n_events
        mean_t4  = float(np.mean(npe_total[mask_t4])) if mask_t4.any() else math.nan

        mask_t20 = npe_total >= T20_THRESH
        frac_t20 = float(mask_t20.sum()) / n_events
        mean_t20 = float(np.mean(npe_total[mask_t20])) if mask_t20.any() else math.nan

        nearest_id = int(common.nearest_top_ids(x_mm, 1)[0])
        results.append({
            "x_mm": x_mm,
            "n_events": n_events,
            "nearest_sipm_id": nearest_id,
            "mean_npe_top_total":   mean_total,
            "mean_npe_top_nearest": mean_nearest,
            "frac_T4":    frac_t4,
            "mean_npe_top_T4":  mean_t4,
            "frac_T20":   frac_t20,
            "mean_npe_top_T20": mean_t20,
        })
        print(
            f"[{fig_id}] x={x_mm:+5d}mm: "
            f"total={mean_total:.1f}  nearest(ID{nearest_id})={mean_nearest:.2f}  "
            f"T4={mean_t4:.1f}({frac_t4:.0%})  T20={mean_t20:.1f}({frac_t20:.0%})"
        )
        del npe_total, npe_nearest

    results_sorted = sorted(results, key=lambda r: r['x_mm'])

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")

    def _make_graph(name, title, y_vals):
        pts = [(r['x_mm'], y) for r, y in zip(results_sorted, y_vals) if math.isfinite(y)]
        if not pts: return None
        g = ROOT.TGraph(len(pts))
        g.SetName(name); g.SetTitle(title)
        for i, (x, y) in enumerate(pts):
            g.SetPoint(i, float(x), float(y))
        return g

    for name, title, vals in [
        (f"g_total_{fig_id}",
         "N_pe Top total (all SiPMs Top);x_{gun} [mm];<N_{pe} Top total>",
         [r['mean_npe_top_total'] for r in results_sorted]),
        (f"g_nearest_{fig_id}",
         "N_pe Top nearest SiPM;x_{gun} [mm];<N_{pe} Top nearest>",
         [r['mean_npe_top_nearest'] for r in results_sorted]),
        (f"g_T4_{fig_id}",
         f"N_pe Top T4 (>={T4_THRESH} PE) mean;x_{{gun}} [mm];<N_{{pe}} Top T4>",
         [r['mean_npe_top_T4'] for r in results_sorted]),
        (f"g_T20_{fig_id}",
         f"N_pe Top T20 (>={T20_THRESH} PE) mean;x_{{gun}} [mm];<N_{{pe}} Top T20>",
         [r['mean_npe_top_T20'] for r in results_sorted]),
    ]:
        g = _make_graph(name, title, vals)
        if g:
            g.Write()
            ROOT.SetOwnership(g, False)

    tf.Close()

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        'x_mm':                  [r['x_mm']               for r in results_sorted],
        'n_events':              [r['n_events']             for r in results_sorted],
        'nearest_sipm_id':       [r['nearest_sipm_id']     for r in results_sorted],
        'mean_npe_top_total':    [r['mean_npe_top_total']  for r in results_sorted],
        'mean_npe_top_nearest':  [r['mean_npe_top_nearest'] for r in results_sorted],
        'frac_T4':               [r['frac_T4']             for r in results_sorted],
        'mean_npe_top_T4':       [r['mean_npe_top_T4']     for r in results_sorted],
        'frac_T20':              [r['frac_T20']            for r in results_sorted],
        'mean_npe_top_T20':      [r['mean_npe_top_T20']   for r in results_sorted],
    })

    # Diagnose T4/T20 saturation
    all_t4_100  = all(r['frac_T4']  >= 1.0 for r in results_sorted)
    all_t20_100 = all(r['frac_T20'] >= 1.0 for r in results_sorted)

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
        "T4_saturated_100pct": all_t4_100,
        "T20_saturated_100pct": all_t20_100,
        "escala": "Y log; X lineal x_gun [mm]",
        "curves": {
            "total":   "media de N_pe_Top_total por evento (suma sobre todos SiPMs Top)",
            "nearest": "media de N_pe del SiPM Top más cercano al muón (nearest_top_ids(x,1))",
            "T4":  f"media de N_pe_Top_total para eventos con N_pe_Top >= {T4_THRESH}",
            "T20": f"media de N_pe_Top_total para eventos con N_pe_Top >= {T20_THRESH}",
        },
        "t4_t20_note": (
            f"T4 y T20 saturadas al 100% en todas las posiciones: "
            f"el Top nunca baja de {T20_THRESH} PE (señal robusta). "
            "Corte no discrimina en el rango simulado. "
            "PENDIENTE CONFIRMACIÓN GERARDO: ¿T4/T20 era corte en N_pe o eficiencia de trigger?"
            if (all_t4_100 and all_t20_100)
            else "T4/T20 varían con la posición; ver CSV."
        ),
        "geometry_note": (
            "Top SiPMs: global_id 16-85 (TOP_LEFT 16-50, TOP_RIGHT 51-85). "
            "Una sola hilera en Y=+30.25mm. nearest_top_ids desde common.py (single source of truth)."
        ),
        "positions": results_sorted,
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "caption_label": (
            f"{material} (EndTop) perfiles N_pe Top vs x_gun (escala log-Y). "
            "Total: todos los SiPMs Top. Nearest: el más cercano al muón. "
            "T4/T20: solo eventos con N_pe_total >= umbral. "
            "PENDIENTE confirmación Gerardo sobre definición T4/T20."
        ),
        "comando": f"python3.12 f4_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    print(f"[{fig_id}] Done → {root_out.name}")
    print(f"[{fig_id}] T4 saturated: {all_t4_100}, T20 saturated: {all_t20_100}")


if __name__ == "__main__":
    generate_f4_sidecar(
        fig_id="fig04", material="EJ-204", dataset="exec07_endtop_2000",
        data_dir=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"),
        sha256_manifest="b67c4f7504c94e546b21f577ca686ee6ca652d46c8ed33b569645347e0cf390b",
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
