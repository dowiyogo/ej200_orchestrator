#!/usr/bin/env python3.12
"""
F6 sidecar generator — EXEC_14.

Redundancia Top: correlación de N_pe entre pares de SiPMs Top cercanos.
Dataset EndTop (EJ-204, exec07_endtop_2000), posiciones cercanas al centro.

Pares analizados (global_id):
  Pair A — (49, 52): X=-32mm y X=+32mm, separación 64mm (par simétrico)
  Pair B — (50, 51): X=-12mm y X=+12mm, separación 24mm (par adyacente al gap)
  Pair C control — (49, 47): X=-32mm y X=-72mm, sep 40mm (2 pasos de 20mm)

Posiciones usadas: x_gun in [-50, 0, +50] mm (fotones van al centro del bar).

Geometría Top:
  ID = 16 + local_idx; X_sipm = -692 + 20*local_idx mm
  ID 47 → local_idx=31 → X=-72mm
  ID 49 → local_idx=33 → X=-32mm
  ID 50 → local_idx=34 → X=-12mm
  ID 51 → local_idx=35 → X=+12mm
  ID 52 → local_idx=36 → X=+32mm
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

# Top SiPM global IDs of interest and their X positions
TOP_IDS = {47: -72, 49: -32, 50: -12, 51: 12, 52: 32}
PAIRS = [
    {"name": "pair_A_49_52", "id1": 49, "id2": 52,
     "label": "IDs 49/52 (X=-32mm / +32mm, sep 64mm)"},
    {"name": "pair_B_50_51", "id1": 50, "id2": 51,
     "label": "IDs 50/51 (X=-12mm / +12mm, sep 24mm, adyacentes gap)"},
    {"name": "pair_C_47_49", "id1": 47, "id2": 49,
     "label": "IDs 47/49 (X=-72mm / -32mm, sep 40mm; control)"},
]
CENTER_POSITIONS = [-50, 0, 50]


def _pearson_r(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3:
        return math.nan
    a_c = a - np.mean(a); b_c = b - np.mean(b)
    denom = math.sqrt(np.dot(a_c, a_c) * np.dot(b_c, b_c))
    return float(np.dot(a_c, b_c) / denom) if denom > 0 else math.nan


def _load_npe_per_sipm(root_path: pathlib.Path, ids_of_interest: list[int]
                        ) -> tuple[dict[int, np.ndarray], int]:
    """
    For each global_id in ids_of_interest, return count of hits per event.
    Returns: ({global_id: npe_array}, n_events)
    """
    with uproot.open(str(root_path)) as f:
        tree = f['sipm_hits']
        arrays = tree.arrays(['event_id', 'global_id'], library='np')

    event_id  = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    del arrays

    n_events = int(np.max(event_id)) + 1
    result = {}
    for gid in ids_of_interest:
        mask = global_id == gid
        ev_g = event_id[mask]
        npe  = np.bincount(ev_g.astype(np.intp), minlength=n_events)
        result[gid] = npe

    del event_id, global_id
    return result, n_events


def generate_f6_sidecar(
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

    ids_needed = list(TOP_IDS.keys())

    # Accumulate N_pe arrays per SiPM across all center positions
    npe_accumulated: dict[int, list[np.ndarray]] = {gid: [] for gid in ids_needed}
    n_events_total = 0
    positions_loaded = []

    for x_mm in CENTER_POSITIONS:
        root_path = data_dir / f"photon_hits_x{x_mm}mm.root"
        verify_root_input(root_path)
        print(f"[{fig_id}] x={x_mm}mm ...", flush=True)
        npe_dict, n_ev = _load_npe_per_sipm(root_path, ids_needed)
        for gid in ids_needed:
            npe_accumulated[gid].append(npe_dict[gid])
        n_events_total += n_ev
        positions_loaded.append(x_mm)
        print(f"[{fig_id}]   n_events={n_ev}, "
              f"ID49 mean={np.mean(npe_dict[49]):.2f}, "
              f"ID52 mean={np.mean(npe_dict[52]):.2f}")

    # Concatenate across positions
    npe_all: dict[int, np.ndarray] = {
        gid: np.concatenate(npe_accumulated[gid]) for gid in ids_needed
    }

    # Compute correlations per pair
    pair_results = []
    for pair in PAIRS:
        id1, id2 = pair["id1"], pair["id2"]
        a = npe_all[id1]; b = npe_all[id2]
        r = _pearson_r(a, b)
        # Corr on events where BOTH have > 0 hits
        both_mask = (a > 0) & (b > 0)
        r_both = _pearson_r(a[both_mask], b[both_mask])
        pair_results.append({
            **pair,
            "x_id1_mm": TOP_IDS[id1],
            "x_id2_mm": TOP_IDS[id2],
            "separation_mm": abs(TOP_IDS[id1] - TOP_IDS[id2]),
            "n_events":        n_events_total,
            "n_events_both_gt0": int(both_mask.sum()),
            "mean_npe_id1":    float(np.mean(a)),
            "mean_npe_id2":    float(np.mean(b)),
            "pearson_r_all":   r,
            "pearson_r_both_gt0": r_both,
        })
        print(
            f"[{fig_id}] {pair['name']}: r={r:.4f}, r(both>0)={r_both:.4f}, "
            f"n_both={both_mask.sum()}"
        )

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")

    for i_pair, pair in enumerate(PAIRS):
        id1, id2 = pair["id1"], pair["id2"]
        a = npe_all[id1]; b = npe_all[id2]
        amax = int(np.percentile(a, 99)) + 2
        bmax = int(np.percentile(b, 99)) + 2
        h2 = ROOT.TH2D(
            pair["name"],
            f"{pair['label']};N_pe ID{id1};N_pe ID{id2}",
            min(amax, 40), 0, amax,
            min(bmax, 40), 0, bmax,
        )
        h2.SetDirectory(ROOT.nullptr)
        for ai, bi in zip(a.tolist(), b.tolist()):
            h2.Fill(float(ai), float(bi))
        h2.Write()

    # TGraph of Pearson r vs pair separation
    g_r = ROOT.TGraph(len(pair_results))
    g_r.SetName(f"g_pearson_{fig_id}")
    g_r.SetTitle("Correlacion N_pe vs separacion pares Top;"
                 "Separacion [mm];Pearson r (ambos>0)")
    for i, pr in enumerate(pair_results):
        g_r.SetPoint(i, float(pr['separation_mm']), pr['pearson_r_both_gt0'])
    g_r.Write()
    ROOT.SetOwnership(g_r, False)
    tf.Close()

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        'pair_name':        [r['name']                 for r in pair_results],
        'id1':              [r['id1']                  for r in pair_results],
        'id2':              [r['id2']                  for r in pair_results],
        'separation_mm':    [r['separation_mm']        for r in pair_results],
        'pearson_r_all':    [r['pearson_r_all']        for r in pair_results],
        'pearson_r_both_gt0': [r['pearson_r_both_gt0'] for r in pair_results],
        'n_events':         [r['n_events']             for r in pair_results],
        'n_events_both_gt0':[r['n_events_both_gt0']   for r in pair_results],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    now = datetime.datetime.now(datetime.UTC).isoformat()
    meta = {
        "fig_id": fig_id, "figura": "F6", "material": material,
        "dataset": dataset, "optica": optica,
        "positions_used": positions_loaded,
        "sha256_manifest_declared": sha256_manifest,
        "top_ids_analyzed": TOP_IDS,
        "pairs": pair_results,
        "n_events_total": n_events_total,
        "geometry_verification": {
            "single_row_Y": True,
            "Y_sipm_mm": 30.25,
            "kBarHalfY_mm": 30.0,
            "kTopHalfY_mm": 0.25,
            "gap_center_mm": "[-12, +12] mm en X",
            "local_idx_formula": "X_sipm = -692 + 20 * local_idx [mm]",
            "id_to_x_mm": {
                "47": -72, "49": -32, "50": -12, "51": 12, "52": 32,
            },
        },
        "escala": "scatter TH2D N_pe(ID_i) vs N_pe(ID_j); Pearson r vs separation",
        "caption_label": (
            f"{material} (EndTop) redundancia Top: correlacion N_pe entre pares. "
            f"Posiciones x_gun = {positions_loaded} mm. "
            "IDs 49/52 (simetrico), 50/51 (adyacente gap), 47/49 (control 40mm)."
        ),
        "time_convention": TIME_CONVENTION,
        "jitter_per_hit_ps": JITTER_PER_HIT_NS * 1000.0,
        "comando": f"python3.12 f6_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    print(f"[{fig_id}] Done → {root_out.name}")


if __name__ == "__main__":
    generate_f6_sidecar(
        fig_id="fig06", material="EJ-204", dataset="exec07_endtop_2000",
        data_dir=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"),
        sha256_manifest="b67c4f7504c94e546b21f577ca686ee6ca652d46c8ed33b569645347e0cf390b",
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
