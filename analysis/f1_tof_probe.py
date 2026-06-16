#!/usr/bin/env python3.12
"""
QA-1c — ToF falsifiability probe for F1 "second component".

Para x=690 y x=0 (control), separa hits de SiPMs End por extremo:
  - END_LEFT  (global_id 0-7)   → extremo en x = -700 mm
  - END_RIGHT (global_id 8-15)  → extremo en x = +700 mm
  - BOTH combined

t_rel = time_ns - min(time_ns per evento, usando TODOS los End SiPMs).

Si en x=+690 mm:
  dist_near (RIGHT) = 700 - 690 =   10 mm → ToF_near ≈  0.036 ns
  dist_far  (LEFT)  = 700 + 690 = 1390 mm → ToF_far  ≈  5.018 ns
  delta_t_pred = 1380 mm / 277 mm/ns ≈ 4.98 ns

Criterio de confirmación:
  |delta_t_medido - delta_t_pred| / delta_t_pred < 0.15
  Y bump desaparece en near-end-only → confirmado ToF geométrico.
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
)

# ── geometry ─────────────────────────────────────────────────────────────────
BAR_HALF_MM   = 700.0    # bar extends ± 700 mm in x
V_EFF_MM_NS   = 277.0    # effective photon velocity in bar [mm/ns]

END_LEFT_IDS  = set(range(0, 8))
END_RIGHT_IDS = set(range(8, 16))

BIN_PS   = 10.0          # 10 ps bins
HI_NS    = 10.0          # 0–10 ns window
BIN_NS   = BIN_PS * 1e-3
N_BINS   = int(round(HI_NS / BIN_NS))  # 1000 bins

BRANCHES = ['event_id', 'global_id', 'time_ns']

TOLERANCE = 0.15   # 15% criterion for ToF confirmation


def _tof_prediction(x_gun_mm: float) -> dict:
    """
    Returns predicted ToF quantities for a given gun position.
    END_RIGHT is at x = +BAR_HALF_MM; END_LEFT at x = -BAR_HALF_MM.
    Near end = RIGHT for positive x_gun; both equal for x_gun = 0.
    """
    dist_right = abs(BAR_HALF_MM - x_gun_mm)   # distance from gun to RIGHT end
    dist_left  = abs(BAR_HALF_MM + x_gun_mm)   # distance from gun to LEFT end
    tof_right  = dist_right / V_EFF_MM_NS
    tof_left   = dist_left  / V_EFF_MM_NS
    delta_pred = abs(tof_left - tof_right)
    near_end   = "RIGHT" if x_gun_mm >= 0 else "LEFT"
    far_end    = "LEFT"  if x_gun_mm >= 0 else "RIGHT"
    return {
        "dist_near_mm":    min(dist_right, dist_left),
        "dist_far_mm":     max(dist_right, dist_left),
        "tof_near_ns":     min(tof_right, tof_left),
        "tof_far_ns":      max(tof_right, tof_left),
        "delta_t_pred_ns": delta_pred,
        "near_end":        near_end,
        "far_end":         far_end,
    }


def _load_end_by_side(root_path: pathlib.Path) -> tuple[dict, int, int]:
    """
    Load End SiPM hits and split by side.
    Returns:
      sides = {
        "left":  t_rel array (END_LEFT only),
        "right": t_rel array (END_RIGHT only),
        "both":  t_rel array (all End SiPMs),
      }
      n_hits_total, n_events
    """
    with uproot.open(str(root_path)) as f:
        tree = f['sipm_hits']
        n_hits_total = int(tree.num_entries)
        arrays = tree.arrays(BRANCHES, library='np')

    event_id  = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    time_ns   = arrays['time_ns'].astype(np.float64)
    del arrays

    n_events = int(np.max(event_id)) + 1

    # End SiPMs: global_id 0-15
    end_mask  = global_id < 16
    ev_end    = event_id[end_mask]
    t_end     = time_ns[end_mask]
    g_end     = global_id[end_mask]

    # t_min per event from ALL End SiPMs (same as f1_sidecar.py)
    t_min = np.full(n_events, np.inf)
    np.minimum.at(t_min, ev_end, t_end)
    invalid = ~np.isfinite(t_min)
    if invalid.any():
        t_min[invalid] = np.nan

    t_min_hit = t_min[ev_end]
    valid     = np.isfinite(t_min_hit)

    t_rel_all = t_end[valid] - t_min_hit[valid]
    g_valid   = g_end[valid]
    t_rel_all = t_rel_all[t_rel_all >= 0.0]
    g_valid   = g_valid[t_end[valid] - t_min_hit[valid] >= 0.0]

    left_mask  = np.isin(g_valid, list(END_LEFT_IDS))
    right_mask = np.isin(g_valid, list(END_RIGHT_IDS))

    sides = {
        "left":  t_rel_all[left_mask][t_rel_all[left_mask] < HI_NS],
        "right": t_rel_all[right_mask][t_rel_all[right_mask] < HI_NS],
        "both":  t_rel_all[t_rel_all < HI_NS],
    }

    del event_id, global_id, time_ns, ev_end, t_end, g_end
    del t_min, t_min_hit, valid, t_rel_all, g_valid
    return sides, n_hits_total, n_events


def _build_th1(t_rel: np.ndarray, name: str, title: str) -> ROOT.TH1D:
    h = ROOT.TH1D(name, title, N_BINS, 0.0, HI_NS)
    h.SetDirectory(ROOT.nullptr)
    for v in t_rel:
        h.Fill(float(v))
    return h


def _peak_stats(h: ROOT.TH1D, lo_ns: float, hi_ns: float) -> dict:
    """
    Mode (bin with max content), onset (first bin > 5% of max), and percentiles
    of the distribution in [lo_ns, hi_ns].
    """
    lo_bin = max(1, h.FindBin(lo_ns + 1e-9))
    hi_bin = min(h.GetNbinsX(), h.FindBin(hi_ns - 1e-9))

    # Mode
    max_c = 0.0
    peak_ns = math.nan
    for b in range(lo_bin, hi_bin + 1):
        c = h.GetBinContent(b)
        if c > max_c:
            max_c = c
            peak_ns = h.GetBinCenter(b)

    # Onset: first bin where content exceeds 5% of max (ignoring first bin to skip the t=0 spike)
    threshold = 0.05 * max_c
    onset_ns = math.nan
    for b in range(lo_bin + 1, hi_bin + 1):   # +1 skips first bin (t=0 spike)
        if h.GetBinContent(b) >= threshold:
            onset_ns = h.GetBinCenter(b)
            break

    # 10th percentile (onset-like, cumulative)
    total = sum(h.GetBinContent(b) for b in range(lo_bin, hi_bin + 1))
    p10_ns = math.nan
    running = 0.0
    if total > 0:
        for b in range(lo_bin, hi_bin + 1):
            running += h.GetBinContent(b)
            if running >= 0.10 * total:
                p10_ns = h.GetBinCenter(b)
                break

    return {"peak_ns": peak_ns, "peak_content": max_c,
            "onset_ns": onset_ns, "p10_ns": p10_ns}


def generate_f1_tof_probe(
    *,
    fig_id: str,
    output_dir: pathlib.Path,
    cases: list[dict],  # each: {material, dataset, root_path, sha256, x_mm}
    optica: str = "sslg4",
) -> None:
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / fig_id

    for ext in ('.root', '.csv', '.meta.json'):
        p = stem.with_suffix(ext)
        if p.exists():
            p.unlink()

    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")

    results     = []
    csv_rows    = []

    for case in cases:
        mat   = case['material']
        x_mm  = case['x_mm']
        path  = pathlib.Path(case['root_path'])
        sha   = case['sha256']
        label = f"{mat}_x{x_mm:+d}"

        verify_root_input(path)
        print(f"[{fig_id}] SHA verify {label} ...", end=' ', flush=True)
        verify_sha(path, sha, label)
        print("OK")

        print(f"[{fig_id}] Loading {label} ...", flush=True)
        sides, n_hits_total, n_events = _load_end_by_side(path)

        pred = _tof_prediction(float(x_mm))

        # Build histograms
        h_both  = _build_th1(sides["both"],
            f"h_both_{label}",
            f"t_rel ALL End;t_{{rel}} [ns];Counts/10 ps  {mat} x={x_mm:+d}mm")
        h_left  = _build_th1(sides["left"],
            f"h_left_{label}",
            f"t_rel END_LEFT (0-7);t_{{rel}} [ns];Counts/10 ps  {mat} x={x_mm:+d}mm")
        h_right = _build_th1(sides["right"],
            f"h_right_{label}",
            f"t_rel END_RIGHT (8-15);t_{{rel}} [ns];Counts/10 ps  {mat} x={x_mm:+d}mm")

        # Write to ROOT file
        h_both.Write(); h_left.Write(); h_right.Write()

        # Peak positions
        # For x=690: near=RIGHT peaks near 0; far=LEFT peaks near 5 ns.
        # Use full range [0,10] — the mode will land naturally.
        # For near-end (expected peak near 0), also search [0, 3] ns to be safe.
        # For far-end (expected peak near 5 ns), also search [2, 9] ns.
        near_end  = pred["near_end"].lower()   # "right" or "left"
        far_end   = pred["far_end"].lower()    # "left" or "right"
        h_near    = h_right if near_end == "right" else h_left
        h_far     = h_left  if far_end  == "left"  else h_right

        stats_near = _peak_stats(h_near, 0.0, HI_NS)
        stats_far  = _peak_stats(h_far,  0.0, HI_NS)
        # Restrict far-end onset search to [2, 9] ns (avoid near-end contamination at t≈0)
        stats_far_hi = _peak_stats(h_far, 2.0, 9.0)

        t_peak_near  = stats_near["peak_ns"]
        t_peak_far   = stats_far["peak_ns"]
        t_onset_near = stats_near["onset_ns"]
        t_onset_far  = stats_far_hi["onset_ns"]  # onset in [2,9] ns = far-end arrival
        t_p10_near   = stats_near["p10_ns"]
        t_p10_far    = stats_far_hi["p10_ns"]    # p10 in [2,9] ns

        # Primary delta_t estimator: onset of far-end (most physical — matches pure ToF)
        delta_onset  = (abs(t_onset_far - t_onset_near)
                        if math.isfinite(t_onset_far) and math.isfinite(t_onset_near)
                        else math.nan)
        # Secondary: mode-based (biased by emission time, but robust)
        delta_mode   = abs(t_peak_far - t_peak_near)

        delta_pred   = pred["delta_t_pred_ns"]
        # Use onset for criterion if available, else mode
        delta_meas   = delta_onset if math.isfinite(delta_onset) else delta_mode
        rel_diff     = abs(delta_meas - delta_pred) / max(delta_pred, 0.1)
        tof_ok       = (rel_diff < TOLERANCE)

        # Check near-end-only: does the second bump (>3 ns) disappear?
        # Count integral of near-end histogram in [3, 9] ns
        bins_hi_near = [h_near.GetBinContent(b)
                        for b in range(h_near.FindBin(3.0), h_near.FindBin(9.0) + 1)]
        integral_near_hi = sum(bins_hi_near)
        total_near       = h_near.GetEntries()
        frac_near_hi     = integral_near_hi / max(total_near, 1.0)

        # Count same range for far-end (should be large fraction)
        bins_hi_far = [h_far.GetBinContent(b)
                       for b in range(h_far.FindBin(3.0), h_far.FindBin(9.0) + 1)]
        integral_far_hi = sum(bins_hi_far)
        total_far       = h_far.GetEntries()
        frac_far_hi     = integral_far_hi / max(total_far, 1.0)

        # Verdict placeholder (computed in two-pass below)
        verdict = "PENDING"

        result = {
            "material":              mat,
            "x_mm":                  x_mm,
            "n_hits_total":          n_hits_total,
            "n_events":              n_events,
            "n_hits_end_left":       int(len(sides["left"])),
            "n_hits_end_right":      int(len(sides["right"])),
            "near_end":              pred["near_end"],
            "far_end":               pred["far_end"],
            "dist_near_mm":          pred["dist_near_mm"],
            "dist_far_mm":           pred["dist_far_mm"],
            "t_peak_near_ns":        t_peak_near,
            "t_peak_far_ns":         t_peak_far,
            "t_onset_near_ns":       t_onset_near,
            "t_onset_far_ns":        t_onset_far,
            "t_p10_near_ns":         t_p10_near,
            "t_p10_far_ns":          t_p10_far,
            "delta_t_onset_ns":      delta_onset,
            "delta_t_mode_ns":       delta_mode,
            "delta_t_meas_ns":       delta_meas,   # primary (onset if available)
            "delta_t_pred_ns":       delta_pred,
            "rel_diff_pct":          rel_diff * 100.0,
            "tof_tolerance_pct":     TOLERANCE * 100.0,
            "frac_near_hi_3to9ns":   frac_near_hi,
            "frac_far_hi_3to9ns":    frac_far_hi,
            "verdict":               verdict,
            "note_estimator":        (
                "delta_t_meas = onset del extremo lejano (primer bin > 5% del max en [2,9] ns) "
                "− onset del extremo cercano. Onset es más físico que el modo porque "
                "no incluye el ensanchamiento por tiempo de emisión."
            ),
        }
        results.append(result)

        print(
            f"[{fig_id}] {label}: near={pred['near_end']}  far={pred['far_end']}\n"
            f"         mode:  t_near={t_peak_near:.3f} ns  t_far={t_peak_far:.3f} ns"
            f"  delta_mode={delta_mode:.3f} ns\n"
            f"         onset: t_near={t_onset_near:.3f} ns  t_far={t_onset_far:.3f} ns"
            f"  delta_onset={delta_onset:.3f} ns\n"
            f"         p10:   t_near={t_p10_near:.3f} ns  t_far={t_p10_far:.3f} ns\n"
            f"         pred={delta_pred:.3f} ns  meas(onset)={delta_meas:.3f} ns"
            f"  rel_diff={rel_diff*100:.1f}%\n"
            f"         frac_near_hi={frac_near_hi:.3f}  frac_far_hi={frac_far_hi:.3f}\n"
            f"         VERDICT: {verdict}"
        )

        csv_rows.append({
            "material":          mat,
            "x_mm":              x_mm,
            "near_end":          pred["near_end"],
            "far_end":           pred["far_end"],
            "t_peak_near_ns":    round(t_peak_near, 4),
            "t_peak_far_ns":     round(t_peak_far, 4),
            "t_onset_far_ns":    round(t_onset_far, 4) if math.isfinite(t_onset_far) else "nan",
            "delta_t_onset_ns":  round(delta_onset, 4) if math.isfinite(delta_onset) else "nan",
            "delta_t_mode_ns":   round(delta_mode, 4),
            "delta_t_pred_ns":   round(delta_pred, 4),
            "rel_diff_pct":      round(rel_diff * 100.0, 1),
            "frac_near_hi":      round(frac_near_hi, 4),
            "frac_far_hi":       round(frac_far_hi, 4),
            "verdict":           verdict,
        })

        del sides, h_both, h_left, h_right

    tf.Close()

    # ── Two-pass verdict ─────────────────────────────────────────────────────
    # Pass 1: build control baseline (x=0 for each material)
    ctrl_frac = {}  # ctrl_frac[material] = frac_near_hi from x=0
    for r in results:
        if r["x_mm"] == 0:
            ctrl_frac[r["material"]] = r["frac_near_hi_3to9ns"]

    # Pass 2: assign verdicts
    BUMP_EXCESS_TOL = 1.30   # near-end frac ≤ 130% of control baseline → bump absent
    for r in results:
        mat  = r["material"]
        x    = r["x_mm"]
        dpred = r["delta_t_pred_ns"]

        if dpred < 0.2:   # control case
            # Expect: both sides have similar frac_hi (symmetric at x=0)
            sym_ratio = (r["frac_near_hi_3to9ns"] /
                         max(r["frac_far_hi_3to9ns"], 0.01))
            sym_ok    = 0.85 < sym_ratio < 1.15
            r["verdict"] = "CONTROL_OK" if sym_ok else "CONTROL_ANOMALY"
            r["verdict_detail"] = (
                f"sym_ratio={sym_ratio:.2f} "
                f"(near={r['frac_near_hi_3to9ns']:.3f}  far={r['frac_far_hi_3to9ns']:.3f})"
            )
        else:   # asymmetric position — apply ToF criterion
            tof_ok = r["rel_diff_pct"] < TOLERANCE * 100.0
            ctrl_baseline = ctrl_frac.get(mat, 0.30)
            bump_ratio    = r["frac_near_hi_3to9ns"] / max(ctrl_baseline, 0.01)
            bump_absent   = bump_ratio <= BUMP_EXCESS_TOL
            r["verdict"] = ("ToF_CONFIRMADO" if (tof_ok and bump_absent)
                            else "ToF_NO_CONFIRMADO")
            r["verdict_detail"] = (
                f"onset_rel_diff={r['rel_diff_pct']:.1f}%<{TOLERANCE*100:.0f}%={'OK' if tof_ok else 'FAIL'}  "
                f"bump_ratio={bump_ratio:.2f}<=1.30={'OK' if bump_absent else 'FAIL'}"
            )
            r["ctrl_baseline_frac_near"] = ctrl_baseline
            r["bump_ratio_vs_ctrl"]      = bump_ratio

    # Update csv_rows with verdicts (post-assignment)
    for i, csv_row in enumerate(csv_rows):
        csv_row["verdict"]        = results[i]["verdict"]
        csv_row["verdict_detail"] = results[i].get("verdict_detail", "")

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        k: [row[k] for row in csv_rows]
        for k in csv_rows[0].keys()
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    overall_verdict = (
        "ToF_CONFIRMADO" if all(
            r["verdict"] in ("ToF_CONFIRMADO", "CONTROL_OK") for r in results
        ) else "REVISAR"
    )
    print("\n[verdicts]")
    for r in results:
        print(f"  {r['material']} x={r['x_mm']:+4d}: {r['verdict']}  "
              f"({r.get('verdict_detail','')})")

    now = datetime.datetime.now(datetime.UTC).isoformat()
    meta = {
        "fig_id":         fig_id,
        "figura":         "QA-1c",
        "descripcion":    (
            "Prueba de falsificabilidad ToF: separa hits End por extremo (LEFT/RIGHT). "
            "Si el pico del extremo lejano aparece a delta_t ≈ 1380mm/277mm/ns ≈ 4.98 ns "
            "Y desaparece en el extremo cercano solo, la segunda componente de F1 es ToF "
            "geométrico, NO modo lento de centelleo."
        ),
        "geometry": {
            "bar_half_mm":    BAR_HALF_MM,
            "v_eff_mm_ns":    V_EFF_MM_NS,
            "end_left_ids":   list(range(0, 8)),
            "end_right_ids":  list(range(8, 16)),
            "end_right_pos_mm": +BAR_HALF_MM,
            "end_left_pos_mm":  -BAR_HALF_MM,
        },
        "histogram_params": {
            "bin_ps": BIN_PS,
            "hi_ns":  HI_NS,
            "n_bins": N_BINS,
        },
        "tof_tolerance_pct":   TOLERANCE * 100.0,
        "overall_verdict":     overall_verdict,
        "cases":               results,
        "interpretation": (
            "ToF_CONFIRMADO → la segunda componente en posiciones asimétricas es "
            "time-of-flight del extremo lejano (~4.98 ns a x=±690 mm), NO modo de emisión. "
            "El ajuste 1-modo es el ajuste físico; el 2-modos se etiqueta "
            "'incluye contribución ToF, no separación de modos de emisión'."
            if overall_verdict == "ToF_CONFIRMADO"
            else "REVISAR: algún caso no satisface el criterio del 15%. Ver detalle en 'cases'."
        ),
        "optica":     optica,
        "timestamp":  now,
        "exec_tag":   "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)

    print(f"\n[{fig_id}] OVERALL VERDICT: {overall_verdict}")
    print(f"[{fig_id}] Done → {root_out.name}")


if __name__ == "__main__":
    DATA_204 = pathlib.Path("/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614")
    DATA_230 = pathlib.Path("/home/reriosto/SHiP/t0minidaq/endonly_mylar_230")
    OUTPUT   = pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs")

    generate_f1_tof_probe(
        fig_id="fig01_tof_probe",
        output_dir=OUTPUT,
        cases=[
            {
                "material": "EJ-204", "x_mm": 690,
                "root_path": DATA_204 / "photon_hits_x690mm.root",
                "sha256": "61435a5c9e69350aed0ca07612483bb3c78b879b1b5edd3fe7215e2e67037d55",
                "dataset": "endonly_mylar_20260614",
            },
            {
                "material": "EJ-204", "x_mm": 0,
                "root_path": DATA_204 / "photon_hits_x0mm.root",
                "sha256": "63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
                "dataset": "endonly_mylar_20260614",
            },
            {
                "material": "EJ-230", "x_mm": 690,
                "root_path": DATA_230 / "photon_hits_x690mm.root",
                "sha256": "c67b636a1e0cf916b2a51dc7b10a9cda4784cfb8807cb2c7aded2bcf06eaad97",
                "dataset": "endonly_mylar_230",
            },
            {
                "material": "EJ-230", "x_mm": 0,
                "root_path": DATA_230 / "photon_hits_x0mm.root",
                "sha256": "0d30b51f40c9ffe46463fc1627cc87faca8900342340c176998e14069050c930",
                "dataset": "endonly_mylar_230",
            },
        ],
    )
