#!/usr/bin/env python3.12
"""
F3 sidecar generator — EXEC_14.

Event display reconstruido desde branches x_mm/y_mm/z_mm/face_type/global_id.
Elige el evento con más hits en End+Top como representativo.
Escribe TNtuple + TH2D (vista superior bar: x_mm vs z_mm) + TH2D (vista lateral: x_mm vs y_mm).
"""
from __future__ import annotations

import os as _os
_os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
_os.environ.setdefault('OMP_NUM_THREADS', '1')

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

BRANCHES_F3 = ['event_id', 'global_id', 'face_type', 'x_mm', 'y_mm', 'z_mm', 'time_ns']


def generate_f3_sidecar(
    *,
    fig_id: str,
    material: str,
    dataset: str,
    root_input: pathlib.Path,
    sha256_expected: str,
    x_mm_gun: int,
    output_dir: pathlib.Path,
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
    print(f"[{fig_id}] SHA verify ...", end=' ', flush=True)
    verify_sha(root_input, sha256_expected, fig_id)
    print("OK")

    print(f"[{fig_id}] Loading hits ...", flush=True)
    with uproot.open(str(root_input)) as f:
        tree = f['sipm_hits']
        n_hits_total = int(tree.num_entries)
        arrays = tree.arrays(BRANCHES_F3, library='np')

    event_id  = arrays['event_id'].astype(np.int64)
    global_id = arrays['global_id'].astype(np.int64)
    face_type = arrays['face_type'].astype(np.int32)
    x_mm_h    = arrays['x_mm'].astype(np.float32)
    y_mm_h    = arrays['y_mm'].astype(np.float32)
    z_mm_h    = arrays['z_mm'].astype(np.float32)
    time_ns   = arrays['time_ns'].astype(np.float32)
    del arrays

    n_events = int(np.max(event_id)) + 1

    # Choose event with most hits (most photons → richest display)
    hits_per_event = np.bincount(event_id.astype(np.intp), minlength=n_events)
    chosen_ev = int(np.argmax(hits_per_event))
    n_hits_ev = int(hits_per_event[chosen_ev])
    print(f"[{fig_id}] n_events={n_events}, chosen event={chosen_ev} ({n_hits_ev} hits)")

    mask = event_id == chosen_ev
    x_ev = x_mm_h[mask]; y_ev = y_mm_h[mask]; z_ev = z_mm_h[mask]
    t_ev = time_ns[mask]; g_ev = global_id[mask]; f_ev = face_type[mask]
    del event_id, global_id, face_type, x_mm_h, y_mm_h, z_mm_h, time_ns

    end_mask = g_ev < 16
    top_mask = g_ev >= 16
    n_end = int(end_mask.sum()); n_top = int(top_mask.sum())

    # ── ROOT sidecar ─────────────────────────────────────────────────────────
    root_out = stem.with_suffix('.root')
    tf = ROOT.TFile(str(root_out), "RECREATE")

    # TNtuple of event hits
    nt = ROOT.TNtuple(
        f"nt_event_{fig_id}", f"Event display {fig_id} ev={chosen_ev}",
        "x_mm:y_mm:z_mm:time_ns:global_id:face_type"
    )
    nt.SetDirectory(tf)
    for i in range(len(x_ev)):
        nt.Fill(float(x_ev[i]), float(y_ev[i]), float(z_ev[i]),
                float(t_ev[i]), float(g_ev[i]), float(f_ev[i]))

    # TH2D: top view (bar cross-section): x_mm vs z_mm
    h_xz = ROOT.TH2D(f"h_xz_{fig_id}",
                     f"Top view x vs z (gun x={x_mm_gun}mm);x_hit [mm];z_hit [mm]",
                     80, -50.0, 50.0, 80, -50.0, 50.0)
    h_xz.SetDirectory(ROOT.nullptr)
    for i in range(len(x_ev)):
        h_xz.Fill(float(x_ev[i]), float(z_ev[i]))

    # TH2D: side view: position along bar vs y_mm
    h_xy = ROOT.TH2D(f"h_xy_{fig_id}",
                     f"Side view bar-x vs y (gun x={x_mm_gun}mm);x_gun [mm];y_hit [mm]",
                     80, -720.0, 720.0, 60, -35.0, 35.0)
    h_xy.SetDirectory(ROOT.nullptr)
    for i in range(len(x_ev)):
        h_xy.Fill(float(x_mm_gun), float(y_ev[i]))

    nt.Write()
    h_xz.Write(); h_xy.Write()
    tf.Close()

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(stem.with_suffix('.csv'), {
        'x_mm':      x_ev.tolist(), 'y_mm': y_ev.tolist(), 'z_mm': z_ev.tolist(),
        'time_ns':   t_ev.tolist(), 'global_id': g_ev.tolist(), 'face_type': f_ev.tolist(),
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    now = datetime.datetime.now(datetime.UTC).isoformat()
    meta = {
        "fig_id": fig_id, "figura": "F3", "material": material,
        "dataset": dataset, "optica": optica,
        "x_mm_gun": x_mm_gun,
        "root_input": str(root_input), "sha256_input": sha256_expected,
        "n_hits_total": n_hits_total, "n_events": n_events,
        "chosen_event_id": chosen_ev, "n_hits_chosen_event": n_hits_ev,
        "n_hits_end": n_end, "n_hits_top": n_top,
        "geometry_note": (
            "Top SiPMs: Y=+30.25mm (una sola hilera); "
            "End SiPMs: 16 canales (0-15); "
            "kBarHalfY=30mm + kTopHalfY=0.25mm; "
            "gap centro -12 a +12 mm en X."
        ),
        "escala": "scatter plot 2D: x vs z (top view) y x vs y (side view)",
        "caption_label": (
            f"{material} evento display reconstruido desde hits (NO re-simulacion). "
            f"x_gun={x_mm_gun}mm. Evento seleccionado: max N_hits={n_hits_ev}."
        ),
        "comando": f"python3.12 f3_sidecar.py fig_id={fig_id}",
        "timestamp": now, "exec_tag": "EXEC_14",
    }
    write_meta(stem.with_suffix('.meta.json'), meta)
    del x_ev, y_ev, z_ev, t_ev, g_ev, f_ev
    print(f"[{fig_id}] Done → {root_out.name}")


if __name__ == "__main__":
    generate_f3_sidecar(
        fig_id="fig03a", material="EJ-204", dataset="endonly_mylar_20260614",
        root_input=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614/photon_hits_x0mm.root"),
        sha256_expected="63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
        x_mm_gun=0, output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
    )
