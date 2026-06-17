#!/usr/bin/env python3.12
"""
F7 sidecar — EXEC_14 (QA-3b).

Order-statistic arrival time <t_n> per Top channel.
Physics: gap between observed RMS and statistical floor (sqrt(n)*tau_d/<Npe>)
isolates optical path dispersion.

Source logic: exec12b_tn_dispersion.py (read-only; imported, not duplicated).
Dataset: EJ-204 EndTop (exec07_endtop_2000).
Positions: x=-690 mm (truncated edge cluster, main figure) + x=0 mm (control).

Sidecar outputs:
  outputs/fig07.root  — TGraphErrors per (channel, position) + stat-floor TGraphs
  outputs/fig07.csv   — flat table: x_mm, ch_id, ch_pos_mm, n, mean_tn, rms, floor
  outputs/fig07.meta.json
  outputs/fig07.pdf   — 3x3 panel figure for x=-690 mm (primary)
  outputs/fig07_x0.pdf — same for x=0 mm (control; not in main beamer slide)
"""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")

# ── Import exec12b logic (read-only; DO NOT copy) ────────────────────────────
_E12B_DIR = pathlib.Path(
    "/home/reriosto/SHiP/ej200_endonly/analysis/exec07"
)
sys.path.insert(0, str(_E12B_DIR.parent))   # .../analysis → enables exec07.common
sys.path.insert(0, str(_E12B_DIR))          # .../exec07   → enables direct import
import exec12b_tn_dispersion as _e12b       # noqa: E402
import common as _c                         # noqa: E402

OUTPUT_DIR  = pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs")
DATA_DIR    = pathlib.Path(
    "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"
)
FIG_ID      = "fig07"
POSITIONS   = (-690, 0)       # primary (truncated edge) + control
_N_MAX      = _e12b.N_MAX     # 20
_DPI        = 150
_CH_COLOR   = "#228B22"       # forest green (Top channels, consistent with exec12b)


# ── SHA-256 helper ────────────────────────────────────────────────────────────

def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


# ── Per-position computation ──────────────────────────────────────────────────

def _compute_position(x_mm: int, *, data_dir: pathlib.Path, tau_d_ns: float) -> dict:
    """Return per-channel stats for the 9 nearest Top channels at x_mm."""
    root_path = _c.expected_file(data_dir, x_mm)
    if not root_path.exists():
        raise FileNotFoundError(f"EndTop ROOT not found: {root_path}")

    event_id, global_id, time_ns = _e12b._load_arrays(data_dir, x_mm)
    ch_ids = _e12b._top_9_neighbors(x_mm)
    top_pos = list(_c.TOP_POSITIONS_MM)

    channels = []
    for ch_id in ch_ids:
        ns, means, rmss, mean_npe = _e12b._per_channel_stats(
            event_id, global_id, time_ns, ch_id
        )
        floor = [float(np.sqrt(n)) * tau_d_ns / mean_npe for n in ns]
        channels.append(
            dict(
                ch_id=ch_id,
                ch_pos_mm=top_pos[ch_id - 16],
                ns=ns,
                means_ns=means,
                rmss_ns=rmss,
                floor_ns=floor,
                mean_npe=mean_npe,
            )
        )
        print(
            f"  x={x_mm:+5d}mm ID{ch_id:3d}"
            f" xw={top_pos[ch_id-16]:+.0f}mm"
            f" Npe={mean_npe:.1f}"
            f" n_pts={len(ns)}"
        )
    return dict(x_mm=x_mm, sha=_sha256(root_path), channels=channels)


# ── ROOT sidecar ──────────────────────────────────────────────────────────────

def _write_root(data: list[dict], out_path: pathlib.Path) -> None:
    tf = ROOT.TFile(str(out_path), "RECREATE")
    for pos_data in data:
        x_mm = pos_data["x_mm"]
        x_tag = f"x{'m' if x_mm < 0 else 'p'}{abs(x_mm)}"
        for ch in pos_data["channels"]:
            ch_id = ch["ch_id"]
            ns_arr = np.array(ch["ns"], dtype=np.float64)
            means  = np.array(ch["means_ns"], dtype=np.float64)
            rmss   = np.array(ch["rmss_ns"], dtype=np.float64)
            floor  = np.array(ch["floor_ns"], dtype=np.float64)
            n_pts  = len(ns_arr)

            if n_pts == 0:
                continue

            # TGraphErrors: <t_n> ± RMS
            g = ROOT.TGraphErrors(n_pts)
            g.SetName(f"g_tn_id{ch_id}_{x_tag}")
            g.SetTitle(
                f"<t_n> +/- RMS | ID {ch_id} | x={x_mm:+d} mm;"
                "n (photon index);<t_n> [ns]"
            )
            for i, (ni, m, r) in enumerate(zip(ns_arr, means, rmss)):
                g.SetPoint(i, ni, m)
                g.SetPointError(i, 0.0, r)
            g.Write()
            ROOT.SetOwnership(g, False)

            # TGraph: statistical floor
            gf = ROOT.TGraph(n_pts)
            gf.SetName(f"g_floor_id{ch_id}_{x_tag}")
            gf.SetTitle(
                f"stat floor sqrt(n)*tau_d/<Npe> | ID {ch_id} | x={x_mm:+d} mm;"
                "n;<sigma>_stat [ns]"
            )
            for i, (ni, fl) in enumerate(zip(ns_arr, floor)):
                gf.SetPoint(i, ni, fl)
            gf.Write()
            ROOT.SetOwnership(gf, False)

    tf.Close()


# ── Matplotlib figure (3×3 grid per position) ────────────────────────────────

def _make_figure(
    pos_data: dict, out_pdf: pathlib.Path, *, material: str, tau_d_ns: float
) -> None:
    x_mm = pos_data["x_mm"]
    channels = pos_data["channels"]
    n_ch = len(channels)      # 9
    n_rows, n_cols = 3, 3

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 8), dpi=_DPI)
    fig.subplots_adjust(hspace=0.55, wspace=0.38,
                        left=0.07, right=0.98, top=0.88, bottom=0.08)
    fig.suptitle(
        f"F7 — Order-statistic arrival time $\\langle t_n\\rangle$ per Top channel"
        f"  |  $x_\\mathrm{{gun}}={x_mm:+d}$ mm ({material} EndTop)\n"
        r"Error bars = RMS (event-to-event). "
        r"Dashed band: stat.\ floor $\sqrt{n}\,\tau_d/\langle N_{pe}\rangle$, "
        f"$\\tau_d={tau_d_ns}$\\,ns",
        fontsize=8, y=0.97,
    )

    for i, ax in enumerate(axes.flat):
        if i >= n_ch:
            ax.set_visible(False)
            continue

        ch = channels[i]
        ns_arr   = np.array(ch["ns"])
        means    = np.array(ch["means_ns"])
        rmss     = np.array(ch["rmss_ns"])
        floor    = np.array(ch["floor_ns"])
        mean_npe = ch["mean_npe"]

        if len(ns_arr) == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=7, color="gray")
            ax.set_title(f"ID {ch['ch_id']} | $x_w$={ch['ch_pos_mm']:+.0f}mm",
                         fontsize=6.5, pad=2)
            continue

        ax.errorbar(
            ns_arr, means, yerr=rmss,
            fmt="o-", color=_CH_COLOR, markersize=2.5, linewidth=1.0,
            elinewidth=0.8, capsize=2,
            label=r"$\langle t_n\rangle\pm\mathrm{RMS}$",
        )
        ax.fill_between(
            ns_arr, means - floor, means + floor,
            alpha=0.18, color="gray",
            label=r"stat. floor",
        )
        ax.plot(ns_arr, means - floor, color="gray", lw=0.7, ls="--", alpha=0.7)
        ax.plot(ns_arr, means + floor, color="gray", lw=0.7, ls="--", alpha=0.7)

        ax.set_xlim(0.5, _N_MAX + 0.5)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(5))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1))
        ax.tick_params(axis="both", labelsize=6)
        ax.set_title(
            f"ID {ch['ch_id']} | $x_w$={ch['ch_pos_mm']:+.0f} mm",
            fontsize=6.5, pad=2,
        )
        ax.annotate(
            f"$\\langle N_{{pe}}\\rangle={mean_npe:.1f}$",
            xy=(0.97, 0.03), xycoords="axes fraction",
            ha="right", va="bottom", fontsize=5.5,
        )

    for ax in axes[n_rows - 1, :]:
        ax.set_xlabel("$n$ (photon index)", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\langle t_n\rangle$ [ns]", fontsize=7)

    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure → {out_pdf.name}")


# ── CSV writer ────────────────────────────────────────────────────────────────

def _write_csv(data: list[dict], out_path: pathlib.Path) -> None:
    rows = ["x_mm,ch_id,ch_pos_mm,n,mean_tn_ns,rms_tn_ns,stat_floor_ns,mean_npe"]
    for pos in data:
        x_mm = pos["x_mm"]
        for ch in pos["channels"]:
            for n, m, r, fl in zip(
                ch["ns"], ch["means_ns"], ch["rmss_ns"], ch["floor_ns"]
            ):
                rows.append(
                    f"{x_mm},{ch['ch_id']},{ch['ch_pos_mm']:.1f},"
                    f"{n},{m:.6f},{r:.6f},{fl:.6f},{ch['mean_npe']:.4f}"
                )
    out_path.write_text("\n".join(rows) + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def generate_f7_sidecar(
    *,
    fig_id: str,
    material: str,
    dataset: str,
    data_dir: pathlib.Path,
    output_dir: pathlib.Path,
    tau_d_ns: float,
    positions: tuple[int, ...],
    sha256_inputs: dict[int, str],
    exec_tag: str = "EXEC_14",
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / fig_id

    print(f"[{fig_id}] Computing order-statistic arrival times …")
    all_data: list[dict] = []
    for x_mm in positions:
        print(f"  x = {x_mm:+d} mm")
        pos_data = _compute_position(x_mm, data_dir=data_dir, tau_d_ns=tau_d_ns)
        all_data.append(pos_data)

    # ── ROOT ──────────────────────────────────────────────────────────────────
    root_out = stem.with_suffix(".root")
    if root_out.exists():
        root_out.unlink()
    _write_root(all_data, root_out)
    print(f"[{fig_id}] ROOT → {root_out.name} ({root_out.stat().st_size/1024:.1f} kB)")

    # ── CSV ───────────────────────────────────────────────────────────────────
    csv_out = stem.with_suffix(".csv")
    _write_csv(all_data, csv_out)
    print(f"[{fig_id}] CSV  → {csv_out.name}")

    # ── Figures ───────────────────────────────────────────────────────────────
    primary_pos = positions[0]
    primary = next(d for d in all_data if d["x_mm"] == primary_pos)
    _make_figure(primary, stem.with_suffix(".pdf"), material=material, tau_d_ns=tau_d_ns)
    for pd in all_data[1:]:
        x_tag = f"xp{pd['x_mm']}" if pd["x_mm"] >= 0 else f"xm{abs(pd['x_mm'])}"
        _make_figure(
            pd,
            output_dir / f"{fig_id}_{x_tag}.pdf",
            material=material,
            tau_d_ns=tau_d_ns,
        )

    # ── meta.json ─────────────────────────────────────────────────────────────
    top_pos = list(_c.TOP_POSITIONS_MM)
    ids_primary = _e12b._top_9_neighbors(primary_pos)
    sha_by_x = {d["x_mm"]: d["sha"] for d in all_data}

    channels_by_pos: dict[str, dict] = {}
    for x_mm in positions:
        ids = _e12b._top_9_neighbors(x_mm)
        key = f"channels_x{'m' if x_mm < 0 else 'p'}{abs(x_mm)}"
        entry: dict = {
            "ch_ids": ids,
            "ch_positions_mm": [top_pos[i - 16] for i in ids],
        }
        if x_mm == positions[0]:
            entry["truncation_note"] = (
                f"At x={x_mm:+d} mm: nearest_idx=0, cluster cannot extend left. "
                "Shifted to 9 leftmost Top channels (IDs 16-24, "
                "x_w=-692...-532 mm). This is expected edge behaviour."
            )
        channels_by_pos[key] = entry

    meta = {
        "fig_id": fig_id,
        "figura": "F7",
        "material": material,
        "dataset": dataset,
        "description": (
            "Order-statistic arrival time <t_n> per Top channel. "
            "9 nearest Top channels (nearest±4 by position index). "
            "Gap between RMS and stat floor isolates optical path dispersion."
        ),
        "positions_analyzed": list(positions),
        "n_max": _N_MAX,
        "min_entries_per_point": _e12b.MIN_ENTRIES,
        "tau_d_ns": tau_d_ns,
        "stat_floor_formula": "sigma_stat(t_n) = sqrt(n) * tau_d / mean_Npe",
        **channels_by_pos,
        "physics_note": (
            f"The statistical floor sqrt(n)*tau_d/<Npe> is the irreducible "
            f"order-statistics limit for pure exponential emission (tau_d={tau_d_ns} ns). "
            "Observed RMS > floor: the excess is optical path length dispersion "
            "(same effect that broadens sigma_t toward edges and shifts v_eff onset "
            "vs nominal). tau_d here is the EMISSION constant (datasheet), "
            "not the fitted arrival-profile tau_f."
        ),
        "sha256_inputs": {
            f"x{x_mm:+d}mm": sha256_inputs.get(x_mm, sha_by_x.get(x_mm, ""))
            for x_mm in positions
        },
        "escala": "linear Y, linear X (n photon index)",
        "caption_label": (
            f"F7 ({material} EndTop): mean arrival time of the n-th photon per Top channel. "
            "Error bars = RMS. Dashed band = statistical floor sqrt(n)*tau_d/<Npe>. "
            "Gap between RMS and floor = optical path dispersion. "
            f"Primary position x={primary_pos:+d} mm (truncated cluster at left edge, IDs 16-24)."
        ),
        "comando": f"python3.12 f7_sidecar.py (fig_id={fig_id})",
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "exec_tag": exec_tag,
        "source_logic": "exec12b_tn_dispersion.py (read-only import; not duplicated)",
    }
    json_out = stem.with_suffix(".meta.json")
    json_out.write_text(json.dumps(meta, indent=2))
    print(f"[{fig_id}] meta → {json_out.name}")
    print(f"[{fig_id}] Done.")


if __name__ == "__main__":
    generate_f7_sidecar(
        fig_id="fig07",
        material="EJ-204",
        dataset="exec07_endtop_2000",
        data_dir=pathlib.Path("/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000"),
        output_dir=pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs"),
        tau_d_ns=float(_c.TAU_D_NS),
        positions=(-690, 0),
        sha256_inputs={},
        exec_tag="EXEC_14",
    )
