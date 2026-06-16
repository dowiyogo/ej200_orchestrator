#!/usr/bin/env python3.12
"""
F2 new sidecar — EXEC_14 QA-3c TRABAJO A.

sigma_t(N_pe) per muon position — timing resolution vs photon yield.
One point per position (31 positions), two series: weak end and mean end.
Model: sigma_t [ps] = sqrt(a^2 / N_pe + b^2)  (fit by weighted linearization).

Sources (both read-only):
  sigma_t  ← t0minidaq/…/analysis/sigma_t_sum4.csv  (curve_fit pre-computed)
  N_pe     ← EndOnly ROOT files via uproot (global_id < 16, split by end)

Outputs per material (fig02a = EJ-204, fig02b = EJ-230):
  outputs/fig02X.pdf   — sigma_t vs N_pe with fit + N_pe* marker
  outputs/fig02X.png
  outputs/fig02X.root  — TGraphErrors (weak + mean) + TF1 fit model
  outputs/fig02X.csv   — x_mm, npe_L, npe_R, npe_weak, npe_mean, sigma_t_ps, sigma_t_err_ps, …
  outputs/fig02X.meta.json
"""
from __future__ import annotations

import csv as _csv
import datetime
import json
import math
import pathlib
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import uproot
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from engine import write_meta, write_csv_arrays  # noqa: E402

OUTPUT_DIR = pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs")
_DPI = 150

DATASETS = [
    dict(
        fig_id="fig02a",
        material="EJ-204",
        data_dir=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614"
        ),
        sigma_csv=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614"
            "/analysis/sigma_t_sum4.csv"
        ),
        col="#1f78b4",
    ),
    dict(
        fig_id="fig02b",
        material="EJ-230",
        data_dir=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_230"
        ),
        sigma_csv=pathlib.Path(
            "/home/reriosto/SHiP/t0minidaq/endonly_mylar_230"
            "/analysis/sigma_t_sum4.csv"
        ),
        col="#e31a1c",
    ),
]


# ── Read pre-computed sigma_t ─────────────────────────────────────────────────

def _read_sigma_csv(path: pathlib.Path) -> dict[int, dict]:
    """sigma_t_sum4.csv → {x_mm_int: {sigma_t_ps, sigma_t_err_ps, …}}"""
    result: dict[int, dict] = {}
    with open(path, newline='') as f:
        for row in _csv.DictReader(f):
            x = int(float(row['x_mm']))
            result[x] = {
                'sigma_t_ps':       float(row['sigma_single_ps']),
                'sigma_t_err_ps':   float(row['sigma_single_error_ps']),
                'n_events':         int(float(row['n_events'])),
                'n_triggered':      int(float(row['n_triggered_lr'])),
                'trigger_efficiency': float(row['trigger_efficiency']),
            }
    return result


# ── N_pe per end (uproot, read-only) ─────────────────────────────────────────

def _compute_npe(root_path: pathlib.Path, n_events: int) -> tuple[float, float]:
    """Return (npe_L, npe_R) = mean photon hits per event for END_LEFT and END_RIGHT."""
    with uproot.open(str(root_path)) as f:
        arr = f['sipm_hits'].arrays(['event_id', 'global_id'], library='np')
    eid = arr['event_id'].astype(np.int64)
    gid = arr['global_id'].astype(np.int64)

    end = gid < 16
    eid, gid = eid[end], gid[end]

    lm = gid < 8       # END_LEFT  IDs 0-7
    rm = gid >= 8      # END_RIGHT IDs 8-15

    npe_L = float(np.mean(np.bincount(eid[lm], minlength=n_events).astype(float)))
    npe_R = float(np.mean(np.bincount(eid[rm], minlength=n_events).astype(float)))
    return npe_L, npe_R


# ── Weighted linearization fit ────────────────────────────────────────────────

def _fit_model(
    npe: np.ndarray,
    sigma: np.ndarray,
    sigma_err: np.ndarray,
) -> dict:
    """
    Fit sigma_t = sqrt(a^2/npe + b^2) via linearization.

    Substitution: y = sigma^2, u = 1/npe → y = A*u + B, A=a^2, B=b^2.
    Weighted by 1/dy^2 where dy = 2*sigma*sigma_err (error propagation).
    """
    y  = sigma ** 2
    dy = 2.0 * sigma * sigma_err
    u  = 1.0 / npe
    w  = 1.0 / dy ** 2

    X = np.column_stack([u, np.ones_like(u)])
    XtW  = X.T * w
    XtWX = XtW @ X
    XtWy = XtW @ y

    coeffs = np.linalg.solve(XtWX, XtWy)
    A_fit, B_fit = float(coeffs[0]), float(coeffs[1])

    y_pred   = A_fit * u + B_fit
    chi2     = float(np.sum(w * (y - y_pred) ** 2))
    ndf      = len(y) - 2
    chi2_ndf = chi2 / ndf if ndf > 0 else -1.0

    # Covariance matrix for parameter errors
    try:
        cov   = np.linalg.inv(XtWX)
        a_err = math.sqrt(abs(cov[0, 0])) / (2.0 * math.sqrt(max(A_fit, 1e-12)))
        b_err = math.sqrt(abs(cov[1, 1])) / (2.0 * math.sqrt(max(B_fit, 1e-12)))
    except np.linalg.LinAlgError:
        a_err = b_err = math.nan

    a_ps     = math.sqrt(max(A_fit, 0.0))
    b_ps     = math.sqrt(max(B_fit, 0.0))
    npe_star = A_fit / B_fit if B_fit > 1e-12 else math.inf

    return {
        'a_ps':     a_ps,
        'a_err_ps': a_err,
        'b_ps':     b_ps,
        'b_err_ps': b_err,
        'npe_star': npe_star,
        'chi2_ndf': chi2_ndf,
        'n_points': len(y),
        'ndf':      ndf,
    }


# ── Matplotlib figure ─────────────────────────────────────────────────────────

def _make_figure(
    *,
    fig_id: str,
    material: str,
    col: str,
    npe_weak: np.ndarray,
    npe_mean: np.ndarray,
    sigma: np.ndarray,
    sigma_err: np.ndarray,
    fit: dict,
) -> pathlib.Path:
    a_ps     = fit['a_ps']
    b_ps     = fit['b_ps']
    npe_star = fit['npe_star']
    chi2_ndf = fit['chi2_ndf']

    fig, ax = plt.subplots(figsize=(7.2, 5.5), dpi=_DPI)

    # ── Weak-end series (main) ────────────────────────────────────────────────
    sort_w = np.argsort(npe_weak)
    ax.errorbar(
        npe_weak[sort_w], sigma[sort_w], yerr=sigma_err[sort_w],
        fmt='o', color=col, ms=5.5, capsize=3, lw=1.2, zorder=5,
        label=r'weak end: $\min(N_{pe}^L,\,N_{pe}^R)$',
    )

    # ── Mean series (secondary) ───────────────────────────────────────────────
    sort_m = np.argsort(npe_mean)
    ax.errorbar(
        npe_mean[sort_m], sigma[sort_m], yerr=sigma_err[sort_m],
        fmt='s', color=col, ms=3.5, capsize=2, lw=0.9, alpha=0.45,
        ls='--', zorder=4,
        label=r'mean: $(N_{pe}^L + N_{pe}^R)/2$',
    )

    # ── Fit curve ─────────────────────────────────────────────────────────────
    npe_lo  = max(float(npe_weak.min()) * 0.75, 1.0)
    npe_hi  = float(npe_mean.max()) * 1.30
    npe_fit = np.logspace(math.log10(npe_lo), math.log10(npe_hi), 400)
    sig_fit = np.sqrt(a_ps ** 2 / npe_fit + b_ps ** 2)
    ax.plot(
        npe_fit, sig_fit, color='#33a02c', lw=2.0, zorder=6,
        label=(
            r'fit: $\sqrt{a^2/N_{pe}+b^2}$'
            f'\n$a$={a_ps:.1f} ps,  $b$={b_ps:.1f} ps'
            f'\n$\\chi^2$/ndf={chi2_ndf:.2f}'
        ),
    )

    # ── N_pe* marker ─────────────────────────────────────────────────────────
    if math.isfinite(npe_star) and npe_lo < npe_star < npe_hi:
        ax.axvline(
            npe_star, color='#6a3d9a', ls='--', lw=1.3, alpha=0.85, zorder=3,
            label=fr'$N_{{pe}}^*={npe_star:.0f}$ (Poisson = floor)',
        )
        # Text at fixed y-fraction (blended transform)
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(
            npe_star * 1.04, 0.97,
            fr'$N_{{pe}}^*$',
            transform=trans, fontsize=8, color='#6a3d9a', va='top',
        )

    ax.set(
        xscale='log',
        xlabel=r'$N_{pe}$ per end',
        ylabel=r'$\sigma_t = \sigma(\Delta T_{LR})/\sqrt{2}$ [ps]',
        title=(
            r'F2: timing resolution vs photon yield' '\n'
            r'where more light stops helping'
            f'  ({material})'
        ),
    )
    ax.legend(fontsize=7.5, loc='upper right')
    ax.grid(True, which='both', alpha=0.25, linewidth=0.5)

    out = OUTPUT_DIR / f'{fig_id}.pdf'
    png = OUTPUT_DIR / f'{fig_id}.png'
    fig.savefig(str(out), format='pdf', bbox_inches='tight')
    fig.savefig(str(png), format='png', dpi=96, bbox_inches='tight')
    plt.close(fig)
    print(f'  [{fig_id}] → {out.name} ({out.stat().st_size / 1024:.1f} kB)')
    return out


# ── ROOT sidecar ──────────────────────────────────────────────────────────────

def _write_root(
    fig_id: str,
    npe_weak: np.ndarray,
    npe_mean: np.ndarray,
    sigma: np.ndarray,
    sigma_err: np.ndarray,
    fit: dict,
) -> None:
    tf = ROOT.TFile(str(OUTPUT_DIR / f'{fig_id}.root'), 'RECREATE')
    n = len(npe_weak)

    sort_w = np.argsort(npe_weak)
    g_weak = ROOT.TGraphErrors(n)
    g_weak.SetName(f'g_weak_{fig_id}')
    g_weak.SetTitle(
        f'#sigma_t vs N_{{pe}} weak end ({fig_id});'
        'N_{pe} weak = min(N_{pe}^{L}, N_{pe}^{R});#sigma_t [ps]'
    )
    for i, idx in enumerate(sort_w):
        g_weak.SetPoint(i, float(npe_weak[idx]), float(sigma[idx]))
        g_weak.SetPointError(i, 0.0, float(sigma_err[idx]))
    g_weak.Write()
    ROOT.SetOwnership(g_weak, False)

    sort_m = np.argsort(npe_mean)
    g_mean = ROOT.TGraphErrors(n)
    g_mean.SetName(f'g_mean_{fig_id}')
    g_mean.SetTitle(
        f'#sigma_t vs N_{{pe}} mean ({fig_id});'
        'N_{pe} mean = (N_{pe}^{L}+N_{pe}^{R})/2;#sigma_t [ps]'
    )
    for i, idx in enumerate(sort_m):
        g_mean.SetPoint(i, float(npe_mean[idx]), float(sigma[idx]))
        g_mean.SetPointError(i, 0.0, float(sigma_err[idx]))
    g_mean.Write()
    ROOT.SetOwnership(g_mean, False)

    npe_lo = max(float(npe_weak.min()) * 0.75, 1.0)
    npe_hi = float(npe_mean.max()) * 1.30
    f_model = ROOT.TF1(
        f'f_model_{fig_id}',
        'sqrt([0]*[0]/x + [1]*[1])',
        npe_lo, npe_hi,
    )
    f_model.SetParameter(0, fit['a_ps'])
    f_model.SetParameter(1, fit['b_ps'])
    f_model.SetParName(0, 'a_ps')
    f_model.SetParName(1, 'b_ps')
    f_model.Write()
    ROOT.SetOwnership(f_model, False)

    tf.Close()


# ── Main per-dataset function ─────────────────────────────────────────────────

def generate_one(cfg: dict) -> dict:
    fig_id   = cfg['fig_id']
    material = cfg['material']
    data_dir = cfg['data_dir']
    col      = cfg['col']

    print(f'\n[{fig_id}] {material}')

    sigma_data = _read_sigma_csv(cfg['sigma_csv'])

    # ── Compute N_pe for each available ROOT file ─────────────────────────────
    root_files = sorted(data_dir.glob('photon_hits_x*.root'))
    rows = []

    for rp in root_files:
        x_str = rp.stem[len('photon_hits_x'):-len('mm')]
        try:
            x_mm = int(x_str)
        except ValueError:
            continue
        if x_mm not in sigma_data:
            continue
        sd = sigma_data[x_mm]
        npe_L, npe_R = _compute_npe(rp, sd['n_events'])
        rows.append({
            'x_mm':               x_mm,
            'npe_L':              npe_L,
            'npe_R':              npe_R,
            'npe_weak':           min(npe_L, npe_R),
            'npe_mean':           (npe_L + npe_R) * 0.5,
            'sigma_t_ps':         sd['sigma_t_ps'],
            'sigma_t_err_ps':     sd['sigma_t_err_ps'],
            'n_events':           sd['n_events'],
            'n_triggered':        sd['n_triggered'],
            'trigger_efficiency': sd['trigger_efficiency'],
        })
        print(
            f'  x={x_mm:+5d}mm  '
            f'npe_L={npe_L:6.1f}  npe_R={npe_R:6.1f}  '
            f'weak={min(npe_L, npe_R):6.1f}  '
            f'σt={sd["sigma_t_ps"]:.1f} ps'
        )

    rows.sort(key=lambda r: r['x_mm'])

    npe_weak_arr = np.array([r['npe_weak']       for r in rows])
    npe_mean_arr = np.array([r['npe_mean']        for r in rows])
    sigma_arr    = np.array([r['sigma_t_ps']      for r in rows])
    sigma_err    = np.array([r['sigma_t_err_ps']  for r in rows])

    # ── Fit ──────────────────────────────────────────────────────────────────
    fit = _fit_model(npe_weak_arr, sigma_arr, sigma_err)
    print(
        f'  Fit: a={fit["a_ps"]:.2f}±{fit["a_err_ps"]:.2f} ps  '
        f'b={fit["b_ps"]:.2f}±{fit["b_err_ps"]:.2f} ps  '
        f'N_pe*={fit["npe_star"]:.1f}  '
        f'χ²/ndf={fit["chi2_ndf"]:.2f} ({fit["n_points"]} pts)'
    )

    # ── Figure ───────────────────────────────────────────────────────────────
    _make_figure(
        fig_id=fig_id, material=material, col=col,
        npe_weak=npe_weak_arr, npe_mean=npe_mean_arr,
        sigma=sigma_arr, sigma_err=sigma_err, fit=fit,
    )

    # ── ROOT ─────────────────────────────────────────────────────────────────
    _write_root(fig_id, npe_weak_arr, npe_mean_arr, sigma_arr, sigma_err, fit)
    print(f'  [{fig_id}] ROOT → {fig_id}.root')

    # ── CSV ──────────────────────────────────────────────────────────────────
    write_csv_arrays(OUTPUT_DIR / f'{fig_id}.csv', {
        'x_mm':               [r['x_mm']               for r in rows],
        'npe_L':              [r['npe_L']               for r in rows],
        'npe_R':              [r['npe_R']               for r in rows],
        'npe_weak':           [r['npe_weak']            for r in rows],
        'npe_mean':           [r['npe_mean']            for r in rows],
        'sigma_t_ps':         [r['sigma_t_ps']          for r in rows],
        'sigma_t_err_ps':     [r['sigma_t_err_ps']      for r in rows],
        'n_events':           [r['n_events']             for r in rows],
        'n_triggered':        [r['n_triggered']          for r in rows],
        'trigger_efficiency': [r['trigger_efficiency']   for r in rows],
    })

    # ── meta.json ─────────────────────────────────────────────────────────────
    npe_star_val = fit['npe_star'] if math.isfinite(fit['npe_star']) else None
    write_meta(OUTPUT_DIR / f'{fig_id}.meta.json', {
        'fig_id':     fig_id,
        'figura':     'F2',
        'material':   material,
        'description': (
            'sigma_t(N_pe) per muon position (31 positions). '
            'sigma_t = sigma(DeltaT_LR)/sqrt(2) from pre-computed sigma_t_sum4.csv. '
            'N_pe per end computed from EndOnly ROOT files (uproot, read-only). '
            'Fit: sigma_t [ps] = sqrt(a^2/N_pe + b^2) by weighted linearization.'
        ),
        'n_positions': len(rows),
        'fit': {
            'model':    'sigma_t [ps] = sqrt(a^2 / N_pe + b^2)',
            'series':   'npe_weak = min(npe_L, npe_R)',
            'a_ps':     fit['a_ps'],
            'a_err_ps': fit['a_err_ps'],
            'b_ps':     fit['b_ps'],
            'b_err_ps': fit['b_err_ps'],
            'npe_star': npe_star_val,
            'chi2_ndf': fit['chi2_ndf'],
            'n_points': fit['n_points'],
            'ndf':      fit['ndf'],
            'method':   'weighted linear regression on sigma_t^2 = a^2/N_pe + b^2',
        },
        'sigma_t_source': str(cfg['sigma_csv']),
        'npe_source':     (
            'EndOnly ROOT files (uproot, read-only); '
            'np.bincount of event_id per END_LEFT (ID 0-7) and END_RIGHT (ID 8-15)'
        ),
        'end_left_ids':  list(range(0, 8)),
        'end_right_ids': list(range(8, 16)),
        'physics_note': (
            'Floor b [ps] = optical path dispersion contribution to sigma_t. '
            'Same effect as RMS-floor gap in F7 order statistics. '
            'N_pe* = a^2/b^2: above this threshold, more photons do not improve sigma_t.'
        ),
        'caption_label': (
            f'F2 ({material}): σt = σ(ΔT_LR)/√2 vs N_pe per end (31 positions). '
            f'Fit σt = √(a²/N_pe + b²): a={fit["a_ps"]:.0f} ps (Poisson), '
            f'b={fit["b_ps"]:.0f} ps (optical floor). '
            + (f'N_pe* = {npe_star_val:.0f}.' if npe_star_val else 'N_pe* undefined.')
        ),
        'timestamp': datetime.datetime.now(datetime.UTC).isoformat(),
        'exec_tag':  'EXEC_14',
        'qa_tag':    'QA-3c TRABAJO A',
    })
    print(f'  [{fig_id}] meta.json written')
    return fit


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}
    for cfg in DATASETS:
        fit = generate_one(cfg)
        summary[cfg['material']] = fit

    print('\n' + '=' * 60)
    print('F2 new sidecar — fit summary')
    print('=' * 60)
    for mat, fit in summary.items():
        nps = fit['npe_star']
        print(
            f'  {mat:8s}  a={fit["a_ps"]:6.1f}±{fit["a_err_ps"]:.1f} ps  '
            f'b={fit["b_ps"]:6.1f}±{fit["b_err_ps"]:.1f} ps  '
            f'N_pe*={nps:.1f}  χ²/ndf={fit["chi2_ndf"]:.2f}'
        )


if __name__ == '__main__':
    main()
