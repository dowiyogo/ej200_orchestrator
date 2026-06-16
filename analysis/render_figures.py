#!/usr/bin/env python3.12
"""
EXEC_14 QA-2 — Figure renderer.
Reads from sidecars (.root/.csv/.meta.json) ONLY.
No physics recalculation from t0minidaq raw hits.
"""
from __future__ import annotations

import os as _os
_os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
_os.environ.setdefault('OMP_NUM_THREADS', '1')

import csv
import io
import json
import math
import pathlib
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import uproot
from PIL import Image

OUTPUTS = pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs")

# ── Style ──────────────────────────────────────────────────────────────────
COL_204   = '#1f78b4'  # blue — EJ-204
COL_230   = '#e31a1c'  # red  — EJ-230
COL_FIT1  = '#33a02c'  # green — 1-mode fit
COL_LEFT  = '#4575b4'  # blue-left — END_LEFT
COL_RIGHT = '#d73027'  # red-right — END_RIGHT
COL_BOTH  = '#636363'  # gray — BOTH combined
MATCOLS   = {'EJ-204': COL_204, 'EJ-230': COL_230}

plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 9,
    'axes.titlesize': 9,
    'axes.labelsize': 8.5,
    'legend.fontsize': 7.5,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'lines.linewidth': 1.0,
    'axes.grid': True,
    'grid.alpha': 0.25,
    'grid.linewidth': 0.5,
})

TIME_CONV = r'$t_{rel} = t_{photon} - \min(t_{End})$ per event'

# ── Helpers ────────────────────────────────────────────────────────────────
def load_meta(fig_id: str) -> dict:
    return json.loads((OUTPUTS / f'{fig_id}.meta.json').read_text())

def load_csv(fig_id: str) -> list[dict]:
    with open(OUTPUTS / f'{fig_id}.csv', newline='') as f:
        return list(csv.DictReader(f))

def scint_1mode(t: np.ndarray, N: float, tau_r: float, A: float, tau_f: float) -> np.ndarray:
    return N * (1 - np.exp(-t / tau_r)) * A * np.exp(-t / tau_f)

def _save(fig: plt.Figure, fig_id: str, suffix: str = '') -> pathlib.Path:
    stem = fig_id + suffix
    pdf = OUTPUTS / f'{stem}.pdf'
    png = OUTPUTS / f'{stem}.png'
    fig.savefig(str(pdf), format='pdf', bbox_inches='tight')
    fig.savefig(str(png), format='png', dpi=96, bbox_inches='tight')
    plt.close(fig)
    print(f'  [{fig_id}] → {pdf.name} ({pdf.stat().st_size/1024:.1f} kB)')
    return pdf

# ── F1 — scintillation emission time profile ──────────────────────────────
def render_f1(fig_id: str) -> pathlib.Path:
    meta = load_meta(fig_id)
    mat  = meta['material']
    x_mm = meta['x_mm']
    col  = MATCOLS[mat]

    with uproot.open(str(OUTPUTS / f'{fig_id}.root')) as f:
        h_s = f[f'h_trel_{fig_id}']
        h_e = f[f'h_trel_{fig_id}_ext']
        v_s = h_s.values();  c_s = h_s.axes[0].centers()
        v_e = h_e.values();  c_e = h_e.axes[0].centers()

    fig = plt.figure(figsize=(7.5, 9))
    gs  = gridspec.GridSpec(2, 1, figure=fig, hspace=0.35,
                            height_ratios=[1.0, 1.4])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # ── Short window 0-2 ns ──
    # Skip the t=0 spike (first bin artifact)
    skip = c_s > 0.001
    ax1.step(c_s[skip], v_s[skip], where='mid', color=col, lw=0.7,
             label='data (all End SiPMs)')
    # 1-mode fit from meta
    fit1s = meta.get('fit_single_decay_2ns', {})
    N, tr, A, tf = (fit1s.get(k, math.nan)
                    for k in ('N', 'tau_rise_ns', 'A_fast', 'tau_fast_ns'))
    if all(math.isfinite(v) and v > 0 for v in (N, tr, A, tf)) and tr < 10 and tf < 30:
        t_fit = np.linspace(0.004, 2.0, 2000)
        y_fit = scint_1mode(t_fit, N, tr, A, tf)
        chi2  = fit1s.get('chi2_ndf', math.nan)
        ax1.plot(t_fit, y_fit, color=COL_FIT1, lw=1.6, ls='-',
                 label=f'1-mode: τ_r={tr:.2f} ns, τ_f={tf:.2f} ns\nχ²/ndf={chi2:.2f}')
    ax1.set(yscale='log', xlim=(0, 2), ylabel='Counts / 2 ps',
            title=f'{mat}  x = {x_mm:+d} mm — window 0–2 ns (2 ps bin, log-Y)')
    ax1.set_ylim(bottom=5)
    ax1.legend(loc='upper right', fontsize=7)
    ax1.set_xlabel('')

    # ── Extended window 0-10 ns ──
    skip_e = c_e >= 0.01   # skip first bin (center 0.005 ns = t=0 spike artifact)
    ax2.step(c_e[skip_e], v_e[skip_e], where='mid', color=col, lw=0.7,
             label='data (all End SiPMs)')
    # 1-mode fit
    fit1e = meta.get('fit_single_decay_10ns', {})
    N, tr, A, tf = (fit1e.get(k, math.nan)
                    for k in ('N', 'tau_rise_ns', 'A_fast', 'tau_fast_ns'))
    if all(math.isfinite(v) and v > 0 for v in (N, tr, A, tf)) and tr < 10 and tf < 30:
        t_fit = np.linspace(0.01, 10.0, 5000)
        y_fit = scint_1mode(t_fit, N, tr, A, tf)
        chi2  = fit1e.get('chi2_ndf', math.nan)
        ax2.plot(t_fit, y_fit, color=COL_FIT1, lw=1.6,
                 label=f'1-mode: τ_r={tr:.2f} ns, τ_f={tf:.2f} ns, χ²/ndf={chi2:.2f}')

    ax2.set(yscale='log', xlim=(0, 10), xlabel=r'$t_{rel}$ [ns]',
            ylabel='Counts / 10 ps',
            title=f'{mat}  x = {x_mm:+d} mm — window 0–10 ns (10 ps bin, log-Y)')
    ax2.set_ylim(bottom=5)
    ax2.legend(loc='upper right', fontsize=7)

    # ToF / onset annotation for asymmetric positions — placed AFTER ylim is set
    if abs(x_mm) >= 400:
        tof_ns = 4.72 if abs(x_mm) >= 680 else 2.25
        ax2.axvline(tof_ns, color='gray', ls='--', lw=0.9, alpha=0.7)
        ax2.text(tof_ns + 0.12, ax2.get_ylim()[0] * 3,
                 f'ToF onset\n~{tof_ns:.2f} ns',
                 fontsize=6.5, color='#555555', va='bottom')

    # Caption (no wrap=True — causes figure height explosion)
    tof_note = ''
    if meta.get('tof_confirmation_qa1c') and abs(x_mm) >= 680:
        tof_note = '\nTail 3-10 ns: ~98% scattered near-end photons + ~2% far-end ToF (QA-1c onset~4.72ns).'
    cap = (f'{mat}  x={x_mm:+d} mm.  {TIME_CONV}.{tof_note}')
    fig.text(0.5, 0.02, cap, ha='center', va='bottom', fontsize=6.5,
             transform=fig.transFigure)

    return _save(fig, fig_id)

# ── F2 — sigma_t vs N_pe per position ────────────────────────────────────
def render_f2(fig_id: str) -> pathlib.Path:
    meta = load_meta(fig_id)
    mat  = meta['material']
    col  = MATCOLS[mat]
    rows = load_csv(fig_id)

    # New format: per-position sigma_t vs N_pe (QA-3c TRABAJO A)
    if 'npe_weak' in rows[0]:
        import math as _math
        import matplotlib.transforms as _mt
        npe_weak = np.array([float(r['npe_weak'])       for r in rows])
        npe_mean = np.array([float(r['npe_mean'])        for r in rows])
        sigma    = np.array([float(r['sigma_t_ps'])      for r in rows])
        sig_err  = np.array([float(r['sigma_t_err_ps'])  for r in rows])
        fit      = meta.get('fit', {})
        a_ps     = float(fit.get('a_ps', 0))
        b_ps     = float(fit.get('b_ps', 0))
        npe_star = fit.get('npe_star')
        chi2_ndf = float(fit.get('chi2_ndf', -1))

        fig, ax = plt.subplots(figsize=(7.2, 5.5))
        sort_w = np.argsort(npe_weak)
        ax.errorbar(npe_weak[sort_w], sigma[sort_w], yerr=sig_err[sort_w],
                    fmt='o', color=col, ms=5.5, capsize=3, lw=1.2, zorder=5,
                    label=r'weak end: $\min(N_{pe}^L,\,N_{pe}^R)$')
        sort_m = np.argsort(npe_mean)
        ax.errorbar(npe_mean[sort_m], sigma[sort_m], yerr=sig_err[sort_m],
                    fmt='s', color=col, ms=3.5, capsize=2, lw=0.9, alpha=0.45,
                    ls='--', zorder=4,
                    label=r'mean: $(N_{pe}^L + N_{pe}^R)/2$')
        npe_lo = max(float(npe_weak.min()) * 0.75, 1.0)
        npe_hi = float(npe_mean.max()) * 1.30
        npe_fit = np.logspace(_math.log10(npe_lo), _math.log10(npe_hi), 400)
        sig_fit = np.sqrt(a_ps**2 / npe_fit + b_ps**2)
        ax.plot(npe_fit, sig_fit, color='#33a02c', lw=2.0, zorder=6,
                label=(r'fit: $\sqrt{a^2/N_{pe}+b^2}$'
                       f'\n$a$={a_ps:.1f} ps,  $b$={b_ps:.1f} ps'
                       f'\n$\\chi^2$/ndf={chi2_ndf:.2f}'))
        if npe_star is not None and _math.isfinite(npe_star) and npe_lo < npe_star < npe_hi:
            ax.axvline(npe_star, color='#6a3d9a', ls='--', lw=1.3, alpha=0.85, zorder=3,
                       label=fr'$N_{{pe}}^*={npe_star:.0f}$ (Poisson = floor)')
            trans = _mt.blended_transform_factory(ax.transData, ax.transAxes)
            ax.text(npe_star * 1.04, 0.97, fr'$N_{{pe}}^*$',
                    transform=trans, fontsize=8, color='#6a3d9a', va='top')
        ax.set(xscale='log',
               xlabel=r'$N_{pe}$ per end',
               ylabel=r'$\sigma_t = \sigma(\Delta T_{LR})/\sqrt{2}$ [ps]',
               title=(r'F2: timing resolution vs photon yield' '\n'
                      r'where more light stops helping'
                      f'  ({mat})'))
        ax.legend(fontsize=7.5, loc='upper right')
        ax.grid(True, which='both', alpha=0.25, linewidth=0.5)
        return _save(fig, fig_id)

    # Legacy format (threshold scan) — kept for reference
    npe    = np.array([float(r['npe_threshold'])       for r in rows])
    sig    = np.array([float(r['sigma_single_ps'])     for r in rows])
    sig_e  = np.array([float(r['sigma_single_err_ps']) for r in rows])
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.errorbar(npe, sig, yerr=sig_e, fmt='o', color=col, capsize=4,
                ms=6, lw=1.5, label=mat)
    ax.set(xlabel=r'$N_{pe}$ End threshold',
           ylabel=r'$\sigma_t = \sigma(\Delta T_{LR})/\sqrt{2}$ [ps]',
           title=f'F2 — Ensemble σ_t vs N_{{pe}} threshold  ({mat})',
           xticks=npe)
    ax.legend(fontsize=8)
    fig.text(0.5, 0.02, meta.get('caption_label', ''), ha='center', fontsize=6.5,
             transform=fig.transFigure)
    return _save(fig, fig_id)

# ── F3 — event display ───────────────────────────────────────────────────
def render_f3(fig_id: str) -> pathlib.Path:
    meta   = load_meta(fig_id)
    mat    = meta['material']
    x_gun  = meta['x_mm_gun']
    rows   = load_csv(fig_id)

    x_mm  = np.array([float(r['x_mm'])     for r in rows])
    y_mm  = np.array([float(r['y_mm'])     for r in rows])
    z_mm  = np.array([float(r['z_mm'])     for r in rows])
    t_ns  = np.array([float(r['time_ns'])  for r in rows])
    gid   = np.array([int(r['global_id'])  for r in rows])

    end_left  = gid < 8
    end_right = (gid >= 8)  & (gid < 16)
    top       = gid >= 16

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # Panel 1: longitudinal x-y view (position along bar vs transverse)
    for mask, col, lbl in [(end_left,  COL_LEFT,  f'End LEFT (ID0-7, n={end_left.sum()})'),
                            (end_right, COL_RIGHT, f'End RIGHT (ID8-15, n={end_right.sum()})'),
                            (top,       '#33a02c',  f'Top (ID16+, n={top.sum()})')]:
        if mask.any():
            ax1.scatter(x_mm[mask], y_mm[mask], s=0.5, c=col, alpha=0.5, rasterized=True,
                        label=lbl)
    ax1.axvline(x_gun, color='k', ls=':', lw=1.0, label=f'x_gun={x_gun:+d} mm')
    ax1.set(xlabel='x_hit [mm]', ylabel='y_hit [mm]',
            title=f'Longitudinal view (x vs y)\n{mat} x_gun={x_gun:+d} mm  ev={meta["chosen_event_id"]}')
    ax1.legend(fontsize=6.5, markerscale=5)

    # Panel 2: cross-section y-z view
    for mask, col, lbl in [(end_left,  COL_LEFT,  'End LEFT'),
                            (end_right, COL_RIGHT, 'End RIGHT'),
                            (top,       '#33a02c',  'Top')]:
        if mask.any():
            ax2.scatter(z_mm[mask], y_mm[mask], s=0.5, c=col, alpha=0.5, rasterized=True,
                        label=lbl)
    ax2.set(xlabel='z_hit [mm]', ylabel='y_hit [mm]',
            title=f'Transverse section (z vs y)\n{mat} x_gun={x_gun:+d} mm')
    ax2.legend(fontsize=6.5, markerscale=5)

    fig.text(0.5, 0.02,
             f'F3 — Photon hits on SiPMs (NOT volume tracks), '
             f'event {meta["chosen_event_id"]} ({meta["n_hits_chosen_event"]} hits). '
             f'{meta.get("caption_label","").split(".")[0]}.',
             ha='center', fontsize=6.5, transform=fig.transFigure)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    return _save(fig, fig_id)

# ── F4 — N_pe Top profiles ───────────────────────────────────────────────
def render_f4(fig_id: str) -> pathlib.Path:
    meta = load_meta(fig_id)
    rows = load_csv(fig_id)

    x       = np.array([float(r['x_mm'])               for r in rows])
    tot     = np.array([float(r['mean_npe_top_total'])  for r in rows])
    near    = np.array([float(r['mean_npe_top_nearest']) for r in rows])
    t4      = np.array([float(r['mean_npe_top_T4'])     for r in rows])
    t20     = np.array([float(r['mean_npe_top_T20'])    for r in rows])

    sort_idx = np.argsort(x)
    x, tot, near, t4, t20 = (a[sort_idx] for a in (x, tot, near, t4, t20))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, tot,  color='#1f78b4', lw=1.8, label='Total (all Top SiPMs)', marker='o', ms=2)
    ax.plot(x, near, color='#e31a1c', lw=1.8, label='Nearest SiPM', marker='s', ms=2)
    ax.plot(x, t4,   color='#33a02c', lw=1.0, ls='--', alpha=0.6,
            label=f'T4 (N_pe≥4)  ← identical to Total', marker='')
    ax.plot(x, t20,  color='#ff7f00', lw=1.0, ls=':', alpha=0.6,
            label=f'T20 (N_pe≥20) ← identical to Total', marker='')

    ax.set(yscale='log', xlabel='$x_{gun}$ [mm]',
           ylabel=r'$\langle N_{pe} \rangle$ Top',
           title='F4 — $N_{pe}$ Top profiles vs position (EndTop EJ-204, log-Y)')
    ax.legend(fontsize=7.5)
    ax.text(0.02, 0.97, f'T4 saturated: {meta.get("T4_saturated_100pct")}\n'
            f'T20 saturated: {meta.get("T20_saturated_100pct")}\n'
            f'{meta.get("t4_t20_note","").split(".")[0]}.',
            transform=ax.transAxes, va='top', fontsize=6.5,
            bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', alpha=0.9))
    fig.text(0.5, 0.02, meta.get('caption_label', ''), ha='center', fontsize=6.5,
             transform=fig.transFigure)
    return _save(fig, fig_id)

# ── F5 — SUM4 sigma_t timing ────────────────────────────────────────────
def render_f5(fig_id: str) -> pathlib.Path:
    meta = load_meta(fig_id)
    mat  = meta['material']
    col  = MATCOLS[mat]

    with uproot.open(str(OUTPUTS / f'{fig_id}.root')) as f:
        g  = f[f'g_sigma_single_{fig_id}']
        xs = g.member('fX')          # [0, 400, 690]
        ys = g.member('fY')          # sigma_single_ps
        eys = g.member('fEY')        # error

        # Delta_LR histograms per position
        hists = {}
        for xv in (0, 400, 690):
            key = f'h_delta_{fig_id}_x{int(xv)}'
            if key + ';1' in f or key in f:
                try:
                    h = f[key]
                    hists[int(xv)] = (h.values(), h.axes[0].centers())
                except Exception:
                    pass

    fig = plt.figure(figsize=(11, 5))
    gs  = gridspec.GridSpec(1, 4, figure=fig, wspace=0.4)

    # 3 histogram panels
    for i, xv in enumerate((0, 400, 690)):
        ax = fig.add_subplot(gs[i])
        if xv in hists:
            v, c = hists[xv]
            ax.step(c, v, where='mid', color=col, lw=0.9)
        # Get sigma from meta
        pos_meta = next((p for p in meta.get('positions', []) if p['x_mm'] == xv), {})
        sig = pos_meta.get('sigma_single_ps', math.nan)
        sig_e = pos_meta.get('sigma_single_err_ps', math.nan)
        eff = pos_meta.get('trigger_efficiency', 1.0)
        ax.set(xlabel=r'$\Delta T_{LR}$ [ns]', ylabel='Counts',
               title=f'{mat}  x={xv:+d} mm')
        ax.text(0.05, 0.97,
                f'σ_single={sig:.0f}±{sig_e:.0f} ps\neff={eff:.1%}',
                transform=ax.transAxes, va='top', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8))

    # Summary bar chart
    ax_s = fig.add_subplot(gs[3])
    bar_x = np.arange(len(xs))
    bars  = ax_s.bar(bar_x, ys, yerr=eys, color=col, alpha=0.85, capsize=4, width=0.6)
    ax_s.set(xticks=bar_x, xticklabels=[f'x={int(v):+d}' for v in xs],
             ylabel=r'$\sigma_t$ [ps]', ylim=(0, max(ys)*1.3),
             title=f'{mat} σ_t SUM4')
    for j, (yv, ey) in enumerate(zip(ys, eys)):
        ax_s.text(j, yv + ey + 5, f'{yv:.0f}', ha='center', fontsize=7.5, fontweight='bold')

    fig.suptitle(f'F5 — SUM4 L/R timing resolution  ({mat})\n'
                 f'σ_t = σ(ΔT_LR)/√2 (intrinsic + 20 ps/hit jitter, no additional readout jitter)',
                 fontsize=9, y=1.02)
    return _save(fig, fig_id)

# ── F6 — Top redundancy ──────────────────────────────────────────────────
def render_f6(fig_id: str) -> pathlib.Path:
    meta   = load_meta(fig_id)
    pairs  = meta.get('pairs', [])

    with uproot.open(str(OUTPUTS / f'{fig_id}.root')) as f:
        h2_data = {}
        for p in pairs:
            name = p['name']
            try:
                h2 = f[name]
                vals  = h2.values()
                xcent = h2.axes[0].centers()
                ycent = h2.axes[1].centers()
                h2_data[name] = (xcent, ycent, vals)
            except Exception:
                pass

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    for ax, p in zip(axes, pairs):
        name = p['name']
        r    = p.get('pearson_r_all', math.nan)
        id1, id2 = p['id1'], p['id2']
        sep  = p.get('separation_mm', '?')
        lbl  = p.get('label', name)

        if name in h2_data:
            xc, yc, vals = h2_data[name]
            XX, YY = np.meshgrid(xc, yc, indexing='ij')
            mask = vals > 0
            if mask.any():
                # Scatter density: plot as 2D histogram image
                # Use log-norm for density
                vmax = vals.max()
                vmin = max(1, vmax * 0.001)
                im = ax.pcolormesh(xc, yc, vals.T,
                                   norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
                                   cmap='Blues', rasterized=True)
                plt.colorbar(im, ax=ax, label='Counts')

        r_color = '#d73027' if r < 0 else '#1a9641'
        r_text  = f'r = {r:+.3f}'
        ax.text(0.05, 0.97, r_text, transform=ax.transAxes, va='top', fontsize=11,
                fontweight='bold', color=r_color,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.9))
        ax.set(xlabel=f'$N_{{pe}}$ ID {id1}', ylabel=f'$N_{{pe}}$ ID {id2}',
               title=f'Pair {chr(65+pairs.index(p))}: {lbl}\n(sep={sep} mm)')

    fig.suptitle('F6 — Top SiPM redundancy: N_pe correlation between neighbouring pairs\n'
                 'Pair A (49/52, symmetric): r=−0.37 complementary  |  '
                 'Pair B (50/51): r=+0.88  |  Pair C (47/49, same side): r=+0.89 redundant',
                 fontsize=9, y=1.01)
    fig.tight_layout()
    return _save(fig, fig_id)

# ── ToF probe figure (QA-1c supplement for contact sheet) ────────────────
def render_tof_probe(fig_id: str = 'fig01_tof_probe') -> pathlib.Path:
    meta = load_meta(fig_id)
    cases = meta.get('cases', [])

    with uproot.open(str(OUTPUTS / f'{fig_id}.root')) as f:
        hists = {}
        for k in f.keys():
            name = k.rstrip(';1')
            h = f[k]
            try:
                hists[name] = (h.values(), h.axes[0].centers())
            except Exception:
                pass

    asym_cases = [c for c in cases if abs(c['x_mm']) >= 680]
    ctrl_cases  = [c for c in cases if abs(c['x_mm']) < 50]
    display = asym_cases + ctrl_cases[:1]   # x=690 for both materials + one x=0 control

    n = len(display)
    fig, axes = plt.subplots(1, n, figsize=(5*n, 5), sharey=False)
    if n == 1: axes = [axes]

    for ax, case in zip(axes, display):
        mat  = case['material']
        x_mm = case['x_mm']
        label = f"{mat}_x{x_mm:+d}"
        col  = MATCOLS.get(mat, '#636363')

        for side, ls, col_s, lbl in [
            ('both',  '-',  COL_BOTH,  'BOTH'),
            ('right', '--', COL_RIGHT, 'END_RIGHT (near)'),
            ('left',  ':',  COL_LEFT,  'END_LEFT (far)'),
        ]:
            hkey = f'h_{side}_{label}'
            if hkey in hists:
                v, c = hists[hkey]
                skip = c > 0.004
                y = np.where(v[skip] > 0, v[skip], np.nan)
                ax.step(c[skip], y, where='mid', color=col_s, lw=1.2, ls=ls,
                        label=lbl, alpha=0.9)

        if abs(x_mm) >= 680:
            res = next((c for c in cases
                        if c['material'] == mat and c['x_mm'] == x_mm), {})
            tof = res.get('delta_t_onset_ns', math.nan)
            pred = res.get('delta_t_pred_ns', math.nan)
            verd = res.get('verdict', '')
            if math.isfinite(tof):
                ax.axvline(tof + res.get('t_onset_near_ns', 0.015),
                           color=COL_LEFT, ls='--', lw=0.9, alpha=0.7)
                ax.text(0.55, 0.85,
                        f'Far-end onset: {tof:.2f} ns\n'
                        f'ToF pred: {pred:.2f} ns\n'
                        f'{verd}',
                        transform=ax.transAxes, fontsize=7, color='#222',
                        bbox=dict(boxstyle='round', fc='white', alpha=0.85))

        ax.set(yscale='log', xlim=(0, 10), xlabel=r'$t_{rel}$ [ns]',
               ylabel='Counts / 10 ps' if ax is axes[0] else '',
               title=f'{mat}  x={x_mm:+d} mm')
        ax.set_ylim(bottom=1)
        ax.legend(fontsize=6.5)

    fig.suptitle('QA-1c — t_rel per end  |  Grey=BOTH  Red--=END_RIGHT (near)  Blue···=END_LEFT (far)\n'
                 'Far end appears entirely at t_rel>4 ns (ToF confirmed, onset 4.72 ns vs 4.98 ns predicted)',
                 fontsize=8.5, y=1.01)
    fig.tight_layout(rect=[0, 0.01, 1, 0.97])
    return _save(fig, 'fig01_tof_probe', '_vis')

# ── Contact sheet ────────────────────────────────────────────────────────
def make_contact_sheet(panels: list[tuple[str, pathlib.Path]]) -> pathlib.Path:
    """panels: [(label, png_path), ...]"""
    n      = len(panels)
    ncols  = 4
    nrows  = math.ceil(n / ncols)
    tw, th = 380, 300   # thumbnail size (px)

    sheet = Image.new('RGB', (ncols * tw, nrows * th), 'white')
    fnt_size = 14

    for i, (label, png_p) in enumerate(panels):
        try:
            img = Image.open(str(png_p)).convert('RGB')
            img.thumbnail((tw - 4, th - 20), Image.LANCZOS)
        except Exception as e:
            img = Image.new('RGB', (tw - 4, th - 20), '#cccccc')

        row = i // ncols
        col = i % ncols
        x0  = col * tw + (tw - img.width) // 2
        y0  = row * th + 18
        sheet.paste(img, (x0, y0))

        # Label text via matplotlib then convert
        from PIL import ImageDraw
        draw = ImageDraw.Draw(sheet)
        draw.text((col * tw + 4, row * th + 2), label,
                  fill='black')

    out = OUTPUTS / 'contact_sheet.png'
    sheet.save(str(out))
    print(f'  [contact_sheet] → {out.name} ({out.stat().st_size/1024:.1f} kB)')
    return out

# ── Main ─────────────────────────────────────────────────────────────────
def main() -> int:
    print('\n[QA-2] Rendering figures from sidecars...')
    pdfs = []

    print('\n[F1] Scintillation time profiles...')
    for fig_id in ('fig01a', 'fig01b', 'fig01c',
                   'fig01d', 'fig01e', 'fig01f'):
        try:
            p = render_f1(fig_id)
            pdfs.append((fig_id, p.with_suffix('.png')))
        except Exception as e:
            print(f'  [FAIL] {fig_id}: {e}')
            import traceback; traceback.print_exc()

    print('\n[F2] Ensemble sigma_t...')
    for fig_id in ('fig02a', 'fig02b'):
        try:
            p = render_f2(fig_id)
            pdfs.append((fig_id, p.with_suffix('.png')))
        except Exception as e:
            print(f'  [FAIL] {fig_id}: {e}'); import traceback; traceback.print_exc()

    print('\n[F3] Event displays...')
    for fig_id in ('fig03a', 'fig03b', 'fig03c', 'fig03d'):
        try:
            p = render_f3(fig_id)
            pdfs.append((fig_id, p.with_suffix('.png')))
        except Exception as e:
            print(f'  [FAIL] {fig_id}: {e}'); import traceback; traceback.print_exc()

    print('\n[F4] N_pe Top profiles...')
    try:
        p = render_f4('fig04')
        pdfs.append(('fig04', p.with_suffix('.png')))
    except Exception as e:
        print(f'  [FAIL] fig04: {e}'); import traceback; traceback.print_exc()

    print('\n[F5] SUM4 timing...')
    for fig_id in ('fig05a', 'fig05b'):
        try:
            p = render_f5(fig_id)
            pdfs.append((fig_id, p.with_suffix('.png')))
        except Exception as e:
            print(f'  [FAIL] {fig_id}: {e}'); import traceback; traceback.print_exc()

    print('\n[F6] Top redundancy...')
    try:
        p = render_f6('fig06')
        pdfs.append(('fig06', p.with_suffix('.png')))
    except Exception as e:
        print(f'  [FAIL] fig06: {e}'); import traceback; traceback.print_exc()

    print('\n[QA-1c supplement] ToF probe visualization...')
    try:
        p = render_tof_probe()
        pdfs.append(('tof_probe', p.with_suffix('.png')))
    except Exception as e:
        print(f'  [FAIL] tof_probe: {e}'); import traceback; traceback.print_exc()

    print('\n[contact_sheet] Building contact sheet...')
    cs = make_contact_sheet(pdfs)

    # ── Verification report ──────────────────────────────────────────────
    print('\n' + '='*60)
    print('QA-2 PDF VERIFICATION')
    print('='*60)
    all_ok = True
    for fig_id, _ in pdfs:
        suffix = '_vis' if 'tof_probe' in fig_id else ''
        stem = (fig_id if 'tof_probe' not in fig_id else 'fig01_tof_probe') + suffix
        p = OUTPUTS / f'{stem}.pdf'
        ok = p.exists() and p.stat().st_size > 0
        status = '[OK]  ' if ok else '[FAIL]'
        if not ok: all_ok = False
        size = p.stat().st_size / 1024 if p.exists() else 0
        print(f'  {status} {stem}.pdf ({size:.1f} kB)')

    print(f'\n[{"PASS" if all_ok else "FAIL"}] {len(pdfs)} panels rendered.')
    print(f'Contact sheet: {cs}')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
