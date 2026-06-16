#!/usr/bin/env python3.12
"""
EXEC_14 manifest checker — QA-1 validation.

Verifica que:
  1. Cada root_input del CSV existe y SHA-256 coincide
  2. Cada sidecar (.root + .csv + .meta.json) existe en outputs/
  3. meta.json tiene n_events en rango [1500, 3000] (donde aplica)
  4. No hay sidecars huérfanos en outputs/ no declarados en el manifiesto
  5. No hay fig_ids duplicados en el manifiesto

Uso:
  python3.12 manifest_checker.py [--manifest path] [--outputs path]
"""
from __future__ import annotations

import csv
import hashlib
import json
import pathlib
import sys
import argparse
import glob


MANIFEST_PATH = pathlib.Path("/home/reriosto/SHiP/orchestrator/beamer_manifest.csv")
OUTPUTS_PATH  = pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs")

# These fig_ids use composite SHA (multiple files) or "all-31" input → skip SHA check
# fig05a/b root_input has literal commas in {0,400,690} → CSV column shift
# Add them to skip lists; their sidecars are checked by extension-pattern match below.
SKIP_SHA_FIGS = {"fig02a", "fig02b", "fig04", "fig05a", "fig05b", "fig06", "fig01_tof_probe"}

# fig05a/b have literal commas in root_input → CSV column shift → wrong salida_* values
# Standard sidecar check skipped; direct check block handles these.
SKIP_SIDECAR_CHECK_FIGS = {"fig05a", "fig05b", "fig01_tof_probe"}

# These fig_ids have no n_events in meta.json at top level (F2, F3, F4, F5 use nested)
SKIP_NEVENTS_FIGS = {"fig02a", "fig02b", "fig03a", "fig03b", "fig03c", "fig03d",
                     "fig04", "fig05a", "fig05b", "fig06", "fig01_tof_probe"}


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def check_manifest(manifest_path: pathlib.Path, outputs_path: pathlib.Path) -> int:
    errors = []
    warnings = []

    # ── Load manifest ─────────────────────────────────────────────────────────
    if not manifest_path.exists():
        print(f"[FATAL] Manifest not found: {manifest_path}")
        return 1

    rows = []
    with open(manifest_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    print(f"Manifest: {len(rows)} rows")

    # ── Check 4: No duplicate fig_ids ────────────────────────────────────────
    seen_ids = {}
    for row in rows:
        fid = row['fig_id']
        if fid in seen_ids:
            errors.append(f"DUPLICATE fig_id: {fid} appears at rows {seen_ids[fid]} and current")
        else:
            seen_ids[fid] = True
    if len(seen_ids) == len(rows):
        print(f"[OK] No duplicate fig_ids ({len(rows)} unique)")

    declared_figs = set(seen_ids.keys())

    # ── Check 1: root_input exists + SHA matches ──────────────────────────────
    print("\n--- SHA check ---")
    for row in rows:
        fid = row['fig_id']
        if fid in SKIP_SHA_FIGS:
            print(f"  [SKIP] {fid}: composite/all-31 input (SHA not single-file)")
            continue
        root_input = pathlib.Path(row['root_input'])
        expected_sha = row['sha256_input'].strip()
        if not root_input.exists():
            errors.append(f"{fid}: root_input not found: {root_input}")
            continue
        actual_sha = sha256_file(root_input)
        if actual_sha == expected_sha:
            print(f"  [OK]  {fid}: {root_input.name}")
        else:
            errors.append(
                f"{fid}: SHA mismatch for {root_input.name}\n"
                f"    expected: {expected_sha}\n"
                f"    actual:   {actual_sha}"
            )

    # ── Check 2: sidecars exist ───────────────────────────────────────────────
    print("\n--- Sidecar existence check ---")
    for row in rows:
        fid = row['fig_id']
        if fid in SKIP_SIDECAR_CHECK_FIGS:
            print(f"  [SKIP] {fid}: CSV column shift (comma in root_input) — see direct check below")
            continue
        for ext_key in ('salida_sidecar', 'salida_csv', 'salida_meta_json'):
            rel = row.get(ext_key, '').strip()
            if not rel:
                continue
            p = pathlib.Path("/home/reriosto/SHiP/orchestrator") / rel
            if p.exists():
                size_kb = p.stat().st_size / 1024
                print(f"  [OK]  {fid} {ext_key}: {p.name} ({size_kb:.1f} kB)")
            else:
                errors.append(f"{fid}: missing sidecar {ext_key}: {p}")

    # ── Check 3: n_events in meta.json ───────────────────────────────────────
    print("\n--- meta.json n_events check ---")
    for row in rows:
        fid = row['fig_id']
        meta_rel = row.get('salida_meta_json', '').strip()
        if not meta_rel:
            continue
        meta_p = pathlib.Path("/home/reriosto/SHiP/orchestrator") / meta_rel
        # Skip if CSV column shift gave us a non-JSON path (fig05a/b commas in root_input)
        if not meta_rel.endswith('.meta.json'):
            warnings.append(f"{fid}: salida_meta_json column has wrong value '{meta_rel}' "
                            "(CSV comma in root_input field caused column shift)")
            continue
        if not meta_p.exists():
            continue  # Already caught in sidecar check
        with open(meta_p) as f:
            meta = json.load(f)

        if fid in SKIP_NEVENTS_FIGS:
            # Check nested (F5 has positions list; F3 has n_events at top level)
            if "n_events" in meta:
                n = meta["n_events"]
                ok = 1500 <= n <= 3000
                status = "[OK] " if ok else "[ERR]"
                print(f"  {status} {fid}: n_events={n}")
                if not ok:
                    errors.append(f"{fid}: n_events={n} outside [1500,3000]")
            elif "positions" in meta:
                for pos in meta["positions"]:
                    n = pos.get("n_events", -1)
                    x = pos.get("x_mm", "?")
                    ok = 1500 <= n <= 3000
                    status = "[OK] " if ok else "[ERR]"
                    print(f"  {status} {fid} x={x}: n_events={n}")
                    if not ok:
                        errors.append(f"{fid} x={x}: n_events={n} outside [1500,3000]")
            elif "n_events_total" in meta:
                n = meta["n_events_total"]
                print(f"  [INFO] {fid}: n_events_total={n} (pooled F6)")
            elif "n_events_total_all_positions" in meta:
                n = meta["n_events_total_all_positions"]
                print(f"  [INFO] {fid}: n_events_total_all_positions={n} (pooled F2)")
            elif "cases" in meta:
                n = sum(c.get("n_events", 0) for c in meta["cases"])
                print(f"  [INFO] {fid}: n_events via cases={n} ({len(meta['cases'])} cases, QA probe)")
            else:
                warnings.append(f"{fid}: no n_events in meta.json")
        else:
            n = meta.get("n_events", -1)
            ok = 1500 <= n <= 3000
            status = "[OK] " if ok else "[ERR]"
            print(f"  {status} {fid}: n_events={n}")
            if not ok:
                errors.append(f"{fid}: n_events={n} outside [1500,3000]")

    # ── Check 2b: direct existence check for fig05a/b + tof_probe (CSV column/comma workaround) ──
    print("\n--- fig05a/b + fig01_tof_probe direct sidecar check ---")
    for fid in ("fig05a", "fig05b", "fig01_tof_probe"):
        for ext in (".root", ".csv", ".meta.json"):
            p = outputs_path / f"{fid}{ext}"
            if p.exists():
                print(f"  [OK]  {fid}: {p.name} ({p.stat().st_size/1024:.1f} kB)")
            else:
                errors.append(f"{fid}: missing sidecar {p.name}")
    # Also check n_events for fig05a/b
    for fid in ("fig05a", "fig05b"):
        meta_p = outputs_path / f"{fid}.meta.json"
        if meta_p.exists():
            with open(meta_p) as mf:
                meta = json.load(mf)
            for pos in meta.get("positions", []):
                n = pos.get("n_events", -1)
                x = pos.get("x_mm", "?")
                ok = 1500 <= n <= 3000
                print(f"  {'[OK] ' if ok else '[ERR]'} {fid} x={x}: n_events={n}")
                if not ok:
                    errors.append(f"{fid} x={x}: n_events={n} outside [1500,3000]")
    # Check n_events for fig01_tof_probe (uses 'cases' list in meta)
    tof_meta_p = outputs_path / "fig01_tof_probe.meta.json"
    if tof_meta_p.exists():
        with open(tof_meta_p) as mf:
            tof_meta = json.load(mf)
        for case in tof_meta.get("cases", []):
            n = case.get("n_events", -1)
            mat = case.get("material", "?")
            x   = case.get("x_mm", "?")
            ok  = 1500 <= n <= 3000
            print(f"  {'[OK] ' if ok else '[ERR]'} fig01_tof_probe {mat} x={x}: n_events={n}")
            if not ok:
                errors.append(f"fig01_tof_probe {mat} x={x}: n_events={n} outside [1500,3000]")
        ov = tof_meta.get("overall_verdict", "?")
        if ov == "ToF_CONFIRMADO":
            print(f"  [OK]  fig01_tof_probe overall_verdict={ov}")
        else:
            warnings.append(f"fig01_tof_probe overall_verdict={ov} (not ToF_CONFIRMADO)")

    # ── Check 5: orphan sidecars ──────────────────────────────────────────────
    print("\n--- Orphan sidecar check ---")
    outputs_roots = set(p.stem for p in outputs_path.glob("*.root")
                        if not p.stem.startswith('_'))
    declared_roots = set(row['fig_id'] for row in rows)
    orphans = outputs_roots - declared_roots
    if orphans:
        for o in sorted(orphans):
            warnings.append(f"Orphan sidecar in outputs/: {o}.root (not in manifest)")
    else:
        print(f"  [OK]  No orphan .root sidecars in outputs/")

    # ── Report ────────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"MANIFEST CHECK RESULT")
    print('='*60)
    if warnings:
        print(f"\nWARNINGS ({len(warnings)}):")
        for w in warnings:
            print(f"  [WARN] {w}")
    if errors:
        print(f"\nERRORS ({len(errors)}):")
        for e in errors:
            print(f"  [ERR]  {e}")
        print(f"\n[FAIL] {len(errors)} errors, {len(warnings)} warnings")
        return 1
    else:
        print(f"\n[PASS] All checks passed. {len(warnings)} warnings.")
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="EXEC_14 manifest checker")
    parser.add_argument("--manifest", default=str(MANIFEST_PATH))
    parser.add_argument("--outputs",  default=str(OUTPUTS_PATH))
    args = parser.parse_args()
    sys.exit(check_manifest(pathlib.Path(args.manifest), pathlib.Path(args.outputs)))
