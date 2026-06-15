#!/usr/bin/env python3.12
"""
EXEC_14 batch runner — QA-1 sidecar generation.

Genera todos los 16 paneles del beamer_manifest.csv.
Orden: fig01a-f → fig02a-b → fig03a-d → fig04 → fig05a-b → fig06.
"""
from __future__ import annotations

import os as _os
_os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
_os.environ.setdefault('OMP_NUM_THREADS', '1')

import pathlib
import sys
import traceback

sys.path.insert(0, str(pathlib.Path(__file__).parent))

import f1_sidecar, f2_sidecar, f3_sidecar, f4_sidecar, f5_sidecar, f6_sidecar

OUTPUT   = pathlib.Path("/home/reriosto/SHiP/orchestrator/outputs")
DATA_204 = pathlib.Path("/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614")
DATA_230 = pathlib.Path("/home/reriosto/SHiP/t0minidaq/endonly_mylar_230")
DATA_TOP = pathlib.Path("/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000")

SHA_204 = {
    0:   "63cd8bbee2be58999c55deb10dc6b4fb236fd8c4dcbf9d818ebc2725a99cfdb6",
    400: "403016634e3192b6a7d349260eb06b18ed4c222d2114d2128efded7a7e3479ec",
    690: "61435a5c9e69350aed0ca07612483bb3c78b879b1b5edd3fe7215e2e67037d55",
}
SHA_230 = {
    0:   "0d30b51f40c9ffe46463fc1627cc87faca8900342340c176998e14069050c930",
    400: "1054527b17854828f663fecc31456723157cda7999785bae652ba51017fc10b6",
    690: "c67b636a1e0cf916b2a51dc7b10a9cda4784cfb8807cb2c7aded2bcf06eaad97",
}
SHA_TOP_MANIFEST = "b67c4f7504c94e546b21f577ca686ee6ca652d46c8ed33b569645347e0cf390b"


def _run(label, fn, **kw):
    print(f"\n{'='*60}")
    print(f"START {label}")
    print('='*60)
    try:
        fn(**kw)
        print(f"[OK] {label}")
        return True
    except Exception:
        traceback.print_exc()
        print(f"[FAIL] {label}")
        return False


def main():
    results = {}

    # ── F1 — scintillation emission time profiles ─────────────────────────────
    for fig_id, material, x_mm, sha_dict, ddir, dataset in [
        ("fig01a", "EJ-204", 0,   SHA_204, DATA_204, "endonly_mylar_20260614"),
        ("fig01b", "EJ-204", 400, SHA_204, DATA_204, "endonly_mylar_20260614"),
        ("fig01c", "EJ-204", 690, SHA_204, DATA_204, "endonly_mylar_20260614"),
        ("fig01d", "EJ-230", 0,   SHA_230, DATA_230, "endonly_mylar_230"),
        ("fig01e", "EJ-230", 400, SHA_230, DATA_230, "endonly_mylar_230"),
        ("fig01f", "EJ-230", 690, SHA_230, DATA_230, "endonly_mylar_230"),
    ]:
        results[fig_id] = _run(fig_id, f1_sidecar.generate_f1_sidecar,
            fig_id=fig_id, material=material, dataset=dataset,
            root_input=ddir / f"photon_hits_x{x_mm}mm.root",
            sha256_expected=sha_dict[x_mm], x_mm=x_mm, output_dir=OUTPUT)

    # ── F2 — sigma_t(N_pe) ensemble ──────────────────────────────────────────
    results["fig02a"] = _run("fig02a", f2_sidecar.generate_f2_sidecar,
        fig_id="fig02a", material="EJ-204", dataset="endonly_mylar_20260614",
        data_dir=DATA_204,
        sha256_manifest="8b8c93099eeae7d9e9e3e3a25cd3400890fbcdeee5e041bda5890c5fae77a8b2",
        output_dir=OUTPUT)

    results["fig02b"] = _run("fig02b", f2_sidecar.generate_f2_sidecar,
        fig_id="fig02b", material="EJ-230", dataset="endonly_mylar_230",
        data_dir=DATA_230,
        sha256_manifest="cb9d3a7abbb40a06dd70c14b624f8bb4f193503725169afb8c1884588774dcaf",
        output_dir=OUTPUT)

    # ── F3 — event displays ───────────────────────────────────────────────────
    for fig_id, material, x_mm, sha_dict, ddir, dataset in [
        ("fig03a", "EJ-204", 0,   SHA_204, DATA_204, "endonly_mylar_20260614"),
        ("fig03b", "EJ-204", 690, SHA_204, DATA_204, "endonly_mylar_20260614"),
        ("fig03c", "EJ-230", 0,   SHA_230, DATA_230, "endonly_mylar_230"),
        ("fig03d", "EJ-230", 690, SHA_230, DATA_230, "endonly_mylar_230"),
    ]:
        results[fig_id] = _run(fig_id, f3_sidecar.generate_f3_sidecar,
            fig_id=fig_id, material=material, dataset=dataset,
            root_input=ddir / f"photon_hits_x{x_mm}mm.root",
            sha256_expected=sha_dict[x_mm], x_mm_gun=x_mm, output_dir=OUTPUT)

    # ── F4 — Top T4 vs T20 profiles ──────────────────────────────────────────
    results["fig04"] = _run("fig04", f4_sidecar.generate_f4_sidecar,
        fig_id="fig04", material="EJ-204", dataset="exec07_endtop_2000",
        data_dir=DATA_TOP, sha256_manifest=SHA_TOP_MANIFEST, output_dir=OUTPUT)

    # ── F5 — SUM4 sigma_t three positions ────────────────────────────────────
    results["fig05a"] = _run("fig05a", f5_sidecar.generate_f5_sidecar,
        fig_id="fig05a", material="EJ-204", dataset="endonly_mylar_20260614",
        data_dir=DATA_204, sha256_per_position=SHA_204, output_dir=OUTPUT)

    results["fig05b"] = _run("fig05b", f5_sidecar.generate_f5_sidecar,
        fig_id="fig05b", material="EJ-230", dataset="endonly_mylar_230",
        data_dir=DATA_230, sha256_per_position=SHA_230, output_dir=OUTPUT)

    # ── F6 — Top redundancy ───────────────────────────────────────────────────
    results["fig06"] = _run("fig06", f6_sidecar.generate_f6_sidecar,
        fig_id="fig06", material="EJ-204", dataset="exec07_endtop_2000",
        data_dir=DATA_TOP, sha256_manifest=SHA_TOP_MANIFEST, output_dir=OUTPUT)

    # Summary
    print(f"\n{'='*60}")
    print("BATCH SUMMARY")
    print('='*60)
    n_ok = sum(1 for v in results.values() if v)
    n_fail = sum(1 for v in results.values() if not v)
    for fig_id, ok in results.items():
        print(f"  {'OK  ' if ok else 'FAIL'} {fig_id}")
    print(f"\n{n_ok}/{len(results)} sidecars OK, {n_fail} failed.")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
