"""
Dataset registry — read-only metadata for registered simulation datasets.

Each entry is a dict with mandatory keys:
  material     str   scintillator name
  tau_d_ns     float decay time constant (ns), from MAC file
  tau_r_ns     float rise time constant (ns), from MAC file
  readout      str   "ENDONLY" or "ENDTOP"
  data_dir     Path  path to ROOT files (t0minidaq/ tree, read-only)
  n_channels   int   number of SiPM channels in the geometry
  opsc_code    str   SSLG4 optical scintillator product code
  provenance   dict  source citations for material/tau_d/readout

Sources must be cited as file paths within the repo so auditors can verify.
"""
from __future__ import annotations

import pathlib

_T0 = pathlib.Path("/home/reriosto/SHiP/t0minidaq")
_EJ200_SRC = pathlib.Path("/home/reriosto/SHiP/ej200/src/external/SSLG4")

DATASETS: dict[str, dict] = {
    "exec07_endtop_2000": {
        "material":   "EJ-204",
        "tau_d_ns":   1.8,
        "tau_r_ns":   0.5,
        "readout":    "ENDTOP",
        "data_dir":   _T0 / "sslg4" / "exec07_endtop_2000",
        "n_channels": 86,
        "opsc_code":  "opsc-102",
        "provenance": {
            "material": "EJ-204 (ELJEN standard fast scintillator)",
            "tau_d_ns": "datasheet + common.py:29 TAU_D_NS=1.8",
            "readout":  "max_gid=85 verified in Stage A audit",
        },
    },

    "ej230_endtop": {
        "material":   "EJ-230",
        "tau_d_ns":   1.5,
        "tau_r_ns":   0.5,
        "readout":    "ENDTOP",
        "data_dir":   _T0 / "results_ej230" / "data",
        "n_channels": 86,
        "opsc_code":  "opsc-106",
        "provenance": {
            "material": (
                "EJ-230/Pilot U2/BC-420 — "
                + str(_EJ200_SRC / "src/OrganicScintillatorFactory.cc")
                + ":36 comment"
            ),
            "tau_d_ns": (
                "1.5 ns — "
                + str(_EJ200_SRC / "macros/oscnt/opsc-106.mac")
                + ":12 SCINTILLATIONTIMECONSTANT1"
            ),
            "readout": (
                "EndTop confirmed: max_gid=85, TopSiPMPV:16-85 present "
                "in t0minidaq/results_ej230/logs/run_x0mm.log"
            ),
        },
    },
}


def get(name: str) -> dict:
    """Return dataset metadata dict; raises KeyError if not registered."""
    return DATASETS[name]
