#!/usr/bin/env python
"""Convert a KDE-Sobol design CSV into a suite YAML for generate_slurm.py.

The design CSV columns
    Sigma_gas, Sigma_star, H_star, Omega, qshear
are mapped to TIGRESS-NCR athinput parameters
    problem/surf, problem/SurfS, problem/zstar, problem/Omega, problem/qshear
and rho_dm is held fixed at the R8 fiducial (0.0064 Msun/pc^3).

Each CSV row becomes one entry in the `models:` list in the output YAML, with
`suffix` = "rowNNNN" so the resulting job IDs are unique.

Usage:
    python project/scripts/csv_to_slurm_yaml.py project/output/design_Sgas10.0_n0032.csv
    python project/scripts/csv_to_slurm_yaml.py <csv> --base R8_8pc --output <yml>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]

# CSV column → TIGRESS-NCR athinput key
PARAM_MAP: dict[str, str] = {
    "Sigma_gas":  "problem/surf",
    "Sigma_star": "problem/SurfS",
    "H_star":     "problem/zstar",
    "Omega":      "problem/Omega",
    "qshear":     "problem/qshear",
}

# Per-base domain defaults (matching athinput.<base> boxes/resolution).
# `dx` in pc; `box_size` in pc; `grid_size` = ranks per axis (subdomain layout).
DOMAIN_DEFAULTS: dict[str, dict] = {
    "R8_8pc":    {"dx": 8,  "box_size": [1024, 1024, 4096], "grid_size": [32, 32, 32]},
    "R8s_4pc":   {"dx": 4,  "box_size": [1024, 1024, 4096], "grid_size": [64, 64, 64]},
    "LGR4_4pc":  {"dx": 4,  "box_size": [512,  512,  4096], "grid_size": [32, 32, 32]},
    "LGR8_8pc":  {"dx": 8,  "box_size": [1024, 1024, 4096], "grid_size": [32, 32, 32]},
}

# Fixed parameters applied to every model. See project/doc/paper1-suite/sections/tigress_mapping.tex
# for the physical mapping (rhodm, Z_gas, Z_dust). The other keys are R8_8pc
# TIGRESS-PHANGS-suite defaults matching Kim et al. runs.
FIXED_PARAMS: dict = {
    # Physics defaults
    "problem/beta":       10,       # plasma beta
    "problem/rhodm":      0.0064,   # M_sun/pc^3, R8 fiducial
    "problem/Z_gas":      1.0,      # solar metallicity
    "problem/Z_dust":     1.0,      # solar dust-to-gas
    # Radiation-pressure defaults
    "radps/eps_extinct":  1.0e-8,
    "radps/xymaxPP":      1024,
    # Output cadence
    "output2/dt":         10.0,
    "output4/dt":         50.0,
}

# Omega in the CSV is in km/s/kpc (PHANGS convention). TIGRESS-NCR uses
# code units where Omega has units of km/s/pc (i.e., velocity/length with
# length in pc). Convert by dividing by 1000.
OMEGA_KMS_PC_PER_KPC = 1.0e-3


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", type=Path,
                   help="Design CSV produced by run_sampling.py")
    p.add_argument("--base", default="R8_8pc",
                   help="TIGRESS-NCR base model / athinput stem (default: R8_8pc)")
    p.add_argument("--output", "-o", type=Path, default=None,
                   help="Output YAML path (default: derived from CSV name)")
    p.add_argument("--machine", type=Path, default=None,
                   help="Optional machine YAML path; if given, its contents are"
                        " merged into the output so a single YAML fully"
                        " specifies the suite. Otherwise pass the machine YAML"
                        " as a second argument to generate_slurm.py.")
    p.add_argument("--decimals", type=int, default=4,
                   help="Rounding precision for physical values in YAML output"
                        " (default: 4)")
    return p.parse_args()


def _round(x: float, ndigits: int) -> float:
    """Round to ndigits and cast to a plain Python float."""
    return float(round(float(x), ndigits))


def build_config(csv_path: Path, base: str, decimals: int,
                 machine_yaml: Path | None) -> dict:
    df = pd.read_csv(csv_path)
    missing = [c for c in PARAM_MAP if c not in df.columns]
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    n_rows = len(df)
    row_width = max(4, len(str(n_rows - 1)))  # 0000..NNNN

    # Per-row varying parameters go in `extra_overrides`. Row index (suffix)
    # keeps job IDs compact (`_row0000` … `_rowNNNN`) so the varying numbers
    # don't clutter the SLURM script names or run directories.
    models = []
    for i, row in df.iterrows():
        varying: dict = {}
        for col, key in PARAM_MAP.items():
            val = float(row[col])
            if col == "Omega":
                # CSV: km/s/kpc (PHANGS convention); code: km/s/pc.
                val *= OMEGA_KMS_PC_PER_KPC
                varying[key] = _round(val, decimals + 3)  # extra digits post-shift
            else:
                varying[key] = _round(val, decimals)
        models.append({
            "base":            base,
            "suffix":          f"row{i:0{row_width}d}",
            "extra_overrides": varying,
        })

    config: dict = {
        "suite": {
            "run_base": "/scratch/gpfs/changgoo/TIGRESS-PHANGS",
            "source_csv": str(csv_path.relative_to(ROOT))
                          if csv_path.is_absolute() else str(csv_path),
        },
    }
    if base in DOMAIN_DEFAULTS:
        config["domain"] = dict(DOMAIN_DEFAULTS[base])
    config["fixed_overrides"] = dict(FIXED_PARAMS)
    config["models"] = models

    if machine_yaml is not None:
        with open(machine_yaml) as f:
            machine = yaml.safe_load(f) or {}
        # generate_slurm.py already knows how to merge machine YAML on top; we
        # bundle it in for convenience so a single YAML is self-contained.
        # Users can still keep them separate if preferred.
        _deep_update(config, machine)

    return config


def _deep_update(base: dict, override: dict) -> None:
    """Recursively merge override into base in place."""
    for k, v in override.items():
        if k in base and isinstance(base[k], dict) and isinstance(v, dict):
            _deep_update(base[k], v)
        else:
            base[k] = v


def default_output_path(csv_path: Path) -> Path:
    stem = csv_path.stem  # e.g. design_Sgas10.0_n0032
    return ROOT / "project" / "suites" / f"{stem}.yml"


def main() -> None:
    args = parse_args()
    if not args.csv.exists():
        print(f"ERROR: CSV not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    out_path = args.output or default_output_path(args.csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    config = build_config(args.csv, args.base, args.decimals, args.machine)

    # Custom dumper: keep top-level mappings block-style but render
    # short numeric lists (box_size, grid_size) in flow style so the file
    # matches the hand-authored convention.
    class _FlowList(list):
        pass

    def _flow_list_representer(dumper, data):
        return dumper.represent_sequence("tag:yaml.org,2002:seq",
                                         data, flow_style=True)

    yaml.SafeDumper.add_representer(_FlowList, _flow_list_representer)

    if "domain" in config:
        d = config["domain"]
        for k in ("box_size", "grid_size"):
            if k in d and isinstance(d[k], list):
                d[k] = _FlowList(d[k])

    with open(out_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False, default_flow_style=False,
                       width=100, indent=2)

    n_models = len(config["models"])
    print(f"Wrote {n_models} models → {out_path}")


if __name__ == "__main__":
    main()
