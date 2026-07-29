#!/usr/bin/env python
"""Generate an augmentation block for an existing KDE-Sobol design.

The submitted suite was drawn from the *unbroadened* PHANGS prior, so its
samples are confined to the observed distribution. This script adds NEW samples
drawn from an *augmented* prior (selected marginals broadened by a per-field dex
amount; default H_star and Omega by 0.3 dex per side, ~x2), appended to the
existing suite as new rows. The existing simulations remain valid; the new block
extends coverage into the broadened tails.

The block uses an independent Sobol seed so its points do not overlap the
existing design, and row suffixes continue from the existing count
(row0032, row0033, ...), so csv_to_slurm_yaml.py emits only the new SLURM
scripts. See project/doc/prior-augmentation.md for the rationale.

Usage:
    python project/scripts/augment_design.py project/output/design_Sgas10.0_n0032.csv
    python project/scripts/augment_design.py <design.csv> \
        --n-augment 32 --augment "H_star=0.3,Omega=0.3" --seed 43
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from prfm.phangs import load_configured_phangs  # noqa: E402
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig  # noqa: E402

OUTPUT_DIR = ROOT / "project" / "output"
CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"

DESIGN_FIELDS = ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]
DEFAULT_AUGMENT = {"H_star": 0.3, "Omega": 0.3}


# --------------------------------------------------------------------------
# Pure helpers (unit-tested; no PHANGS data required)
# --------------------------------------------------------------------------
def parse_augment(spec: str) -> dict[str, float]:
    """Parse "H_star=0.3,Omega=0.3" into {"H_star": 0.3, "Omega": 0.3}."""
    out: dict[str, float] = {}
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        key, val = part.split("=")
        out[key.strip()] = float(val)
    return out


def parse_sigma_gas_from_stem(stem: str) -> float | None:
    """Extract the target Sigma_gas from a stem like 'design_Sgas10.0_n0032'."""
    m = re.search(r"Sgas([0-9.]+)", stem)
    return float(m.group(1)) if m else None


def add_block_suffixes(block: pd.DataFrame, offset: int, width: int = 4) -> pd.DataFrame:
    """Attach a `suffix` column row{offset+i} so job IDs continue the suite."""
    block = block.reset_index(drop=True).copy()
    block["suffix"] = [f"row{offset + i:0{width}d}" for i in range(len(block))]
    return block


# --------------------------------------------------------------------------
# Data-dependent generation
# --------------------------------------------------------------------------
def select_reference(
    target_sigma_gas: float, delta: float, config_path: Path, aux_size: int
):
    """Load PHANGS and return (designer, band-selected reference table)."""
    result = load_configured_phangs(str(config_path), base_dir=str(ROOT))
    cfg = SamplingConfig(
        target_sigma_gas=target_sigma_gas,
        delta_sigma_gas=delta,
        kde_aux_sample_size=aux_size,
    )
    designer = PHANGSSamplingDesigner(result["table"], config=cfg)
    designer = designer.with_config(target_sigma_gas=target_sigma_gas)
    return designer, designer.select_reference_pixels()


def generate_block(
    designer, reference, n_augment: int, augment_dex: dict[str, float], seed: int
) -> pd.DataFrame:
    """Draw n_augment samples from the augmented prior (independent seed)."""
    return designer.synthesize_kde_sobol(
        reference, n_samples=n_augment, seed=seed, augment_dex=augment_dex
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("design_csv", type=Path, help="Existing physical design CSV")
    p.add_argument("--n-augment", type=int, default=32, help="New samples (default 32)")
    p.add_argument(
        "--augment",
        default="H_star=0.3,Omega=0.3",
        help="Per-field dex broadening (default 'H_star=0.3,Omega=0.3')",
    )
    p.add_argument(
        "--seed", type=int, default=43,
        help="Sobol seed for the block; must differ from the base design (default 43)",
    )
    p.add_argument("--sigma-gas", type=float, default=None,
                   help="Target Sigma_gas (default: parsed from CSV name)")
    p.add_argument("--delta", type=float, default=0.3, help="Band half-width (dex)")
    p.add_argument("--config", type=Path, default=CONFIG_PATH)
    p.add_argument("--aux-size", type=int, default=100_000)
    p.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.design_csv.exists():
        print(f"ERROR: design CSV not found: {args.design_csv}", file=sys.stderr)
        sys.exit(1)

    existing = pd.read_csv(args.design_csv)
    stem = args.design_csv.stem
    sigma_gas = args.sigma_gas or parse_sigma_gas_from_stem(stem)
    if sigma_gas is None:
        print("ERROR: could not determine target Sigma_gas; pass --sigma-gas",
              file=sys.stderr)
        sys.exit(1)

    augment_dex = parse_augment(args.augment)
    if args.seed == 42:
        print("WARNING: seed 42 matches the base design seed; use a different "
              "seed so the block does not overlap the existing points.",
              file=sys.stderr)

    print(f"Target Sigma_gas = {sigma_gas} (delta={args.delta}); "
          f"augment = {augment_dex}; seed = {args.seed}")
    designer, reference = select_reference(
        sigma_gas, args.delta, args.config, args.aux_size
    )
    print(f"  reference pixels: {len(reference)}")

    block = generate_block(designer, reference, args.n_augment, augment_dex, args.seed)
    block = add_block_suffixes(block[DESIGN_FIELDS], offset=len(existing))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    block_path = args.output_dir / f"augment_{stem}_n{args.n_augment:04d}.csv"
    block.to_csv(block_path, index=False)

    # Combined CSV (existing rows + block) for plotting / record-keeping.
    existing_out = existing.copy()
    existing_out["suffix"] = [f"row{i:04d}" for i in range(len(existing_out))]
    combined = pd.concat([existing_out, block], ignore_index=True)
    combined_path = args.output_dir / f"{stem}_plus_aug{args.n_augment:04d}.csv"
    combined.to_csv(combined_path, index=False)

    print(f"Wrote {len(block)} augmentation rows -> {block_path}")
    print(f"Wrote combined {len(combined)}-row design -> {combined_path}")


if __name__ == "__main__":
    main()
