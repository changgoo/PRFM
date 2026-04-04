"""
Download the PHANGS megatable v4.0 public release.

Usage
-----
# Download all annulus files (default, ~10 MB):
python scripts/download_phangs.py

# Download a specific aperture type:
python scripts/download_phangs.py --aperture hexagon

# Download all three aperture types (~200 MB total):
python scripts/download_phangs.py --aperture all

# Download a subset of galaxies:
python scripts/download_phangs.py --galaxies NGC0628 NGC1097 IC1954

# Force re-download even if files already exist:
python scripts/download_phangs.py --force

# Custom output directory:
python scripts/download_phangs.py --dest data/my_phangs
"""

import argparse
from pathlib import Path

# Add repo root to path so we can import prfm without installing
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from prfm.phangs import download, list_galaxies


def main():
    parser = argparse.ArgumentParser(
        description="Download PHANGS megatable ECSV files (Sun+2022 v4.0)"
    )
    parser.add_argument(
        "--dest",
        default="data/phangs_megatable",
        help="Destination directory (default: data/phangs_megatable)",
    )
    parser.add_argument(
        "--aperture",
        default="annulus",
        choices=["annulus", "gauss", "hexagon", "all"],
        help="Aperture type to download (default: annulus)",
    )
    parser.add_argument(
        "--galaxies",
        nargs="+",
        metavar="GALAXY",
        default=None,
        help="Subset of galaxy names (default: all 90)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if files already exist",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Print available galaxy names and exit",
    )
    args = parser.parse_args()

    if args.list:
        for g in list_galaxies():
            print(g)
        return

    dest = Path(args.dest)
    print(f"Downloading to: {dest.resolve()}")

    paths = download(
        dest_dir=dest,
        aperture=args.aperture,
        galaxies=args.galaxies,
        force=args.force,
    )

    print(f"Done. {len(paths)} file(s) in {dest}")


if __name__ == "__main__":
    main()
