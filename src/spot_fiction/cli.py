# spot-fiction — synthesize transcript density channels for MERSCOPE/VIZGEN datasets
# Copyright (C) 2024  Laurent Guerard <laurent.guerard@unibas.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""CLI entry point for spot-fiction."""

import argparse
from pathlib import Path

from .core import Z_PLANES, load_manifest, load_transcripts, process_z, update_manifest


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="spot-fiction",
        description=(
            "Synthesize a transcript density channel for a MERSCOPE region. "
            "Reads detected_transcripts.csv and the affine transform from "
            "images/micron_to_mosaic_pixel_transform.csv, renders Gaussian-blurred "
            "spot counts, and writes mosaic_<NAME>_z{0-6}.tif to images/."
        ),
    )
    p.add_argument(
        "data_dir",
        type=Path,
        help="MERSCOPE region directory (contains detected_transcripts.csv and images/)",
    )
    p.add_argument(
        "--sigma",
        type=float,
        default=10.0,
        metavar="PX",
        help="Gaussian sigma in pixels (default: 10 ≈ 1.08 µm at 0.108 µm/px)",
    )
    p.add_argument(
        "--genes",
        nargs="*",
        default=None,
        metavar="GENE",
        help="Include only these gene names (default: all genes)",
    )
    p.add_argument(
        "--name",
        default="Transcripts",
        metavar="NAME",
        help="Stain name used in output filenames (default: Transcripts)",
    )
    p.add_argument(
        "--z",
        type=int,
        default=None,
        metavar="Z",
        help=f"Process only this z-plane index 0-{Z_PLANES - 1} (default: all)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    data_dir: Path = args.data_dir.resolve()
    img_dir = data_dir / "images"

    if not (data_dir / "detected_transcripts.csv").exists():
        raise FileNotFoundError(f"detected_transcripts.csv not found in {data_dir}")
    if not img_dir.is_dir():
        raise FileNotFoundError(f"images/ directory not found in {data_dir}")

    manifest = load_manifest(img_dir)
    height = manifest["mosaic_height_pixels"]
    width = manifest["mosaic_width_pixels"]

    z_list = [args.z] if args.z is not None else list(range(Z_PLANES))

    print(f"data_dir : {data_dir}")
    print(f"sigma    : {args.sigma} px ({args.sigma * 0.108:.2f} µm)")
    print(f"name     : {args.name}")
    print(f"z-planes : {z_list}")
    print(f"image    : {height} × {width} px")
    if args.genes:
        print(f"genes    : {args.genes}")

    df = load_transcripts(data_dir, img_dir, args.genes)

    for z in z_list:
        z_df = df[df["pz"] == z]
        out_path = img_dir / f"mosaic_{args.name}_z{z}.tif"
        process_z(z_df, z, args.sigma, out_path, height, width)

    update_manifest(img_dir, args.name, z_list)
    print("Done.")
