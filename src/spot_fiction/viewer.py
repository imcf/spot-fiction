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
"""Napari viewer: density map + transcript overlay."""

import argparse
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import tifffile
import zarr

from .core import Z_PLANES, load_manifest, load_transform

# Palette for gene coloring (cycles if more genes than colors)
_COLORS = [
    "cyan", "magenta", "yellow", "lime", "orange",
    "red", "deepskyblue", "hotpink", "chartreuse", "gold",
]


def _load_pyramid_lazy(tif_path: Path) -> list[da.Array]:
    """Return pyramid levels as dask arrays (lazy — tiles loaded on demand)."""
    store = tifffile.imread(tif_path, aszarr=True)
    zobj = zarr.open(store, mode="r")

    if isinstance(zobj, zarr.Group):
        # Pyramidal: levels stored as '0', '1', '2', ...
        keys = sorted(zobj.keys(), key=int)
        return [da.from_zarr(store, component=k) for k in keys]
    else:
        # Single level
        return [da.from_zarr(store)]


def _load_transcripts(
    data_dir: Path,
    img_dir: Path,
    z: int,
    genes: list[str] | None,
    max_points: int,
) -> pd.DataFrame:
    """
    Load detected_transcripts.csv for a single z-plane.
    Returns DataFrame with px, py, gene columns.
    Subsamples to max_points if needed.
    """
    usecols = ["global_x", "global_y", "global_z", "gene"]
    print("Loading transcripts…", flush=True)
    df = pd.read_csv(data_dir / "detected_transcripts.csv", usecols=usecols)

    # Filter z-plane first (cheap before transform)
    df = df[np.round(df["global_z"]).clip(0, Z_PLANES - 1).astype(int) == z]

    if genes is not None:
        df = df[df["gene"].isin(set(genes))]

    mat = load_transform(img_dir)
    manifest = load_manifest(img_dir)
    height, width = manifest["mosaic_height_pixels"], manifest["mosaic_width_pixels"]

    mx, my = df["global_x"].to_numpy(), df["global_y"].to_numpy()
    df = df.copy()
    df["px"] = np.round(mat[0, 0] * mx + mat[0, 1] * my + mat[0, 2]).astype(np.int32)
    df["py"] = np.round(mat[1, 0] * mx + mat[1, 1] * my + mat[1, 2]).astype(np.int32)
    df = df[(df.px >= 0) & (df.px < width) & (df.py >= 0) & (df.py < height)]

    print(f"  {len(df):,} transcripts at z={z}", flush=True)

    if len(df) > max_points:
        df = df.sample(max_points, random_state=42)
        print(f"  subsampled to {max_points:,}", flush=True)

    return df[["px", "py", "gene"]].reset_index(drop=True)


def launch(
    data_dir: Path,
    z: int = 3,
    stain: str = "Transcripts",
    genes: list[str] | None = None,
    max_points: int = 500_000,
) -> None:
    """
    Open napari with the density image and transcript point overlay.

    Parameters:
        data_dir:   MERSCOPE region directory.
        z:          Z-plane index (0-6).
        stain:      Stain name used in mosaic filename (default: Transcripts).
        genes:      Gene names to show as points (None = all, subsampled).
        max_points: Max transcripts to display (subsampled randomly if exceeded).
    """
    import napari

    img_dir = data_dir / "images"
    manifest = load_manifest(img_dir)
    px_um: float = manifest.get("microns_per_pixel", 0.108)
    scale = (px_um, px_um)  # (row, col) scale in µm

    # ── Image ──────────────────────────────────────────────────────────────
    tif_path = img_dir / f"mosaic_{stain}_z{z}.tif"
    if not tif_path.exists():
        raise FileNotFoundError(
            f"{tif_path} not found — run spot-fiction first to generate it"
        )

    print(f"Opening {tif_path.name}…", flush=True)
    pyramid = _load_pyramid_lazy(tif_path)
    print(f"  {len(pyramid)} resolution levels", flush=True)

    # ── Transcripts ─────────────────────────────────────────────────────────
    df = _load_transcripts(data_dir, img_dir, z, genes, max_points)

    # napari points are (row, col) = (y, x) = (py, px)
    all_coords = df[["py", "px"]].to_numpy(dtype=float)
    gene_list = df["gene"].to_numpy()
    unique_genes = sorted(df["gene"].unique())

    # ── Viewer ──────────────────────────────────────────────────────────────
    viewer = napari.Viewer(title=f"spot-fiction · {stain} z={z}")

    viewer.add_image(
        pyramid,
        multiscale=True,
        name=stain,
        colormap="inferno",
        scale=scale,
        blending="additive",
    )

    if len(unique_genes) == 1 or genes is None:
        # Single layer — all points same color
        viewer.add_points(
            all_coords,
            name="transcripts",
            size=2,
            face_color="cyan",
            border_color="transparent",
            opacity=0.6,
            scale=scale,
        )
    else:
        # One layer per gene for independent toggle + distinct colors
        for i, gene in enumerate(unique_genes):
            mask = gene_list == gene
            color = _COLORS[i % len(_COLORS)]
            viewer.add_points(
                all_coords[mask],
                name=gene,
                size=2,
                face_color=color,
                border_color="transparent",
                opacity=0.7,
                scale=scale,
            )

    print(
        f"napari ready — {len(unique_genes)} gene(s), {len(df):,} points",
        flush=True,
    )
    napari.run()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="spot-fiction-view",
        description="Open napari with the density map and transcript overlay.",
    )
    p.add_argument("data_dir", type=Path, help="MERSCOPE region directory")
    p.add_argument("--z", type=int, default=3, help="Z-plane index (default: 3)")
    p.add_argument(
        "--stain", default="Transcripts",
        help="Stain name in mosaic filename (default: Transcripts)",
    )
    p.add_argument(
        "--genes", nargs="*", default=None, metavar="GENE",
        help="Show only these genes as colored layers (default: all, subsampled)",
    )
    p.add_argument(
        "--max-points", type=int, default=500_000, metavar="N",
        help="Max transcript points to display (default: 500000)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    launch(
        data_dir=args.data_dir.resolve(),
        z=args.z,
        stain=args.stain,
        genes=args.genes,
        max_points=args.max_points,
    )
