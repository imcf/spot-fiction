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
"""Core processing: transcript CSV → Gaussian density → uint16 BigTIFF."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter

Z_PLANES = 7
STRIP_ROWS = 1024  # ~285 MB float32 per strip with padding


def load_transform(img_dir: Path) -> np.ndarray:
    """Return 3×3 affine matrix (micron → mosaic pixel)."""
    mat = np.loadtxt(img_dir / "micron_to_mosaic_pixel_transform.csv")
    if mat.shape != (3, 3):
        raise ValueError(f"Expected 3×3 transform matrix, got {mat.shape}")
    return mat


def load_manifest(img_dir: Path) -> dict:
    with open(img_dir / "manifest.json") as f:
        return json.load(f)


def load_transcripts(
    data_dir: Path,
    img_dir: Path,
    genes: list[str] | None = None,
) -> pd.DataFrame:
    """
    Load detected_transcripts.csv, apply affine transform, return DataFrame
    with integer pixel columns px, py and z-plane index pz.

    Parameters:
        data_dir: MERSCOPE region directory containing detected_transcripts.csv.
        img_dir:  images/ subdirectory containing the transform CSV.
        genes:    Optional list of gene names to include (None = all).
    """
    manifest = load_manifest(img_dir)
    height = manifest["mosaic_height_pixels"]
    width = manifest["mosaic_width_pixels"]

    print("Loading transcripts…", flush=True)
    usecols = ["global_x", "global_y", "global_z"]
    if genes is not None:
        usecols.append("gene")
    df = pd.read_csv(data_dir / "detected_transcripts.csv", usecols=usecols)
    print(f"  {len(df):,} transcripts", flush=True)

    if genes is not None:
        df = df[df["gene"].isin(set(genes))]
        print(f"  {len(df):,} after gene filter", flush=True)

    mat = load_transform(img_dir)
    mx, my = df["global_x"].to_numpy(), df["global_y"].to_numpy()
    df["px"] = np.round(mat[0, 0] * mx + mat[0, 1] * my + mat[0, 2]).astype(np.int32)
    df["py"] = np.round(mat[1, 0] * mx + mat[1, 1] * my + mat[1, 2]).astype(np.int32)
    df["pz"] = np.round(df["global_z"]).clip(0, Z_PLANES - 1).astype(np.int8)

    n_before = len(df)
    df = df[
        (df.px >= 0) & (df.px < width) & (df.py >= 0) & (df.py < height)
    ]
    dropped = n_before - len(df)
    if dropped:
        print(f"  {dropped:,} out-of-bounds transcripts dropped", flush=True)

    return df[["px", "py", "pz"]].reset_index(drop=True)


def process_z(
    z_df: pd.DataFrame,
    z: int,
    sigma: float,
    out_path: Path,
    height: int,
    width: int,
) -> None:
    """
    Render one z-plane to a uint16 BigTIFF.

    Two-pass strategy to avoid holding the full image in RAM:
      Pass 1 — fill a float32 numpy memmap with Gaussian-blurred spot counts.
      Pass 2 — normalize to uint16, write via tifffile.memmap (single-page BigTIFF).

    Peak disk: ~22 GB (float32 temp) + ~11 GB (uint16 output) per plane.
    The float32 temp file is deleted after each plane.
    """
    px_arr = z_df["px"].to_numpy()
    py_arr = z_df["py"].to_numpy()
    print(f"z={z}: {len(px_arr):,} transcripts → {out_path.name}", flush=True)

    pad = int(np.ceil(3 * sigma))
    tmp_path = out_path.with_suffix(".density.bin")

    # ── Pass 1: float32 density ────────────────────────────────────────────
    tmp = np.memmap(tmp_path, dtype=np.float32, mode="w+", shape=(height, width))

    for r0 in range(0, height, STRIP_ROWS):
        r1 = min(r0 + STRIP_ROWS, height)
        h = r1 - r0

        mask = (py_arr >= r0 - pad) & (py_arr < r1 + pad)
        strip = np.zeros((h + 2 * pad, width), dtype=np.float32)
        if mask.any():
            lpy = py_arr[mask] - r0 + pad
            np.add.at(strip, (lpy, px_arr[mask]), 1.0)
        strip = gaussian_filter(strip, sigma=sigma)
        tmp[r0:r1] = strip[pad: pad + h]

        print(f"  density {r1 / height * 100:4.0f}%", end="\r", flush=True)

    tmp.flush()
    max_val = float(tmp.max())
    scale = 65535.0 / max_val if max_val > 0 else 1.0
    print(f"\n  max={max_val:.5f}", flush=True)

    # ── Pass 2: uint16 BigTIFF ─────────────────────────────────────────────
    out = tifffile.memmap(out_path, shape=(height, width), dtype=np.uint16, bigtiff=True)

    for r0 in range(0, height, STRIP_ROWS):
        r1 = min(r0 + STRIP_ROWS, height)
        out[r0:r1] = (tmp[r0:r1] * scale).astype(np.uint16)
        print(f"  writing  {r1 / height * 100:4.0f}%", end="\r", flush=True)

    out.flush()
    del out
    size_gb = out_path.stat().st_size / 1e9
    print(f"\n  saved {out_path.name} ({size_gb:.1f} GB)", flush=True)

    del tmp
    tmp_path.unlink(missing_ok=True)


def update_manifest(img_dir: Path, name: str, z_list: list[int]) -> None:
    """Append new stain entries to manifest.json mosaic_files list."""
    manifest = load_manifest(img_dir)
    existing = {e["file_name"] for e in manifest["mosaic_files"]}
    added = 0
    for z in sorted(z_list):
        fname = f"mosaic_{name}_z{z}.tif"
        if fname not in existing:
            manifest["mosaic_files"].append({"stain": name, "z": z, "file_name": fname})
            added += 1
    with open(img_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=4)
    print(f"manifest.json: {added} entries added for stain '{name}'", flush=True)
