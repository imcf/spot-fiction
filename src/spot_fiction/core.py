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
"""Core processing: transcript CSV → Gaussian density → pyramidal OME-TIFF."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter

Z_PLANES = 7
STRIP_ROWS = 1024       # rows processed per strip during density + pyramid build
TILE_SIZE = 512         # tile dimensions for the output OME-TIFF
RAM_THRESHOLD = 2 << 30  # 2 GB: pyramid levels larger than this use a temp memmap


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


def _downsample_strip(src: np.ndarray, r0: int, r1: int, pw: int) -> np.ndarray:
    """2× block-mean downsample a horizontal strip (rows r0:r1) of src."""
    pw_even = pw - pw % 2
    strip = src[r0:r1, :pw_even].astype(np.float32)
    h, w = strip.shape
    h2, w2 = h // 2, w // 2
    return (
        strip.reshape(h2, 2, w2, 2)
        .mean(axis=(1, 3))
        .clip(0, 65535)
        .astype(np.uint16)
    )


def _build_pyramid(
    level0: np.ndarray,
    tmp_dir: Path,
) -> tuple[list[np.ndarray], list[Path]]:
    """
    Build all sub-resolution levels from a uint16 level-0 array/memmap.

    Levels larger than RAM_THRESHOLD bytes are backed by temp memmaps so
    they don't consume physical RAM. Returns (levels_list, tmp_paths).
    """
    levels: list[np.ndarray] = [level0]
    tmp_paths: list[Path] = []
    prev = level0

    while True:
        ph, pw = prev.shape
        if ph < 256 and pw < 256:
            break

        ph2, pw2 = ph // 2, pw // 2
        byte_size = ph2 * pw2 * 2  # uint16

        if byte_size > RAM_THRESHOLD:
            tmp_path = tmp_dir / f"tmp_pyramid_l{len(levels)}.bin"
            level_n: np.ndarray = np.memmap(
                tmp_path, dtype=np.uint16, mode="w+", shape=(ph2, pw2)
            )
            tmp_paths.append(tmp_path)
        else:
            level_n = np.empty((ph2, pw2), dtype=np.uint16)

        ph_even = ph - ph % 2
        for r0 in range(0, ph_even, STRIP_ROWS * 2):
            r1 = min(r0 + STRIP_ROWS * 2, ph_even)
            block = _downsample_strip(prev, r0, r1, pw)
            level_n[r0 // 2: r1 // 2, :] = block

        if isinstance(level_n, np.memmap):
            level_n.flush()

        levels.append(level_n)
        prev = level_n
        print(
            f"  pyramid l{len(levels) - 1}: {ph2}×{pw2}", end="\r", flush=True
        )

    n = len(levels) - 1
    print(f"\n  {n} sub-resolution level{'s' if n != 1 else ''}", flush=True)
    return levels, tmp_paths


def _write_ome_pyramid(
    out_path: Path,
    levels: list[np.ndarray],
    channel_name: str,
    pixel_size_um: float,
) -> None:
    """Write pyramid levels as a tiled, LZW-compressed OME-TIFF."""
    n_sublevels = len(levels) - 1
    metadata = {
        "axes": "YX",
        "Channel": {"Name": channel_name},
        "PhysicalSizeX": pixel_size_um,
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": pixel_size_um,
        "PhysicalSizeYUnit": "µm",
    }
    write_opts = dict(
        tile=(TILE_SIZE, TILE_SIZE),
        compression="lzw",
        photometric="minisblack",
    )
    print(f"  writing OME-TIFF ({n_sublevels + 1} levels)…", flush=True)
    with tifffile.TiffWriter(out_path, bigtiff=True, ome=True) as tif:
        tif.write(levels[0], subifds=n_sublevels, metadata=metadata, **write_opts)
        for i, level in enumerate(levels[1:], 1):
            tif.write(level, subfiletype=1, metadata=None, **write_opts)
            print(f"  level {i}/{n_sublevels}", end="\r", flush=True)
    print(flush=True)


def process_z(
    z_df: pd.DataFrame,
    z: int,
    sigma: float,
    out_path: Path,
    height: int,
    width: int,
    channel_name: str = "Transcripts",
    pixel_size_um: float = 0.108,
) -> None:
    """
    Render one z-plane to a pyramidal OME-TIFF.

    Strategy (all large arrays are disk-backed memmaps):
      Pass 1 — float32 Gaussian density memmap (~22 GB).
      Pass 2 — normalize → uint16 level-0 memmap (~11 GB).
      Pass 3 — build sub-resolution pyramid; large levels use temp memmaps.
      Pass 4 — write all levels as tiled, LZW-compressed OME-TIFF.
      Cleanup — all temp files deleted.

    Peak concurrent disk: ~33 GB (float32 + uint16 level 0) during pass 2.
    Final OME-TIFF is ~15 GB (pyramid adds ~33% over the full-res plane).
    """
    px_arr = z_df["px"].to_numpy()
    py_arr = z_df["py"].to_numpy()
    print(f"z={z}: {len(px_arr):,} transcripts → {out_path.name}", flush=True)

    pad = int(np.ceil(3 * sigma))
    density_path = out_path.with_suffix(".density.bin")
    level0_path = out_path.with_suffix(".l0.bin")

    # ── Pass 1: float32 Gaussian density ──────────────────────────────────
    density = np.memmap(density_path, dtype=np.float32, mode="w+", shape=(height, width))

    for r0 in range(0, height, STRIP_ROWS):
        r1 = min(r0 + STRIP_ROWS, height)
        h = r1 - r0
        mask = (py_arr >= r0 - pad) & (py_arr < r1 + pad)
        strip = np.zeros((h + 2 * pad, width), dtype=np.float32)
        if mask.any():
            lpy = py_arr[mask] - r0 + pad
            np.add.at(strip, (lpy, px_arr[mask]), 1.0)
        strip = gaussian_filter(strip, sigma=sigma)
        density[r0:r1] = strip[pad: pad + h]
        print(f"  density {r1 / height * 100:4.0f}%", end="\r", flush=True)

    density.flush()
    max_val = float(density.max())
    scale = 65535.0 / max_val if max_val > 0 else 1.0
    print(f"\n  max={max_val:.5f}", flush=True)

    # ── Pass 2: normalize → uint16 level-0 memmap ─────────────────────────
    level0 = np.memmap(level0_path, dtype=np.uint16, mode="w+", shape=(height, width))

    for r0 in range(0, height, STRIP_ROWS):
        r1 = min(r0 + STRIP_ROWS, height)
        level0[r0:r1] = (density[r0:r1] * scale).astype(np.uint16)
        print(f"  normalize {r1 / height * 100:4.0f}%", end="\r", flush=True)

    level0.flush()
    del density
    density_path.unlink(missing_ok=True)
    print(flush=True)

    # ── Pass 3: build pyramid ──────────────────────────────────────────────
    levels, pyr_tmp_paths = _build_pyramid(level0, out_path.parent)

    # ── Pass 4: write OME-TIFF ─────────────────────────────────────────────
    _write_ome_pyramid(out_path, levels, channel_name, pixel_size_um)

    size_gb = out_path.stat().st_size / 1e9
    print(f"  saved {out_path.name} ({size_gb:.1f} GB)", flush=True)

    # ── Cleanup ────────────────────────────────────────────────────────────
    del level0
    level0_path.unlink(missing_ok=True)
    for p in pyr_tmp_paths:
        p.unlink(missing_ok=True)


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
