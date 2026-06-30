# spot-fiction

> Synthesize fake fluorescence channels from MERSCOPE transcript coordinates.

Renders the spatial distribution of detected transcripts as a Gaussian-blurred
density image, producing `mosaic_<NAME>_z{0-6}.tif` files that drop directly
into any MERSCOPE region folder alongside the real DAPI and PolyT channels.

## Requirements

- [pixi](https://pixi.sh) ≥ 0.71
- Linux x86-64 (other platforms: add to `platforms` in `pixi.toml`)
- Disk space: ~22 GB temp + ~11 GB output **per z-plane**
  (~77 GB total for all 7 planes)

## Installation

```bash
git clone git@github.com:imcf/spot-fiction.git
cd spot-fiction
pixi install
```

## Input and Expected Output

### Input

A standard MERSCOPE region directory, e.g. `region_0/`:

```
region_0/
├── detected_transcripts.csv   # required — global_x/y/z + gene columns
└── images/
    ├── manifest.json                        # image dimensions + stain list
    ├── micron_to_mosaic_pixel_transform.csv # affine transform (3×3)
    ├── mosaic_DAPI_z0.tif … z6.tif
    └── mosaic_PolyT_z0.tif … z6.tif
```

### Output

One BigTIFF per z-plane written into `images/`:

```
images/mosaic_Transcripts_z0.tif … z6.tif
```

- Format: single-page uint16 BigTIFF, identical dimensions to DAPI/PolyT
- Values: Gaussian-blurred transcript counts, normalized to full uint16 range
- `manifest.json` is updated automatically with the new stain entries

## Usage

```bash
# All z-planes, all genes, default sigma (10 px ≈ 1.08 µm)
pixi run spot-fiction /path/to/region_0

# Single z-plane — useful for testing before committing 77 GB
pixi run spot-fiction /path/to/region_0 --z 3

# Custom Gaussian blur radius
pixi run spot-fiction /path/to/region_0 --sigma 5

# Subset of genes → separate named channel
pixi run spot-fiction /path/to/region_0 \
    --genes Pomc Rbfox1 Slc17a7 \
    --name ExcNeuron_markers

# All options
pixi run spot-fiction --help
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `data_dir` | *(required)* | Path to MERSCOPE region directory |
| `--sigma PX` | `10` | Gaussian blur radius in pixels (0.108 µm/px) |
| `--genes GENE …` | all | Restrict to specific gene names |
| `--name NAME` | `Transcripts` | Stain name used in output filenames |
| `--z Z` | all (0–6) | Process only one z-plane |

## Changelog

```bash
pixi run changelog   # writes CHANGELOG.md via git-cliff
```

## How It Works

1. **Load** — reads `detected_transcripts.csv` (all columns: `global_x`,
   `global_y`, `global_z`, optionally `gene`).

2. **Transform** — applies the 3×3 affine matrix from
   `micron_to_mosaic_pixel_transform.csv` to convert micron coordinates to
   integer mosaic pixel positions.

3. **Z-plane assignment** — `round(global_z).clip(0, 6)` maps the continuous
   z coordinate to one of the 7 optical sections.

4. **Density rendering** — for each z-plane, transcripts are placed as unit
   impulses into a float32 image processed in horizontal strips (1024 rows
   at a time with 3σ overlap padding). `scipy.ndimage.gaussian_filter`
   smooths each strip; the padded borders prevent edge artefacts.

5. **Normalization** — the full float32 density map is written to a temporary
   numpy memmap (`~22 GB`), the global maximum is computed, and the result
   is linearly scaled to `[0, 65535]` uint16.

6. **Output** — a single-page uint16 BigTIFF is created via `tifffile.memmap`,
   matching the format of the original mosaic files. The temp file is deleted.

7. **Manifest** — `images/manifest.json` is updated with the new stain entries
   so downstream tools (e.g. [SOPA](https://github.com/gustaveroussy/sopa),
   [napari](https://napari.org)) discover the channel automatically.

## License

GNU General Public License v3.0 or later — see [LICENSE](LICENSE).
