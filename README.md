# spot-fiction

> Synthesize fake fluorescence channels from MERSCOPE transcript coordinates.

Renders the spatial distribution of detected transcripts as a Gaussian-blurred
density image, producing `mosaic_<NAME>_z{0-6}.tif` files that drop directly
into any MERSCOPE region folder alongside the real DAPI and PolyT channels.

## Requirements

- [pixi](https://pixi.sh) ≥ 0.71
- Linux x86-64 (other platforms: add to `platforms` in `pixi.toml`)
- Disk space: ~33 GB peak temp + ~15 GB output **per z-plane**
  (~105 GB total for all 7 planes; pyramid adds ~33% over flat)

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

- Format: pyramidal OME-TIFF (tiled 512×512, LZW, multi-resolution)
- Dtype: uint16, normalized to full [0, 65535] range per z-plane
- Sub-resolutions: 2× block-mean downsample at each level, down to < 256 px
- OME-XML encodes channel name and physical pixel size (µm)
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

5. **Normalization** — the full float32 density map (~22 GB temp memmap) is
   linearly scaled to `[0, 65535]` uint16 into a second temp memmap (~11 GB).
   The float32 temp is deleted at this point.

6. **Pyramid** — sub-resolution levels are built by 2× block-mean
   downsampling. Levels > 2 GB use temp numpy memmaps; smaller levels are
   held in RAM. Downsampling stops when both dimensions drop below 256 px
   (typically ~10 levels for this dataset).

7. **Output** — all levels are written as a single tiled (512×512), LZW-
   compressed pyramidal OME-TIFF via tifffile. OME-XML embeds the channel
   name and physical pixel size. Final file ~15 GB per z-plane. All temp
   files are deleted.

8. **Manifest** — `images/manifest.json` is updated with the new stain entries
   so downstream tools (e.g. [SOPA](https://github.com/gustaveroussy/sopa),
   [napari](https://napari.org)) discover the channel automatically.

## License

GNU General Public License v3.0 or later — see [LICENSE](LICENSE).
