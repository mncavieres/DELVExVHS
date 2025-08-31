# Batch-fitting ESO X-shooter spectra with RVSpecFit

This package provides a command-line script that scans a folder of X-shooter FITS files, processes **one arm per run** (UVB, VIS, or NIR), and appends the fitted results to a single CSV catalog—**updating the file after every spectrum** so runs are safe to interrupt and resume.

---

## What you get

* A single CSV catalog with one row per processed spectrum.
* Works arm-by-arm: `--arm UVB|VIS|NIR` (default: `NIR`).
* Reasonable **default wavelength windows** (in **nm**):

  * **UVB:** 305–555
  * **VIS:** 550–1020
  * **NIR:** 1024–2480
    You can override per run with `--wave-min-nm / --wave-max-nm`.
* Per-arm **setup** (`--setup`) and **config** (`--config`) so you can match instrument settings and template libraries.
* **Resumable** processing: optional `--skip-existing` avoids re-adding rows for the same `(OBJECT, ARM)`.

---

## Quick start

1. **Install dependencies**

```bash
# minimal python deps
pip install numpy astropy
# install RVSpecFit (follow its instructions: pip/conda/from source)
```

2. **Prepare RVSpecFit templates per arm (one-time per setup)**
   For each arm you plan to process, you must build:

* an **interpolated template grid**,
* the **n-D interpolator**, and
* **CCF templates**,

and you must use the **same `--setup` string** you’ll pass to the script.

> RVSpecFit template tools expect wavelengths in **Ångström**. Convert nm→Å by ×10.

Example shell plan (adjust paths, resolution model, and options):

```bash
# 0) One-time: create PHOENIX database (reused by all arms)
rvs_read_grid \
  --prefix /path/to/PHOENIX/v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/ \
  --templdb /tmpl/xshooter/files.db

# --- UVB (setup: xshooter_uvb), 305–555 nm => 3050–5550 Å ---
rvs_make_interpol --setup xshooter_uvb --lambda0 3050 --lambda1 5550 \
  --resol_func '...' --step 0.5 \
  --templdb /tmpl/xshooter/files.db --oprefix /tmpl/xshooter/ \
  --templprefix /path/to/PHOENIX/... --wavefile WAVE_PHOENIX-ACES-AGSS-COND-2011.fits --air
rvs_make_nd --prefix /tmpl/xshooter/ --setup xshooter_uvb
rvs_make_ccf --setup xshooter_uvb --lambda0 3050 --lambda1 5550 \
  --every 30 --vsinis 0,10,300 --prefix /tmpl/xshooter/ --step 0.5

# --- VIS (setup: xshooter_vis), 550–1020 nm => 5500–10200 Å ---
rvs_make_interpol --setup xshooter_vis --lambda0 5500 --lambda1 10200 \
  --resol_func '...' --step 0.5 \
  --templdb /tmpl/xshooter/files.db --oprefix /tmpl/xshooter/ \
  --templprefix /path/to/PHOENIX/... --wavefile WAVE_PHOENIX-ACES-AGSS-COND-2011.fits --air
rvs_make_nd --prefix /tmpl/xshooter/ --setup xshooter_vis
rvs_make_ccf --setup xshooter_vis --lambda0 5500 --lambda1 10200 \
  --every 30 --vsinis 0,10,300 --prefix /tmpl/xshooter/ --step 0.5

# --- NIR (setup: xshooter_nir), 1024–2480 nm => 10240–24800 Å ---
rvs_make_interpol --setup xshooter_nir --lambda0 10240 --lambda1 24800 \
  --resol_func '...' --step 0.5 \
  --templdb /tmpl/xshooter/files.db --oprefix /tmpl/xshooter/ \
  --templprefix /path/to/PHOENIX/... --wavefile WAVE_PHOENIX-ACES-AGSS-COND-2011.fits --air
rvs_make_nd --prefix /tmpl/xshooter/ --setup xshooter_nir
rvs_make_ccf --setup xshooter_nir --lambda0 10240 --lambda1 24800 \
  --every 30 --vsinis 0,10,300 --prefix /tmpl/xshooter/ --step 0.5
```

3. **Create a config YAML per arm**
   Point RVSpecFit to the template library you just built and set your search ranges. Minimal example:

```yaml
# nir_config.yaml
template_lib: "/tmpl/xshooter/"   # directory used above by the template tools
min_vel: -1500
max_vel: 1500
min_vel_step: 0.2
vel_step0: 5
min_vsini: 0.01
max_vsini: 500
second_minimizer: 1
```

4. **Run the script** (assumes it’s saved as `make_catalog_xsh.py`)

* **NIR** (defaults to 1024–2480 nm):

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm NIR --config /configs/nir_config.yaml --setup xshooter_nir \
  --recursive --skip-existing
```

* **VIS** (defaults to 550–1020 nm):

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm VIS --config /configs/vis_config.yaml --setup xshooter_vis
```

* **UVB** (defaults to 305–555 nm):

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm UVB --config /configs/uvb_config.yaml --setup xshooter_uvb
```

Override the window if needed:

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm VIS --config /configs/vis_config.yaml --setup xshooter_vis \
  --wave-min-nm 560 --wave-max-nm 1000
```

---

## Command reference

```
python make_catalog_xsh.py INPUT_DIR OUTPUT_CSV [options]
```

**Core options**

* `--arm {UVB,VIS,NIR}`: choose the X-shooter arm to process (default: `NIR`)
* `--config PATH`: RVSpecFit YAML config for **this arm/setup**
* `--setup STRING`: setup name used when building templates for this arm (default: `xshooter`)
* `--wave-min-nm FLOAT`, `--wave-max-nm FLOAT`: wavelength range in **nm** (override arm defaults)

**Convenience**

* `--recursive`: search subdirectories for FITS
* `--skip-existing`: don’t re-add rows already present for the same `(OBJECT, ARM)`

---

## Input data expectations

* **Primary header** must contain:

  * `OBJECT` — target identifier (becomes the first CSV column)
  * `HIERARCH ESO SEQ ARM` — should be `UVB`, `VIS`, or `NIR` (used to select files)
* **HDU 1** (table) must provide:

  * `WAVE` (in **nm**)
  * `FLUX_REDUCED`
  * `ERR_REDUCED`

Files with missing headers/columns are skipped and a warning is printed; a minimal “failure” row is still added to the CSV for traceability.

---

## What the CSV contains

Columns (in order):

```
OBJECT, arm,
teff, logg, feh, alpha,
vsini,
vel, vel_err, vel_skewness, vel_kurtosis,
teff_err, logg_err, feh_err, alpha_err,
minimize_success, bad_hessian,
chisq, logl, npix, chisq_red,
filename
```

Notes:

* `arm` is the arm processed for that row.
* `chisq_red` is computed as `chisq / npix` when available.
* On exceptions, a row is still written with flags set to `False` and blanks elsewhere.

You can safely **append multiple arms to the same CSV**; use `--skip-existing` to avoid duplicate `(OBJECT, ARM)` entries.

---

## Tips & troubleshooting

* **Template/CCF per arm is mandatory.** Each arm has different λ-coverage and resolution, so build templates and CCFs **per arm** and keep the **same `--setup`** across template tools, your YAML config, and the script run.
* **Units:** The script’s wavelength range is in **nm**; RVSpecFit template tools (`rvs_make_interpol`, `rvs_make_ccf`) use **Å**.
* **Resolution model:** Choose a realistic `--resol_func` for each arm (constant R or wavelength-dependent).
* **Continuum & masks:** If you see large χ² or unrealistically small errors, check your config (variance weighting, telluric/bad-pixel masks, continuum model).
* **Multiple exposures per object:** If you want one row per exposure regardless of object name, don’t use `--skip-existing`; the `filename` column preserves provenance.

---

## How the script works (at a glance)

For each FITS in `INPUT_DIR` that matches the chosen arm:

1. Read `WAVE`/`FLUX_REDUCED`/`ERR_REDUCED` from HDU 1.
2. Restrict to the selected `[wave_min_nm, wave_max_nm]` (nm), convert to Å.
3. Build an in-memory spectrum with your `--setup`.
4. Run RVSpecFit to estimate RV and stellar parameters for that spectrum.
5. Write a row to `OUTPUT_CSV` immediately (checkpoint).
