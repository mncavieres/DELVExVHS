# README — Batch-fitting X-shooter spectra with `rvspecfit`

This script walks a folder of FITS files, **selects one X-shooter arm** (UVB/VIS/NIR), runs your `rvspecfit` pipeline on each spectrum, and **appends results to a CSV catalog** after every file (safe to interrupt and resume).

---

## Features

* Process **one arm per run**: `--arm UVB|VIS|NIR` (default: `NIR`).
* Arm-specific **default wavelength ranges** (nm):

  * **UVB:** 305–555
  * **VIS:** 550–1020
  * **NIR:** 1024–2480
    Override with `--wave-min-nm` / `--wave-max-nm`.
* Uses your `rvspecfit` flow: `fitter_ccf` → `vel_fit.process`.
* **Setup string** for `SpecData` via `--setup` (default: `xshooter`).
* Optional **per-arm config** via `--config` (YAML).
* CSV includes `OBJECT` (first column) and `arm`, plus parameters, errors, flags.
* **Resumable**: with `--skip-existing`, it skips rows that already exist for the same `(OBJECT, ARM)`.

---

## Requirements

* Python 3.9+ recommended
* Packages:

  * `rvspecfit`
  * `astropy`
  * `numpy`

Install (example):

```bash
pip install astropy numpy
# Install rvspecfit per its instructions (e.g., pip, conda, or from source).
```

---

## Expected FITS format

* Primary header must include:

  * `OBJECT` — target identifier (used as the first CSV column)
  * `HIERARCH ESO SEQ ARM` — one of `UVB`, `VIS`, `NIR` (used to filter)
* HDU 1 (table) must contain columns:

  * `WAVE` (nm)
  * `FLUX_REDUCED`
  * `ERR_REDUCED`

---

## Usage


```bash
python make_catalog_xsh.py /path/to/fits_folder /path/to/catalog.csv [OPTIONS]
```

### Common options

* `--arm {UVB,VIS,NIR}`: which arm to process (default `NIR`)
* `--config /path/to/config.yaml`: rvspecfit config for this arm (optional)
* `--setup xshooter`: SpecData setup string (default `xshooter`)
* `--recursive`: search subdirectories for FITS
* `--skip-existing`: don’t re-process rows already present for the same `(OBJECT, ARM)`
* `--wave-min-nm  ...` / `--wave-max-nm ...`: override arm defaults (nm)

---

## Examples

### NIR (defaults to 1024–2480 nm)

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm NIR \
  --config /configs/nir_config.yaml \
  --setup xshooter \
  --recursive \
  --skip-existing
```

### VIS (defaults updated to 550–1020 nm)

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm VIS \
  --config /configs/vis_config.yaml \
  --setup xshooter_vis
```

Override VIS wavelength window:

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm VIS \
  --config /configs/vis_config.yaml \
  --setup xshooter_vis \
  --wave-min-nm 560 --wave-max-nm 1000
```

### UVB (defaults updated to 305–555 nm)

```bash
python make_catalog_xsh.py /data/xshooter /data/catalog.csv \
  --arm UVB \
  --config /configs/uvb_config.yaml \
  --setup xshooter_uvb
```

### Building a single catalog across arms

You can append results from different arms into the **same** CSV. The catalog includes an `arm` column, and `--skip-existing` compares `(OBJECT, ARM)` so multiple arms per object are allowed:

```bash
# UVB pass
python make_catalog_xsh.py /data/uvb /data/catalog.csv --arm UVB --config /configs/uvb.yaml --setup xshooter_uvb --skip-existing

# VIS pass
python make_catalog_xsh.py /data/vis /data/catalog.csv --arm VIS --config /configs/vis.yaml --setup xshooter_vis --skip-existing

# NIR pass
python make_catalog_xsh.py /data/nir /data/catalog.csv --arm NIR --config /configs/nir.yaml --setup xshooter --skip-existing
```

---

## What the script writes

Each processed spectrum yields one CSV row with these columns:

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

* `OBJECT` comes from the primary header.
* `arm` is the arm being processed for that file.
* `minimize_success`, `bad_hessian` come from the fit result (booleans).
* `chisq_red` is computed as `chisq / npix` when available.
* On **failure**, a minimal row is still written with flags set to `False` and blanks elsewhere (keeping a complete ledger of attempts).

---

## How it filters files

Only files where the **primary header** has:

```
HIERARCH ESO SEQ ARM == <selected --arm>
```

are processed; others are skipped.

---

## Notes & tips

* **Setup/config per arm:** Choose appropriate `--setup` and `--config` for each arm’s resolution and preprocessing.
* **Wavelengths:** Arm defaults are sensible starting points; refine with `--wave-*` if your reduction masks specific regions.
* **Resuming:** Use `--skip-existing` to avoid re-processing already cataloged `(OBJECT, ARM)` pairs.
* **Provenance:** The `filename` column helps trace back results.

---

## Internals (fit flow)

For each file:

1. Read HDU 1 table (`WAVE`, `FLUX_REDUCED`, `ERR_REDUCED`).
2. Restrict to `[wave_min_nm, wave_max_nm]` (nm).
3. Convert wavelength to Å (×10) and construct `spec_fit.SpecData(setup, ...)`.
4. `fitter_ccf.fit(...)` to get initial parameters (carry `best_vsini` if present).
5. `vel_fit.process(...)` runs the maximum-likelihood fit with `options={'npoly': 15}`.
6. Flatten results → CSV row; append immediately (checkpoint).

---

## Troubleshooting

* **Missing columns** (`WAVE`, `FLUX_REDUCED`, `ERR_REDUCED`) or headers (`OBJECT`, `HIERARCH ESO SEQ ARM`): the file will be skipped with a warning; a minimal failure row is still written.
* **Weirdly small uncertainties / large χ²:** ensure correct variance, masks, and LSF handling in your `rvspecfit` config.
* **Multiple spectra per object:** If you need per-exposure rows instead of per-object, remove `--skip-existing` or change the resumable key to include `filename`.


# Template preparation per arm (UVB / VIS / NIR)

rvspecfit requires you to prepare template files **per spectral configuration** (i.e., per arm), and the **`--setup` string must match** what you use in your code (the `SpecData(setup, ...)` value and the `--setup` CLI option). Each arm (UVB, VIS, NIR) should have its **own setup name**, its **own interpolator**, and its **own CCF** files.

> Summary of steps (run once per arm & setup):
>
> 1. `rvs_read_grid` (create PHOENIX DB; usually once overall)
> 2. `rvs_make_interpol` (interpolated templates at your arm’s λ-range & resolution)
> 3. `rvs_make_nd` (N-D interpolator)
> 4. `rvs_make_ccf` (CCF templates)

These are the **official rvspecfit tools**; see the docs for details and extra options. ([rvspecfit.readthedocs.io][1])

## Units and ranges

* `--lambda0/--lambda1` for **template tools are in Ångström** (not nm).
  Use ×10 to convert from the defaults you pass to the script:

  * UVB default: **305–555 nm → 3050–5550 Å**
  * VIS default: **550–1020 nm → 5500–10200 Å**
  * NIR default: **1024–2480 nm → 10240–24800 Å**

## Example: build templates for each arm

Pick a **root directory** where you’ll store template outputs, e.g. `/tmpl/xshooter/`. You can reuse the same PHOENIX database for all arms, but run `make_interpol`, `make_nd`, and `make_ccf` **separately per arm** with distinct `--setup` names.

### 0) One-time: create PHOENIX database

```bash
rvs_read_grid \
  --prefix /path/to/PHOENIX/v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/ \
  --templdb /tmpl/xshooter/files.db
```

(Do this once; reuse `files.db` below.) ([rvspecfit.readthedocs.io][1])

---

### UVB (setup example: `xshooter_uvb`, 3050–5550 Å)

```bash
# Interpolated spectra for UVB
rvs_make_interpol \
  --setup xshooter_uvb \
  --lambda0 3050 --lambda1 5550 \
  --resol_func '...' \
  --step 0.5 \
  --templdb /tmpl/xshooter/files.db \
  --oprefix /tmpl/xshooter/ \
  --templprefix /path/to/PHOENIX/v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/ \
  --wavefile WAVE_PHOENIX-ACES-AGSS-COND-2011.fits \
  --air \
  --revision v2020x

# Build the n-D interpolator
rvs_make_nd \
  --prefix /tmpl/xshooter/ \
  --setup xshooter_uvb \
  --revision v2020x

# Build the CCF templates
rvs_make_ccf \
  --setup xshooter_uvb \
  --lambda0 3050 --lambda1 5550 \
  --every 30 \
  --vsinis 0,10,300 \
  --prefix /tmpl/xshooter/ \
  --step 0.5 \
  --revision v2020x
```

---

### VIS (setup example: `xshooter_vis`, 5500–10200 Å)

```bash
rvs_make_interpol \
  --setup xshooter_vis \
  --lambda0 5500 --lambda1 10200 \
  --resol_func '...' \
  --step 0.5 \
  --templdb /tmpl/xshooter/files.db \
  --oprefix /tmpl/xshooter/ \
  --templprefix /path/to/PHOENIX/... \
  --wavefile WAVE_PHOENIX-ACES-AGSS-COND-2011.fits \
  --air \
  --revision v2020x

rvs_make_nd \
  --prefix /tmpl/xshooter/ \
  --setup xshooter_vis \
  --revision v2020x

rvs_make_ccf \
  --setup xshooter_vis \
  --lambda0 5500 --lambda1 10200 \
  --every 30 \
  --vsinis 0,10,300 \
  --prefix /tmpl/xshooter/ \
  --step 0.5 \
  --revision v2020x
```

---

### NIR (setup example: `xshooter_nir`, 10240–24800 Å)

```bash
rvs_make_interpol \
  --setup xshooter_nir \
  --lambda0 10240 --lambda1 24800 \
  --resol_func '...' \
  --step 0.5 \
  --templdb /tmpl/xshooter/files.db \
  --oprefix /tmpl/xshooter/ \
  --templprefix /path/to/PHOENIX/... \
  --wavefile WAVE_PHOENIX-ACES-AGSS-COND-2011.fits \
  --air \
  --revision v2020x

rvs_make_nd \
  --prefix /tmpl/xshooter/ \
  --setup xshooter_nir \
  --revision v2020x

rvs_make_ccf \
  --setup xshooter_nir \
  --lambda0 10240 --lambda1 24800 \
  --every 30 \
  --vsinis 0,10,300 \
  --prefix /tmpl/xshooter/ \
  --step 0.5 \
  --revision v2020x
```

> Tip: choose a realistic `--resol_func` for each arm (e.g., a constant resolving power or a wavelength-dependent model). See the rvspecfit docs for details and options. ([rvspecfit.readthedocs.io][1])

## Config files per arm

Create a **separate config YAML** for each arm (or a single YAML that points to different `template_lib` dirs depending on your `--setup`). The key bit is that it must point at the directory where the above tools wrote templates:

```yaml
# nir_config.yaml
template_lib: "/tmpl/xshooter/"   # root where xshooter_nir outputs live
min_vel: -1500
max_vel: 1500
min_vel_step: 0.2
vel_step0: 5
min_vsini: 0.01
max_vsini: 500
second_minimizer: 1
```

Use the matching setup when running the script, e.g.:

```bash
# UVB pass
python make_catalog_xsh.py /data/uvb /data/catalog.csv \
  --arm UVB --setup xshooter_uvb --config uvb_config.yaml

# VIS pass
python make_catalog_xsh.py /data/vis /data/catalog.csv \
  --arm VIS --setup xshooter_vis --config vis_config.yaml

# NIR pass
python make_catalog_xsh.py /data/nir /data/catalog.csv \
  --arm NIR --setup xshooter_nir --config nir_config.yaml
```

## Why separate builds per arm?

Each **arm** has its **own wavelength range and resolution**, so you need **arm-specific template grids** and **CCF templates** that match those characteristics. Using **the same `--setup` string** consistently (in template tools, in your config, and in `SpecData`) ensures the fitter loads the correct interpolator and CCF for that arm. ([rvspecfit.readthedocs.io][1])

---

**Reference:** rvspecfit docs (template preparation pipeline and commands). ([rvspecfit.readthedocs.io][1])

[1]: https://rvspecfit.readthedocs.io/ "RVSpecFit: Automated Spectroscopic Pipeline — rvspecfit  documentation"
