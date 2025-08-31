#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv
import math
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table


from rvspecfit import fitter_ccf, vel_fit, spec_fit, utils


CATALOG_COLUMNS = [
    "OBJECT",                 # first column (from primary header)
    "teff", "logg", "feh", "alpha",
    "vsini",
    "vel", "vel_err", "vel_skewness", "vel_kurtosis",
    "teff_err", "logg_err", "feh_err", "alpha_err",
    "minimize_success", "bad_hessian",
    "chisq", "logl", "npix", "chisq_red",
    "filename",               # provenance is handy
]


def _to_float(x):
    try:
        if x is None:
            return math.nan
        return float(np.asarray(x))
    except Exception:
        return math.nan

def _safe_header_get(header, key, default=""):
    try:
        return str(header.get(key, default))
    except Exception:
        return str(default)

def _normalize_key_get(d, k):
    """Access dict where keys may be np.str_ etc., using a plain string key."""
    if not isinstance(d, dict):
        return None
    for key in d.keys():
        if str(key) == k:
            return d[key]
    return None

def _ensure_catalog(csv_path):
    csv_path = Path(csv_path)
    if not csv_path.exists():
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CATALOG_COLUMNS)
            writer.writeheader()

def _iter_fits_files(input_dir, recursive=False):
    exts = ("*.fits", "*.fit", "*.fits.fz", "*.fz")
    if recursive:
        for pat in exts:
            yield from Path(input_dir).rglob(pat)
    else:
        for pat in exts:
            yield from Path(input_dir).glob(pat)


def fit_one_file_with_rvspecfit(fits_path, config_path=None,
                                wave_min_nm=1200.0, wave_max_nm=2478.0):
    """
    Reproduces your code:
    - Table.read(..., hdu=1) with columns WAVE, FLUX_REDUCED, ERR_REDUCED
    - restrict wavelength to [wave_min_nm, wave_max_nm] in nm
    - convert wavelength to Angstrom (×10)
    - fitter_ccf.fit to get initial guess
    - vel_fit.process to do the fit
    Returns the result dict from rvspecfit (as in your example).
    """
    # Config (optional)
    config = utils.read_config(config_path) if config_path else None

    # Read spectral table from HDU 1
    tab = Table.read(str(fits_path), hdu=1)
    wavelength = tab["WAVE"]              # likely in nm
    spec = tab["FLUX_REDUCED"]
    espec = tab["ERR_REDUCED"]

    # Mask to your range (in nm)
    wv = getattr(wavelength, "value", wavelength)  # support Quantity or ndarray
    mask = (wv > wave_min_nm) & (wv < wave_max_nm)
    wavelength = wavelength[mask]
    spec = spec[mask]
    espec = espec[mask]

    # Build SpecData (convert nm→Å by ×10)
    specdata = [spec_fit.SpecData(
        "xshooter",
        wavelength * 10.0,
        spec,
        espec
    )]

    # Initial guess from CCF
    res0 = fitter_ccf.fit(specdata, config)
    param0 = dict(res0["best_par"])  # copy

    # If present, carry over vsini
    if res0.get("best_vsini") is not None:
        param0["vsini"] = res0["best_vsini"]

    # Start off with higher feh, as in your example
    #param0["feh"] = 0.5

    options = {"npoly": 15}

    # Main fit
    res1 = vel_fit.process(
        specdata,
        param0,
        fixParam=[],           
        config=config,
        options=options
    )

    return res1

def _extract_row(result, object_name, fname):
    param = result.get("param", {}) or {}
    perr  = result.get("param_err", {}) or {}

    chisq = result.get("chisq")
    npix = None
    if isinstance(result.get("npix_array"), (list, tuple)) and len(result["npix_array"]) > 0:
        npix = result["npix_array"][0]
    if npix is None:
        npix = result.get("npix")

    chisq_red = math.nan
    if chisq is not None and npix not in (None, 0):
        chisq_red = _to_float(chisq) / _to_float(npix)

    return {
        "OBJECT": object_name,
        "teff": _to_float(_normalize_key_get(param, "teff")),
        "logg": _to_float(_normalize_key_get(param, "logg")),
        "feh": _to_float(_normalize_key_get(param, "feh")),
        "alpha": _to_float(_normalize_key_get(param, "alpha")),
        "vsini": _to_float(result.get("vsini")),
        "vel": _to_float(result.get("vel")),
        "vel_err": _to_float(result.get("vel_err")),
        "vel_skewness": _to_float(result.get("vel_skewness")),
        "vel_kurtosis": _to_float(result.get("vel_kurtosis")),
        "teff_err": _to_float(_normalize_key_get(perr, "teff")),
        "logg_err": _to_float(_normalize_key_get(perr, "logg")),
        "feh_err": _to_float(_normalize_key_get(perr, "feh")),
        "alpha_err": _to_float(_normalize_key_get(perr, "alpha")),
        "minimize_success": bool(result.get("minimize_success", False)),
        "bad_hessian": bool(result.get("bad_hessian", False)),
        "chisq": _to_float(chisq),
        "logl": _to_float(result.get("logl")),
        "npix": int(_to_float(npix)) if not math.isnan(_to_float(npix)) else "",
        "chisq_red": _to_float(chisq_red),
        "filename": str(fname),
    }


def process_folder(input_dir, output_csv, config_path=None, recursive=False,
                   wave_min_nm=1200.0, wave_max_nm=2478.0, skip_existing=False):
    """
    Walk input_dir, filter to NIR arm spectra, fit, and append to CSV.
    """
    input_dir = Path(input_dir)
    output_csv = Path(output_csv)
    _ensure_catalog(output_csv)

    # optional resumability by OBJECT
    existing = set()
    if skip_existing:
        with output_csv.open("r", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                existing.add(r.get("OBJECT", ""))

    with output_csv.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_COLUMNS)

        for fp in _iter_fits_files(input_dir, recursive=recursive):
            try:
                # Primary header for filter + OBJECT
                hdr0 = fits.getheader(fp, 0)
                arm = _safe_header_get(hdr0, "HIERARCH ESO SEQ ARM", default="")
                if arm != "NIR":
                    continue  # skip non-NIR

                obj = _safe_header_get(hdr0, "OBJECT", default=Path(fp).stem)
                if skip_existing and obj in existing:
                    continue

                # Run your fitter
                result = fit_one_file_with_rvspecfit(
                    fp, config_path=config_path,
                    wave_min_nm=wave_min_nm, wave_max_nm=wave_max_nm
                )

                # Write row immediately (checkpoint)
                row = _extract_row(result, obj, Path(fp).name)
                writer.writerow(row)
                f.flush()
                existing.add(obj)

            except Exception as e:
                # Record a failure row so progress is still tracked
                row = {k: "" for k in CATALOG_COLUMNS}
                row.update({
                    "OBJECT": Path(fp).stem,
                    "minimize_success": False,
                    "bad_hessian": False,
                    "filename": Path(fp).name,
                })
                writer.writerow(row)
                f.flush()
                print(f"[WARN] Failed on {fp}: {e}", file=sys.stderr)

def main():
    p = argparse.ArgumentParser(description="Batch-fit NIR X-shooter spectra with rvspecfit and update a CSV catalog.")
    p.add_argument("input_dir", help="Folder containing FITS files")
    p.add_argument("output_csv", help="Path to the catalog CSV to create/update")
    p.add_argument("--config", help="Path to rvspecfit config.yaml (optional)")
    p.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    p.add_argument("--skip-existing", action="store_true", help="Skip OBJECTs already present in the CSV")
    p.add_argument("--wave-min-nm", type=float, default=1200.0, help="Lower wavelength bound in nm (default 1200)")
    p.add_argument("--wave-max-nm", type=float, default=2478.0, help="Upper wavelength bound in nm (default 2478)")
    args = p.parse_args()

    process_folder(
        args.input_dir,
        args.output_csv,
        config_path=args.config,
        recursive=args.recursive,
        wave_min_nm=args.wave_min_nm,
        wave_max_nm=args.wave_max_nm,
        skip_existing=args.skip_existing,
    )

if __name__ == "__main__":
    main()
