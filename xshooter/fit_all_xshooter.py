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
    "OBJECT", "arm",                 # keep OBJECT as first column, add ARM next
    "teff", "logg", "feh", "alpha",
    "vsini",
    "vel", "vel_err", "vel_skewness", "vel_kurtosis",
    "teff_err", "logg_err", "feh_err", "alpha_err",
    "minimize_success", "bad_hessian",
    "chisq", "logl", "npix", "chisq_red",
    "filename",
]

# Sensible X-shooter arm defaults (nm). Can be overridden by CLI args.
DEFAULT_RANGES_NM = {
    "UVB": (300.0, 559.0),
    "VIS": (559.0, 1024.0),
    "NIR": (1024.0, 2480.0),
}


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
    p = Path(input_dir)
    if recursive:
        for pat in exts:
            yield from p.rglob(pat)
    else:
        for pat in exts:
            yield from p.glob(pat)


def fit_one_file_with_rvspecfit(
    fits_path,
    config_path=None,
    setup="xshooter",
    wave_min_nm=None,
    wave_max_nm=None,
):
    """
    Reads HDU=1 table (WAVE, FLUX_REDUCED, ERR_REDUCED), restricts to [min,max] nm,
    converts nm->Å, runs fitter_ccf for init (carries best_vsini), then vel_fit.process.
    """
    config = utils.read_config(config_path) if config_path else None

    tab = Table.read(str(fits_path), hdu=1)
    wavelength = tab["WAVE"]
    spec = tab["FLUX_REDUCED"]
    espec = tab["ERR_REDUCED"]

    wv = getattr(wavelength, "value", wavelength)
    if wave_min_nm is not None and wave_max_nm is not None:
        mask = (wv > wave_min_nm) & (wv < wave_max_nm)
        wavelength = wavelength[mask]
        spec = spec[mask]
        espec = espec[mask]

    specdata = [spec_fit.SpecData(
        setup,
        wavelength * 10.0,  # nm -> Å
        spec,
        espec
    )]

    res0 = fitter_ccf.fit(specdata, config)
    param0 = dict(res0["best_par"])

    if res0.get("best_vsini") is not None:
        param0["vsini"] = res0["best_vsini"]

    options = {"npoly": 15}

    res1 = vel_fit.process(
        specdata,
        param0,
        fixParam=[],
        config=config,
        options=options
    )
    return res1

def _extract_row(result, object_name, fname, arm_str):
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
        "arm": arm_str,
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

def process_folder(
    input_dir,
    output_csv,
    arm="NIR",
    config_path=None,
    setup="xshooter",
    recursive=False,
    wave_min_nm=None,
    wave_max_nm=None,
    skip_existing=False
):
    input_dir = Path(input_dir)
    output_csv = Path(output_csv)
    arm = arm.upper()
    if arm not in ("UVB", "VIS", "NIR"):
        raise ValueError(f"--arm must be one of UVB, VIS, NIR (got {arm})")

    # fill wavelength defaults if not provided
    if wave_min_nm is None or wave_max_nm is None:
        wave_min_nm, wave_max_nm = DEFAULT_RANGES_NM[arm]

    if wave_min_nm >= wave_max_nm:
        raise ValueError("wave-min-nm must be < wave-max-nm")

    _ensure_catalog(output_csv)

    # resumability: key on (OBJECT, ARM) so same object can appear in multiple arms
    existing = set()
    if skip_existing:
        with output_csv.open("r", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                key = (r.get("OBJECT", ""), r.get("arm", "").upper())
                existing.add(key)

    with output_csv.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_COLUMNS)

        for fp in _iter_fits_files(input_dir, recursive=recursive):
            try:
                hdr0 = fits.getheader(fp, 0)
                hdr_arm = _safe_header_get(hdr0, "HIERARCH ESO SEQ ARM", default="").upper()
                if hdr_arm != arm:
                    continue  # only process selected arm

                obj = _safe_header_get(hdr0, "OBJECT", default=Path(fp).stem)
                key = (obj, arm)
                if skip_existing and key in existing:
                    continue

                result = fit_one_file_with_rvspecfit(
                    fp,
                    config_path=config_path,
                    setup=setup,
                    wave_min_nm=wave_min_nm,
                    wave_max_nm=wave_max_nm,
                )

                row = _extract_row(result, obj, Path(fp).name, arm)
                writer.writerow(row)
                f.flush()
                existing.add(key)

            except Exception as e:
                # still write a trace row on failure
                try:
                    hdr_arm = fits.getheader(fp, 0).get("HIERARCH ESO SEQ ARM", arm)
                    obj = fits.getheader(fp, 0).get("OBJECT", Path(fp).stem)
                except Exception:
                    hdr_arm = arm
                    obj = Path(fp).stem
                row = {k: "" for k in CATALOG_COLUMNS}
                row.update({
                    "OBJECT": obj,
                    "arm": str(hdr_arm),
                    "minimize_success": False,
                    "bad_hessian": False,
                    "filename": Path(fp).name,
                })
                writer.writerow(row)
                f.flush()
                print(f"[WARN] Failed on {fp}: {e}", file=sys.stderr)

def main():
    p = argparse.ArgumentParser(
        description="Batch-fit X-shooter spectra (UVB/VIS/NIR) with rvspecfit and update a CSV catalog."
    )
    p.add_argument("input_dir", help="Folder containing FITS files")
    p.add_argument("output_csv", help="Path to the catalog CSV to create/update")
    p.add_argument("--arm", default="NIR", choices=["UVB", "VIS", "NIR"],
                   help="Which arm to process (default: NIR)")
    p.add_argument("--config", help="Path to rvspecfit config.yaml for this arm (optional)")
    p.add_argument("--setup", default="xshooter",
                   help='SpecData setup string for this arm (default: "xshooter")')
    p.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip rows whose (OBJECT, ARM) already exist in the CSV")
    p.add_argument("--wave-min-nm", type=float, default=None,
                   help="Lower wavelength bound in nm (override arm default)")
    p.add_argument("--wave-max-nm", type=float, default=None,
                   help="Upper wavelength bound in nm (override arm default)")
    args = p.parse_args()

    process_folder(
        args.input_dir,
        args.output_csv,
        arm=args.arm,
        config_path=args.config,
        setup=args.setup,
        recursive=args.recursive,
        wave_min_nm=args.wave_min_nm,
        wave_max_nm=args.wave_max_nm,
        skip_existing=args.skip_existing,
    )

if __name__ == "__main__":
    main()
