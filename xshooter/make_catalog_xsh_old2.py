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
from tqdm import tqdm

CATALOG_COLUMNS = [
    "OBJECT", "arm",
    "teff", "logg", "feh", "alpha",
    "vsini",
    "vel", "vel_err", "vel_skewness", "vel_kurtosis",
    "teff_err", "logg_err", "feh_err", "alpha_err",
    "minimize_success", "bad_hessian",
    "chisq", "logl", "npix", "chisq_red",
    "filename",
]

# Arm defaults (nm)
DEFAULT_RANGES_NM = {
    "UVB": (300.0, 559.0),
    "VIS": (552.0, 1019.0),
    "NIR": (1024.0, 2480.0),
}

# ---------- Helpers ----------
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

def _compute_chi2_red(result):
    """Prefer array totals if present; else fall back to chisq/npix."""
    chi2 = result.get("chisq")
    npix = result.get("npix")
    chi2_arr = result.get("chisq_array")
    npix_arr = result.get("npix_array")
    if chi2_arr is not None:
        try:
            chi2 = float(np.sum(np.asarray(chi2_arr)))
        except Exception:
            pass
    if npix_arr is not None:
        try:
            npix = int(np.sum(np.asarray(npix_arr)))
        except Exception:
            pass
    if chi2 is None or npix in (None, 0):
        return math.nan, _to_float(result.get("chisq")), ("" if npix is None else int(_to_float(npix)))
    return (float(chi2) / float(npix)), float(chi2), int(npix)

def fit_one_file_with_rvspecfit(
    fits_path,
    config_path=None,
    setup="xshooter",
    wave_min_nm=None,
    wave_max_nm=None,
    hdu=1,
    arm="NIR",
):
    """
    Reads HDU table (WAVE, FLUX_REDUCED, ERR_REDUCED), optional nm window,
    converts nm->Å, runs fitter_ccf (carry best_vsini), then vel_fit.process.
    """
    config = utils.read_config(config_path) if config_path else None

    tab = Table.read(str(fits_path), hdu=hdu)
    wavelength = tab["WAVE"]             # expected in nm (Quantity or ndarray)
    spec = np.asarray(tab["FLUX_REDUCED"], dtype=float)
    espec = np.asarray(tab["ERR_REDUCED"], dtype=float)

    # nm window
    wv = getattr(wavelength, "value", wavelength)
    if wave_min_nm is not None and wave_max_nm is not None:
        mask_nm = (wv > wave_min_nm) & (wv < wave_max_nm)
    else:
        mask_nm = np.ones_like(np.asarray(wv), dtype=bool)

    # ---- universal error inflation ----
    # Always de-weight non-positive uncertainties for ALL arms
    espec = np.where(espec <= 0, 1e10, espec)

    # ---- VIS-specific telluric O2 band (nm) ----
    if arm.upper() == "VIS":
        espec = np.where((wv > 757.5) & (wv < 767.5), 1e10, espec)

    # Keep only finite values (no espec>0 check needed since we inflated <=0)
    finite = np.isfinite(wv) & np.isfinite(spec) & np.isfinite(espec)
    mask = mask_nm & finite

    wave_A = (wavelength[mask] * 10.0)  # nm -> Å
    flux = spec[mask]
    eflux = espec[mask]

    # Build SpecData
    specdata = [spec_fit.SpecData(setup, wave_A, flux, eflux)]

    # Initial guess
    res0 = fitter_ccf.fit(specdata, config)
    param0 = dict(res0["best_par"])
    if res0.get("best_vsini") is not None:
        param0["vsini"] = res0["best_vsini"]

    options = {"npoly": 15}
    res1 = vel_fit.process(specdata, param0, fixParam=[], config=config, options=options)
    return res1


def _extract_row(result, object_name, fname, arm_str):
    param = result.get("param", {}) or {}
    perr  = result.get("param_err", {}) or {}

    chi2_red, chi2_val, npix_val = _compute_chi2_red(result)

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
        "chisq": _to_float(chi2_val),
        "logl": _to_float(result.get("logl")),
        "npix": npix_val if npix_val != "" else "",
        "chisq_red": _to_float(chi2_red),
        "filename": str(fname),
    }

def process_folder(
    input_dir,
    output_csv,
    arm="NIR",
    config_path=None,
    setup=None,            # <- allow None so we can default per-arm
    recursive=False,
    wave_min_nm=None,
    wave_max_nm=None,
    skip_existing=False,
    object_key="OBJECT",
    hdu=1,
):
    input_dir = Path(input_dir)
    output_csv = Path(output_csv)
    arm = arm.upper()
    if arm not in ("UVB", "VIS", "NIR"):
        raise ValueError(f"--arm must be one of UVB, VIS, NIR (got {arm})")

    # ----- default setup per arm -----
    if setup is None:
        setup = "vis_xshooter" if arm == "VIS" else "xshooter"

    # fill wavelength defaults if not provided
    if wave_min_nm is None or wave_max_nm is None:
        wave_min_nm, wave_max_nm = DEFAULT_RANGES_NM[arm]
    if wave_min_nm >= wave_max_nm:
        raise ValueError("wave-min-nm must be < wave-max-nm")

    _ensure_catalog(output_csv)

    # resumability: key on (OBJECT, ARM)
    existing = set()
    if skip_existing:
        with output_csv.open("r", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                key = (r.get("OBJECT", ""), r.get("arm", "").upper())
                existing.add(key)

    n_processed = 0
    n_skipped = 0
    n_failed = 0

    with output_csv.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_COLUMNS)

        for fp in tqdm(_iter_fits_files(input_dir, recursive=recursive), 
                       desc="Processing files",
                         total=len(list(_iter_fits_files(input_dir, recursive=recursive)))):
            try:
                hdr0 = fits.getheader(fp, 0)
                hdr_arm = _safe_header_get(hdr0, "HIERARCH ESO SEQ ARM", default="").upper()
                if hdr_arm != arm:
                    n_skipped += 1
                    continue  # only process selected arm

                # OBJECT with fallback
                obj = _safe_header_get(hdr0, object_key, default="").strip()
                if not obj:
                    obj = _safe_header_get(hdr0, "HIERARCH ESO OBS TARG NAME", default=Path(fp).stem).strip()

                key = (obj, arm)
                if skip_existing and key in existing:
                    n_skipped += 1
                    continue

                result = fit_one_file_with_rvspecfit(
                    fp,
                    config_path=config_path,
                    setup=setup,
                    wave_min_nm=wave_min_nm,
                    wave_max_nm=wave_max_nm,
                    hdu=hdu,
                    arm=arm,  # pass through for VIS-specific masking
                )

                row = _extract_row(result, obj, Path(fp).name, arm)
                writer.writerow(row)
                f.flush()
                existing.add(key)
                n_processed += 1

            except Exception as e:
                # Still write a trace row on failure
                try:
                    obj_try = fits.getheader(fp, 0).get(object_key, Path(fp).stem)
                except Exception:
                    obj_try = Path(fp).stem
                row = {k: "" for k in CATALOG_COLUMNS}
                row.update({
                    "OBJECT": str(obj_try),
                    "arm": arm,
                    "minimize_success": False,
                    "bad_hessian": False,
                    "filename": Path(fp).name,
                })
                try:
                    writer.writerow(row)
                    f.flush()
                except Exception:
                    pass
                n_failed += 1
                print(f"[WARN] Failed on {fp}: {e}", file=sys.stderr)

    print(f"[INFO] Done. processed={n_processed}, skipped={n_skipped}, failed={n_failed}", file=sys.stderr)

def main():
    p = argparse.ArgumentParser(
        description="Batch-fit X-shooter spectra (UVB/VIS/NIR) with rvspecfit and update a CSV catalog."
    )
    p.add_argument("input_dir", help="Folder containing FITS files")
    p.add_argument("output_csv", help="Path to the catalog CSV to create/update")
    p.add_argument("--arm", default="NIR", choices=["UVB", "VIS", "NIR"],
                   help="Which arm to process (default: NIR)")
    p.add_argument("--config", help="Path to rvspecfit config.yaml for this arm (optional)")
    # default is arm-dependent; if user omits --setup, VIS -> vis_xshooter, else xshooter
    p.add_argument("--setup", default=None,
                   help='SpecData setup string; default is "vis_xshooter" for VIS, otherwise "xshooter"')
    p.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip rows whose (OBJECT, ARM) already exist in the CSV")
    p.add_argument("--wave-min-nm", type=float, default=None,
                   help="Lower wavelength bound in nm (override arm default)")
    p.add_argument("--wave-max-nm", type=float, default=None,
                   help="Upper wavelength bound in nm (override arm default)")
    p.add_argument("--object-key", default="OBJECT",
                   help='Primary header keyword for target name (default: "OBJECT")')
    p.add_argument("--hdu", type=int, default=1,
                   help="HDU index of the spectral table (default: 1)")
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
        object_key=args.object_key,
        hdu=args.hdu,
    )

if __name__ == "__main__":
    main()
