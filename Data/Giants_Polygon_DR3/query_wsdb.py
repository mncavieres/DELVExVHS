# pip install sqlutilpy astropy tqdm
import os, json
import sqlutilpy
from astropy.table import Table
from collections import Counter
from tqdm import tqdm

# credentials
with open('/Users/mncavieres/Documents/2024-1/Delve/credentials_wsdb') as f:
    user, password = f.read().split(', ')

HOST = 'wsdb.ast.cam.ac.uk'
DB   = 'wsdb'

# table names might change to delve dr2 for backwards compatibility in the other scripts
TBL_DELVE = 'delve_dr3.main'   # DELVE DR3
TBL_VHS   = 'vhs_1603.des'     # VHS
TBL_GAIA  = 'gaia_source'      # Gaia DR3

# Polygon vertices in (g_i, i_K) space
RADIUS_ARCSEC = 1.0
PX = [1.63, 1.08, 1.54, 2.00, 2.64, 2.04]  # polygon x = g_i
PY = [2.43, 2.17, 2.61, 2.95, 3.28, 2.65]  # polygon y = i_K
MIN_X, MAX_X = min(PX), max(PX)
MIN_Y, MAX_Y = min(PY), max(PY)

# We will sweep all hpix_32 IDs in [0, 12288)
# (each covers ~3.36 sq.deg, so 12288 covers full sky
NSIDE = 32
NPIX  = 12 * NSIDE * NSIDE  # 12,288
HP_START = 0
HP_END   = NPIX

# path setup
OUTDIR = "/Users/mncavieres/Documents/2024-1/Delve/Data/Giants_Polygon_DR3"
LEDGER = os.path.join(OUTDIR, "_done_hpix32.json")
ERRLOG = os.path.join(OUTDIR, "_errors.log")
os.makedirs(OUTDIR, exist_ok=True)

# ---------- one-time metadata: alias ONLY duplicate column names ----------
def _parse_schema_table(qualified):
    return qualified.split('.', 1) if '.' in qualified else (None, qualified)

def _resolve_schema_for(table_name):
    q = """
    SELECT table_schema, COUNT(*) AS ncols
    FROM information_schema.columns
    WHERE table_name = %s
    GROUP BY table_schema
    ORDER BY ncols DESC
    LIMIT 1
    """
    r = sqlutilpy.get(q, db=DB, host=HOST, user=user, password=password,
                      asDict=True, params=(table_name,))
    return r['table_schema'][0] if r and 'table_schema' in r and len(r['table_schema']) else None

def _list_columns(schema, table):
    if schema:
        q = """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s AND table_name = %s
        ORDER BY ordinal_position
        """
        r = sqlutilpy.get(q, db=DB, host=HOST, user=user, password=password,
                          asDict=True, params=(schema, table))
    else:
        q = """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_name = %s
        ORDER BY ordinal_position
        """
        r = sqlutilpy.get(q, db=DB, host=HOST, user=user, password=password,
                          asDict=True, params=(table,))
    return list(r['column_name'])

d_schema, d_table = _parse_schema_table(TBL_DELVE)
v_schema, v_table = _parse_schema_table(TBL_VHS)
g_schema, g_table = _parse_schema_table(TBL_GAIA)
if g_schema is None:
    g_schema = _resolve_schema_for(g_table)

cols_d = _list_columns(d_schema, d_table)
cols_v = _list_columns(v_schema, v_table)
cols_g = _list_columns(g_schema, g_table)

name_counts = Counter(cols_d + cols_v + cols_g)
SUFFIX = {'d': 'delve', 'v': 'vhs', 'g': 'gaia'}

def _alias_name(col, alias_letter):
    return col if name_counts[col] == 1 else f"{col}_{SUFFIX[alias_letter]}"

def _sel(alias_sql, alias_letter, cols):
    return ",\n    ".join(f'{alias_sql}."{c}" AS "{_alias_name(c, alias_letter)}"' for c in cols)

sel_d     = _sel('dfull', 'd', cols_d)
sel_vfull = _sel('vfull', 'v', cols_v)
sel_gfull = _sel('gfull', 'g', cols_g)

tbl_d = f'{d_schema}.{d_table}'
tbl_v = f'{v_schema}.{v_table}'
tbl_g = f'{g_schema+"." if g_schema else ""}{g_table}'

# main SQL query template
SQL_TEMPLATE = f"""
WITH params AS (
  SELECT %s/3600.0 AS rdeg
),
d0 AS MATERIALIZED (
  SELECT
    ctid AS d_ctid,
    ra, dec,
    wavg_mag_psf_g, wavg_mag_psf_i
  FROM {tbl_d} d
  WHERE d.hpix_32 = %s
    AND d.ext_coadd = 1
    AND d.s_extractor_flags_g = 0
    AND d.s_extractor_flags_i = 0
    AND d.s_extractor_flags_r = 0
    AND d.s_extractor_flags_z = 0
    AND (d.wavg_mag_psf_g - d.wavg_mag_psf_i) > 1.2
    AND (d.wavg_mag_psf_g - d.wavg_mag_psf_i) < 2.6
),
pairs AS (
  SELECT
    d0.d_ctid,
    vbest.v_ctid,
    gbest.source_id,
    q3c_dist(d0.ra, d0.dec, vbest.ra, vbest.dec) AS dist_d_v_deg,
    q3c_dist(d0.ra, d0.dec, gbest.ra, gbest.dec) AS dist_d_g_deg,
    (d0.wavg_mag_psf_g - d0.wavg_mag_psf_i) AS g_i,
    (d0.wavg_mag_psf_i - vbest.kapercormag3) AS i_K
  FROM d0
  JOIN LATERAL (
    SELECT
      v.ctid AS v_ctid,
      v.ra, v.dec,
      v.kapercormag3
    FROM {tbl_v} v, params p
    WHERE q3c_join(d0.ra, d0.dec, v.ra, v.dec, p.rdeg)
      AND v.mergedclass = -1
      AND v.jerrorbit = 0
      AND v.kerrorbit = 0
      AND v.prim = 1
      AND v.javconf > 95
      AND v.kavconf > 95
      AND (v.jx < 8800 OR v.jy < 12300)
      AND v.japercormag4 BETWEEN 12 AND 18
      AND v.kapercormag4 BETWEEN 11 AND 19
    ORDER BY q3c_dist(d0.ra, d0.dec, v.ra, v.dec)
    LIMIT 1
  ) AS vbest ON TRUE
  JOIN LATERAL (
    SELECT
      g.source_id,
      g.ra, g.dec
    FROM {tbl_g} g, params p
    WHERE q3c_join(d0.ra, d0.dec, g.ra, g.dec, p.rdeg)
      AND g.parallax < 0.1
    ORDER BY q3c_dist(d0.ra, d0.dec, g.ra, g.dec)
    LIMIT 1
  ) AS gbest ON TRUE
),
poly AS (
  SELECT
    ARRAY[{", ".join(str(x) for x in PX)}]::float8[] AS px,
    ARRAY[{", ".join(str(y) for y in PY)}]::float8[] AS py
),
filtered AS (
  SELECT p.*
  FROM pairs p, poly
  WHERE p.g_i BETWEEN {MIN_X} AND {MAX_X}
    AND p.i_K BETWEEN {MIN_Y} AND {MAX_Y}
    AND (
      SELECT (SUM(
               CASE
                 WHEN ((py[i] > p.i_K) <> (py[j] > p.i_K))
                      AND (p.g_i < ((px[j]-px[i]) * (p.i_K - py[i])
                                    / NULLIF((py[j]-py[i]), 0) + px[i]))
                 THEN 1 ELSE 0 END
             ) %% 2)           -- NOTE: '%%' to escape modulo for sqlutilpy
      FROM generate_subscripts(px,1) AS i,
           LATERAL (SELECT CASE WHEN i < array_length(px,1) THEN i+1 ELSE 1 END AS j) AS nxt
    ) = 1
)
SELECT
  {sel_d},
  {sel_vfull},
  {sel_gfull},
  f.g_i, f.i_K, f.dist_d_v_deg, f.dist_d_g_deg
FROM filtered f
JOIN {tbl_d} dfull ON dfull.ctid = f.d_ctid
JOIN {tbl_v} vfull ON vfull.ctid = f.v_ctid
JOIN {tbl_g} gfull ON gfull.source_id = f.source_id;
"""

# resume
done = set()
if os.path.exists(LEDGER):
    try:
        done = set(json.load(open(LEDGER)))
    except Exception:
        done = set()

# main loop
pbar = tqdm(range(HP_START, HP_END), desc="hpix_32", unit="pix")
for pix in pbar:
    spix = str(pix)
    if spix in done:
        continue
    outfile = os.path.join(OUTDIR, f"xmatch_hpix32_{pix}.fits")
    if os.path.exists(outfile):
        done.add(spix)
        json.dump(sorted(done, key=int), open(LEDGER, "w"))
        continue

    try:
        # params: (radius_arcsec, hpix_32)
        res = sqlutilpy.get(SQL_TEMPLATE,
                            db=DB, host=HOST, user=user, password=password,
                            asDict=True, params=(RADIUS_ARCSEC, pix))
        if not res or len(res.keys()) == 0:
            Table().write(outfile, overwrite=True)
        else:
            Table(res).write(outfile, overwrite=True)
    except Exception as e:
        with open(ERRLOG, "a") as fh:
            fh.write(f"hpix_32={pix}: {e}\n")
        continue

    done.add(spix)
    json.dump(sorted(done, key=int), open(LEDGER, "w"))
