"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-09-16 with heavy help from Co-Pilot. It has written almost all of this script, so we'll see how that goes.



"""

import sys


sys.path.append('../')

##!/usr/bin/env python3
#"""
#Crossmatch Hunt & Reffert (2024) open cluster catalogue (J/A+A/686/A42)
#with Gentile Fusillo et al. (2021) white dwarf catalogue (J/MNRAS/508/3877)
#using Gaia source_ids, filtering by membership probability.

#Pure Astropy Table workflow — preserves units and metadata.
#"""

#from astroquery.vizier import Vizier
#from astropy.table import Table, join, vstack
#import numpy as np

## -----------------------------
## Config
## -----------------------------
#CAT_HR24 = "J/A+A/686/A42/members"       # Hunt & Reffert 2024 open cluster catalogue
#CAT_GF21 = "J/MNRAS/508/3877/maincat"    # Gentile Fusillo et al. 2021 Gaia EDR3 WD catalogue

#PMEM_THRESHOLD = 0.7  # keep only members with Pmem >= 0.7
#OUT_MATCH = "crossmatch_HR2024xGF2021.fits"

## -----------------------------
## Helpers
## -----------------------------
#def vizier_to_table(catalog, row_limit=-1):
    #"""
    #Download all tables from a VizieR catalog and vertically stack them
    #into a single Astropy Table, adding a __table column for provenance.
    #"""
    #Vizier.ROW_LIMIT = row_limit
    #Vizier.TIMEOUT = 1200
    #v = Vizier(columns=["**"])
    #tables = v.get_catalogs(catalog)  # TableList

    #if len(tables) == 0:
        #raise RuntimeError(f"No tables returned for {catalog}")

    #for idx, t in enumerate(tables):
        ##t["__table"] = t.meta.get("name", f"{catalog}_table{idx+1}")
        #print(catalog,idx, t.meta.get("name"))
    #if len(tables) > 1:
        #combined = vstack(tables, metadata_conflicts="silent")
    #else:
        #combined = tables[0]

    #return combined

#def find_gaia_id_column(tbl):
    #candidates = ["source_id", "GaiaDR3", "GaiaEDR3", "Source", "Gaia"]
    #for col in candidates:
        #if col in tbl.colnames:
            #return col
    #for col in tbl.colnames:
        #if "source" in col.lower() and "id" in col.lower():
            #return col
    #raise KeyError("No Gaia source_id-like column found.")

#def find_probability_column(tbl):
    #candidates = ["Pmem", "PMemb", "prob", "MemberProb", "Prob", "pmem"]
    #for col in candidates:
        #if col in tbl.colnames:
            #return col
    #for col in tbl.colnames:
        #if "prob" in col.lower():
            #return col
    #raise KeyError("No membership probability column found.")

## -----------------------------
## Load catalogs
## -----------------------------
#print("Downloading Hunt & Reffert 2024...")
#hr = vizier_to_table(CAT_HR24, row_limit=-1)

#print("Downloading Gentile Fusillo 2021...")
#gf = vizier_to_table(CAT_GF21, row_limit=-1)

#gf.pprint()
#print(len(gf))
#hr.pprint()
#print(len(hr))
## -----------------------------
## Identify Gaia ID columns
## -----------------------------
#hr_id_col = find_gaia_id_column(hr)
#gf_id_col = find_gaia_id_column(gf)

## -----------------------------
## Show HR24 columns for inspection
## -----------------------------
#print("\n--- Hunt & Reffert 2024 columns ---")
#for col in hr.colnames:
    #print(col)
#print("\nFirst 5 rows of HR24:")
#print(hr[:5])

## -----------------------------
## Filter HR24 by membership probability
## -----------------------------
#prob_col = find_probability_column(hr)
#print(f"\nFiltering HR24 members with {prob_col} >= {PMEM_THRESHOLD}")
#mask = np.array(hr[prob_col], dtype=float) >= PMEM_THRESHOLD
#hr = hr[mask]

## -----------------------------
## Join on Gaia ID
## -----------------------------
#match = join(
    #hr, gf,
    #keys_left=hr_id_col,
    #keys_right=gf_id_col,
    #join_type="inner",
    #table_names=["HR24", "GF21"],
    #uniq_col_name="{col_name}_{table_name}"
#)

## -----------------------------
## Save output
## -----------------------------
#match.write(OUT_MATCH, overwrite=True)

## -----------------------------
## Report
## -----------------------------
#print(f"\nRows in Hunt & Reffert after filtering: {len(hr)}")
#print(f"Rows in Gentile Fusillo raw: {len(gf)}")
#print(f"Cross-matched rows: {len(match)}")
#print(f"Output saved to: {OUT_MATCH}")

##!/usr/bin/env python3
#"""
#Crossmatch Hunt & Reffert (2024) open cluster members table
#(J/A+A/686/A42/members)
#with Gentile Fusillo et al. (2021) main white dwarf catalogue
#(J/MNRAS/508/3877/maincat)
#using Gaia source_ids, filtering by membership probability.

#Pure Astropy Table workflow — preserves units and metadata.
#"""

#from astroquery.vizier import Vizier
#from astropy.table import Table, join
#import numpy as np

## -----------------------------
## Config
## -----------------------------
#CAT_HR24_MEMBERS = "J/A+A/686/A42/members"     # star-level members
#CAT_GF21_MAIN    = "J/MNRAS/508/3877/maincat"  # main WD catalogue

#PMEM_THRESHOLD = 0.7  # keep only members with Pmem >= 0.7
#OUT_MATCH = "crossmatch_HR2024members_x_GF2021maincat.fits"

## -----------------------------
## Helpers
## -----------------------------
#def vizier_single_table(catalog):
    #"""
    #Download a single table from VizieR and return as an Astropy Table.
    #"""
    #Vizier.ROW_LIMIT = -1
    #Vizier.TIMEOUT = 1200
    #v = Vizier(columns=["**"],row_limit=-1)
    ##tables = v.get_catalogs(catalog)
    #print('catalog',catalog)
    #tables = v.query_constraints(catalog=catalog)
    #if len(tables) == 0:
        #raise RuntimeError(f"No tables returned for {catalog}")
    #if len(tables) > 1:
        #print(f"Warning: more than one table returned for {catalog}, using the first.")
    #return tables[0]

#def find_gaia_id_column(tbl):
    #candidates = ["source_id", "GaiaDR3", "GaiaEDR3", "Source", "Gaia"]
    #for col in candidates:
        #if col in tbl.colnames:
            #return col
    #for col in tbl.colnames:
        #if "source" in col.lower() and "id" in col.lower():
            #return col
    #raise KeyError("No Gaia source_id-like column found.")

#def find_probability_column(tbl):
    #candidates = ["Pmem", "PMemb", "prob", "MemberProb", "Prob", "pmem"]
    #for col in candidates:
        #if col in tbl.colnames:
            #return col
    #for col in tbl.colnames:
        #if "prob" in col.lower():
            #return col
    #raise KeyError("No membership probability column found.")

## -----------------------------
## Load specific sub-tables
## -----------------------------
#print("Downloading Hunt & Reffert 2024 members table...")
#hr = vizier_single_table(CAT_HR24_MEMBERS)

#print("Downloading Gentile Fusillo 2021 maincat table...")
#gf = vizier_single_table(CAT_GF21_MAIN)

## -----------------------------
## Identify Gaia ID columns
## -----------------------------
#hr_id_col = find_gaia_id_column(hr)
#gf_id_col = find_gaia_id_column(gf)

## -----------------------------
## Show HR24 columns for inspection
## -----------------------------
#print("\n--- Hunt & Reffert 2024 members table columns ---")
#for col in hr.colnames:
    #print(col)
#print("\nFirst 5 rows of HR24 members table:")
#print(hr[:5])

## -----------------------------
## Filter HR24 by membership probability
## -----------------------------
#prob_col = find_probability_column(hr)
#print(f"\nFiltering HR24 members with {prob_col} >= {PMEM_THRESHOLD}")
#mask = np.array(hr[prob_col], dtype=float) >= PMEM_THRESHOLD
#hr = hr[mask]

## -----------------------------
## Join on Gaia ID
## -----------------------------
#match = join(
    #hr, gf,
    #keys_left=hr_id_col,
    #keys_right=gf_id_col,
    #join_type="inner",
    #table_names=["HR24", "GF21"],
    #uniq_col_name="{col_name}_{table_name}"
#)

## -----------------------------
## Save output
## -----------------------------
#match.write(OUT_MATCH, overwrite=True)

## -----------------------------
## Report
## -----------------------------
#print(f"\nRows in Hunt & Reffert members after filtering: {len(hr)}")
#print(f"Rows in Gentile Fusillo maincat: {len(gf)}")
#print(f"Cross-matched rows: {len(match)}")
#print(f"Output saved to: {OUT_MATCH}")


##I'm actually going to try this as an ADQL script in the Gaia portal per Co-pilot's recommendation. Never mind it was delusional and those tables aren't mirrored in the Gaia portal.

from astroquery.utils.tap.core import TapPlus
import numpy as np
import astropy.table as tbl

def normalize_scaled_integers(tbl, verbose=True):
    """
    Convert integer columns with units like mmag, mas, km/s into floats
    and apply the correct scale factor (e.g. mmag -> mag).
    """
    # Map of unit string (lowercase) to scale factor and new unit
    scale_map = {
        'mmag': (1/1000.0, 'mag'),       # millimagnitudes -> magnitudes
        'mas': (1.0, 'mas'),             # milliarcseconds (keep as float)
        'mas/yr': (1.0, 'mas/yr'),       # proper motions
        'km/s': (1.0, 'km/s'),           # radial velocities
    }

    fixed = []
    for name in tbl.colnames:
        col = tbl[name]
        unit = str(getattr(col, 'unit', '')).lower().strip()

        if col.dtype.kind in ('i', 'u') and unit in scale_map:
            scale, new_unit = scale_map[unit]
            tbl[name] = col.astype('float64') * scale
            try:
                tbl[name].unit = new_unit
            except Exception:
                pass
            fixed.append((name, unit, new_unit, scale))

    if verbose:
        if fixed:
            print("Converted integer columns to floats with scaling:")
            for name, old_u, new_u, scale in fixed:
                print(f" - {name}: {old_u} -> {new_u} (×{scale})")
        else:
            print("No integer columns with matching units found.")
    return
# Example usage:
# result = job.get_results()

# Connect to the VizieR TAP service
tapvizier = TapPlus(url="http://tapvizier.cds.unistra.fr/TAPVizieR/tap")

# ADQL query: server-side join
adql_query = """
SELECT
    hr.GaiaDR3,
    hr.ID
FROM
    "J/A+A/686/A42/members" AS hr
JOIN
    "J/MNRAS/508/3877/maincat" AS gf
    ON hr.GaiaDR3 = gf.GaiaEDR3
JOIN
    "J/A+A/686/A42/clusters" AS cl
    ON hr.ID = cl.ID
WHERE
    hr.Prob >=0.5
"""

#adql_query=""" SELECT TOP 1 * FROM  "J/A+A/686/A42/members" """
#adql_query=""" SELECT TOP 1 * FROM  "J/A+A/686/A42/clusters" """

#adql_query=""" SELECT TOP 1 * FROM  "J/MNRAS/508/3877/maincat" """


# Launch the job asynchronously (good for large queries)
job = tapvizier.launch_job_async(adql_query)
#job = tapvizier.launch_job(adql_query).get_results()
#result=job
cross_ids = job.get_results()
cross_ids.pprint()

from astroquery.vizier import Vizier

# Remove row limit
Vizier.ROW_LIMIT = -1
Vizier.TIMEOUT = 1200
# Convert IDs to Python list
id_list = [str(sid) for sid in cross_ids['GaiaDR3']]
cluster_list = [str(sid) for sid in cross_ids['ID']]

# Helper: query a catalogue by Gaia DR3 ID list
def fetch_catalog_params(catalog_id, id_list, id_col='GaiaDR3'):
    v = Vizier(columns=["**"])
    # Build constraint dict: {column_name: [list_of_values]}
    constraints = {id_col: id_list}
    tables = v.query_constraints(catalog=catalog_id, **constraints)
    return tables[0] if tables else None

#from astroquery.vizier import Vizier
import math

Vizier.ROW_LIMIT = -1
Vizier.TIMEOUT = 300  # increase timeout to 5 minutes

def fetch_catalog_params_chunked(catalog_id, id_list, id_col='source_id', chunk_size=500):
    Vizier.TIMEOUT = 1200  # increase timeout to 5 minutes
    v = Vizier(columns=["**"])
    Vizier.TIMEOUT = 1200  # increase timeout to 5 minutes

    results = []
    n_chunks = math.ceil(len(id_list) / chunk_size)
    for i in range(0, len(id_list), chunk_size):
        chunk = id_list[i:i+chunk_size]
        constraints = {id_col: chunk}
        tables = v.query_constraints(catalog=catalog_id, **constraints)
        if tables:
            results.append(tables[0])
        print(f"Fetched chunk {i//chunk_size+1}/{n_chunks} ({len(chunk)} IDs)")
    if results:
        from astropy.table import vstack
        return vstack(results, join_type='exact')
    return None


# Fetch HR24 members parameters
print('querying HR24 for cluster members info')
#hr24_params = fetch_catalog_params("J/A+A/686/A42/members", id_list)
hr24_params = fetch_catalog_params_chunked("J/A+A/686/A42/members", id_list,id_col='GaiaDR3')

# Fetch GF21 maincat parameters
print('querying GF21 for WD info')

#gf21_params = fetch_catalog_params("J/MNRAS/508/3877/maincat", id_list,id_col='GaiaDR3')
gf21_params = fetch_catalog_params_chunked("J/MNRAS/508/3877/maincat", id_list,id_col='GaiaEDR3')
print(gf21_params)
print('querying  HR 24 for cluster info') 
#hr24_clusters=fetch_catalog_params("J/A+A/686/A42/clusters", cluster_list,id_col='ID')
hr24_clusters=fetch_catalog_params_chunked("J/A+A/686/A42/clusters", cluster_list,id_col='ID')

gf_hr=tbl.join(hr24_params, gf21_params, keys_left='GaiaDR3', keys_right='GaiaEDR3')
gf_hr_cluster=tbl.join(gf_hr,hr24_clusters,keys='ID')
#result.write("HR24members_GF21maincat_crossmatch.fits", overwrite=True)

#normalize_scaled_integers(result)

#result['Gmag','Gmag_2'].pprint()
## Save locally in your preferred format
#result.write("HR24members_GF21maincat_crossmatch.csv", overwrite=True)
gf_hr_cluster.write("HR24members_GF21maincat_crossmatch.csv", overwrite=True)

print(f"Downloaded {len(result)} matched rows")
