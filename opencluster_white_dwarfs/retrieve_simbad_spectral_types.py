"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-07-31 Moved over for the open clusters on 2025-09-16

Created with the help of Copilot





"""


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column, vstack, join, MaskedColumn
#import scipy.interpolate as scinterp
import time
#import pyvo
#from astroquery.vizier import Vizier

#from astroquery.gaia import Gaia
#from astroquery.xmatch import XMatch


sys.path.append('../')
sys.path.append('/Users/BenKaiser/Desktop/radial_velocity_calculations/')
import cal_params as cp
#import spec_plot_tools as spt

plt.rc('lines',linewidth=0.5)

#input_name='gf21_GaiaeDR3_faintWDs_veryinterestingwds.csv'
#input_name='gf21_GaiaeDR3_faintWDs_prettyinterestingwds.csv'
#input_name='gf21_GaiaeDR3_faintWDs_gaiaadded_veryinterestingwds.csv'
#input_name='gf21_GaiaeDR3_faintWDs_gaiaadded_prettyinterestingwds.csv'
input_name='HR24members_GF21maincat_crossmatch.csv'
input_table=Table.read(input_name)
#input_table=Table.read('gf21_GaiaeDR3_faintWDs_prettyinterestingwds.csv')
                       
                       
import re
from astroquery.simbad import Simbad
import re
from astropy.table import Table, vstack
from astroquery.simbad import Simbad

default_chunk_size=1000


def query_simbad_name_and_sptype(source_ids, chunk_size=default_chunk_size):
    """
    Batch‐query SIMBAD for each Gaia DR3 ID’s canonical name and spectral type.

    Parameters
    ----------
    source_ids : array‐like of int
        Gaia DR3 source IDs.
    chunk_size : int
        Number of IDs per SIMBAD batch.

    Returns
    -------
    Astropy Table
        Columns: source_id, simbad_name, spectral_type
    """
    simbad = Simbad()
    simbad.reset_votable_fields()
    simbad.add_votable_fields('ids', 'main_id', 'sptype','sp_qual')

    tables = []
    for start in range(0, len(source_ids), chunk_size):
        batch = source_ids[start:start+chunk_size]
        names = [f"Gaia DR3 {sid}" for sid in batch]
        result = simbad.query_objects(names)
        if result is None:
            continue

        #sids, simbad_names, sptypes = [], [], []
        sids, simbad_names, sptypes,spquals = [], [], [],[]
        print('Results')
        result.pprint()
        print(result.colnames)
        for row in result:
            # parse out the Gaia DR3 source_id from IDS field
            ids_str = (row['ids'].decode()
                       if isinstance(row['ids'], bytes)
                       else str(row['ids']))
            m = re.search(r'Gaia DR3\s+(\d+)', ids_str)
            #sid = int(m.group(1)) if m else None
            sid = int(m.group(1)) if m else 0

            # get SIMBAD’s preferred name
            name = (row['main_id'].decode()
                    if isinstance(row['main_id'], bytes)
                    else row['main_id'])

            # get spectral type (may be None)
            sp = row['sp_type']
            sp = sp.decode().strip() if isinstance(sp, bytes) else sp
            
            spqual = row['sp_qual']
            spqual = spqual.decode().strip() if isinstance(spqual, bytes) else spqual

            sids.append(sid)
            simbad_names.append(name)
            sptypes.append(sp)
            spquals.append(spqual)

        #tables.append(
            #Table([sids, simbad_names, sptypes],
                  #names=('source_id', 'simbad_name', 'simbad_sp_type'))
        #)
        tables.append(
            Table([sids, simbad_names, sptypes,spquals],
                  names=('source_id', 'simbad_name', 'simbad_sp_type','simbad_sp_qual'),dtype=('int64', 'U50', 'U20','U20'))
        )

    if tables:
        return vstack(tables)
    else:
        # empty table with correct column names & types
        return Table(names=('source_id', 'simbad_name', 'spectral_type'),
                     dtype=('int64', 'U50', 'U20'))


# ---------------- Example Usage ----------------

# 1. Load your existing Astropy table (must contain 'source_id' column)
# tbl = Table.read('your_input_table.fits')

# 2. Extract Gaia DR3 IDs
# source_ids = [int(x) for x in tbl['source_id'] if x is not None]

# 3. Query SIMBAD
# simbad_tbl = query_simbad_name_and_sptype(source_ids)

# 4. Merge back into your table
# merged = tbl.join(simbad_tbl, keys='source_id', join_type='left')

# 5. Inspect or save
# print(merged['source_id', 'simbad_name', 'spectral_type'])
# merged.write('with_simbad_names_and_sptypes.fits', overwrite=True)

import re
from astropy.table import Table, vstack
from astroquery.simbad import Simbad

def query_simbad_all_gaia_ids(source_ids, chunk_size=default_chunk_size):
    """
    For each Gaia DR3 ID in `source_ids`, query SIMBAD and extract:
      – Gaia DR2 ID (if present)
      – Gaia eDR3 ID (if present)
      – Gaia DR3 ID  (should match your input)
      – SIMBAD main name
      – Spectral type

    Returns an Astropy Table with columns:
      ['input_dr3','gaia_dr2','gaia_eDR3','simbad_name','spectral_type']
    """
    # prepare the SIMBAD query
    simbad = Simbad()
    simbad.reset_votable_fields()
    simbad.add_votable_fields('ids', 'main_id', 'sp_type')

    # regexes to find each Gaia version
    pat_dr2  = re.compile(r'Gaia\s*DR2\s+(\d+)')
    pat_edr3 = re.compile(r'Gaia\s*eDR3\s+(\d+)')
    pat_dr3  = re.compile(r'Gaia\s*DR3\s+(\d+)')

    all_chunks = []
    for start in range(0, len(source_ids), chunk_size):
        batch = source_ids[start:start+chunk_size]
        names = [f"Gaia DR3 {sid}" for sid in batch]

        result = simbad.query_objects(names)
        result.pprint()
        print(result.colnames)
        if result is None:
            continue

        rows = []
        for row in result:
            # decode the IDS field
            raw_ids = row['ids']
            ids_str = raw_ids.decode() if isinstance(raw_ids, bytes) else str(raw_ids)

            # extract each Gaia ID if present
            m2  = pat_dr2.search(ids_str)
            me3 = pat_edr3.search(ids_str)
            m3  = pat_dr3.search(ids_str)

            dr2_id   = int(m2.group(1)) if m2 else None
            edr3_id  = int(me3.group(1)) if me3 else None
            input_id = int(m3.group(1)) if m3 else None

            # SIMBAD name and spectral type
            raw_name = row['main_id']
            name     = raw_name.decode() if isinstance(raw_name, bytes) else str(raw_name)

            raw_sp   = row['main_id']
            sp_type  = raw_sp.decode().strip() if isinstance(raw_sp, bytes) else str(raw_sp)
            

            rows.append((input_id, dr2_id, edr3_id, name, sp_type))
        chunk_tbl = Table(rows,
                          names=('input_dr3','gaia_dr2','gaia_eDR3',
                                 'simbad_name','spectral_type'))
        all_chunks.append(chunk_tbl)

    if all_chunks:

        return vstack(all_chunks)
    else:
        # empty table with correct schema
        return Table(names=('input_dr3','gaia_dr2','gaia_eDR3',
                            'simbad_name','spectral_type'),
                     dtype=('i8','i8','i8','U50','U20'))


# Example usage:

# 1. Suppose you have an Astropy table `tbl` with a column 'GaiaEDR3' of integers:
#    source_ids = [int(x) for x in tbl['GaiaEDR3']]

# 2. Call the function:
#    simbad_ids = query_simbad_all_gaia_ids(source_ids)

# 3. Merge back on your input DR3:
#    merged = tbl.join(simbad_ids,
#                      left_on='GaiaEDR3',
#                      right_on='input_dr3',
#                      join_type='left')

# 4. Now `merged` has columns ['gaia_dr2','gaia_eDR3','simbad_name','spectral_type'] you can inspect or save.

unique_vals, unique_start=np.unique(input_table['GaiaDR3'],return_index=True)
input_table=input_table[unique_start]

simbad_matches=query_simbad_name_and_sptype(input_table['GaiaDR3'])
print(simbad_matches.colnames)
#simbad_matches=query_simbad_all_gaia_ids(input_table['GaiaEDR3'])
#simbad_matches['source_id']=simbad_matches['source_id'].astype('int64')
#print(simbad_matches['source_id'].dtype)
print(input_table['GaiaDR3'].dtype)
merged=join(input_table,simbad_matches, keys_left='GaiaDR3',keys_right='source_id', join_type='left')
#merged.remove_column('source_id_2')
#unique_inds=np.unique(merged['source_id'])
#merged_trim=merged[unique_inds]
#simbad_matches.pprint()

merged.pprint()

output_name=input_name.split('.')[0]+'_simbadadded.csv'
merged.write(output_name)
#merged_trim.write(output_name)

