"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-07-29 with the assistance of Microsoft's Copilot

This should get the white dwarfs in the Gentile Fusillo Catalogue for Gaia eDR3 that are sufficiently intrinsically faint for lithium to be present in their atmospheres and then it should get the Panstarrs photometry for all of those objects too. It should probably first filter on quality filters.

Ideally it would then save all of the faint white dwarfs with their associated panstarrs photometry... I suppose within this catalogue I could then have it somehow pull out and highlight the DZs with Li and then I could double check the panstarrs color space that those all fall inside.



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
from astropy.table import Table, Column, vstack, join
#import scipy.interpolate as scinterp
import time
import pyvo
from astroquery.vizier import Vizier

from astroquery.gaia import Gaia
from astroquery.xmatch import XMatch


sys.path.append('../')



GentileFusillo_catname="J/MNRAS/508/3877"
#gfviz=Vizier(catalog=GentileFusillo_catname,columns=['GaiaDR2'])
#print(Vizier(catalog=GentileFusillo_catname).get_catalog_metadata())

#from astroquery.vizier import Vizier
#from astroquery.gaia   import Gaia
from astroquery.mast   import Catalogs
import astropy.table    as tbl

gabslimit=14.8 #Li detection limit I used for the el-Badry selection based on LHS 2534
output_name='gf21_GaiaeDR3_faintWDs.csv'

Gaia.login(credentials_file='../Gaia_credentials.txt')

# 1. Query Gentile Fusillo (Gaia DR3) for Gmag > 15
Vizier.ROW_LIMIT = -1
Vizier.TIMEOUT   = 120
viz = Vizier(
    columns        = ["source_id", "Gmag"],
    column_filters = {"Gmag": "["+str(gabslimit)+",]"},
)
#gf = viz.get_catalogs("J/MNRAS/508/3877")
gf=Vizier.query_constraints(catalog=GentileFusillo_catname, Gmag='>'+str(gabslimit), RPlx='>10')[0]
gf.pprint()
# 2. Fetch Gaia→PS1 neighbor mapping via TAP
#ids_list = ",".join(map(str, gf["GaiaEDR3"]))
#adql = f"""
#SELECT source_id, original_ext_source_id
#FROM gaiadr3.panstarrs1_best_neighbour
#WHERE source_id IN ({ids_list})
#"""
#job   = Gaia.launch_job_async(adql)
#neigh = job.get_results()

ids=np.unique(gf["GaiaEDR3"])
chunks = np.array_split(ids, 10)  # ten ~equal pieces
results = []

for chunk in chunks:
    id_list = ",".join(map(str, chunk))
    #adql = f"""
      #SELECT source_id, original_ext_source_id
      #FROM gaiadr3.panstarrs1_best_neighbour
      #WHERE source_id IN ({id_list})
    #"""
    adql = f"""
      SELECT *
      FROM gaiadr3.panstarrs1_best_neighbour AS ps1
      WHERE ps1.source_id IN ({id_list})
      JOIN gaiadr2.panstarrs1_original_valid as ps1photo
      ON ps1.original_ext_source_id = ps1photo.obj_id
      
    """
    job = Gaia.launch_job_async(adql)
    results.append(job.get_results())

from astropy.table import vstack
neigh = vstack(results)

# 3. Merge Gentile Fusillo & neighbor tables on source_id
#gf_neigh = tbl.join(gf, neigh, keys="source_id")
gf_neigh = tbl.join(gf, neigh, keys_left='GaiaEDR3', keys_right="source_id")

## 4. Retrieve Pan-STARRS photometry by PS1 objID via MAST
#ps1_ids = list(set(gf_neigh["original_ext_source_id"]))
#ps1     = Catalogs.query_criteria(catalog="Panstarrs", objID=ps1_ids)

# 5. Final join: GF+Gaia info ↔ PS1 photometry
#final = tbl.join(
    #gf_neigh,
    #ps1,
    #left_on  = "original_ext_source_id",
    #right_on = "objID"
#)

final=gf_neigh

print(f"Total matches: {len(final)}")
print(final[:5])

final.write(output_name)
