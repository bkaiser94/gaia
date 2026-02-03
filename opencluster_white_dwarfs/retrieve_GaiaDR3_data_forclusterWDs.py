"""
Adapation by Ben Kaiser (UNC-Chapel Hill) 2025-12-18 that he made for the interesting WDs search. Basically, this will take the existing table's Gaia DR3 IDs and find them in the Gaia DR3 table to get the actual Gaia DR3 data.

Created by Ben Kaiser (UNC-Chapel Hill) 2025-08-01 with the assistance of Microsoft's Copilot, spun out of retrieve_faint_wds.py

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



#GentileFusillo_catname="J/MNRAS/508/3877"
#gfviz=Vizier(catalog=GentileFusillo_catname,columns=['GaiaDR2'])
#print(Vizier(catalog=GentileFusillo_catname).get_catalog_metadata())

#from astroquery.vizier import Vizier
from astroquery.gaia   import Gaia
from astroquery.mast   import Catalogs
import astropy.table    as tbl

gabslimit=14.8 #Li detection limit I used for the el-Badry selection based on LHS 2534
#output_name='gf21_GaiaeDR3_faintWDs_gaiaadded.csv'
input_name='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded_wdagesadded_agediffsig.csv'
output_name=input_name.split('.')[0]+'_gaiaadded.csv'

#Gaia.login(credentials_file='../Gaia_credentials.txt')

# 1. Query Gentile Fusillo (Gaia DR3) for Gmag > 15
Vizier.ROW_LIMIT = -1
Vizier.TIMEOUT   = 120

#gf = viz.get_catalogs("J/MNRAS/508/3877")
#It's the allcaps GMAG that is the absolute Gmag while Gmag is the apparent G magnitude. Yep, thanks a lot, Nicola.
#gf=Vizier.query_constraints(catalog=GentileFusillo_catname, GMAG='>'+str(gabslimit), RPlx='>10')[0]
#gf.pprint()
# 2. Fetch Gaia→PS1 neighbor mapping via TAP
#ids_list = ",".join(map(str, gf["GaiaEDR3"]))
#adql = f"""
#SELECT source_id, original_ext_source_id
#FROM gaiadr3.panstarrs1_best_neighbour
#WHERE source_id IN ({ids_list})
#"""
#job   = Gaia.launch_job_async(adql)
#neigh = job.get_results()

#gf=Table.read('gf21_GaiaeDR3_faintWDs.csv')
gf=Table.read(input_name)


ids=np.unique(gf["GaiaEDR3"])
print('n_ids:',ids.shape)
chunks = np.array_split(ids, 4)  # ten ~equal pieces
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
      FROM gaiadr3.gaia_source as gaia
    WHERE gaia.source_id IN ({id_list})
    """
    job = Gaia.launch_job_async(adql)
    results.append(job.get_results())

from astropy.table import vstack
new_gaia = vstack(results)



#for colname in new_gaia.colnames:
    #print(colname)
# 3. Merge Gentile Fusillo & neighbor tables on source_id
#gf_neigh = tbl.join(gf, neigh, keys="source_id")
gf_gaia = tbl.join(gf, new_gaia, keys_left='GaiaDR3', keys_right="source_id")

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

final=gf_gaia

print(f"Total matches: {len(final)}")
print(final[:5])

final.write(output_name)
