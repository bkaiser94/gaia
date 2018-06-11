"""
Created by Ben Kaiser (UNC- Chapel Hill) 2018-06-08

Take several input parameters (such as being located within 100 pc), and query the Gaia database for DR2 data.
Then return the output to a text file that can be read by other programs, such as one to make a color-magnitude
diagram.

Based on the tutorial examples available at https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
"""
import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table

"""
SELECT TOP 500 source_id,ra,ra_error,dec,dec_error,parallax,parallax_error,phot_g_mean_mag,bp_rp,radial_velocity,radial_velocity_error,phot_variable_flag,teff_val,a_g_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',0,0,360))=1    AND  (parallax>=0 AND phot_g_mean_mag<=6)
"""

output_name = 'top500_nearby_gaia.csv'

distance_limit = 100 #pc

parallax_min = 1./distance_limit #arcseconds of parallax
parallax_min= parallax_min*1e-3 #milliarcseconds of parallax

parallax_over_error =5. #max uncertainty allowed in data.


#search_statement= "SELECT ALL source_id,ra,ra_error,dec,dec_error,parallax,parallax_error,phot_g_mean_mag,bp_rp,radial_velocity,radial_velocity_error,phot_variable_flag,teff_val,a_g_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',179.99997916666666,0,359.9999583333333,180))=1    AND  (parallax>=1e-5 AND parallax_over_error>=5)"

start_string = "SELECT TOP 500 * FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',179.99997916666666,0,359.9999583333333,180))=1    AND  "

condition_string = "(parallax>=" +str(parallax_min) + " AND parallax_over_error>=" + str(parallax_over_error) + ")"

search_statement= start_string + condition_string

#search_statement= "SELECT TOP 100 source_id,ra,ra_error,dec,dec_error,parallax,parallax_error,phot_g_mean_mag,bp_rp,radial_velocity,radial_velocity_error,phot_variable_flag,teff_val,a_g_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',179.99997916666666,0,359.9999583333333,180))=1    AND  (parallax>=1e-5 AND parallax_over_error>=5)"



print(search_statement)

######
def coordinate_search():
    coordinate = coord.SkyCoord(ra = 280, dec = -60, unit=(u.degree, u.degree), frame = 'icrs')
    width = 0.1*u.deg
    height = 0.1*u.deg
    r = Gaia.query_object_async(coordinate =coordinate, width =width, height = height)
    r.pprint()
    
    
def cone_search():
    coordinate = coord.SkyCoord(ra = 280, dec = -60, unit=(u.degree, u.degree), frame = 'icrs')
    radius = 0.1*u.deg
    j= Gaia.cone_search_async(coordinate, radius)
    r=j.get_results()
    r.pprint()
    
def get_public_tables(tap_plus = True):
    if tap_plus:
        tables= Gaia.load_tables(only_names = True)
        for table in tables:
            print(table.get_qualified_name())
    else:
        tables = Gaia.load_tables()
        for table in tables:
            print(table.get_qualified_name())
            
def load_table(table_name):
    """
    Re-open a previously saved gaia table and mess with it
    """
    table = Gaia.load_table(table_name)
    print(table)
    #to inspect columns
    for column in (table.get_columns()):
        print(column.get_name())
    
def synchronous_query(file_dump = False):
    """
    Results are not stored server-side. No more than 2000 rows are allowed in this method. For more than 2000 rows, asynchronous queries must be used.
    """
    if file_dump:
        job = Gaia.launch_job("select top 100 \
            solution_id, ref_epoch, ra_dec_corr, astrometric_n_obs_al, matched_observations, duplicated_source, phot_variable_flag \
            from gaiadr1.gaia_source order by source_id", dump_to_file = True) #original from example
        print(job)
        print('next thing')
        r= job.get_results()
        print(r)
    else:
        job = Gaia.launch_job("select top 100 \
            solution_id, ref_epoch, ra_dec_corr, astrometric_n_obs_al, matched_observations, duplicated_source, phot_variable_flag \
            from gaiadr1.gaia_source order by source_id") #original from example
        #job = Gaia.launch_job("select top 100 \
            #solution_id, ref_epoch, ra_dec_corr, astrometric_n_obs_al, matched_observations, duplicated_source, phot_variable_flag \
            #from gaiadr2.gaia_source order by source_id")
        print(job)
        r = job.get_results()
        print(r['solution_id'])
        print(r)
        
        
def synchronous_query_uptab():
    """
    I'm not totally clear on how this works since the xml file they're using isn't actually revealed.
    Perhaps it's got fields corresponding to those for the search, but who's to know?
    """
    upload_resource= 'my_table.xml'
    j = Gaia.launch_job(query = 'select * from tap_upload.table_test', upload_resource= upload_resource, \
        upload_table_name = 'table_test', verbose =True)
    r= j.get_results
    r.pprint()
    
    


def asynchronous_query():
    #job = Gaia.launch_job_async("select top 100 * from gaiadr2.gaia_source order by source_id")
    #job = Gaia.launch_job_async("select top 100 where parallax <= 1e-5 from gaiadr2.gaia_source order by source_id")
    #print(search_statement)
    job = Gaia.launch_job_async(search_statement)
    r= job.get_results()
   
    #print(r['solution_id'])
    r.pprint() #the pretty printing method for astropy tables
    return r
#############3
#coordinate_search()

#cone_search()

#get_public_tables()
#get_public_tables(tap_plus= False)

#synchronous_query()

#synchronous_query_uptab()


output_table = asynchronous_query()
output_table.write(output_name, format = 'ascii.csv')



