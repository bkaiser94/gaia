"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-09-06


Take a Gaia dataset (that includes all of the bells and whistles) and generate galactic kinematic quantities if they
aren't known and save them into the file in addition to the previously available values. If they are already known,
just go ahead and make the plots that show memberships and stuff.

This is probably not going to work correctly for awhile.

THIS SHOULD BE RUN IN PYTHON 3 BECAUSE ASTROPY NEEDS TO BE ON ITS GAME!!!

According to astropy documentation
https://docs.astropy.org/en/stable/api/astropy.coordinates.builtin_frames.GalacticLSR.html#astropy.coordinates.builtin_frames.GalacticLSR.v_bary
GalacticLSR already assumes the Schoenrich et al (2010) barycentric velocity relative to LSR which is the one I have in v_sun_pec, so I'm not going to try to provide it as an input.

"""

from __future__ import print_function
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable, Column
import matplotlib.pyplot as plt
import scipy.stats as scistats
import matplotlib.patches as mp


#import plotting_dicts as pod


###########
#Galactic parameters from Gaia collaboration (D. Katz et al. 2018) "Mapping the Milky Way Disc Kinematics"
sun_height= 27.*u.pc
sun_dist= 8.34*u.kpc
circular_v_at_sun=240.*u.km/u.s
v_sun_pec= [11.1,12.24, 7.25] *u.km/u.s
v_sun_bobylev=[7.90, 11.73, 7.39]*u.km/u.s #Bobylev 2017 value, used by Torres et al. 2019


v_sun_tot= coord.CartesianDifferential([v_sun_pec[0], v_sun_pec[1]+circular_v_at_sun, v_sun_pec[2]])
v_sun_bobylev_diff= coord.CartesianDifferential(v_sun_bobylev)

gc_frame= coord.Galactocentric(galcen_distance=sun_dist,
                               galcen_v_sun=v_sun_tot,
                               z_sun = sun_height)

galcart= coord.Galactic(representation_type= coord.CartesianRepresentation, differential_type=coord.CartesianDifferential)

#altGalacticLSR=coord.GalacticLSR(v_bary=v_sun_bobylev_diff)

altGalacticLSR=coord.GalacticLSR(v_bary=v_sun_bobylev_diff, differential_type='cartesian')


#galLSR_base=coord.GalacticLSR
galLSR_base=altGalacticLSR

############
#disk kinematics
thin_disk_pec= [2.5, -8.9, -1.0]*u.km/u.s
thin_disk_disp= [30.8, 17.4, 14.0] *u.km/u.s


##############3
print('\n\n YOU BETTER BE USING PYTHON3!\n\n\n')

list_color = '#1ca1f2'


#######3error distribution variables
mc_number = 100000
percent_off = 34 #1-sigma equivalent
#############
rv_zero=0.
rv_sigma=100.
#####





target_input='20190829_DZNas.csv'

target_table = Table.read(target_input)

##############################
#Torres et al. 2019 values
thin_disc_vel=[2.52, -1.45, 2.60]
thin_disc_disp=[22.5, 17.4, 16.2]
thick_disc_vel=[-18.58,-30.03, 1.06]
thick_disc_disp=[50.4, 29.3, 33.1]
halo_vel=[-27.31, -92.29, 3.41]
halo_disp=[100.8, 67.4, 66.9]

###############################################
###############################################
###############################################

def get_mc_distribution(value, error):
    error_distribution = np.random.normal(loc= value, scale = error, size = mc_number)
    #plt.hist(error_distribution)
    #plt.axvline(x=value, color='r')
    #plt.axvline(x=np.median(error_distribution), linestyle='--', color='k')
    #plt.show()
    #print('value:', value, 'error:',error, 'median:', np.median(error_distribution))
    return error_distribution


def remove_negative(array, verbose= True):
    output_array = array[np.where(array>0)]
    if (verbose and array.shape[0]-output_array.shape[0] >0):
        print('Removed ' +str(array.shape[0]-output_array.shape[0]) + ' negatives')
    return output_array

def match_sizes(change_array, match_array):
    """
    Intended to keep compatibility with an array that has had negatives removed
    """
    try:
        min_inds = np.nanmin([change_array.shape[0], match_array.shape[0]])
        
        return change_array[:min_inds], match_array[:min_inds]
    except AttributeError:
        #the inputs aren't actually arrays
        return change_array, match_array


def get_errors(distribution, percent_off = percent_off):
    """
    values for the error bars on the plot
    
    Returns
    [- value, + value], so to get the points on the plot where they'd be located you do
    value - [- value] , value + [- value]
    Basically these are the width of the uncertainty range on either side.
    """
    low_bar = np.nanpercentile(distribution, 50-percent_off)
    median = np.nanmedian(distribution)
    high_bar = np.nanpercentile(distribution, 50+percent_off)
    try:
        return np.array([[median-low_bar],[high_bar-median]])
    #except astropy.units.core.UnitsError as error:
    except u.core.UnitsError as error:
        return  np.array([[median.value-low_bar],[high_bar-median.value]])
    
def test_errorbars(group_vel, group_disp, label):
    group_vel=np.array(group_vel)
    group_disp=np.array(group_disp)
    u_dist= get_mc_distribution(group_vel[0], group_disp[0])
    v_dist= get_mc_distribution(group_vel[1], group_disp[1])
    w_dist=get_mc_distribution(group_vel[2], group_disp[2])
    uw_dist= np.sqrt(v_dist**2+w_dist**2)
    
    uw_mean= np.nanmean(uw_dist)
    uw_std= np.std(uw_dist)
    print(label, 'UW', uw_mean, '+/-', uw_std)
    print('quadrature version', np.sqrt(group_vel[0]**2+group_vel[2]**2),np.sqrt(group_disp[0]**2+group_disp[2]**2))
    print(label, 'V MC', np.nanmean(v_dist), '+/-', np.std(v_dist), np.median(v_dist))
    print(label, 'V reported', group_vel[1], group_disp[1])
    #plt.hist(uw_dist)
    #plt.show()
    return [group_vel[1], group_disp[1]], [uw_mean, uw_std]


def generate_toomre_diagram():
    #fig= plt.figure()
    def plot_kinematic_group(group_vel, group_disp, label):
        group_v, group_uw= test_errorbars(group_vel, group_disp, label)
        #plt.errorbar(group_vel[1], np.sqrt(group_vel[0]**2+group_vel[2]**2), xerr= group_disp[1], yerr= np.sqrt(group_disp[0]**2+group_disp[2]**2),label=label)
        plt.errorbar(group_v[0], group_uw[0], xerr= group_v[1], yerr= group_uw[1] ,label=label, marker='s')
        return
    td_v, td_uw=test_errorbars(thin_disc_vel, thin_disc_disp,'thin disk')
    plot_kinematic_group(thin_disc_vel, thin_disc_disp, 'thin disk')
    plot_kinematic_group(thick_disc_vel, thick_disc_disp, 'thick disk')
    plot_kinematic_group(halo_vel, halo_disp, 'halo')
    #generate 3-sigma and 5-sigma ellipse
    sig3_ellipse= mp.Ellipse((td_v[0], td_uw[0]), td_v[1]*3*2, td_uw[1]*3*2, fill=False, edgecolor='k', linestyle='--')
    sig5_ellipse= mp.Ellipse((td_v[0], td_uw[0]), td_v[1]*5*2, td_uw[1]*5*2, fill=False, edgecolor='k', linestyle='--')
    ax=plt.axes()
    ax.add_patch(sig3_ellipse)
    ax.add_patch(sig5_ellipse)
    plt.xlabel('V (km/s)')
    plt.ylabel(r'$(U^2+W^2)^{1/2}$ (km/s)')
    plt.legend(loc='best')
    plt.xlim(-300, 100)
    plt.show()
    

#test_errorbars(thin_disc_vel, thin_disc_disp, 'thin disk')
#test_errorbars(thick_disc_vel, thick_disc_disp, 'thick disk')
#test_errorbars(halo_vel, halo_disp, 'halo')
###############################################
###############################################
###############################################



#astropy says you can create sky coords just using an array, so maybe I should just be doing that.
def get_galactic_coords(row):
    star_coord= coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec= row['pmra']*u.mas/u.yr, pm_dec= row['pmdec']*u.mas/u.yr, radial_velocity=0. *u.km/u.s, distance= 1000./row['parallax'] *u.pc , frame='icrs')
    #print('l', row['l'], 'b', row['b'])
    #print(star_coord.galactic)
    #print(star_coord.galactic.l.value - row['l'], star_coord.galactic.b.value - row['b'])
    #galactic_coords= star_coord.transform_to(coord.Galactocentric)
    galactic_coords= star_coord.transform_to(gc_frame)
    print('\n\n=====')
    print(row['name'])
    print(galactic_coords)
    return galactic_coords

def get_galLSR_coords(row, do_mc=False, vary_rv=False):
    star_coord= coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec= row['pmra']*u.mas/u.yr, pm_dec= row['pmdec']*u.mas/u.yr, radial_velocity=0. *u.km/u.s, distance= 1000./row['parallax'] *u.pc , frame='icrs')
    #star_coord= coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec= row['pmra']/np.cos(row['dec']/180*np.pi)*u.mas/u.yr, pm_dec= row['pmdec']*u.mas/u.yr, radial_velocity=0. *u.km/u.s, distance= (1000./row['parallax']) *u.pc, frame='icrs')
    #galLSR_coords= star_coord.transform_to(coord.GalacticLSR)
    galLSR_coords= star_coord.transform_to(galLSR_base)

    if (do_mc and (not vary_rv)) :
        print('Not varying the RV')
        pmra_dist=get_mc_distribution(row['pmra'], row['pmra_error'])
        pmdec_dist=get_mc_distribution(row['pmdec'], row['pmdec_error'])
        parallax_dist= get_mc_distribution(row['parallax'],row['parallax_error'])
        parallax_dist=remove_negative(parallax_dist)
        pmdec_dist, parallax_dist=match_sizes(pmdec_dist, parallax_dist)
        pmra_dist, parallax_dist= match_sizes(pmra_dist, parallax_dist)
        dist_coord=coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec=pmra_dist*u.mas/u.yr, pm_dec= pmdec_dist*u.mas/u.yr, radial_velocity=0. *u.km/u.s, distance= 1000./parallax_dist *u.pc , frame='icrs')
        #dist_coord=dist_coord.transform_to(coord.GalacticLSR)
        dist_coord=dist_coord.transform_to(galLSR_base)
        #print(np.sum(pmdec_dist-dist_coord.pm_dec.value))
        return galLSR_coords, dist_coord
    elif (do_mc and vary_rv) :
        pmra_dist=get_mc_distribution(row['pmra'], row['pmra_error'])
        pmdec_dist=get_mc_distribution(row['pmdec'], row['pmdec_error'])
        parallax_dist= get_mc_distribution(row['parallax'],row['parallax_error'])
        rv_dist=get_mc_distribution(rv_zero, rv_sigma)
        parallax_dist=remove_negative(parallax_dist)
        pmdec_dist, parallax_dist=match_sizes(pmdec_dist, parallax_dist)
        pmra_dist, parallax_dist= match_sizes(pmra_dist, parallax_dist)
        rv_dist, parallax_dist=match_sizes(rv_dist, parallax_dist)
        dist_coord=coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec=pmra_dist*u.mas/u.yr, pm_dec= pmdec_dist*u.mas/u.yr, radial_velocity=rv_dist *u.km/u.s, distance= 1000./parallax_dist *u.pc , frame='icrs')
        #dist_coord=dist_coord.transform_to(coord.GalacticLSR)
        dist_coord=dist_coord.transform_to(galLSR_base)
        return galLSR_coords, dist_coord
    else:
        return galLSR_coords
        #for pmra, pmdec, parallax in zip(pmra_dist, pmdec_dist, parallax_dist):
            #dist_coord=coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec= pmdec*u.mas/u.yr, radial_velocity=0. *u.km/u.s, distance= 1000./parallax *u.pc , frame='icrs')
    #galLSR_coords= star_coord.transform_to(galcart)
    
    #return galLSR_coords

def gal_vel_curve(l, d):
    """
    From page 912 of Carroll and Ostlie (B.O.B.)
    
    
    """
    l=l/180.*np.pi
    A=14.8*u.km/u.s/u.kpc
    B= -12.4*u.km/u.s/u.kpc
    v_rad= A*d*np.sin(2*l)
    v_tan=A*d*np.cos(2*l)+B*d
    return v_rad, v_tan

#def get_v_phi(gal_coord):
    #v_phi= gal_coord.
    #return v_phi
#for row in target_table:
    #get_galactic_coords(row)
    
#galLSR_coords= get_galLSR_coords(target_table).differential_type(coord.CartesianDifferential)
#galLSR_coords= get_galLSR_coords(target_table)

def plot_values(target_table, plot_vals=['V', 'UW'], do_mc=True, vary_rv=False, color=list_color):
    legend_needed=True
    for row in target_table:
        galLSR_single, galLSR_dist= get_galLSR_coords(row, do_mc=do_mc, vary_rv=vary_rv)
        galLSR_vel_dist=galLSR_dist.velocity
        U_dist=galLSR_vel_dist.d_x
        V_dist= galLSR_vel_dist.d_y
        W_dist= galLSR_vel_dist.d_z
        UW_dist=np.sqrt(U_dist**2. +W_dist**2.)
        V_err= get_errors(V_dist)
        UW_err=get_errors(UW_dist)
        U=galLSR_single.velocity.d_x
        V=galLSR_single.velocity.d_y
        W=galLSR_single.velocity.d_z
        UW=np.sqrt(U**2. +W**2.)
        print('\n')
        print(row['name'], U, V, W)
        print('+/-', get_errors(U_dist), get_errors(V_dist),get_errors(W_dist))
        #plt.hist(UW_dist.value)
        #plt.show()
        if plot_vals[0]=='V':
            xvals=V
            xerr=V_err
            x_dist=V_dist
        if plot_vals[1]=='UW':
            yvals=UW
            yerr=UW_err
            y_dist=UW_dist
        if plot_vals[0]=='U':
            xvals=U
            xerr=get_errors(U_dist)
            x_dist=U_dist
        if plot_vals[1]=='W':
            yvals=W
            yerr=get_errors(W_dist)
            y_dist=W_dist
        if plot_vals[1]=='U':
            yvals=U
            yerr=get_errors(U_dist)
            y_dist=U_dist
        xvals=xvals.value
        yvals=yvals.value
        xerr=xerr
        yerr=yerr
        x_dist=x_dist.value
        y_dist=y_dist.value
        #plt.errorbar(xvals, yvals, xerr=xerr, yerr=yerr, marker='o', color=color)
        if (legend_needed and vary_rv):
            plt.errorbar(np.median(x_dist), np.median(y_dist), xerr=xerr, yerr=yerr, marker='o', color=color, label='Gaussian RV dist. around 0')
            legend_needed=False
        if (legend_needed and not vary_rv):
            plt.errorbar(np.median(x_dist), np.median(y_dist), xerr=xerr, yerr=yerr, marker='o', color=color, label=r'RV$\equiv$0')
            legend_needed=False
        else:
            plt.errorbar(np.median(x_dist), np.median(y_dist), xerr=xerr, yerr=yerr, marker='o', color=color)
        plt.annotate(str(row['name']),xy=(np.median(x_dist), np.median(y_dist)), xycoords='data', xytext=(np.median(x_dist), np.median(y_dist)), textcoords= 'data' , fontsize=6, color =color)
        #plt.annotate(str(row['name']),xy=(xvals, yvals), xycoords='data', xytext=(xvals,yvals), textcoords= 'data' , fontsize=6, color =color)
    plt.xlabel(plot_vals[0])
    plt.ylabel(plot_vals[1])
    return

#plot_values(target_table, plot_vals=['V','W'], color='r', vary_rv=True)
#plot_values(target_table, plot_vals=['V','W'], color=list_color)

#plt.show()

#plot_values(target_table, plot_vals=['V','U'], color='r', vary_rv=True)
#plot_values(target_table, plot_vals=['V','U'], color=list_color)

#plt.show()

#plot_values(target_table, plot_vals=['V','UW'], color='r', vary_rv=True)
plot_values(target_table, plot_vals=['V','UW'], color=list_color)
generate_toomre_diagram()

plt.show()

all_list=[]
uw_err_list=[]
v_err_list=[]

for row in target_table:
    galLSR_single, galLSR_dist= get_galLSR_coords(row, do_mc=True, vary_rv=False)
    
    #output_vel=get_galLSR_coords(row,do_mc=False, vary_rv=False)
    output_vel=galLSR_single.velocity
    #output_row=[row['name'], row['ra'],row['dec'],output_vel.d_x.to(u.km/u.s).value, output_vel.d_y.to(u.km/u.s).value, output_vel.d_z.to(u.km/u.s).value, 1000./row['parallax'], row['parallax']]
    #all_list.append(output_row)
    
    
    galLSR_vel_dist=galLSR_dist.velocity
    U_dist=galLSR_vel_dist.d_x
    V_dist= galLSR_vel_dist.d_y
    W_dist= galLSR_vel_dist.d_z
    
    UW_dist=np.sqrt(U_dist**2. +W_dist**2.)
    V_err= get_errors(V_dist)
    UW_err=get_errors(UW_dist)
    U=galLSR_single.velocity.d_x
    V=galLSR_single.velocity.d_y.value
    W=galLSR_single.velocity.d_z
    UW=np.sqrt(U**2. +W**2.).value
    
    output_row=[row['name'], row['ra'],row['dec'], output_vel.d_y.to(u.km/u.s).value, np.median(UW_dist.value), 1000./row['parallax'], row['parallax'], V_err[0,0], V_err[1,0], UW_err[0,0], UW_err[1,0]]
    all_list.append(output_row)
    
    
    #v_err_list.append([V_err])
    #uw_err_list.append([UW_err])
  
    print('\n\n')
    print(row['name'])
    print('U', 'cat. val:',U, 'MC dist med:',np.median(U_dist.value))
    print('V', 'cat. val:',V, 'MC dist med:',np.median(V_dist.value))
    print('W', 'cat. val:',W, 'MC dist med:',np.median(W_dist.value))
    print('**********\n\n')
    #plt.hist(UW_dist.value)
    #plt.show()
    plt.errorbar(np.median(V_dist.value),np.median(UW_dist.value), xerr=V_err, yerr=UW_err, marker='o', color='r')
    plt.annotate(str(row['name']),xy=(np.median(V_dist.value),np.median(UW_dist.value)), xycoords='data', xytext=(np.median(V_dist.value),np.median(UW_dist.value)), textcoords= 'data' , fontsize=6, color ='r')


all_list_array=np.array(all_list)
#v_err_col= Column(v_err_list)
#uw_err_col=Column(uw_err_list)

#v_err_col.pprint()

#text_names=['name','ra','dec','u', 'v', 'w', 'distance', 'parallax']
text_names=['name','ra','dec','v', 'uw2', 'distance', 'parallax', 'v_err_lo', 'v_err_hi','uw2_err_lo', 'uw2_err_hi']
#text_header=','.join(text_names)
#np.savetxt('20200517_output_velocities_test.csv',all_list_array, delimiter=',')
#np.savetxt('J1644_paper_outputs.csv',all_list, delimiter=',', header=text_header)
output_table=Table(all_list_array, names=text_names)

#output_table.add_column(v_err_col, name='v_err')
#output_table.add_column(uw_err_col, name='uw2_err')
output_name='kinematics_outputs.csv'
print('\n\nSaving', output_name, '\n\n')
output_table.write(output_name)

for row in target_table:
    galLSR_single, galLSR_dist= get_galLSR_coords(row, do_mc=True)
    galLSR_vel_dist=galLSR_dist.velocity
    U_dist=galLSR_vel_dist.d_x
    V_dist= galLSR_vel_dist.d_y
    W_dist= galLSR_vel_dist.d_z
    UW_dist=np.sqrt(U_dist**2. +W_dist**2.)
    V_err= get_errors(V_dist)
    UW_err=get_errors(UW_dist)
    U=galLSR_single.velocity.d_x
    V=galLSR_single.velocity.d_y.value
    W=galLSR_single.velocity.d_z
    UW=np.sqrt(U**2. +W**2.).value
    
    #plt.scatter(galLSR_dist.pm_dec.value, UW_dist.value)
    #plt.title('UW vs. pmdec')
    #plt.show()
    
    #plt.scatter(galLSR_dist.pm_ra_cosdec.value, UW_dist.value)
    #plt.title('UW vs. pmra')
    #plt.show()
    
    #plt.scatter(galLSR_dist.distance.value, UW_dist.value)
    #plt.title('UW vs. distance')
    #plt.show()
    
    #plt.scatter(galLSR_dist.pm_dec.value, V_dist.value)
    #plt.title('V vs. pmdec')
    #plt.show()
    
    #plt.scatter(galLSR_dist.pm_ra_cosdec.value, V_dist.value)
    #plt.title('V vs. pmra')
    #plt.show()
    
    #plt.scatter(galLSR_dist.distance.value, V_dist.value)
    #plt.title('V vs. distance')
    #plt.show()
    
    plt.errorbar(np.median(V_dist.value),np.median(UW_dist.value), xerr=V_err, yerr=UW_err, marker='o', color='g')
    plt.annotate(str(row['name']),xy=(np.median(V_dist.value),np.median(UW_dist.value)), xycoords='data', xytext=(np.median(V_dist.value),np.median(UW_dist.value)), textcoords= 'data' , fontsize=6, color ='g')
    plt.errorbar(V,UW, xerr=V_err, yerr=UW_err, marker='o', color=list_color)
    #plt.hist(V_dist.value)
    #plt.axvline(x=np.median(V_dist.value),linestyle='--', color='k')
    #plt.axvline(x=V, linestyle='--', color=list_color)
    #plt.title(row['name'])
    #plt.show()

galLSR_coords= get_galLSR_coords(target_table)
galLSR_vel= galLSR_coords.velocity
U=galLSR_vel.d_x
V= galLSR_vel.d_y
W=galLSR_vel.d_z
print(target_table['name'])
print('U', U)
print('V', V)
print('W', W)
UW= np.sqrt(U**2. +W**2.)
for row, v, uw in zip(target_table,V, UW):
    plt.annotate(str(row['name']),xy=(v.value, uw.value), xycoords='data', xytext=(v.value, uw.value), textcoords= 'data' , fontsize=6, color =list_color)
plt.xlabel('V (km/s)')
plt.ylabel(r'$(U^2+W^2)^{1/2}$ (km/s)')
plt.title(target_input + ' kinematics')
plt.scatter(v_sun_pec[1], np.sqrt(v_sun_pec[0]**2+v_sun_pec[2]**2), color='orange')
plt.xlim(-550, 150)
plt.ylim(0, 450)
plt.show()
    

galLSR_vel= galLSR_coords.velocity
#galLSR_vel.write(target_input.split('.')[0]+'_velocities.csv')
print('galLSR_coords')
print(galLSR_coords)
print(galLSR_vel)
print(galLSR_coords.velocity)
print(type(galLSR_vel))
U=galLSR_vel.d_x
V= galLSR_vel.d_y
W=galLSR_vel.d_z
UW= np.sqrt(U**2. +W**2.)
v_tot=np.sqrt(U**2. +V**2. + W**2.)


plt.scatter(V, np.sqrt(U**2. +W**2.))
for row, v, uw in zip(target_table,V, UW):
    plt.annotate(str(row['name']),xy=(v.value, uw.value), xycoords='data', xytext=(v.value, uw.value), textcoords= 'data' , fontsize=6, color =list_color)
plt.scatter(v_sun_pec[1], np.sqrt(v_sun_pec[0]**2+v_sun_pec[2]**2), color='yellow')
plt.xlabel('V (km/s)')
plt.ylabel(r'$(U^2+W^2)^{1/2}$ (km/s)')
plt.title(target_input + ' kinematics')
plt.xlim(-550, 150)
plt.ylim(0, 450)
plt.show()



#plt.scatter(target_table['l'], v_tot)
#for row, v  in zip(target_table,v_tot):
    #plt.annotate(str(row['name']),xy=(row['l'],v.value), xycoords='data', xytext=(row['l'], v.value), textcoords= 'data' , fontsize=6, color =list_color)
#plt.ylabel(r'$(U^2+V^2+W^2)^{1/2}$ (km/s)')
#plt.xlabel(r'$l$')
#plt.ylim(0,500)
#plt.show()


l_vals= np.linspace(0,360, 1000)

plt.scatter(target_table['l'], V)
v_rad, v_tan= gal_vel_curve(l_vals, np.max(1./target_table['parallax'])*u.kpc)
plt.plot(l_vals, v_tan, label='v_tan')
plt.plot(l_vals, v_rad, label='v_rad')
for row, v  in zip(target_table,V):
    plt.annotate(str(row['name']),xy=(row['l'],v.value), xycoords='data', xytext=(row['l'], v.value), textcoords= 'data' , fontsize=6, color =list_color)
plt.ylabel(r'$(V)$ (km/s)')
plt.xlabel(r'$l$')
plt.xlim(0,360)
plt.legend()
plt.show()


plt.scatter(V, U, color=list_color)
plt.xlabel('V (km/s)')
plt.ylabel('U (km/s)')
for row, v, o  in zip(target_table, V, U):
    plt.annotate(str(row['name']),xy=(v.value, o.value), xycoords='data', xytext=(v.value, o.value), textcoords= 'data' , fontsize=6, color =list_color)
plt.title(target_input + ' kinematics')
plt.show()


galactic_coords= get_galactic_coords(target_table)

#print(galactic_coords[0])
#plt.scatter(galactic_coords.v_x, galactic_coords.v_y)
#for row, gal_coord  in zip(target_table, galactic_coords):
    #plt.annotate(str(row['name']),xy=(gal_coord.v_x.value, gal_coord.v_y.value), xycoords='data', xytext=(gal_coord.v_x.value, gal_coord.v_y.value), textcoords= 'data' , fontsize=8, color =list_color)
#plt.xlabel('v_x')
#plt.ylabel('v_y')
#plt.show()

plt.scatter(galactic_coords.x.to(u.kpc), galactic_coords.y, label="WD")
plt.scatter(-1*sun_dist, 0., color='orange', label='sun')
#plt.quiver(galactic_coords.x.to(u.kpc), galactic_coords.y, galactic_coords.v_x, galactic_coords.v_y, headwidth=3, width=0.005)
#plt.quiver(galactic_coords.x.to(u.kpc), galactic_coords.y,U, V, headwidth=2, width=0.005, label='proj. pec. velocity (km/s)')
plt.quiver(galactic_coords.x.to(u.kpc), galactic_coords.y,U, -V, headwidth=2, width=0.005, label='proj. pec. velocity (km/s)') #assuming the y-axis will be inverted. The quiver plot doesn't actually correct the arrows for some reason....
for row, gal_coord  in zip(target_table, galactic_coords):
    plt.annotate(str(row['name']),xy=(gal_coord.x.to(u.kpc).value, gal_coord.y.value), xycoords='data', xytext=(gal_coord.x.to(u.kpc).value, gal_coord.y.value), textcoords= 'data' , fontsize=8, color =list_color)
plt.xlabel('X (kpc)')
plt.ylabel('Y (pc)')
plt.legend(loc='best')
plt.gca().invert_yaxis()

plt.title(target_input + ' kinematics')
plt.show()


plt.scatter(galactic_coords.x.to(u.kpc), galactic_coords.z, label="WD")
plt.scatter(-1*sun_dist, sun_height, color='orange',label='sun')
#plt.quiver(galactic_coords.x.to(u.kpc), galactic_coords.z, galactic_coords.v_x, galactic_coords.v_z, headwidth=3, width=0.005)
plt.quiver(galactic_coords.x.to(u.kpc), galactic_coords.z, U, W, headwidth=2, width=0.005, label='proj. pec. velocity (km/s)')
for row, gal_coord  in zip(target_table, galactic_coords):
    plt.annotate(str(row['name']),xy=(gal_coord.x.to(u.kpc).value, gal_coord.z.value), xycoords='data', xytext=(gal_coord.x.to(u.kpc).value, gal_coord.z.value), textcoords= 'data' , fontsize=8, color =list_color)

print('difference:', U-galactic_coords.v_x)
plt.xlabel('X (kpc)')
plt.ylabel('Z (pc)')
plt.title(target_input + ' kinematics')
plt.legend(loc='best')
plt.show()



#plt.scatter(galactic_coords.x.to(u.pc), galactic_coords.z, label="WD")
#plt.scatter(-1*sun_dist*1000., sun_height, color='orange',label='sun')
##plt.quiver(galactic_coords.x.to(u.kpc), galactic_coords.z, galactic_coords.v_x, galactic_coords.v_z, headwidth=3, width=0.005)
#plt.quiver(galactic_coords.x.to(u.pc), galactic_coords.z, U, W, headwidth=2, width=0.005, label='proj. pec. velocity (km/s)')
#for row, gal_coord  in zip(target_table, galactic_coords):
    #plt.annotate(str(row['name']),xy=(gal_coord.x.to(u.pc).value, gal_coord.z.value), xycoords='data', xytext=(gal_coord.x.to(u.pc).value, gal_coord.z.value), textcoords= 'data' , fontsize=8, color =list_color)

#print('difference:', U-galactic_coords.v_x)
#plt.xlabel('X (pc)')
#plt.ylabel('Z (pc)')
#plt.title(target_input + ' kinematics')
#plt.legend(loc='best')
#plt.show()

plt.scatter(galactic_coords.z, galactic_coords.v_z)
for row, gal_coord  in zip(target_table, galactic_coords):
    plt.annotate(str(row['name']),xy=(gal_coord.z.value, gal_coord.v_z.value), xycoords='data', xytext=(gal_coord.z.value, gal_coord.v_z.value), textcoords= 'data' , fontsize=8, color =list_color)
plt.xlabel('z')
plt.ylabel('v_z')
plt.show()
