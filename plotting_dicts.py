"""
This is an externally located file that will contain the dictionaries called by plot_gaia_cmd.py for plotting multiple target lists at one time. It should include the file name to be plotted, the color of the points for that file, and the label that should be associated with those points.

"""

target_lists=['elm_survey_gaia.csv',
    'pre_elms_gaia.csv',
    'BLAPs_gaia.csv',
    'pulsar_companions_gaia.csv',
    'hot_wind_wds_gaia.csv',
    'hv_wds_gaia.csv']

plot_dict={'BLAPs_gaia.csv':{
    "color":'b',
    "label":'BLAP'},
    'pulsar_companions_gaia.csv':{
        'color':'#1ca1f2',
        'label':'PSR'},
    'hot_wind_wds_gaia.csv':{
        'color':'magenta',
        'label':'UHE WDs'},
    'hv_wds_gaia.csv':{
        'color':'g',
        'label':'Hyper Velocity WDs'},
    'pre_elms_gaia.csv':{
        'color':'#66ff33',
        'label':'Pre-ELM'},
    'elm_survey_gaia.csv':{
        'color':'#9966ff',
        'label':'ELM Survey'}}


#absmag:{ colours[0]:
cmd_border={'g':{
    'bp':{'x':[-0.7, 5], 'y':[-3,16]},
    'g':{'x':[-0.5,2], 'y':[-3,16]}}}
