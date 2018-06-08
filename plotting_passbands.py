"""
Created by Ben Kaiser (UNC-Chapel Hill)

Plot the DR2 passbands and the "Revised DR2 passbands" for the GAIA mission.

Requires one to already have downloaded the gaia passband files from the gaia website and update the file paths to get to them.

"""


import numpy as np
import matplotlib.pyplot as plt

dr2_passband_file = 'GaiaDR2_Passbands_ZeroPoints/GaiaDR2_Passbands.dat'
dr2_rev_passband_file = 'GaiaDR2_Revised_Passbands_ZeroPoints/GaiaDR2_RevisedPassbands.dat'


dr2_all = np.genfromtxt(dr2_passband_file).T
rev_all = np.genfromtxt(dr2_rev_passband_file).T
print dr2_all.shape
wavelengths_dr2 = dr2_all[0]
wavelengths_rev = rev_all[0]
print wavelengths_dr2.shape

def good_for_plots_dr2(bandpass, sigma):
    good_vals = np.where(bandpass < 99)
    return wavelengths_dr2[good_vals], bandpass[good_vals], sigma[good_vals]

def good_for_plots_rev(bandpass, sigma):
    good_vals = np.where(bandpass < 99)
    return wavelengths_rev[good_vals], bandpass[good_vals], sigma[good_vals]


Gband_dr2 = dr2_all[1]
Gband_rev = rev_all[1]

Gband_sig_dr2 = dr2_all[2]
Gband_sig_rev = rev_all[2]

BPband_dr2 = dr2_all[3]
BPband_rev = rev_all[3]

BPband_sig_dr2 = dr2_all[4]
BPband_sig_rev = rev_all[4]

RPband_dr2 = dr2_all[5]
RPband_rev = rev_all[5]

RPband_sig_dr2 = dr2_all[6]
RPband_sig_rev = rev_all[6]

dr2_G_tuple = good_for_plots_dr2(Gband_dr2, Gband_sig_dr2)
rev_G_tuple = good_for_plots_rev(Gband_rev, Gband_sig_rev)

dr2_RP_tuple = good_for_plots_dr2(RPband_dr2, RPband_sig_dr2)
rev_RP_tuple = good_for_plots_rev(RPband_rev, RPband_sig_rev)

dr2_BP_tuple = good_for_plots_dr2(BPband_dr2, BPband_sig_dr2)
rev_BP_tuple = good_for_plots_rev(BPband_rev, BPband_sig_rev)

def plot_tuple_error(intuple, label):
    plt.errorbar(intuple[0], intuple[1], intuple[2], label = label)
    return

def plot_tuple(intuple, label, color, linestyle= 'default'):
    plt.plot(intuple[0], intuple[1], label=label, color= color, linestyle= linestyle)
    return

#plt.errorbar(dr2_G_tuple[0], dr2_G_tuple[1], dr2_G_tuple[2], label = "G DR2")
#plt.errorbar(rev_G_tuple[0], rev_G_tuple[1], rev_G_tuple[2], label = "G DR2 revised")
#plt.legend()
#plt.show()

plt.axhline(y=0, color ='k')
plot_tuple_error(dr2_G_tuple, "G DR2")
plot_tuple_error(rev_G_tuple, "G _DR2 revised")
plt.legend()
plt.show()

plt.axhline(y=0, color ='k')
plot_tuple_error(dr2_RP_tuple, "$G_{RP}$ DR2")
plot_tuple_error(rev_RP_tuple, "$G _{RP}$DR2 revised")
plt.legend()
plt.show()

plt.axhline(y=0, color ='k')
plot_tuple_error(dr2_BP_tuple, "$G_{BP}$ DR2")
plot_tuple_error(rev_BP_tuple, "$G _{BP}$ DR2 revised")
plt.legend()
plt.show()

######
plt.axhline(y=0, color ='k')

plot_tuple(dr2_G_tuple, "G DR2", color = 'g', linestyle = '--')
plot_tuple(rev_G_tuple, "G _DR2 revised", color = 'g')

plot_tuple(dr2_RP_tuple, "$G_{RP}$ DR2", color = 'r', linestyle = '--')
plot_tuple(rev_RP_tuple, "$G _{RP}$DR2 revised", color = 'r')

plot_tuple(dr2_BP_tuple, "$G_{BP}$ DR2", color = 'b', linestyle = '--')
plot_tuple(rev_BP_tuple, "$G _{BP}$ DR2 revised", color = 'b')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmissivity')
plt.legend()
plt.show()
