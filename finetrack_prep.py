##!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
finetrack_prep.py
Jesse Wimert (11/2023)

This routine reads calibration data and TEP from an ATL03 file.  
TEP is used to create expected waveform tables.
An h5 file is created which contains the calibration data for 
first photon bias correction and expected 
waveform tables (created from the TEP).

INPUTS:
        -g1 ATL03 granule

OUTPUTS:
        atl03_calibrations_and_wf_tables.h5 


Steps:


Notes:

contents of atl03_calibrations_and_wf_tables.h5:

-atl03 filename						Filename of ATL03 processed
-orientation						Spacecraft orintation
-CAL19/dead_time (:)				CAL19 deadtime look-up table
-CAL19/fpb_corr (:, : , :)			CAL19 correction look-up table
-CAL19/strength (:, :)				CAL19 strength look-up table
-CAL19/width (:, :)					CAL19 width look-up table
-exp_wf_table/binz (:)				expected wf table bins values
-exp_wf_table/mz (:)				expected wf table mean values
-exp_wf_table/sdz (:)				expected wf table sigma values
-exp_wf_table/wf_table1 (:, :, :)   expected wf table using TEP1
-exp_wf_table/wf_table2 (:. :. :)   expected wf table using TEP2
-gtx/dead_time						dead time for beam gtx
-gtx/strong_weak					strong (1) or weak (0) beam label
-gtx/wf_table_select				process wf_table1 (1) or wf_table2 (2)



python finetrack_prep.py -g1 ATL03_20220726163210_05311604_006_02.h5

"""
#
import numpy as np
import h5py
#
import argparse
import math
import os

##
## Function to parse command line arguments
##
 
def parse():
    """ Parse command line input arguments. """

    parser = argparse.ArgumentParser(description='Compute ATL07 sea ice heights')
    parser.add_argument('-g1', action='store', default=None, dest='granule_1',
                        help='Path to granule #1',
                        type=str)
    inps = parser.parse_args()
    return inps

##
## Print welcome message
##

os.system('clear')
print(' ')
print('==============================================')
print(' ')
print('finetrack_prep.py')
print(' ')
print(' ')
print('This routine reads calibration data and TEP from an ATL03 file.')
print('TEP is used to create expected waveform tables.')
print('An h5 file is created which contains the calibration data for ')
print('first photon bias correction and expected waveform tables')
print(' ')
print(' ')
print('Jesse Wimert')
print('Last update: November 2023')
print(' ')
print('==============================================')
print(' ')
print(' ')


# Read command line input parameters
inps = parse()

# Check if input parameters are provided
if inps.granule_1 is None:
    print(' ')
    print('Missing path for granule #1. Use option -G1')
    sys.exit(1)

# Assign input parameters to variables
granule_1 = inps.granule_1


##
## Define speed of light
##


constant_c = 299792458.0

##
## Define beam ground track labels
##


beams = "gt1l gt1r gt2l gt2r gt3l gt3r"
beams = beams.split()

##
## Set waveform table constants
##


lb_bin=-3.5
ub_bin=3.5
bin_delta=0.025

lb_mu=-0.5
ub_mu=0.5
mu_delta=0.01

lb_sig=0.0
ub_sig=1.5
sig_delta=0.01

n_bin = math.ceil((ub_bin - lb_bin)/bin_delta) + 1
n_mu = math.ceil((ub_mu - lb_mu) / mu_delta) + 1
n_sig = math.ceil((ub_sig - lb_sig) / sig_delta) + 1

##
## Set waveform table axes
##


binz = np.zeros(n_bin)
mz = np.zeros(n_mu)
sdz = np.zeros(n_sig)

for ii in np.arange(0,n_bin):
    binz[ii] = lb_bin + ii*bin_delta

for ii in np.arange(0,n_mu):
    mz[ii] = lb_mu + ii*mu_delta

for ii in np.arange(0,n_sig):
    sdz[ii] = lb_sig + ii*sig_delta
sdz[0] = 0.001

##
## Define gaussian distribution
##


def gaussian(x, mu, sig):
    return (
        np.exp(-((x - mu)**2)/(2.0*(sig**2))) / (np.sqrt(2.0 * np.pi) * sig)
    )

###
### Start Routine
###

###
### Open ATL03, save filename
###


hf = h5py.File(granule_1, 'r')
atl03_filename = os.path.basename(granule_1)

#
# Read orientation, set strong and weak beam labels
#


orientation = hf.get('/orbit_info/sc_orient')[0]

if orientation == 0:
    strong_beams = "gt1l gt2l gt3l"
    weak_beams = "gt1r gt2r gt3r"
elif orientation == 1:
    strong_beams = "gt1l gt2l gt3l"
    weak_beams = "gt1r gt2r gt3r"

strong_beams = strong_beams.split()
weak_beams = weak_beams.split()

strong_weak = {}
for beam in beams:
    if beam in strong_beams:
        strong_weak[beam] = 1
    else:
        strong_weak[beam] = 0

###
### Read calibration data
###

#
# Read dead_time tables from ATL03 (seconds)
#
# Note, dead_time values are repeated for each beam of a beam pair 
# (gt1l and gt1r have same dead_time values)
# The 20 entries in the dead_time table correspond to the 16 channels for
# the strong beam and 4 channels for the weak beam.
# To compute average deadtime for each beam, average dead time over the
# activate channels for each ground track
#


deadtime_read = {}

for beam in beams:
    deadtime_read[beam] = hf['ancillary_data/calibrations/dead_time/' + beam + '/dead_time'][:]

#
# Compute dead_time for each beam
# Average and convert to nanoseconds
#


dead_time = {}

for beam in strong_beams:
    dead_time[beam] = (np.sum(deadtime_read[beam][0:16])/16.0)* 1.0e9

for beam in weak_beams:
    dead_time[beam] = (np.sum(deadtime_read[beam][16:20])/4.0)* 1.0e9

#
# Read CAL19 tables from ATL03 
#
# Note, the CAL19 tables are identical over all six beams, only need to
# read one
# 


cal19_dead_time = hf['ancillary_data/calibrations/first_photon_bias/gt1l/dead_time'][:]
cal19_strength = hf['ancillary_data/calibrations/first_photon_bias/gt1l/strength'][:][:]
cal19_width = hf['ancillary_data/calibrations/first_photon_bias/gt1l/width'][:][:]
cal19_corr = hf['ancillary_data/calibrations/first_photon_bias/gt1l/ffb_corr'][:][:][:]

###
### Read TEP data
###

#
# Read tep_range_prim to get bounds of primary pulse
#


tep_range_prim = hf.get('/ancillary_data/tep/tep_range_prim')[:]

#
# Read tep_valid_spot to determine which TEP to use for each beam
#


tep_valid_spot = hf.get('/ancillary_data/tep/tep_valid_spot')[:]

wf_table_select = {}
ii=0
for beam in beams:
    wf_table_select[beam] = tep_valid_spot[ii]
    ii = ii + 1

##
## Read tep_hist and tep_hist_time for both pce1 and pce2
##


tep_pce1_hist = hf['atlas_impulse_response/pce1_spot1/tep_histogram/tep_hist'][:]
tep_pce1_time = hf['atlas_impulse_response/pce1_spot1/tep_histogram/tep_hist_time'][:]

tep_pce2_hist = hf['atlas_impulse_response/pce2_spot3/tep_histogram/tep_hist'][:]
tep_pce2_time = hf['atlas_impulse_response/pce2_spot3/tep_histogram/tep_hist_time'][:]

#
# Close h5 file
#


hf.close()

###
### Adjust TEP histograms so that it may be convolved with gaussian distributions
###

#
# Isolate primary pulse using tep_range_prim (same for both TEPS)
#


i0 = 0
for i in np.arange(0, len(tep_pce1_hist)):
    if (tep_pce1_time[i] > tep_range_prim[0]):
        continue
    i0 = i

i1 = i0
for i in np.arange(i0, len(tep_pce1_hist)):
    if (tep_pce1_time[i] >= tep_range_prim[1]):
        continue
    i1 = i

tep_hist1 = tep_pce1_hist[i0:i1 + 1]
tep_time1 = tep_pce1_time[i0:i1 + 1]

tep_hist2 = tep_pce2_hist[i0:i1 + 1]
tep_time2 = tep_pce2_time[i0:i1 + 1]

#
# Find mean of histograms
#


hist_mean1 = np.sum(tep_hist1*tep_time1)/np.sum(tep_hist1)
hist_mean2 = np.sum(tep_hist2*tep_time2)/np.sum(tep_hist2)

#
# Find tep_hist_time closest to mean
#


i_mean1 = 0
min_diff = abs(tep_time1[0] - hist_mean1)
for i in np.arange(0,len(tep_time1)-1):
    diff = abs(tep_time1[i]-hist_mean1)
    if (diff < min_diff):
        min_diff = diff
        i_mean1 = i

i_mean2 = 0
min_diff = abs(tep_time2[0] - hist_mean2)
for i in np.arange(0,len(tep_time2)-1):
    diff = abs(tep_time2[i]-hist_mean2)
    if (diff < min_diff):
        min_diff = diff
        i_mean2 = i

#
# Center the mean and normalize
#


mean_x1 = tep_time1[i_mean1]
tep_time1_center = tep_time1 - mean_x1
tep_hist1_normal = tep_hist1/np.sum(tep_hist1)

mean_x2 = tep_time2[i_mean2]
tep_time2_center = tep_time2 - mean_x2
tep_hist2_normal = tep_hist2/np.sum(tep_hist2)

#
# Convert from time to meters
#


tep_meter1 = tep_time1_center * (constant_c/2.0)
tep_meter2 = tep_time2_center * (constant_c/2.0)

#
# Squish TEP
#


tep_hist1_squish = np.zeros(len(tep_meter1))
tep_hist1_final = np.zeros(len(tep_meter1))

tep_hist2_squish = np.zeros(len(tep_meter2))
tep_hist2_final = np.zeros(len(tep_meter2))

for i in np.arange(0, i_mean1/2 - 1):
    tep_hist1_squish[int(i_mean1+i)] = tep_hist1_normal[int(i_mean1 + i*2)]
    tep_hist1_squish[int(i_mean1-i)] = tep_hist1_normal[int(i_mean1 - i*2)]
tep_hist1_squish[0:int(i_mean1 - i_mean1/2 -2)] = 0.0
tep_hist1_squish[int(i_mean1 + i_mean1/2 -1):len(tep_time1_center)-1] = 0.0

for i in np.arange(0, i_mean2/2 - 1):
    tep_hist2_squish[int(i_mean2+i)] = tep_hist2_normal[int(i_mean2 + i*2)]
    tep_hist2_squish[int(i_mean2-i)] = tep_hist2_normal[int(i_mean2 - i*2)]
tep_hist2_squish[0:int(i_mean2 - i_mean2/2 -2)] = 0.0
tep_hist2_squish[int(i_mean2 + i_mean2/2 -1):len(tep_time2_center)-1] = 0.0

tep_hist1_squish = np.where(tep_hist1_squish > 0.0, tep_hist1_squish, 0.0)
tep_hist2_squish = np.where(tep_hist2_squish > 0.0, tep_hist2_squish, 0.0)

#
# shift 2 bins to the left
# (IS THIS CORRECT?)
#


tep_hist1_final[2:len(tep_time1_center)-1] = tep_hist1_squish[0:len(tep_time1_center)-3]
tep_hist2_final[2:len(tep_time2_center)-1] = tep_hist2_squish[0:len(tep_time2_center)-3]

#
# normalize
#


tep_hist1_final_normal = tep_hist1_final/np.sum(tep_hist1_final)
tep_hist2_final_normal = tep_hist2_final/np.sum(tep_hist2_final)

#
# Reverse order of histogram and flip signs to switch from range space to height space
#

tx_hts1 = np.zeros(len(tep_meter1))
tx_pulse1 = np.zeros(len(tep_hist1_final))

tx_hts2 = np.zeros(len(tep_meter2))
tx_pulse2 = np.zeros(len(tep_hist2_final))

for ii in np.arange(0,len(tep_hist1_final)):
    jj = (len(tep_hist1_final)-1 ) - ii
    tx_hts1[jj] = tep_meter1[ii]
    tx_pulse1[jj] = tep_hist1_final[ii]

for ii in np.arange(0,len(tep_hist2_final)):
    jj = (len(tep_hist2_final)-1 ) - ii
    tx_hts2[jj] = tep_meter2[ii]
    tx_pulse2[jj] = tep_hist2_final[ii]

tx_hts1 = -1.0 * tx_hts1
tx_hts2 = -1.0 * tx_hts2

#
# Compute centroid and shift
#

hist_ctr1 = np.sum(tx_pulse1 * tx_hts1)/np.sum(tx_hts1)
hist_ctr2 = np.sum(tx_pulse2 * tx_hts2)/np.sum(tx_hts2)


tx_hts1 = tx_hts1 - hist_ctr1
tx_hts2 = tx_hts2 - hist_ctr2

##
## Interpolate to binz values - TEP1
##

#
# subset interpolation arguments
#


lb=0
rb=n_bin

for ii in np.arange(0,n_bin):
    if (binz[ii] >= tx_hts1[0]):
        continue
    lb = ii

for ii in reversed(np.arange(0,n_bin)):
    if (binz[ii] <= tx_hts1[len(tep_hist1_final)-1]):
        continue
    rb = ii

#
# interpolate
#


tep_binz1 = np.zeros(n_bin)
tep_binz1[lb:rb] = np.interp(binz[lb:rb], tx_hts1, tx_pulse1)

##
## Interpolate to binz values - TEP2
##

#
# subset interpolation arguments
#


lb=0
rb=n_bin

for ii in np.arange(0,n_bin):
    if (binz[ii] >= tx_hts2[0]):
        continue
    lb = ii

for ii in reversed(np.arange(0,n_bin)):
    if (binz[ii] <= tx_hts2[len(tep_hist2_final)-1]):
        continue
    rb = ii

#
# interpolate
#


tep_binz2 = np.zeros(n_bin)
tep_binz2[lb:rb] = np.interp(binz[lb:rb], tx_hts2, tx_pulse2)

###
### Create expected waveform tables
###


wf_table1 = np.zeros([n_mu, n_sig, n_bin])
wf_table2 = np.zeros([n_mu, n_sig, n_bin])
temp_waveform = np.zeros(n_bin)
conv_waveform1 = np.zeros(n_bin)
conv_waveform2 = np.zeros(n_bin)

#
# fill bins of both expected waveform tables
#


for i_sig in np.arange(0, n_sig):
    for i_mu in np.arange(0, n_mu):
    #
    # compute gaussian distribution with values of mu and sd
    #
        temp_waveform = gaussian(binz, mz[i_mu], sdz[i_sig])
    #
    # compute convolved waveform of gaussian and TEP
    #
        conv_waveform1 = np.convolve(temp_waveform, tep_binz1, 'same')
        conv_waveform2 = np.convolve(temp_waveform, tep_binz2, 'same')
    #
    # load convolved waveform to expected waveform table
    #
        wf_table1[i_mu, i_sig, :] = conv_waveform1 / np.sum(conv_waveform1)
        wf_table2[i_mu, i_sig, :] = conv_waveform2 / np.sum(conv_waveform2)

###
### write to h5 file
###

hf = h5py.File('atl03_calibrations_and_wf_tables.h5', 'w')
print('1')
dt = h5py.special_dtype(vlen=str)
print('2')
# hf.create_dataset('atl03_filename', data=atl03_filename, dtype = dt)
print('3')
hf.create_dataset('atl03_filename', (1,), data=atl03_filename, dtype = dt)
hf.create_dataset('orientaion', (1,), data=orientation)


g1 = hf.create_group('CAL19')
g1.create_dataset('dead_time', data=cal19_dead_time)
g1.create_dataset('strength', data=cal19_strength)
g1.create_dataset('width', data=cal19_width)
g1.create_dataset('fpb_corr', data=cal19_corr)

g2 = hf.create_group('exp_wf_table')
g2.create_dataset('binz', data=binz)
g2.create_dataset('mz', data=mz)
g2.create_dataset('sdz', data=sdz)
g2.create_dataset('wf_table1', data=wf_table1)
g2.create_dataset('wf_table2', data=wf_table2)
# g2.create_dataset('tep_hist1', data=tep_hist1)
# g2.create_dataset('tep_hist1_final', data=tep_hist1_final)
# g2.create_dataset('tx_hts1', data=tx_hts1)
# g2.create_dataset('tx_pulse1', data=tx_pulse1)
# g2.create_dataset('tep_binz1', data=tep_binz1)



ii=3
for beam in beams:
    temp_str = 'g' + str(ii)
    ii = ii + 1
    temp_str = hf.create_group(beam)
    temp_str.create_dataset('dead_time', (1,),  data=dead_time[beam])
    temp_str.create_dataset('strong_weak', (1,),  data=strong_weak[beam])
    temp_str.create_dataset('wf_table_select', (1,),  data=wf_table_select[beam])

#
# Close h5 file
#


hf.close()

###
### Print summary log
###


print(' ')
print('==============================================')
print(' ')
print('atl03_calibrations_and_wf_tables.h5 Summary:')
print(' ')
print('==============================================')
print(' ')
print(' ')
print('ATL03 file processed:')
print(atl03_filename)
print(' ')
print('Spacecraft Orientation')
print(orientation)
print(' ')
print(' ')
print('Strong Beam ground tracks')
print(' ')
print(strong_beams)
print(' ')
print('Weak Beam ground tracks')
print(' ')
print(weak_beams)
print(' ')
print(' ')
print('Beam gt, strong/weak, wf_table_select, dead_time')
print(' ')
for beam in beams:
	print(beam, strong_weak[beam], wf_table_select[beam], dead_time[beam])
print(' ')
print(' ')
print('Size of CAL19 tables')
print(' ')
print('dead_time',cal19_dead_time.shape)
print('strength', cal19_strength.shape)
print('width',cal19_width.shape)
print('fppb_corr', cal19_corr.shape)
print(' ')
print('Size of wf_table1')
print(wf_table1.shape)
print('Size of wf_table2')
print(wf_table2.shape)
print(' ')
print('n sigma bins')
print(n_sig)
print(' ')
print('n mu bins')
print(n_mu)
print(' ')
print('n hist bins')
print(n_bin)
print(' ')
print('==============================================')
print(' ')
print(' ')

