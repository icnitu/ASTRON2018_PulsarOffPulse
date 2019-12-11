import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import argparse
from collections import Iterable
import sys
import inspect


# ------- Get which txt files I need to read:

def flatten(items):
    """Yield items from any nested iterable; see Reference."""
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x


current_dir = os.path.dirname(os.path.realpath(__file__)) # can just use ./

input_folder = "0034_cgf" # this should be inputted when running later

try:
	input_folder_name = glob(input_folder)[0]
except IndexError:
	print "\nERROR: The input folder doesn't exist!\n"
	sys.exit()


# this file:
filename = inspect.getframeinfo(inspect.currentframe()).filename
# the directory this file is in:
current_folder = os.path.dirname(os.path.abspath(filename))+'/'

path = input_folder_name + '/'

# count the number of observations/folders in input_folder:
nobs = 0

for _, dirnames, _ in os.walk(path):
  # ^ this idiom means "we won't be using this value"
	if len(dirnames) > 0:
		list_obs = dirnames
	nobs += len(dirnames)


list_obs = np.sort(list_obs)
# list_obs is now the list of observations, ordered

analysis_file = "AnalysisResults.txt"
analysis_files = []
for iobs in range(nobs):
	analysis_files.append(glob(path+list_obs[iobs]+'/'+analysis_file))

# analysis_files is a list of the full path +name to the txt file "AnalysisResults.txt"

analysis_files = list(flatten(analysis_files))
#print analysis_files


# ----- Read from txt files:
# read lines: 
# Non-averaged
#9	scintill_bw = ... chans = ... MHz
#10	scintill_bw2 = ... chans = ... MHz
#11	scintill_time = ... chans = ... min
# Averaged 1/4
#25	scintill_bw = ... chans = ... MHz
#26	scintill_bw2 = ... chans = ... MHz
#27	scintill_time = ... chans = ... min
# Averaged 1/2
#41	scintill_bw = ... chans = ... MHz
#42	scintill_bw2 = ... chans = ... MHz
#43	scintill_time = ... chans = ... min

# Note: scintill_bw is used for eta;
# 	if one Gaussian fitted and not two, scintill_bw = scintill_bw2


list_sbw = []
list_sbw2 = []
list_st = []

list_sbw_avg4 = []
list_sbw2_avg4 = []
list_st_avg4 = []

list_sbw_avg2 = []
list_sbw2_avg2 = []
list_st_avg2 = []

for ifile in analysis_files:
	
	if not os.path.isfile(ifile):
		print("File path {} does not exist. Exiting...".format(ifile))
		sys.exit()

	with open(ifile) as fp:
		iline = 0
		for line in fp:
			if iline == 9:
				list_sbw.append(line)
			if iline == 10:
				list_sbw2.append(line)
			if iline == 11:
				list_st.append(line)

			if iline == 25:
				list_sbw_avg4.append(line)
			if iline == 26:
				list_sbw2_avg4.append(line)
			if iline == 27:
				list_st_avg4.append(line)

			if iline == 41:
				list_sbw_avg2.append(line)
			if iline == 42:
				list_sbw2_avg2.append(line)
			if iline == 43:
				list_st_avg2.append(line)
			
			iline += 1

list_obs_times = [] # just the "01Nov2013" etc
full_obs_times = [] # the "0034-0534-385.22Dec2010.med"

sbw = np.empty(nobs)
sbw2 = np.empty(nobs)
st = np.empty(nobs)

sbw_avg4 = np.empty(nobs)
sbw2_avg4 = np.empty(nobs)
st_avg4 = np.empty(nobs)

sbw_avg2 = np.empty(nobs)
sbw2_avg2 = np.empty(nobs)
st_avg2 = np.empty(nobs)

for i in range(nobs):
	list_obs_times.append(list_obs[i].split('.')[1])
	full_obs_times.append(list_obs[i].split('_')[0])

	dummy = list_sbw[i].split(' = ')[-1]
	sbw[i] = float (dummy.split(' ')[0])
	dummy = list_sbw2[i].split(' = ')[-1]
	sbw2[i] = float (dummy.split(' ')[0])
	dummy = list_st[i].split(' = ')[-1]
	st[i] = float (dummy.split(' ')[0])

	dummy = list_sbw_avg4[i].split(' = ')[-1]
	sbw_avg4[i] = float (dummy.split(' ')[0])
	dummy = list_sbw2_avg4[i].split(' = ')[-1]
	sbw2_avg4[i] = float (dummy.split(' ')[0])
	dummy = list_st_avg4[i].split(' = ')[-1]
	st_avg4[i] = float (dummy.split(' ')[0])

	dummy = list_sbw_avg2[i].split(' = ')[-1]
	sbw_avg2[i] = float (dummy.split(' ')[0])
	dummy = list_sbw2_avg2[i].split(' = ')[-1]
	sbw2_avg2[i] = float (dummy.split(' ')[0])
	dummy = list_st_avg2[i].split(' = ')[-1]
	st_avg2[i] = float (dummy.split(' ')[0])

all_s = np.empty((nobs,9))
all_s[:,0] = sbw; all_s[:,1] = sbw2; all_s[:,2] = st
all_s[:,3] = sbw_avg4; all_s[:,4] = sbw2_avg4; all_s[:,5] = st_avg4
all_s[:,6] = sbw_avg2; all_s[:,7] = sbw2_avg2; all_s[:,8] = st_avg2

headr = "BW (MHz); BW2 (MHz); T (min); BW  (avg 1/4)(MHz); BW2 (avg 1/4)(MHz); T (avg 1/4)(min); BW  (avg 1/2)(MHz); BW2 (avg 1/2)(MHz); T (avg 1/2)(min)" 

all_towrite = [[]]

for i in range(nobs):
	all_towrite[i].append(list_obs_times[i])
	all_towrite[i].append(all_s[i])
	all_towrite.append([])

#print all_towrite
obs_pulsar = full_obs_times[0].split('.')[0]
savefile = obs_pulsar + "_scintillParams.txt"
np.savetxt(savefile, all_s, fmt='%.10f', delimiter=' ', newline='\n', header=headr, comments='# ', encoding=None)

plt.hist(sbw,color='black',bins = 5,range = [0,4])
plt.show()
plt.hist(sbw_avg4,color='green',bins = 5,range = [0,4])
plt.show()
plt.hist(sbw_avg2,color='red',bins = 5,range = [0,4])
plt.show()




