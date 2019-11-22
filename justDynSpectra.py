import psrchive as pc
import numpy as np
import matplotlib.pyplot as plt

from readInputs import *
from plotFrequencyVsPhase import *
from plotPulseProfile import *
from Bandpass import *
from DynSpectrum_corrections import *
from plotDynSpectrum import *
from CorrelationsFULL import *
from printResults import *
from averaging import *




# read file with archive
archive_file_message = "\nInput the file name (including path) containing the archive you would like to use.\n"
archive_file = read_good_file(archive_file_message)
# read folder where to save plots & text files
saved_folder_message = "\nInput the folder name (including path) where you want to save the plots/text files generated.\n"
saved_folder = read_good_folder(saved_folder_message)


print "\n -----> Getting the archive...\n"

archive = pc.Archive_load(archive_file)
print "\n -----> Dedispersing the archive...\n"
archive.dedisperse()
print "\n -----> Pscrunching the archive...\n"
archive.pscrunch()
Nsubint = archive.get_nsubint()
Nchan = archive.get_nchan()	
Nbin = archive.get_nbin()
			
bandwidth = (archive.get_bandwidth())
centreFreq = archive.get_centre_frequency()
freqMin = centreFreq-bandwidth/2.0
freqMax = centreFreq+bandwidth/2.0
chanToMHz = bandwidth/float(Nchan)
		
obstime = archive.integration_length()
binToTime = (obstime/60.)/float(Nsubint) # to minutes

	
archive_copy = archive.clone()
archive_copy.remove_baseline()
data3_initial = archive_copy.get_data()
data3_initial = data3_initial[:,0,:,:] # (Nsubint,Nchan,Nbin)
data_freqphase = data3_initial.mean(0) # (Nchan,Nbin)
data_profile = data_freqphase.mean(0) # (Nbin)
	
			
print "\n -----> Plotting Frequency Vs Phase...\n"
# save data_freqphase to file as numpy array
file_freqphase_name = saved_folder + 'data_FreqPhase'
np.save(file_freqphase_name,data_freqphase)
plot_FreqVsPhase(data_freqphase,freqMin,freqMax,saved_folder)
				
print "\n -----> Plotting Pulse Profile...\n"
# save data_profile to file as numpy array
file_profile_name = saved_folder + 'data_PulseProfile'
np.save(file_profile_name,data_profile)
plot_PulseProfile(data_profile,saved_folder)


	
binOffset_message = "\nInput the bin number at which you want to place the peak of the profile, as an integer in the range [" + str(0) + ", " + str(Nbin) + "), or -1 if you don't want to change it. You can see the current Pulse Profile in " + str(saved_folder) + "PulseProfile.png\n"

binOffset = read_good_int(binOffset_message,-1.,Nbin-1)
phaseOffset = binOffset/float(Nbin)


if binOffset >= 0:
	# offset:
	print "\n -----> Offsetting peak of profile in phase to bin",binOffset,"...\n"
	archive_copy.centre_max_bin(phase_offset=phaseOffset)
	archive.centre_max_bin(phase_offset=phaseOffset)
	data3_initial = archive_copy.get_data()
	data3_initial = data3_initial[:,0,:,:] # (Nsubint,Nchan,Nbin)
	data_freqphase = data3_initial.mean(0) # (Nchan,Nbin)
	data_profile = data_freqphase.mean(0) # (Nbin)

	# redo freqphase and profile:	
			
	print "\n -----> Plotting (offset) Frequency Vs Phase...\n"
	# save data_freqphase to file as numpy array
	file_freqphase_name = saved_folder + 'data_FreqPhase'
	np.save(file_freqphase_name,data_freqphase)
	plot_FreqVsPhase(data_freqphase,freqMin,freqMax,saved_folder)
		
		
	print "\n -----> Plotting (offset) Pulse Profile...\n"
	# save data_profile to file as numpy array
	file_profile_name = saved_folder + 'data_PulseProfile'
	np.save(file_profile_name,data_profile)
	plot_PulseProfile(data_profile,saved_folder)

elif binOffset == -1:
	print "\n -----> No offsetting phase, moving on...\n"	

else:
	print "\n -----> Shouldn't get here, something wrong!\n"	

	

ranges_message = "\n -- You will next need to input the on- and off-pulse ranges. You can see the Pulse Profile in " + str(saved_folder) + "PulseProfile.png\n"

print ranges_message

offBinRange_message = "\nInput the off-pulse phase range, as integer bin values of the form 'min_bin,max_bin' or 'max_bin,min_bin', within the range [" + str(0) + ", " + str(Nbin) + "). These will be used for all the analysis to follow.\n"
		
offBin_min,offBin_max = read_good_binRange(offBinRange_message,0,Nbin-1)
	
	
onBinRange_message = "\nInput the on-pulse phase range, as integer bin values of the form 'min_bin,max_bin' or 'max_bin,min_bin', withing the range [" + str(0) + ", " + str(Nbin) + "). These will be used for all the analysis to follow.\n"
		
onBin_min,onBin_max = read_good_binRange(onBinRange_message,0,Nbin-1)


print "\n -----> Removing baseline...\n"
# Base removal:
archive.remove_baseline()
data3_gain_base = archive.get_data()
data3_gain_base = data3_gain_base[:,0,:,:] # (Nsubint,Nchan,Nbin)


Noff = offBin_max - offBin_min
Non = onBin_max - onBin_min 	
				
data3_copy = data3_gain_base.copy()
				
data_on = data3_copy[:,:,onBin_min:onBin_max]
data_off = data3_copy[:,:,offBin_min:offBin_max]
			
			
data_on = data_on.mean(2) # (Nsubint,Nchan)
data_off = data_off.mean(2)
				

rms_off = get_rms_off(data_off)		

rms_off[rms_off==0] = np.nan
data_on_corr = np.divide(data_on,rms_off)
data_on_corr = data_on_corr*np.sqrt(Non/float(Noff))
data_off_corr = np.divide(data_off,rms_off)


# save data_on(off)_corr to file as DS, as numpy arrays	
file_ds_on_name = saved_folder + 'data_DS_on'
np.save(file_ds_on_name,data_on_corr)

file_ds_off_name = saved_folder + 'data_DS_off'
np.save(file_ds_off_name,data_off_corr)


# non-averaged:

print '\n\n------ Non-averaged ------\n'
				
print "\n -----> Plotting Dynamic Spectrum on-pulse...\n"
plot_DynSpectrum(data_on_corr, freqMin, freqMax, binToTime, 'DynSpectrum_on.png', saved_folder)
print "\n -----> Plotting Dynamic Spectrum off-pulse...\n"
plot_DynSpectrum(data_off_corr, freqMin, freqMax,  binToTime, 'DynSpectrum_off.png', saved_folder)
	
	




