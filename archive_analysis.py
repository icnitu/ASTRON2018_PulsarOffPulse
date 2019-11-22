import psrchive as pc
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import sys

#import readInputs
import plotFrequencyVsPhase
import plotPulseProfile
import Bandpass_multiple
import DynSpectrum_corrections
import plotDynSpectrum
import CorrelationsFULL
import printResults
import averaging



def check_output_folder(input_value):
	
	if os.path.isdir(input_value)==False:
	
		try:
			print "\n--- The folder doesn't exist, creating the folder...\n"
			os.mkdir(input_value)
		
		except OSError:
			print "\n--- The folder path doesn't exist, try again!\n"
			sys.exit()
		
						
	input_value = str(input_value)
	if input_value[-1] != '/':
		input_value = input_value + '/'
		
	return input_value


def archive_analysis(archive_file, output_folder, fphaseOffset = -1, ffreq_fract=0.03, fNbands_number = 1,fdoOnPulse = False, fdoOffPulse = False, fonBin_f=np.empty(2),foffBin_f=np.empty(2),small_fit = "n", small_fit_t = "n", ffinishEarly = False):
	

	output_folder = check_output_folder(output_folder)


	print "\n -----> Getting the archive...\n"
	#print (archive_file)
	archive = pc.Archive_load(archive_file)
	print "\n -----> Dedispersing the archive...\n"
	archive.dedisperse()
	print "\n -----> Pscrunching the archive...\n"
	archive.pscrunch()
	Nsubint = archive.get_nsubint()
	Nchan = archive.get_nchan()	
	Nbin = archive.get_nbin()		

	offBin_min=int(round(foffBin_f[0]*Nbin))
	offBin_max=int(min(Nbin-1,round(foffBin_f[1]*Nbin)))
	
	onBin_min=int(round(fonBin_f[0]*Nbin))
	onBin_max=int(min(Nbin-1,round(fonBin_f[1]*Nbin)))



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
	file_freqphase_name = output_folder + 'data_FreqPhase'
	np.save(file_freqphase_name,data_freqphase)
	plotFrequencyVsPhase.plot_FreqVsPhase(data_freqphase,freqMin,freqMax,output_folder)
				
	print "\n -----> Plotting Pulse Profile...\n"
	# save data_profile to file as numpy array
	file_profile_name = output_folder + 'data_PulseProfile'
	np.save(file_profile_name,data_profile)
	plotPulseProfile.plot_PulseProfile(data_profile,output_folder,onBin_fmin = fonBin_f[0],onBin_fmax = fonBin_f[1], offBin_fmin = foffBin_f[0], offBin_fmax = foffBin_f[1])





	if fphaseOffset:
		# offset:
		print "\n -----> Offsetting peak of profile to phase",fphaseOffset,"...\n"
		archive_copy.centre_max_bin(phase_offset=fphaseOffset)
		archive.centre_max_bin(phase_offset=fphaseOffset)
		data3_initial = archive_copy.get_data()
		data3_initial = data3_initial[:,0,:,:] # (Nsubint,Nchan,Nbin)
		data_freqphase = data3_initial.mean(0) # (Nchan,Nbin)
		data_profile = data_freqphase.mean(0) # (Nbin)

		# redo freqphase and profile:	
			
		print "\n -----> Plotting (offset) Frequency Vs Phase...\n"
		# save data_freqphase to file as numpy array
		file_freqphase_name = output_folder + 'data_FreqPhase'
		np.save(file_freqphase_name,data_freqphase)
		plotFrequencyVsPhase.plot_FreqVsPhase(data_freqphase,freqMin,freqMax,output_folder)
		
		
		print "\n -----> Plotting (offset) Pulse Profile...\n"
		# save data_profile to file as numpy array
		file_profile_name = output_folder + 'data_PulseProfile'
		np.save(file_profile_name,data_profile)
		plotPulseProfile.plot_PulseProfile(data_profile,output_folder,onBin_fmin = fonBin_f[0],onBin_fmax = fonBin_f[1], offBin_fmin = foffBin_f[0], offBin_fmax = foffBin_f[1])

	else:
		print "\n -----> No offsetting phase, moving on...\n"


	if fdoOffPulse:

		# TODO: error handling for bandpass fit
		print "\n -----> Getting Bandpass data and fit...\n"
		weights = archive.get_weights()		
		data4 = archive.get_data()
		data3_copy = data4.copy()
		data3_copy = data3_copy[:,0,:,:] # (Nsubint,Nchan,Nbin)
		
		
		# get bandpass data and fit
		data_off_bandpass,poly_vals_fits = Bandpass_multiple.get_multiple_bandpass_with_fit(data3_copy,weights,offBin_min,offBin_max,Nbands=fNbands_number,poly_order=10)
		# to modify poly_order etc - here!
			

		print "\n -----> Plotting Bandpass with fit...\n"
		# save data_off_bandpass to file as numpy array
		file_bandpass_name = output_folder + 'data_Bandpass'
		np.save(file_bandpass_name,data_off_bandpass)
		Bandpass_multiple.plot_multiple_bandpass_with_fit(data_off_bandpass,poly_vals_fits,freqMin,freqMax,output_folder)
			
		
		if fdoOnPulse:	
			# 
			print "\n -----> Removing baseline...\n"
			# Base removal:
			archive.remove_baseline()
			data3_gain_base = archive.get_data()
			data3_gain_base = data3_gain_base[:,0,:,:] # (Nsubint,Nchan,Nbin)

			print "\n -----> Applying Gain corrections...\n"
			# Gain correction:
			xchan = np.arange(Nchan)
			# have poly_vals_fits
	
			normalGain = poly_vals_fits/poly_vals_fits.mean()
			#print "normalGain =",normalGain				
			for ichan in range(Nchan):
				data3_gain_base[:,ichan,:] = data3_gain_base[:,ichan,:]/normalGain[ichan]		

		
			Noff = offBin_max - offBin_min
			Non = onBin_max - onBin_min 	
				
			data3_copy = data3_gain_base.copy()
				
				
			data_on = data3_copy[:,:,onBin_min:onBin_max]
			data_off = data3_copy[:,:,offBin_min:offBin_max]
			
			
			data_on = data_on.mean(2) # (Nsubint,Nchan)
			data_off = data_off.mean(2)
				
			# TODO: maybe add subint_number as a param?
			subint_number = 1
			data_on_corr,data_off_corr = DynSpectrum_corrections.correct_data_on_off(data_on,data_off,weights,Non,Noff,fract = ffreq_fract,Nsubint_last=subint_number)
			# save data_on(off)_corr to file as DS, as numpy arrays	
			file_ds_on_name = output_folder + 'data_DS_on'
			np.save(file_ds_on_name,data_on_corr)

			file_ds_off_name = output_folder + 'data_DS_off'
			np.save(file_ds_off_name,data_off_corr)

			# TODO: small Gaussian stuff

#			small_message = "\nWould you like to fit Gaussians to just the central part of ACon(f)? Answer with 'y' or 'n'.\n"

#			small_fit = read_good_input(small_message,["y","n"])
			#small_fit = "n"

			if str(small_fit)=="y":
				small_fit = True
			elif str(small_fit)=="n":
				small_fit = False

			small_fit_avg = small_fit

#			small_message_t = "\nWhat about for the ACon(t)? Answer with 'y' or 'n'.\n"

#			small_fit_t = read_good_input(small_message_t,["y","n"])
			
			#small_fit_t = "n"
			if str(small_fit_t)=="y":
				small_fit_t = True
			elif str(small_fit_t)=="n":
				small_fit_t = False

			small_fit_t_avg = small_fit_t

			

			# non-averaged:

			print '\n\n------ Non-averaged ------\n'
				
			print "\n -----> Plotting Dynamic Spectrum on-pulse...\n"
			plotDynSpectrum.plot_DynSpectrum(data_on_corr, freqMin, freqMax, binToTime, 'DynSpectrum_on.png', output_folder)
			print "\n -----> Plotting Dynamic Spectrum off-pulse...\n"
			plotDynSpectrum.plot_DynSpectrum(data_off_corr, freqMin, freqMax,  binToTime, 'DynSpectrum_off.png', output_folder)
	
	
				
			# CORRELATIONS:

			# open file:
			writeFile_name = output_folder + 'AnalysisResults.txt'
			writeFile = open(writeFile_name,"w")
			
			writeFile.write('\nArchive file: %s' %archive_file)
			writeFile.write('\n\nOn-pulse fractional range: [%4.3f, %4.3f]' % (fonBin_f[0], fonBin_f[1]) )
			writeFile.write('\nOff-pulse fractional range: [%4.3f, %4.3f]' % (foffBin_f[0], foffBin_f[1]) )
			writeFile.write('\n\n')

			I_on = data_on_corr.copy()
			I_off = data_off_corr.copy()
			
			if ffinishEarly:
				results_values_small = printResults.print_results_small(I_on, bandwidth, obstime, chanToMHz,binToTime,output_folder,writeFile,small0=small_fit,small0t=small_fit_t)
				scintill_time = results_values_small[2]
				scintill_bw = results_values_small[0]
				Neff=results_values_small[4]
			else:
				results_values = printResults.print_results(I_on,I_off, bandwidth, obstime, chanToMHz,binToTime,output_folder,writeFile,small0=small_fit,small0t=small_fit_t)
	
	
				scintill_time = results_values[8]
				scintill_bw = results_values[6] # first Gaussian fitted in case of 2!
				Neff=results_values[11]
				#print '\ntest: ', scintill_time,scintill_bw, '\n'

			# ------ Averaged by 1/4:

			size_t = averaging.get_size_avg(scintill_time,fraction=1/4.)
			size_f = averaging.get_size_avg(scintill_bw,fraction=1/4.)

			data_on_corr_avg = averaging.average(data_on_corr,size_t,size_f)
			data_off_corr_avg = averaging.average(data_off_corr,size_t,size_f)

			# save data_on(off)_corr to file as DS, as numpy arrays	
			file_ds_on_avg_name = output_folder + 'data_DS_on_avg4'
			np.save(file_ds_on_avg_name,data_on_corr_avg)

			file_ds_off_avg_name = output_folder + 'data_DS_off_avg4'
			np.save(file_ds_off_avg_name,data_off_corr_avg)
	
			# new values for Nchan etc in case of averaging:
			Nchan_avg4 = data_on_corr_avg.shape[1]
			Nsubint_avg4 = data_on_corr_avg.shape[0]
			chanToMHz_avg4 = bandwidth/float(Nchan_avg4)
			binToTime_avg4 = (obstime/60.)/float(Nsubint_avg4) # to minutes

			scintill_bw_forAvg4 = scintill_bw*Nchan_avg4/float(Nchan)
			scintill_time_forAvg4 = scintill_time*Nsubint_avg4/float(Nsubint)

	
			# 		plot avg..
			print '\n\n------ Averaged by 1/4*scintillation parameters ------\n'

			print "\n -----> Plotting Averaged (by 1/4) Dynamic Spectrum on-pulse...\n"
			plotDynSpectrum.plot_DynSpectrum(data_on_corr_avg, freqMin, freqMax, binToTime_avg4, 'DynSpectrum_on_avg4.png', output_folder)
			print "\n -----> Plotting Averaged (by 1/4) Dynamic Spectrum off-pulse...\n"
			plotDynSpectrum.plot_DynSpectrum(data_off_corr_avg, freqMin, freqMax,  binToTime_avg4, 'DynSpectrum_off_avg4.png', output_folder)	
		
			
			I_on_avg = data_on_corr_avg.copy()
			I_off_avg = data_off_corr_avg.copy()
	
			if ffinishEarly:
				results_values_avg_small = printResults.print_results_small(I_on_avg, bandwidth, obstime, chanToMHz_avg4, binToTime_avg4, output_folder, writeFile, is_avg=True, scintill_bw0=scintill_bw_forAvg4, scintill_time0=scintill_time_forAvg4, small0=small_fit_avg, small0t=small_fit_t_avg, size_t_avg = size_t, size_f_avg=size_f)
			else:
				results_values_avg_small = printResults.print_results(I_on_avg,I_off_avg, bandwidth, obstime, chanToMHz_avg4, binToTime_avg4, output_folder, writeFile, is_avg=True, scintill_bw0=scintill_bw_forAvg4, scintill_time0=scintill_time_forAvg4, small0=small_fit_avg, small0t=small_fit_t_avg, size_t_avg = size_t, size_f_avg=size_f)
	

			# ------ Averaged by 1/2:

			size_t2 = averaging.get_size_avg(scintill_time,fraction=1/2.)
			size_f2 = averaging.get_size_avg(scintill_bw,fraction=1/2.)

			data_on_corr_avg2 = averaging.average(data_on_corr,size_t2,size_f2)
			data_off_corr_avg2= averaging.average(data_off_corr,size_t2,size_f2)

			# save data_on(off)_corr to file as DS, as numpy arrays	
			file_ds_on_avg2_name = output_folder + 'data_DS_on_avg2'
			np.save(file_ds_on_avg2_name,data_on_corr_avg2)

			file_ds_off_avg2_name = output_folder + 'data_DS_off_avg2'
			np.save(file_ds_off_avg2_name,data_off_corr_avg2)
	
			# new values for Nchan etc in case of averaging:
			Nchan_avg2 = data_on_corr_avg2.shape[1]
			Nsubint_avg2 = data_on_corr_avg2.shape[0]
			chanToMHz_avg2 = bandwidth/float(Nchan_avg2)
			binToTime_avg2 = (obstime/60.)/float(Nsubint_avg2) # to minutes
			scintill_bw_forAvg2 = scintill_bw*Nchan_avg2/float(Nchan)
			scintill_time_forAvg2 = scintill_time*Nsubint_avg2/float(Nsubint)

	
			# 		plot avg..
			print "\n\n------ Averaged by 1/2*scintillation parameters ------\n"

			print "\n -----> Plotting Averaged (by 1/2) Dynamic Spectrum on-pulse...\n"
			plotDynSpectrum.plot_DynSpectrum(data_on_corr_avg2, freqMin, freqMax, binToTime_avg2, 'DynSpectrum_on_avg2.png', output_folder)
			print "\n -----> Plotting Averaged (by 1/2) Dynamic Spectrum off-pulse...\n"
			plotDynSpectrum.plot_DynSpectrum(data_off_corr_avg2, freqMin, freqMax,  binToTime_avg2, 'DynSpectrum_off_avg2.png', output_folder)	
		

			I_on_avg2 = data_on_corr_avg2.copy()
			I_off_avg2 = data_off_corr_avg2.copy()

			if ffinishEarly:
				results_values_avg2_small = printResults.print_results_small(I_on_avg2, bandwidth, obstime, chanToMHz_avg4, binToTime_avg4, output_folder, writeFile, is_avg=True, scintill_bw0=scintill_bw_forAvg4, scintill_time0=scintill_time_forAvg4, small0=small_fit_avg, small0t=small_fit_t_avg, size_t_avg = size_t, size_f_avg=size_f)
			else:
				results_values_avg2_small = printResults.print_results(I_on_avg2, I_off_avg2, bandwidth, obstime, chanToMHz_avg2, binToTime_avg2, output_folder, writeFile, is_avg=True, scintill_bw0=scintill_bw_forAvg2, scintill_time0=scintill_time_forAvg2, fract=1/2., small0=small_fit_avg, small0t=small_fit_t_avg, size_t_avg = size_t2, size_f_avg=size_f2)
	

	 		writeFile.close()
		



