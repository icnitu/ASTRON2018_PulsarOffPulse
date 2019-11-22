import numpy as np
import matplotlib.pyplot as plt

from Correlations import *
from GaussFit_ACon import *
from Eta import *
from plotCorrelations import *


def CorrelationsFULL(I_on, I_off, bandwidth, obstime, saved_folder, isAVG=False, scintill_bwi=0, scintill_timei=0,f=4, small = False, smallt=False, size_t = 1,size_f = 1):
	
	Nchan = I_on.shape[1]
	Nsubint = I_on.shape[0]
	chanToMHz = bandwidth/float(Nchan)
	binToTime = (obstime/60.)/float(Nsubint) # to minutes

				
	deltaI_on = deltaIx(I_on)
	deltaI_off = deltaIx(I_off)
		
	deltaI_on_NaNto0 = array_NaNto0(deltaI_on)
	deltaI_off_NaNto0 = array_NaNto0(deltaI_off)
			
	print '\n -----> Calculating CrossCorrelation...\n'	
	# --- Crosscorrelation		
	Ngood_CC = Ngood_inCrossCorr(deltaI_on,deltaI_off)
	CCft = CrossCorr_ft(deltaI_on_NaNto0,deltaI_off_NaNto0,Ngood_CC)
	CCt = CrossCorr_t(CCft)
	CCf = CrossCorr_f(CCft)
	
	
	
	print '\n -----> Ploting CrossCorrelation...\n'	

	if isAVG==False:
		# save CCft to file as numpy array
		file_CC_name = saved_folder + 'data_CC'	
		np.save(file_CC_name,CCft)

		plot_CC(CCft, CCt, CCf, binToTime, chanToMHz, 'CrossCorrelation.png', saved_folder)
	if isAVG:
		# save CCft to file as numpy array
		file_CC_avg_name = saved_folder + 'data_CC_avg'	+ str(f)
		np.save(file_CC_avg_name,CCft)
		
		plot_name = 'CrossCorrelation_avg' + str(f) + '.png'
		plot_CC(CCft, CCt, CCf, binToTime, chanToMHz, plot_name, saved_folder)
	
	
	CCt0 = CrossCorr_0(CCt)
	CCf0 = CrossCorr_0(CCf)
	CC00,stdCC00 = CrossCorr_00(CCt0,CCf0)
			
	print '\n -----> Calculating off-pulse AutoCorrelation...\n'	
			
	# --- Autocorrelation off-pulse		
	Ngood_AC_off = Ngood_inCrossCorr(deltaI_off,deltaI_off)
	ACft_off = CrossCorr_ft(deltaI_off_NaNto0,deltaI_off_NaNto0,Ngood_AC_off)
	ACt_off = CrossCorr_t(ACft_off)
	ACf_off = CrossCorr_f(ACft_off)
	
	print '\n -----> Ploting off-pulse AutoCorrelation...\n'	
	
	if isAVG==False:
		# save ACft_off to file as numpy array
		file_ACoff_name = saved_folder + 'data_ACoff'
		np.save(file_ACoff_name,ACft_off)
		
		plot_CC(ACft_off, ACt_off, ACf_off, binToTime, chanToMHz, 'AutoCorrelation_off.png', saved_folder)
	if isAVG:
		# save ACft_off to file as numpy array
		file_ACoff_avg_name = saved_folder + 'data_ACoff_avg' + str(f)
		np.save(file_ACoff_avg_name,ACft_off)
		plot_name = 'AutoCorrelation_off_avg' + str(f) + '.png'
		plot_CC(ACft_off, ACt_off, ACf_off, binToTime, chanToMHz, plot_name, saved_folder)
	
					
	ACt0_off = CrossCorr_0(ACt_off)
	ACf0_off = CrossCorr_0(ACf_off)
	AC00_off,stdAC00_off = CrossCorr_00(ACt0_off,ACf0_off)
		
		
	print '\n -----> Calculating on-pulse AutoCorrelation...\n'	
	
	# --- Autocorrelation on-pulse		
	Ngood_AC_on = Ngood_inCrossCorr(deltaI_on,deltaI_on)
	ACft_on = CrossCorr_ft(deltaI_on_NaNto0,deltaI_on_NaNto0,Ngood_AC_on)
	ACt_on = CrossCorr_t(ACft_on)
	ACf_on = CrossCorr_f(ACft_on)		
			
	ACt0_on = CrossCorr_0(ACt_on)
	ACf0_on = CrossCorr_0(ACf_on)
	AC00_on,stdAC00_on = CrossCorr_00(ACt0_on,ACf0_on)
		
	print '\n -----> Fitting Gaussians to on-pulse AutoCorrelation...\n'	
	
	# only fit one Gaussian to ACt
	if smallt:
		ACt_on_small = ACt_on[(Nsubint/2-1):(3*Nsubint/2)]
		GaussFitParams_t = fitGauss_ACon(ACt_on_small)
	elif smallt==False:	
		GaussFitParams_t = fitGauss_ACon(ACt_on)
		
	# GaussFitParams_t = [a,x0,sigma,yoffset]
	sigma_t = GaussFitParams_t[2]
	
	try:
		# ! fit two Gaussians to frequency plot
		if small:
			ACf_on_small = ACf_on[(Nchan/2-1):(3*Nchan/2)]
			twoGaussFitParams_f = fit2Gauss_ACon(ACf_on_small)
			
		elif small==False:
			twoGaussFitParams_f = fit2Gauss_ACon(ACf_on)
			#twoGaussFitParams_f = fitGauss_ACon(ACf_on)
			
		# twoGaussFitParams_f = [a1,a2,x0,sigma1,sigma2,yoffset]
		print '\n -----> Fitted two Gaussians to ACf...\n'	
		
	
		sigma_f = min(twoGaussFitParams_f[3],twoGaussFitParams_f[4])
		sigma_f2 = max(twoGaussFitParams_f[3],twoGaussFitParams_f[4])
	
	
	except RuntimeError:
		if small:
			twoGaussFitParams_f = fitGauss_ACon(ACf_on_small)
		elif small==False:	
			twoGaussFitParams_f = fitGauss_ACon(ACf_on)
		
		print '\n -----> Fitted one Gaussian to ACf...\n'
		
		sigma_f = twoGaussFitParams_f[2]
		sigma_f2 = sigma_f
	
		
	print '\n -----> Ploting on-pulse AutoCorrelation...\n'	
	
	if isAVG==False:
		# save ACft_on to file as numpy array
		file_ACon_name = saved_folder + 'data_ACon'	
		np.save(file_ACon_name,ACft_on)

		plot_AC_withGauss(ACft_on, ACt_on, ACf_on, binToTime, chanToMHz, 'AutoCorrelation_on.png', saved_folder,GaussFitParams_t, twoGaussFitParams_f, smallFit0=small, smallFit0t=smallt)
	if isAVG:
		# save ACft_on to file as numpy array
		file_ACon_avg_name = saved_folder + 'data_ACon_avg' + str(f)	
		np.save(file_ACon_avg_name,ACft_on)
		plot_name = 'AutoCorrelation_on_avg' + str(f) + '.png'
		plot_AC_withGauss(ACft_on, ACt_on, ACf_on, binToTime, chanToMHz, plot_name, saved_folder,GaussFitParams_t, twoGaussFitParams_f, smallFit0=small, smallFit0t=smallt)
	
	print '\n -----> Calculating scintillation scales...\n'	
	
	# get scintillation scales
	scintill_bw = get_scintill_bw(sigma_f)
	scintill_time = get_scintill_time(sigma_t)
	scintill_bw2 = get_scintill_bw(sigma_f2)	
	
	print '\n -----> Calculating kappa and Neff...\n'	
			
	# get kappa,Neff
	kappa = get_kappa(I_on,size_t0=size_t,size_f0=size_f)
	kappa_off = get_kappa(I_off,size_t0=size_t,size_f0=size_f)
	
	
	
	print '\n -----> Calculating eta and etaMax...\n'	
	
	if isAVG == False:
		# get eta..
		Neff = get_Neff(Nsubint,Nchan,scintill_bw,scintill_time,kappa)
		
				
	if isAVG:
		# keep the old (non-avg) scintill_bw,time
		Neff = get_Neff(Nsubint,Nchan,scintill_bwi,scintill_timei,kappa)
	
	eta, sigma_eta,syst_eta = get_eta(deltaI_on,deltaI_off,CC00,stdCC00,AC00_on,stdAC00_on,Neff)
				
	etaMax,sigma_etaMax,syst_etaMax = get_etaMax(ACft_off,AC00_off,stdAC00_off,AC00_on,stdAC00_on,Neff)
		
	return AC00_on, stdAC00_on, AC00_off, stdAC00_off, CC00, stdCC00, scintill_bw, scintill_bw2, scintill_time, kappa, kappa_off, Neff, eta, sigma_eta, syst_eta, etaMax, sigma_etaMax, syst_etaMax
		

def CorrelationsFULL_small(I_on, bandwidth, obstime, saved_folder, isAVG=False, scintill_bwi=0, scintill_timei=0,f=4, small = False, smallt=False, size_t = 1,size_f = 1):
	
	Nchan = I_on.shape[1]
	Nsubint = I_on.shape[0]
	chanToMHz = bandwidth/float(Nchan)
	binToTime = (obstime/60.)/float(Nsubint) # to minutes

				
	deltaI_on = deltaIx(I_on)
		
	deltaI_on_NaNto0 = array_NaNto0(deltaI_on)
				
	print '\n -----> Calculating on-pulse AutoCorrelation...\n'	
	
	# --- Autocorrelation on-pulse		
	Ngood_AC_on = Ngood_inCrossCorr(deltaI_on,deltaI_on)
	ACft_on = CrossCorr_ft(deltaI_on_NaNto0,deltaI_on_NaNto0,Ngood_AC_on)
	ACt_on = CrossCorr_t(ACft_on)
	ACf_on = CrossCorr_f(ACft_on)		
			
	ACt0_on = CrossCorr_0(ACt_on)
	ACf0_on = CrossCorr_0(ACf_on)
	AC00_on,stdAC00_on = CrossCorr_00(ACt0_on,ACf0_on)
		
	print '\n -----> Fitting Gaussians to on-pulse AutoCorrelation...\n'	
	
	# only fit one Gaussian to ACt
	if smallt:
		ACt_on_small = ACt_on[(Nsubint/2-1):(3*Nsubint/2)]
		GaussFitParams_t = fitGauss_ACon(ACt_on_small)
	elif smallt==False:	
		GaussFitParams_t = fitGauss_ACon(ACt_on)
		
	# GaussFitParams_t = [a,x0,sigma,yoffset]
	sigma_t = GaussFitParams_t[2]
	
	try:
		# ! fit two Gaussians to frequency plot
		if small:
			ACf_on_small = ACf_on[(Nchan/2-1):(3*Nchan/2)]
			twoGaussFitParams_f = fit2Gauss_ACon(ACf_on_small)
			
		elif small==False:
			twoGaussFitParams_f = fit2Gauss_ACon(ACf_on)
			#twoGaussFitParams_f = fitGauss_ACon(ACf_on)
			
		# twoGaussFitParams_f = [a1,a2,x0,sigma1,sigma2,yoffset]
		print '\n -----> Fitted two Gaussians to ACf...\n'	
		
	
		sigma_f = min(twoGaussFitParams_f[3],twoGaussFitParams_f[4])
		sigma_f2 = max(twoGaussFitParams_f[3],twoGaussFitParams_f[4])
	
	
	except RuntimeError:
		if small:
			twoGaussFitParams_f = fitGauss_ACon(ACf_on_small)
		elif small==False:	
			twoGaussFitParams_f = fitGauss_ACon(ACf_on)
		
		print '\n -----> Fitted one Gaussian to ACf...\n'
		
		sigma_f = twoGaussFitParams_f[2]
		sigma_f2 = sigma_f
	
		
	print '\n -----> Ploting on-pulse AutoCorrelation...\n'	
	
	if isAVG==False:
		# save ACft_on to file as numpy array
		file_ACon_name = saved_folder + 'data_ACon'	
		np.save(file_ACon_name,ACft_on)

		plot_AC_withGauss(ACft_on, ACt_on, ACf_on, binToTime, chanToMHz, 'AutoCorrelation_on.png', saved_folder,GaussFitParams_t, twoGaussFitParams_f, smallFit0=small, smallFit0t=smallt)
	if isAVG:
		# save ACft_on to file as numpy array
		file_ACon_avg_name = saved_folder + 'data_ACon_avg' + str(f)	
		np.save(file_ACon_avg_name,ACft_on)
		plot_name = 'AutoCorrelation_on_avg' + str(f) + '.png'
		plot_AC_withGauss(ACft_on, ACt_on, ACf_on, binToTime, chanToMHz, plot_name, saved_folder,GaussFitParams_t, twoGaussFitParams_f, smallFit0=small, smallFit0t=smallt)
	
	print '\n -----> Calculating scintillation scales...\n'	
	
	# get scintillation scales
	scintill_bw = get_scintill_bw(sigma_f)
	scintill_time = get_scintill_time(sigma_t)
	scintill_bw2 = get_scintill_bw(sigma_f2)	
	
	print '\n -----> Calculating kappa and Neff...\n'	
			
	# get kappa,Neff
	kappa = get_kappa(I_on,size_t0=size_t,size_f0=size_f)	


	if isAVG == False:
		# get eta..
		Neff = get_Neff(Nsubint,Nchan,scintill_bw,scintill_time,kappa)
		
				
	if isAVG:
		# keep the old (non-avg) scintill_bw,time
		Neff = get_Neff(Nsubint,Nchan,scintill_bwi,scintill_timei,kappa)

			
	return scintill_bw, scintill_bw2, scintill_time, kappa, Neff
		

		
		
