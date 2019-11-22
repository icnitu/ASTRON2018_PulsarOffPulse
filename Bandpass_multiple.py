import numpy as np
import matplotlib.pyplot as plt
import warnings

from Bandpass import *


def get_each_fit_bandpass(data_off_bandpass,xchans,poly_order_i=8,Nrepeat_i = 10,Naway0_i = 10):
	
#	warnings.warn("badfit", RankWarning)

	Nchan = data_off_bandpass.shape[0]
	xi = xchans.copy()
	yi = data_off_bandpass.copy()
	# clean of NaNs:
	xi = xi[~np.isnan(yi)]
	yi = yi[~np.isnan(yi)]
	try:	

		for irepeat in range(Nrepeat_i):
			Naway = max(float(Naway0_i)/(1.5**irepeat),3.)
			print "\n !!!!!1.0 DEBUG!!!!!\n"
			poly_coeff = np.polyfit(xi,yi,poly_order_i)
			residuals = yi-np.polyval(poly_coeff,xi)
			residuals_mean = residuals.mean()
			residuals_std = residuals.std()
		
			awayMin = residuals_mean - Naway*residuals_std
			awayMax = residuals_mean + Naway*residuals_std
			# 'clean' of points outside [awayMin,awayMax]
			xClean = xi[(residuals>=awayMin) & (residuals<=awayMax)]
			yClean = yi[(residuals>=awayMin) & (residuals<=awayMax)]
			xi = xClean
			yi = yClean
		print "\n !!!!!2.0 DEBUG!!!!!\n"
			
		# now fit final x,y
		poly_coeff_fit = np.polyfit(xi,yi,poly_order_i)	
	
	except RankWarning:
		print "\nWARNING: Polyfit may be poorly conditioned.\n"		

	# return coefficients of the polynomial
	return poly_coeff_fit


#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")


def get_multiple_fits_bandpass(data_off_bandpass,Nbands_mult=1,poly_order_mult=8,Nrepeat_mult=10,Naway0_mult=10):

	# for bandpass from Nbands different bands, which need to be fit separately

	Nchan = data_off_bandpass.shape[0]
	Nchan_subset = int(Nchan/Nbands_mult)

	if int(Nchan/Nbands_mult) != (Nchan/float(Nbands_mult)):
		print "/nError: the bandpass does not have an integer number of bands!/n"

	subsets_bandpass = data_off_bandpass.reshape(Nbands_mult,Nchan_subset)

	poly_vals_fits = np.empty(Nchan)
	
	for iband in range(Nbands_mult):	
		ipoly_coeff_fit_subset = None
		ixchan = np.arange(iband*Nchan_subset,(iband+1)*Nchan_subset)
		ipoly_coeff_fit_subset = get_each_fit_bandpass(subsets_bandpass[iband], ixchan, poly_order_i = poly_order_mult, Nrepeat_i=Nrepeat_mult, Naway0_i=Naway0_mult)
		
		poly_vals_fits[iband*Nchan_subset:(iband+1)*Nchan_subset] = np.polyval(ipoly_coeff_fit_subset,ixchan)


	return poly_vals_fits



def get_multiple_bandpass_with_fit(data3, weights, offBin_min, offBin_max, Nbands=1, poly_order=8, Nrepeat=10, Naway0=10):

	data_off_bandpass = get_initial_bandpass(data3, weights, offBin_min, offBin_max)

	if Nbands==0:
		# no bandpass correction!
		Nchan = data_off_bandpass.shape[0]
		poly_vals_fits = np.ones([Nchan])
	else:
		poly_vals_fits = get_multiple_fits_bandpass(data_off_bandpass, Nbands_mult=Nbands, poly_order_mult=poly_order, Nrepeat_mult=Nrepeat, Naway0_mult=Naway0)
		
	return data_off_bandpass,poly_vals_fits


def plot_multiple_bandpass_with_fit(data_off_bandpass, poly_vals_fits, freqMin,freqMax, saved_folder):
	

	Nchan = data_off_bandpass.shape[0]
	xfreqs = np.linspace(freqMin,freqMax,num=Nchan)
	xchan = np.arange(Nchan)
	
	fig = plt.figure()
	
	frameMain = fig.add_axes((.1,.3,.8,.6))
	# plot data
	plt.plot(xfreqs,data_off_bandpass,'.k')
	# plot fit on top
	plt.plot(xfreqs,poly_vals_fits,'-r')
	plt.xlim(freqMin,freqMax)

	#plt.ylim(data_off_bandpass.min(),data_off_bandpass.max())	
	#plt.ylim(-2000,4000)
	
	frameMain.set_xticklabels([])
	plt.ylabel('Intensity (arb.u.)')
	
	
	# plot residuals below:
	residuals = data_off_bandpass - poly_vals_fits
	residuals_noNaN = residuals[~np.isnan(residuals)]
	residuals_std = residuals_noNaN.std()
	
	frameRes = fig.add_axes((.1,.1,.8,.15))
	plt.plot(xfreqs,residuals,'xk')
	# horizontal line at 0 (fit)
	plt.locator_params(axis='y',tight=True,nbins=6)
	plt.axhline(y=0,c='r')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Resid.')
	plt.ylim((-3)*residuals_std,3*residuals_std)
	plt.xlim(freqMin,freqMax)
	plot_name = saved_folder+'Bandpass_with_fit.png'
	plt.savefig(plot_name, dpi=240)
	

	
