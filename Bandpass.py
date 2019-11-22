import numpy as np
import matplotlib.pyplot as plt


def get_initial_bandpass(data3,weights,offBin_min,offBin_max):
	
	print "\n -----> Applying initial weights corrections to bandpass...\n"
					
	data3_copy = data3.copy() # (Nsubint,Nchan,Nbin)
	Nsubint = data3_copy.shape[0]
	Nchan = data3_copy.shape[1]
	
	data_off_weighted = data3_copy[:,:,offBin_min:offBin_max]
	data_off_weighted = data_off_weighted.mean(2) # (Nsubint,Nchan)
	# turn points with weight = 0 to NaNs
	data_off_weighted[weights == 0] = np.nan
	# delete subint channels that are fully NaNs
	data_off_weighted = data_off_weighted[~np.all(np.isnan(data_off_weighted),axis=1)]

	#data_off_weighted

	data_off_weighted_T = data_off_weighted.transpose()
	
	data_off_bandpass = np.empty(Nchan)

	for ichan in range(Nchan):
		data_off_bandpass[ichan] = np.median((data_off_weighted_T[ichan])[~np.isnan(data_off_weighted_T[ichan])],0) # (Nchan)
	
	#print "\n\n ***********DEBUG \n\n"
	
	
	xy=np.empty([Nsubint,Nchan])
	for ix in range(Nsubint):
		for iy in range(Nchan):
			xy[ix][iy] = ix*10000+iy
	
	#print "/n xy = ", xy
	#print "data_off_weighted = ", xy[np.isnan(data_off_weighted)], 
	return data_off_bandpass

def get_fit_bandpass(data_off_bandpass,poly_order_i=8,Nrepeat_i = 10,Naway0_i = 10):

	Nchan = data_off_bandpass.shape[0]
	xi = np.arange(Nchan)
	yi = data_off_bandpass.copy()
	# clean of NaNs:
	xi = xi[~np.isnan(yi)]
	yi = yi[~np.isnan(yi)]
	
	for irepeat in range(Nrepeat_i):
		Naway = max(float(Naway0_i)/(1.5**irepeat),3.)
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

	# now fit final x,y
	poly_coeff_fit = np.polyfit(xi,yi,poly_order_i)	
	
	# return coefficients of the polynomial
	return poly_coeff_fit
	
def get_bandpass_with_fit(data3,weights,offBin_min,offBin_max,poly_order=8,Nrepeat=10,Naway0=10):

	data_off_bandpass = get_initial_bandpass(data3,weights,offBin_min,offBin_max)
	poly_coeff_fit = get_fit_bandpass(data_off_bandpass,poly_order_i = poly_order,Nrepeat_i = Nrepeat,Naway0_i = Naway0)
	
	return data_off_bandpass,poly_coeff_fit


def plot_bandpass_with_fit(data_off_bandpass,poly_coeff_fit,freqMin,freqMax,saved_folder):

	Nchan = data_off_bandpass.shape[0]
	xfreqs = np.linspace(freqMin,freqMax,num=Nchan)
	xchan = np.arange(Nchan)
	poly_val_fit = np.polyval(poly_coeff_fit,xchan)
	
	fig = plt.figure()
	
	frameMain = fig.add_axes((.1,.3,.8,.6))
	# plot data
	plt.plot(xfreqs,data_off_bandpass,'.k')
	# plot fit on top
	plt.plot(xfreqs,poly_val_fit,'-r')
	plt.xlim(freqMin,freqMax)
	
	frameMain.set_xticklabels([])
	plt.ylabel('Intensity (arb.u.)')
	
	
	# plot residuals below:
	residuals = data_off_bandpass - poly_val_fit
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
	
	
		
	
	
	
	

	
	

