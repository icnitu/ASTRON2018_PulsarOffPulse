import numpy as np




def get_rms_off(data_off):
	
	Nchan = data_off.shape[1]
	rms_off = np.empty(Nchan)
	for ichan in range(Nchan):
		data_off_ichan = data_off[:,ichan].copy()
		data_off_ichan = data_off_ichan[~np.isnan(data_off_ichan)]
		if len(data_off_ichan) >0:
			rms_off[ichan]=data_off_ichan.std()
		else:
			rms_off[ichan] = np.nan
	return rms_off	

# weights, rms, cut freq edges, cut last subint
def do_corrections1(data_ox,weights,Nox,Noff,rms_off,fract1=0.03,Nsubint_last1=0):
	# data_ox is {data_on,data_off}, (Nsubint,Nchan)
	# rms_off = (Nchan)
	
	Nchan = data_ox.shape[1]
	Nsubint = data_ox.shape[0]
	data_ox_corr = data_ox.copy()
	# weight correction:
	print "\n -------> Applying weight corrections...\n"
	data_ox_corr[weights==0] = np.nan
	
	# rms:
	print "\n -------> Applying RMS corrections...\n"
	rms_off[rms_off==0] = np.nan
	data_ox_corr = np.divide(data_ox_corr,rms_off)
	data_ox_corr = data_ox_corr*np.sqrt(Nox/float(Noff))
	
	# cut freq edges:
	print "\n -------> Cutting a fraction of " + str(fract1) + " of frequency edges...\n"
	dNchan = int(Nchan*fract1)
	data_ox_corr[:,0:dNchan] = np.nan
	data_ox_corr[:,(Nchan-dNchan):Nchan] = np.nan
	
	# cut last (1/few) subint channel(s)
	print "\n -------> Cutting " + str(Nsubint_last1) + " of last subintegration channels...\n"
	
	data_ox_corr[(Nsubint-Nsubint_last1):Nsubint,:] = np.nan
	
	return data_ox_corr
	

def correct_data_on_off(data_on,data_off,weights,Non,Noff,fract = 0.03,Nsubint_last=0):
	
	Nchan = data_off.shape[1]
	Nsubint = data_off.shape[0]
	
	
	rms_off = get_rms_off(data_off)
	
	print "\n -----> Applying corrections to on-pulse...\n"
	data_on_corr = do_corrections1(data_on,weights,Non,Noff,rms_off,fract1=fract,Nsubint_last1=Nsubint_last)
	
	print "\n -----> Applying corrections to off-pulse...\n"
	data_off_corr = do_corrections1(data_off,weights,Noff,Noff,rms_off,fract1=fract,Nsubint_last1=Nsubint_last)
	
	return data_on_corr,data_off_corr
	
	
	
	
	
	
	
	
	
