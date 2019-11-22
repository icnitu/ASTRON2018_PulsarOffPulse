import matplotlib.pyplot as plt
import numpy as np


def plot_PulseProfile(data_profile, saved_folder, onBin_fmin=np.empty(1), onBin_fmax = np.empty(1), offBin_fmin = np.empty(1), offBin_fmax = np.empty(1)):
	
	Nbin = data_profile.shape[0]
	xbins = np.arange(Nbin)
	
	xphase = xbins/float(Nbin)

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	plt.plot(xphase,data_profile,'-b')
		
	plt.xlabel('Pulse Phase')
	plt.ylabel('Intensity (arb.u.)')
	
	#plt.xticks(np.arange(0,Nbin,step=100))
	#plt.grid()
	# to see the lines if at edges of graphs:
	if onBin_fmin<=0.001:
		onBin_fmin = 0.002
	if offBin_fmax>=0.999:
		offBin_fmax = 1-0.002
	
	ax.set_xticks(np.arange(1.1,step=0.1))
	ax.set_xticks(np.arange(1.1,step=0.02),minor=True)
	ax.grid(which='minor',alpha = 0.4)
	ax.grid(which='major',alpha = 0.5)
	plt.xlim(0,1)
	plt.axvline(x=onBin_fmin, ymin=data_profile.min(), ymax=data_profile.max(),c='r',label = 'On-pulse')
	plt.axvline(x=onBin_fmax, ymin=data_profile.min(), ymax=data_profile.max(),c='r')
	plt.axvline(x=offBin_fmin, ymin=data_profile.min(), ymax=data_profile.max(),ls='--',c='g',label = 'Off-pulse')
	plt.axvline(x=offBin_fmax, ymin=data_profile.min(), ymax=data_profile.max(),ls='--',c='g')
	plt.legend()

	plot_name = saved_folder+'PulseProfile.png'
	
	plt.savefig(plot_name,dpi = 240)
	
	
