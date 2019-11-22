import matplotlib.pyplot as plt



def plot_FreqVsPhase(data_freqphase,freqMin,freqMax,saved_folder):

	Nbin = data_freqphase.shape[1]
	plt.figure()
	plt.imshow(data_freqphase,aspect = 'auto',interpolation = 'nearest',extent = [0,1,freqMax,freqMin])
	plt.xlabel('Pulse Phase')
	plt.ylabel('Frequency (MHz)')
	
	plot_name = saved_folder+'FrequencyVsPhase.png'
	
	plt.savefig(plot_name,dpi = 240)
	
	
