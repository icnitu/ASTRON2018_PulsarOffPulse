import matplotlib.pyplot as plt
import numpy as np


def plot_DynSpectrum(data_ox, freqMin, freqMax, binToTime, plot_name,saved_folder,plot_Naway = 5):
	
	Nsubint = data_ox.shape[0]
	data_ox_noNaN = data_ox.copy()
	data_ox_noNaN = data_ox_noNaN[~np.isnan(data_ox_noNaN)]
	
	mean = data_ox_noNaN.mean()
	std = data_ox_noNaN.std()
	
	intensMin = mean - plot_Naway*std
	intensMax = mean + plot_Naway*std
	
	data_ox_plot = data_ox.transpose()
	
	plt.figure()
	plt.imshow(data_ox_plot,aspect = 'auto',interpolation = 'nearest',cmap = 'afmhot',extent = [0,Nsubint*binToTime,freqMax,freqMin], vmin=intensMin, vmax=intensMax)
	plt.colorbar()
	plt.xlabel('Time (min)')
	plt.ylabel('Frequency (MHz)')
	plot_path = saved_folder + plot_name
	plt.savefig(plot_path, dpi=240)
	
	
	
