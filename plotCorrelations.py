import matplotlib.pyplot as plt

from GaussFit_ACon import *

def plot_CC(CCft, CCt, CCf, binToTime, chanToMHz, plot_name, saved_folder, DOfitGauss = False, GaussFitParams_t = [], twoGaussFitParams_f = [], smallFit=False,   smallFitt=False):

	Nsubint = (CCft.shape[0]+1)/2
	Nchan = (CCft.shape[1]+1)/2
	# set limits for plots:
	
	# min:
	CCfNoNaN = CCf[~np.isnan(CCf)]
	CCtNoNaN = CCt[~np.isnan(CCt)]
	if CCfNoNaN.min() > 0:
		Minf = CCfNoNaN.min()*0.9
		
	else:
		Minf = CCfNoNaN.min()*1.1
	
	if CCtNoNaN.min() > 0:
		Mint = CCtNoNaN.min()*0.9
	
	else:
		Mint = CCtNoNaN.min()*1.1
	Minft = min(Minf,Mint)
	
	# max:	
	CCt_no0 = np.delete(CCt,Nsubint-1)
	CCf_no0 = np.delete(CCf,Nchan-1)
	
	CCfNoNaN_no0 = CCf_no0[~np.isnan(CCf_no0)]
	CCtNoNaN_no0 = CCt_no0[~np.isnan(CCt_no0)]
	
	Maxf = CCfNoNaN_no0.max()*1.1
	Maxt = CCtNoNaN_no0.max()*1.1
	Maxft = max(Maxf,Maxt)
		
	
	# plot:
	CCft_plot = CCft.transpose()
		
	fig = plt.figure()
	# 2D CC(f,t)
	frame1 = fig.add_axes((.32,.1,.52,.6))
	plt.imshow(CCft_plot,aspect = 'auto',interpolation = 'nearest',extent = [(-1)*(Nsubint-1)*binToTime,(Nsubint-1)*binToTime,(Nchan-1)*chanToMHz,(-1)*(Nchan-1)*chanToMHz],vmin = Minft,vmax = Maxft)
	#plt.colorbar()
	frame1.set_yticklabels([])
	plt.xlabel('Time Lag (min)')
	
	
	# 1D CC(t)
	frame2 = fig.add_axes((.32,.73,.52,.23))
	timeLags = np.arange((-1)*(Nsubint-1),Nsubint)
	timeLags = timeLags*binToTime
	plt.plot(timeLags,CCt,'-k')
	if DOfitGauss:
	
		if smallFitt:
			GaussFit_t = Gauss(np.arange(Nsubint+1),*GaussFitParams_t)
			timeLags_small=timeLags[(Nsubint/2-1):(3*Nsubint/2)]
			plt.plot(timeLags_small, GaussFit_t,'-r')
		elif smallFitt==False:
			GaussFit_t = Gauss(np.arange(2*Nsubint-1),*GaussFitParams_t)
			
			plt.plot(timeLags,GaussFit_t,'-r')
			
	
	frame2.set_xticklabels([])
	#plt.yticks(np.linspace(CCt.min(),CCt.max(),num=3))
	plt.locator_params(axis='y', nbins=5)
	plt.ticklabel_format(style = 'sci',axis='y',scilimits=(1,2))
	plt.xlim((-1)*(Nsubint-1)*binToTime,(Nsubint-1)*binToTime)
	plt.ylim(Mint,Maxt)
	plt.ylabel('I (arb.u.)')
	
	
	# 1D CC(f), on the side
	frame3 = fig.add_axes((.1,.1,.2,.6))
	freqLags = np.arange((-1)*(Nchan-1),Nchan)
	freqLags = freqLags*chanToMHz
	
	plt.plot(CCf,freqLags,'-k')
	
	if DOfitGauss:
		
		try:
			if smallFit:
				twoGaussFit_f = twoGauss(np.arange(Nchan+1),*twoGaussFitParams_f)
			elif smallFit==False:
				twoGaussFit_f = twoGauss(np.arange(2*Nchan-1),*twoGaussFitParams_f)
			# twoGaussFit_f = [a1,a2,x0,sigma1,sigma2,yoffset]
			firstGaussFitParams_f = [twoGaussFitParams_f[0],twoGaussFitParams_f[2],twoGaussFitParams_f[3],twoGaussFitParams_f[5]]
			secondGaussFitParams_f = [twoGaussFitParams_f[1],twoGaussFitParams_f[2],twoGaussFitParams_f[4],twoGaussFitParams_f[5]]
			
			if smallFit:
				firstGaussFit_f = Gauss(np.arange(Nchan+1),*firstGaussFitParams_f)
				secondGaussFit_f = Gauss(np.arange(Nchan+1),*secondGaussFitParams_f)
				freqLags_small = freqLags[(Nchan/2-1):(3*Nchan/2)]
				plt.plot(firstGaussFit_f,freqLags_small,'-c')
				plt.plot(secondGaussFit_f,freqLags_small,'-g')
				plt.plot(twoGaussFit_f,freqLags_small,'-r')
				
			elif smallFit==False:
				firstGaussFit_f = Gauss(np.arange(2*Nchan-1),*firstGaussFitParams_f)
				secondGaussFit_f = Gauss(np.arange(2*Nchan-1),*secondGaussFitParams_f)
				plt.plot(firstGaussFit_f,freqLags,'-c')
				plt.plot(secondGaussFit_f,freqLags,'-g')
				plt.plot(twoGaussFit_f,freqLags,'-r')
		

		except TypeError:
			if smallFit:
				GaussFit_f = Gauss(np.arange(Nchan+1),*twoGaussFitParams_f)
				freqLags_small=freqLags[(Nchan/2-1):(3*Nchan/2)]
				plt.plot(GaussFit_f,freqLags_small,'-r')
			elif smallFit==False:
				GaussFit_f = Gauss(np.arange(2*Nchan-1),*twoGaussFitParams_f)
			
				plt.plot(GaussFit_f,freqLags,'-r')
		
			
	plt.xlim(Maxf,Minf)
	plt.ylim((Nchan-1)*chanToMHz,((-1)*(Nchan-1))*chanToMHz)
	#plt.xticks(np.linspace(CCf.max(),CCf.min(),num=3))
	plt.locator_params(axis='x', nbins=5)
	plt.ticklabel_format(style = 'sci',axis='x',scilimits=(1,2))
	plt.xlabel('I (arb.u.)') 
	plt.ylabel('Frequency Lag (MHz)')
	
	
	# colourbar:
	cbar = fig.add_axes((.85,.1,.01,.6))
	axp = frame1.imshow(CCft_plot,aspect = 'auto',interpolation = 'nearest',extent = [(-1)*(Nsubint-1)*binToTime,(Nsubint-1)*binToTime,(Nchan-1)*chanToMHz,(-1)*(Nchan-1)*chanToMHz],vmin = Minft,vmax = Maxft)
	plt.colorbar(axp, cax = cbar)

	plot_path = saved_folder + plot_name
	plt.savefig(plot_path,dpi = 240)
	
	

def plot_AC_withGauss(ACft, ACt, ACf, binToTime, chanToMHz, plot_name, saved_folder,gaussFitParams_t, twogaussFitParams_f,smallFit0=False, smallFit0t=False):

	plot_CC(ACft, ACt, ACf, binToTime, chanToMHz, plot_name, saved_folder, DOfitGauss = True, GaussFitParams_t = gaussFitParams_t, twoGaussFitParams_f = twogaussFitParams_f,smallFit=smallFit0, smallFitt=smallFit0t)
	

	
	
	
	
	
	
	
	
	
	
	
