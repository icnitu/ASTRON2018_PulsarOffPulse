import numpy as np


def get_kappa(I_on,size_t0 = 1,size_f0 = 1):

	Nsubint = I_on.shape[0]
	Nchan = I_on.shape[1]
	I_scintles = I_on[I_on > (3.5/np.sqrt(float(size_t0*size_f0)))] # more than 3*std
	# I_scintles flattened
	Nk = I_scintles.size
	#print '\nNk = ',Nk, '\n'
	k = Nk/float(Nsubint*Nchan)

	return k

def get_Neff(Nsubint,Nchan,scintill_bw,scintill_time,k):
	
	Nt = 1 + k*(Nsubint/float(scintill_time))
	Nf = 1 + k*(Nchan/float(scintill_bw))
	
	Neff = Nt*Nf
	return Neff
	

def get_eta(deltaI_on,deltaI_off,CC00,stdCC00,AC00,stdAC00,Neff):
	
	eta = CC00/AC00
	U = deltaI_off - eta*deltaI_on
	sigmaU = U[~np.isnan(U)].std()
	
	sigma_eta = sigmaU/(np.sqrt(Neff*AC00))
	
	syst_eta = np.sqrt((stdCC00/AC00)**2 + (stdAC00*CC00/(AC00**2))**2)
	# syst_eta is due to the std of AC00 etc 
	# 	when calculating from both t & f slices
	
	return eta, sigma_eta,syst_eta
	
def get_etaMax(ACft_off,AC00_off,stdAC00_off,AC00,stdAC00,Neff):

	if AC00_off < 0:
		etaMax = np.nan
		syst_etaMax = np.nan
		print "\n Warning: Autocorrelation for off pulse at (0,0) is negative, so etaMax cannot be calculated.\n"
	else:
		etaMax = np.sqrt(AC00_off/AC00)
		syst_etaMax = etaMax*np.sqrt((stdAC00/(2.*AC00))**2 + (stdAC00_off/(2.*AC00_off))**2)
		
	
	Nsubint = (ACft_off.shape[0]+1)/2
	Nchan = (ACft_off.shape[1]+1)/2
	
	sigma_off = np.sqrt(ACft_off[Nsubint-1,Nchan-1] - AC00_off)
	
	sigma_etaMax = sigma_off/np.sqrt(Neff*AC00)
	
	return etaMax,sigma_etaMax,syst_etaMax
	
	
