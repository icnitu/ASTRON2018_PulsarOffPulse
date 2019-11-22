import numpy as np
from scipy.optimize import curve_fit


#def Gauss(x,a,x0,sigma,yoffset):
#	return yoffset+np.fabs(a)*np.exp(-((x-x0)**2)/(2*(sigma**2)))
def Gauss(x,a,x0,sigma):
	return np.fabs(a)*np.exp(-((x-x0)**2)/(2*(sigma**2)))

def fitGauss_ACon(ACx_on):
	
	Nx = ACx_on.shape[0]
	Nx0 = (Nx-1)/2
	
	xfit = np.arange(Nx)
	ACxfit = ACx_on.copy()
	
	# delete centre (noise) - point
	xfit = np.delete(xfit,Nx0)
	ACxfit = np.delete(ACxfit,Nx0)
	
	# check for NaNs
	xfit = xfit[~np.isnan(ACxfit)]
	ACxfit = ACxfit[~np.isnan(ACxfit)]
	
	# guess initial parameters
	a_guess = max(ACxfit)
	#yoffset_guess = min(ACxfit)
	yoffset_guess = 0
	x0_guess = int(Nx/2)
	sigma_guess = np.sqrt(np.abs(np.sum((ACxfit-yoffset_guess)*((xfit-x0_guess)**2))/float(np.sum(ACxfit-yoffset_guess))))
	
	
	#fitParams,fitCov = curve_fit(Gauss,xfit,ACxfit,p0=[a_guess,x0_guess,sigma_guess,yoffset_guess])
	fitParams,fitCov = curve_fit(Gauss,xfit,ACxfit,p0=[a_guess,x0_guess,sigma_guess])
	# fitParams = [a,x0,sigma,yoffset]
	fitParams[0] = np.fabs(fitParams[0])
	fitParams[2] = np.fabs(fitParams[2]) # +ve sigma	
	
	return fitParams
	
#def twoGauss(x,a1,a2,x0,sigma1,sigma2,yoffset):
	
#	return yoffset + np.fabs(a1) * np.exp(-(x - x0)**2 / (2 * sigma1**2)) + np.fabs(a2) * np.exp(-(x - x0)**2 / (2 * sigma2**2))

def twoGauss(x,a1,a2,x0,sigma1,sigma2):
	
	return np.fabs(a1) * np.exp(-(x - x0)**2 / (2 * sigma1**2)) + np.fabs(a2) * np.exp(-(x - x0)**2 / (2 * sigma2**2))


def fit2Gauss_ACon(ACx_on):

	Nx = ACx_on.shape[0]
	Nx0 = (Nx-1)/2
	
	xfit = np.arange(Nx)
	ACxfit = ACx_on.copy()
	
	# delete centre (noise) - point
	xfit = np.delete(xfit,Nx0)
	ACxfit = np.delete(ACxfit,Nx0)
	
	# check for NaNs
	xfit = xfit[~np.isnan(ACxfit)]
	ACxfit = ACxfit[~np.isnan(ACxfit)]
		
	# guess initial parameters
	a1_guess = max(ACxfit)
	a2_guess = a1_guess/2.
	#yoffset_guess = min(ACxfit)
	yoffset_guess = 0
	x0_guess = int(Nx/2)
	sigma1_guess = np.sqrt(np.abs(np.sum((ACxfit-yoffset_guess)*((xfit-x0_guess)**2))/float(np.sum(ACxfit-yoffset_guess))))
	sigma2_guess = sigma1_guess*3.
	
	#fitParams,fitCov = curve_fit(twoGauss,xfit,ACxfit,p0=[a1_guess,a2_guess,x0_guess,sigma1_guess,sigma2_guess,yoffset_guess])
	fitParams,fitCov = curve_fit(twoGauss,xfit,ACxfit,p0=[a1_guess,a2_guess,x0_guess,sigma1_guess,sigma2_guess])
	# fitParams = [a1,a2,x0,sigma1,sigma2,yoffset]
	fitParams[0] = np.fabs(fitParams[0])
	fitParams[1] = np.fabs(fitParams[1])
	fitParams[3] = np.fabs(fitParams[3]) # +ve sigma1	
	fitParams[4] = np.fabs(fitParams[4]) # +ve sigma2	
	
	return fitParams




def get_scintill_bw(sigma_f):
	# HWHM
	sigma_f = np.fabs(sigma_f)
	scintill_bw = sigma_f*np.sqrt(2*np.log(2.))
	
	return scintill_bw

def get_scintill_time(sigma_t):
	# HW(1/e)M
	sigma_t = np.fabs(sigma_t)
	scintill_time = sigma_t*np.sqrt(2.)
	
	return scintill_time




