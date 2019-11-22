import numpy as np
from scipy import signal




def deltaIx(Ix):
	meanIx = np.mean(Ix[~np.isnan(Ix)])
	deltaIx= Ix - meanIx
	return deltaIx

def array_NaNto0(array):
	array_NaNto0 = array.copy()
	array_NaNto0[np.isnan(array_NaNto0)]=0.
	
	return array_NaNto0
	

def Ngood_inCrossCorr(deltaI_on,deltaI_off):
	# number of not NaN pixels in each pixel of CC
	bool_on = (~np.isnan(deltaI_on))*1.
	bool_off = (~np.isnan(deltaI_off))*1.
	Ngood = signal.correlate2d(bool_on,bool_off) * 1.
	Ngood [Ngood==0] = np.nan
	
	return Ngood
	
	
# these can be used for AC_on,AC_off,CC
def CrossCorr_ft(deltaI_on_NaNto0,deltaI_off_NaNto0,Ngood):

	# Ngood.shape == CC.shape
	Nsubint = deltaI_on_NaNto0.shape[0]
	Nchan = deltaI_on_NaNto0.shape[1]
	
	
	CCft = signal.correlate2d(deltaI_on_NaNto0,deltaI_off_NaNto0)
	CCft = np.divide(CCft,np.sqrt(float(Nsubint*Nchan)*Ngood))
	
	return CCft
	

def CrossCorr_t(CCft):
	
	Nchan = (CCft.shape[1]+1)/2
	CCt = CCft[:,Nchan-1]
	return CCt

def CrossCorr_f(CCft):
	
	Nsubint = (CCft.shape[0]+1)/2
	CCf = CCft[Nsubint-1,:]
	return CCf

def CrossCorr_0(CCx):
	
	Nx = (CCx.shape[0]+1)/2
	# ignore central point as it's probably noise
	# try to fit a -ve parabola to the 4 points around
	# 	as an approx to Gaussian
	
	# fit a parabola to (x,y):
	y = np.empty(4)
	y[0] = CCx[Nx-3]; y[1] = CCx[Nx-2]; y[2] = CCx[Nx]; y[3] = CCx[Nx+1]
	x = np.array([Nx-3,Nx-2,Nx,Nx+1])
	
	# check for NaNs
	y = y[~np.isnan(y)]
	x = x[~np.isnan(y)]
	
	Nfit = len(y)
	if Nfit == 0:
		print " \n Warning: all 4 points around centre are NaNs!\n"
		return np.nan
		
	# fit 2nd order polynomial
	poly2_coeff = np.polyfit(x,y,2)
	# p[0]* x**2 + p[1]*x + p[2]
	if poly2_coeff[0] < 0:
		# good, -ve parabola
		# max of parabola
		CCx0 =  poly2_coeff[2] - ((poly2_coeff[1])**2)/(4.* poly2_coeff[0])
	else:
		CCx0 = y.mean() # fit horizontal line
		
	return CCx0

def CrossCorr_00(CCt0,CCf0):
	if (~np.isnan(CCt0)) and (~np.isnan(CCf0)):
		CC00 = (CCf0+CCt0)/2.
		stdCC00 = np.sqrt((CCt0-CC00)**2 + (CCf0-CC00)**2)
	elif (np.isnan(CCt0)) and (~np.isnan(CCf0)):
		CC00 = CCf0
		stdCC00 = np.sqrt((CCf0-CC00)**2)
	elif (~np.isnan(CCt0)) and (np.isnan(CCf0)):
		CC00 = CCt0
		stdCC00 = np.sqrt((CCt0-CC00)**2)
	else:
		CC00 = np.nan
		stdCC00 = np.nan
	
	return CC00,stdCC00 
		

def CrossCorr_00_full(CCft):
	
	deltaI_on_NaNto0 = array_NaNto0(deltaI_on)
	deltaI_off_NaNto0 = array_NaNto0(deltaI_off)
					
	Ngood = Ngood_inCrossCorr(deltaI_on,deltaI_off)
	CCft = CrossCorr_ft(deltaI_on_NaNto0,deltaI_off_NaNto0,Ngood)
	
	CCt = CrossCorr_t(CCft)
	CCf = CrossCorr_f(CCft)
	CCt0 = CrossCorr_0(CCt)
	CCf0 = CrossCorr_0(CCf)
	CC00,stdCC00 = CrossCorr_00(CCt0,CCf0)
	
	return CC00,stdCC00

