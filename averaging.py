import numpy as np


def get_size_avg(scintill_size,fraction=1/4.):
	
	if (int(scintill_size*fraction))%2 == 0:
		size = int(scintill_size*fraction)+1
	else:
		size = int(scintill_size*fraction)

	return size

def average(data,size_t,size_f):

# TODO: say max of size_i for which can avg...
# ...is such that i/size_i >~5 ?

	# size_t and size_f need to be odd
	# make sure they're int
	size_t = int(size_t)
	size_f = int(size_f)

	Nt = data.shape[0]
	Nf = data.shape[1]
	
	# make sure they're odd
	if size_t%2 == 0:
		size_t = size_t+1
	if size_f%2 == 0:
		size_f = size_f+1
	
	Nmax_t = int(Nt/size_t)
	Nmax_f = int(Nf/size_f)
	
	if Nmax_t < 5:
		size_t = int(Nt/5.)
		if size_t%2 == 0:
			size_t = size_t-1
		Nmax_t = 5
		print " !!! WARNING: the time averaging length too large; used a fifth of the total time interval instead."
		
	if Nmax_f < 5:
		size_f = int(Nf/5.)
		if size_f%2 == 0:
			size_f = size_f-1
		Nmax_f = 5	
		print " !!! WARNING: the frequency averaging length too large; used a fifth of the total frequency interval instead."


	half_t = (size_t-1)/2
	half_f = (size_f-1)/2	
	
	
	data_avg = np.empty([(Nmax_t+1),(Nmax_f+1)])
	tavg = 0
	favg = 0
	
	for it in range(half_t,half_t+Nmax_t*size_t,size_t):
		range_t = data[it-half_t:it+half_t+1,:]
		# region 1
		favg = 0
		for ifr in range(half_f,half_f+Nmax_f*size_f,size_f):
			range_ft = range_t[:,ifr-half_f:ifr+half_f+1]
			range_ft_noNaN = range_ft[~np.isnan(range_ft)]
			Ni = range_ft_noNaN.size # how many 'usable' pixels
			if Ni ==0:
				data_avg[tavg][favg] = np.nan
			else:			
				data_avg[tavg][favg] = np.sum(range_ft_noNaN)*np.sqrt(1./float(size_t*size_f*Ni))			
		
			favg = favg+1
		# region 2
		range_ft = range_t[:,Nmax_f*size_f:Nf]
		range_ft_noNaN = range_ft[~np.isnan(range_ft)]
		Ni = range_ft_noNaN.size # how many 'usable' pixels
		if Ni ==0:
			data_avg[tavg][favg] = np.nan
		else:			
			data_avg[tavg][favg] = np.sum(range_ft_noNaN)*np.sqrt(1./float(size_t*size_f*Ni))	
		
		tavg = tavg+1	
	
	range_t = data[Nmax_t*size_t:Nt,:]
	# region 3 
	favg = 0
	for ifr in range(half_f,half_f+Nmax_f*size_f,size_f):
		range_ft = range_t[:,ifr-half_f:ifr+half_f+1]
		range_ft_noNaN = range_ft[~np.isnan(range_ft)]
		Ni = range_ft_noNaN.size # how many 'usable' pixels
		if Ni ==0:
			data_avg[tavg][favg] = np.nan
		else:
			data_avg[tavg][favg] = np.sum(range_ft_noNaN)*np.sqrt(1./float(size_t*size_f*Ni))			
		
		favg = favg+1
	# region 4
	range_ft = range_t[:,Nmax_f*size_f:Nf]
	range_ft_noNaN = range_ft[~np.isnan(range_ft)]
	Ni = range_ft_noNaN.size # how many 'usable' pixels
	if Ni ==0:
		data_avg[tavg][favg] = np.nan
	else:			
		data_avg[tavg][favg] = np.sum(range_ft_noNaN)*np.sqrt(1./float(size_t*size_f*Ni))	

	return data_avg
	
	
	
		
