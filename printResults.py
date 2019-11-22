from CorrelationsFULL import *



def print_results(I_on,I_off, bandwidth, obstime, chanToMHz, binToTime, saved_folder, writeFile, is_avg = False, scintill_bw0=0,scintill_time0=0,fract=1/4., small0=False, vsmall0=False, small0t=False, vsmall0t=False, size_t_avg = 1,size_f_avg=1):

	if is_avg==False:
		writeFile.write('\n------ Non-averaged ------\n\n' )
		AC00_on, stdAC00_on, AC00_off, stdAC00_off, CC00, stdCC00, scintill_bw, scintill_bw2, scintill_time, kappa, kappa_off, Neff, eta, sigma_eta, syst_eta, etaMax, sigma_etaMax, syst_etaMax = CorrelationsFULL(I_on, I_off, bandwidth, obstime, saved_folder, small=small0, vsmall = vsmall0, smallt=small0t, vsmallt=vsmall0t)
		
	if is_avg:
		writeFile.write('\n------ Averaged ------\n' )
		strfract = '\n------ 1/'+str(int(1/fract))+' ------\n\n'
		writeFile.write(strfract)
		AC00_on, stdAC00_on, AC00_off, stdAC00_off, CC00, stdCC00, scintill_bw, scintill_bw2, scintill_time, kappa, kappa_off, Neff, eta, sigma_eta, syst_eta, etaMax, sigma_etaMax, syst_etaMax = CorrelationsFULL(I_on, I_off, bandwidth, obstime, saved_folder, isAVG=True, scintill_bwi=scintill_bw0, scintill_timei=scintill_time0, f=int(1/fract), small=small0, vsmall=vsmall0, smallt=small0t, vsmallt=vsmall0t, size_t = size_t_avg, size_f = size_f_avg)
	
	
	
	print '\n -----> Printing results...\n'
	
	# WRITE IN FILE:
		
	scintill_bw_write = 'scintill_bw = ' + str(scintill_bw) + ' chans = ' + str(scintill_bw*np.fabs(chanToMHz)) + ' MHz\n'
	print scintill_bw_write
	writeFile.write(scintill_bw_write)
	
	scintill_bw2_write = 'scintill_bw2 = ' + str(scintill_bw2) + ' chans = ' + str(scintill_bw2*np.fabs(chanToMHz)) + ' MHz\n'
	print scintill_bw2_write
	writeFile.write(scintill_bw2_write)
		
	scintill_time_write = 'scintill_time = ' + str(scintill_time) + ' chans = ' + str(scintill_time*binToTime) + ' min\n'
	print scintill_time_write
	writeFile.write(scintill_time_write)
		
	kappa_write = 'kappa = ' + str(kappa) + ' (within 3.5*sigma)\n'
	print kappa_write
	writeFile.write(kappa_write)
	
	kappa_off_write = 'kappa_off = ' + str(kappa_off) + ' (within 3.5*sigma)\n'
	print kappa_off_write
	writeFile.write(kappa_off_write)
		
	Neff_write = 'Neff = ' + str(Neff) + '\n'
	print Neff_write
	writeFile.write(Neff_write)
		
	AC00_write = 'AC00 = ' + str(AC00_on)  + ' +/- ' + str(stdAC00_on) + '\n'
	print AC00_write
	writeFile.write(AC00_write)
		
	AC00_off_write = 'AC00_off = ' + str(AC00_off)  + ' +/- ' + str(stdAC00_off) + '\n'
	print AC00_off_write
	writeFile.write(AC00_off_write)

	CC00_write = 'CC00 = ' + str(CC00)  + ' +/- ' + str(stdCC00) + '\n'
	print CC00_write
	writeFile.write(CC00_write)
	
	eta_write = 'eta = ' + str(eta) + ' +/- ' + str(sigma_eta) + ' +/- ' + str(syst_eta) + ' (syst)' + '\n'
	print eta_write
	writeFile.write(eta_write)
		
	etaMax_write = 'etaMax = ' + str(etaMax) + ' +/- ' + str(sigma_etaMax) + ' +/- ' + str(syst_etaMax) + ' (syst)' + '\n'
	print etaMax_write
	writeFile.write(etaMax_write)
	
	return AC00_on, stdAC00_on, AC00_off, stdAC00_off, CC00, stdCC00, scintill_bw, scintill_bw2, scintill_time, kappa, kappa_off,Neff, eta, sigma_eta, syst_eta, etaMax, sigma_etaMax, syst_etaMax



def print_results_small(I_on, bandwidth, obstime, chanToMHz, binToTime, saved_folder, writeFile, is_avg = False, scintill_bw0=0,scintill_time0=0,fract=1/4., small0=False, vsmall0=False, small0t=False, vsmall0t=False,size_t_avg = 1,size_f_avg=1):

	if is_avg==False:
		writeFile.write('\n------ Non-averaged ------\n\n' )
		scintill_bw, scintill_bw2, scintill_time, kappa, Neff = CorrelationsFULL_small(I_on, bandwidth, obstime, saved_folder, small=small0, vsmall=vsmall0, smallt=small0t,vsmallt=vsmall0t)
		
	if is_avg:
		writeFile.write('\n------ Averaged ------\n' )
		strfract = '\n------ 1/'+str(int(1/fract))+' ------\n\n'
		writeFile.write(strfract)
		scintill_bw, scintill_bw2, scintill_time, kappa, Neff = CorrelationsFULL_small(I_on, bandwidth, obstime, saved_folder, isAVG=True, scintill_bwi=scintill_bw0, scintill_timei=scintill_time0, f=int(1/fract), small=small0, vsmall=vsmall0, smallt=small0t,vsmallt=vsmall0t,size_t = size_t_avg, size_f = size_f_avg)
	
	
	
	print '\n -----> Printing results...\n'
	
	# WRITE IN FILE:
		
	scintill_bw_write = 'scintill_bw = ' + str(scintill_bw) + ' chans = ' + str(scintill_bw*np.fabs(chanToMHz)) + ' MHz\n'
	print scintill_bw_write
	writeFile.write(scintill_bw_write)
	
	scintill_bw2_write = 'scintill_bw2 = ' + str(scintill_bw2) + ' chans = ' + str(scintill_bw2*np.fabs(chanToMHz)) + ' MHz\n'
	print scintill_bw2_write
	writeFile.write(scintill_bw2_write)
		
	scintill_time_write = 'scintill_time = ' + str(scintill_time) + ' chans = ' + str(scintill_time*binToTime) + ' min\n'
	print scintill_time_write
	writeFile.write(scintill_time_write)
		
	kappa_write = 'kappa = ' + str(kappa) + ' (within 3.5*sigma)\n'
	print kappa_write
	writeFile.write(kappa_write)
			
	Neff_write = 'Neff = ' + str(Neff) + '\n'
	print Neff_write
	writeFile.write(Neff_write)


	return scintill_bw, scintill_bw2, scintill_time, kappa, Neff


	
