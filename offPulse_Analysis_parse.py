import psrchive as pc
import numpy as np
import matplotlib.pyplot as plt
#from glob2 import glob
from glob import glob
import os
import argparse
from collections import Iterable
import sys
import inspect

import archive_analysis


#  Need: PATH
#	FILES NAMES/ */ etc
#  Use glob2 


def read_txt(txt_files_name):

	txt_files = open(txt_files_name,'r')
	ispath=True
	data_files = []
	while True:
		line=txt_files.readline()
		if not line: break
		
		# allow comments with '#'
		if line[0]!='#' and line!='\n':
			if ispath:
				path = line.strip()
				ispath = False
			else:
				data_files.append(line.strip())

	txt_files.close()
	return path,data_files

def flatten(items):
    """Yield items from any nested iterable; see Reference."""
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x


def fraction(string):
		
	value=float(string)

	if value <0 or value >1:
		message = "%r is not a fractional value between [0,1)." %string
		raise argparse.ArgumentTypeError(message)

	return value

def fractionHalf(string):
		
	value=float(string)

	if value <0 or value >0.5:
		message = "%r is not a fractional value between [0,0.5)." %string
		raise argparse.ArgumentTypeError(message)

	return value


def posint(string):
		
	value=int(string)

	if value <=0:
		message = "%r is not a positive integer." %string
		raise argparse.ArgumentTypeError(message)

	return value

def setValue(X):	

	try:
		x=X[0]
	except:
		x=X
	
	return x

def setOnValue(onpulse_array):

	try:
		onnBin_fmin = min(onpulse_array[0],onpulse_array[1]) 
		onnBin_fmax = max(onpulse_array[0],onpulse_array[1])
		doOnnPulse = True
	except:
		print "\nSomething wrong with on pulse range!\n"

	return doOnnPulse, onnBin_fmin, onnBin_fmax

def setOffValue(offpulse_array):
	doOfPulse = False
	# TODO: update this
	try:
		ofBin_fmin = min(offpulse_array[0], offpulse_array[1]) 
		ofBin_fmax = max(offpulse_array[0], offpulse_array[1])
		doOfPulse=True
		#print "\nThis is the full off-pulse analysis.\n"
	except:
		ofBin_fmin=None
		ofBin_fmax=None
		print "\nThis analysis will stop at giving the Pulse Profile.\n"	
	
	# TODO change this to appropriate message!

	return doOfPulse, ofBin_fmin, ofBin_fmax



# ---------- MAIN

# ^^^^^^^^^^ ADD ARGS

current_dir = os.path.dirname(os.path.realpath(__file__)) # can just use ./

# TODO: plot the limits of on\off pulse on pulse profile



# TODO: change prog to the final name of file
parser = argparse.ArgumentParser(prog='python offPulse_Analysis_parse.py')

# file (.txt) to contain path and files
# Note: works with *

parser.add_argument('files_list', metavar = 'text_file', 
			help='Text document containing path and data files names')
# OPTIONAL ARGS:

parser.add_argument('-s', '--savef', dest = 'output_folder',metavar = 'output_folder', 	
			help='Folder containing all the output files', default=current_dir)

parser.add_argument('-on', '--onpulse', nargs=2, type=fraction, metavar = 'X', default = [0,1-1e-5],
			help='On-pulse range as a sequence of min and max fractions of total phase. Default is the full observation. Default is the full signal.')

parser.add_argument('-off','--offpulse', nargs=2, type=fraction, metavar = 'X', default = [0,1-1e-5],
			help='Off-pulse range as a sequence of min and max fractions of total phase. Default is the full signal.')

parser.add_argument('-t', '--translate', nargs=1, type=fraction, metavar = 'X', 
			help='Translate the peak to a specified position as a fraction of total phase.')

parser.add_argument('-zf', '--zapfreq', nargs=1, type=fractionHalf, metavar = 'X', default=0.0,
			help='Zap this fraction of frequency channels from the edges. Default is 0 (i.e. 0%%).')

parser.add_argument('-b', '--nband', nargs=1, type=posint, metavar = 'X', default=1,
			help='Number of bands in the bandpass. Default is 1.')

parser.add_argument('-cgf', '--centralGfreq', action='store_const', const='y', default='n',
			help='Only attempt a Gaussian fit to the central part of the frequency AutoCorrelation function. Necessary sometimes because of high noise.')

parser.add_argument('-cgt', '--centralGtime', action='store_const', const='y', default='n',
			help='Only attempt a Gaussian fit to the central part of the time AutoCorrelation function. Necessary sometimes because of high noise.')

parser.add_argument('-f', '--finish', action='store_const', const=True, default = False,
			help='Finish the analysis at scintillation bandwidth and timescale.' )


# TODO: make sure zapping is done BEFORE translating..	
# 	more complicated then it seems!!

# TODO: warning: scintill_bw value used is always the first one printed!


parser_args = parser.parse_args()


# ^^^^^^^^^^ CHECK INPUTS

# Check if input file exists
try:
	txt_files_name = glob(parser_args.files_list)[0]
except IndexError:
	print "\nERROR: The file doesn't exist!\n"
	sys.exit()


# Translation of peak	
phaseOffset = setValue(parser_args.translate)

# Zapping freq channels
freq_fract = setValue(parser_args.zapfreq)

# Multiple bands in bandpass
Nbands_number = setValue(parser_args.nband)

# On/Off bin fraction ranges
doOnPulse, onBin_fmin, onBin_fmax = setOnValue(parser_args.onpulse)
doOffPulse, offBin_fmin, offBin_fmax = setOffValue(parser_args.offpulse)
			
# Only central Gaussian fit
centralGf = setValue(parser_args.centralGfreq)
centralGt = setValue(parser_args.centralGtime)

# Finish?
finishEarly = setValue(parser_args.finish)
#print 'finishEarly = ',finishEarly
if finishEarly:
	print "\nThis analysis will stop at giving the scintillation parameters.\n"
else:
	print "\nThis is the full off-pulse analysis.\n"


# ^^^^^^^^^^ READ IN DATA FILES. CREATE OUTPUT FOLDER

filename = inspect.getframeinfo(inspect.currentframe()).filename
current_folder = os.path.dirname(os.path.abspath(filename))+'/'

saved_folder = parser_args.output_folder

path,data_files = read_txt(txt_files_name)

if path[0]=='~':
	path=os.path.expanduser("~")+path[1:]
if saved_folder[0]=='~':
	saved_folder=os.path.expanduser("~")+saved_folder[1:]
if saved_folder[0]=='.' and saved_folder[1]=='/':
	saved_folder=current_folder+saved_folder[2:]
if saved_folder[-1]!='/':
	saved_folder = saved_folder+'/'


#print 'saved_folder = ',saved_folder

archive_files=[]
data_files_actual=[]
for ifile in data_files:
	ifile_actual = glob(path+ifile)
	
	for jifile in ifile_actual:
		archive_files.append(jifile)
		data_files_actual.append(jifile.split('/')[-1])  


archive_files = list(flatten(archive_files))
data_files_actual = list(flatten(data_files_actual))

# TODO: replace '----->' with logging


# ^^^^^^^^^^ ANALYSIS
for ifile in range(len(archive_files)):
	output_folder = saved_folder+data_files_actual[ifile]+'_offPulseAnalysis'
	print ifile, '. ',data_files_actual[ifile]

	if doOffPulse:
		archive_analysis.archive_analysis(archive_files[ifile], output_folder, fphaseOffset = phaseOffset, ffreq_fract=freq_fract, fNbands_number = Nbands_number,fdoOnPulse = doOnPulse, fdoOffPulse = doOffPulse, fonBin_f=np.array([onBin_fmin,onBin_fmax]), foffBin_f=np.array([offBin_fmin,offBin_fmax]),small_fit = centralGf, small_fit_t = centralGt,ffinishEarly=finishEarly)
	else:
		archive_analysis.archive_analysis(archive_files[ifile], output_folder, fphaseOffset = phaseOffset, ffreq_fract=freq_fract, fNbands_number = Nbands_number,fdoOnPulse = doOnPulse, fdoOffPulse = doOffPulse, fonBin_f=np.array([onBin_fmin,onBin_fmax]),ffinishEarly=finishEarly)

		
		
# TODO: make a different folder inside the saved_folder for each set of obs
# TODO: save .npy data files as .txt as well





