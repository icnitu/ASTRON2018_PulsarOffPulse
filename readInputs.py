import os
import numpy as np



def read_good_input(message,list_good_inputs):

	bad_input = True
	input_value = 'Dummy'
	while bad_input:
		try:
			input_value = raw_input(message)
			index = list_good_inputs.index(input_value.lower())
			bad_input = False
		except ValueError:
			print "\n--- Answer not appropriate, please try again!\n"
	return input_value
	

def read_good_file(message):

	bad_input = True
	input_value = 'Dummy'
	while bad_input:
		try:
			archive_file = raw_input(message)
			file_try = open(archive_file,"r")
			file_try.close()
			bad_input = False
		except IOError:
			print "\n--- The file doesn't exist/cannot be opened, try again!\n"
	return archive_file

def read_good_folder(message):
	
	bad_input = True
	input_value = 'Dummy'
	while bad_input:
		input_value = raw_input(message)
		if os.path.isdir(input_value):
			bad_input = False
		else:
		
			print "\n--- The folder doesn't exist, creating the folder...\n"
			try:
				os.mkdir(input_value)
				bad_input = False
			except OSError:
				print "\n--- The folder path doesn't exist, try again!\n"
			
						
	input_value = str(input_value)
	if input_value[-1] != '/':
		input_value = input_value + '/'
		
	return input_value
	
			
def read_good_list_commas(message,list_good_inputs):
	
	bad_input = True
	input_values = []
	while bad_input:
		
		input_string = str(raw_input(message))
		# ignore any spaces around:
		input_values = [value.strip() for value in input_string.split(',')]
		N_values = len(input_values)
		count_good = 0	
		
		if (input_values[0] == 'none' or input_values[0] == 'all') and N_values > 1:
			# if 'none' or 'all', no other plots
			print "\n--- Answer not appropriate, please try again!\n"
		elif (input_values[0] == 'none' or input_values[0] == 'all') and N_values == 1:
			bad_input = False
		
		else:			

			for i in range(N_values):		
				try: 
					index = list_good_inputs.index(input_values[i].lower())
					count_good = count_good + 1
				except ValueError:
					print "\n--- Answer not appropriate, please try again!\n"
					break
	
		# if all values in list were good
		if count_good == N_values:
			bad_input = False
	
	return input_values
				
		
def read_good_binRange(message,binRange_min,binRange_max):
	
	bad_input = True
	input_values = []
	bin_min = np.nan
	bin_max = np.nan
	while bad_input:
		
		input_string = raw_input(message)
		try:
			input_values = map(int, input_string.split(','))
			
			if len(input_values) != 2:
				print "\n--- Answer not appropriate - exactly two numbers needed, please try again!\n"
			else:
				bin_min = min(input_values)
				bin_max = max(input_values)
				if bin_min < binRange_min:
					print "\n--- Answer not appropriate - out of range, please try again!\n"
				elif bin_max > binRange_max:
					print "\n--- Answer not appropriate - out of range, please try again!\n"
				elif bin_min == bin_max:
					print "\n--- Answer not appropriate - range 0, please try again!\n"
				else:
					bad_input = False
			
		except ValueError:
			print "\n--- Answer not appropriate, please try again!\n"
			
			
	return bin_min,bin_max
		
		
def read_good_float(message,min_value,max_value):
		
	bad_input = True
	input_value = 0.
	
	while bad_input:
		try:
			input_value = raw_input(message)
			input_float = float(input_value)
			if input_float >= min_value and input_float <=max_value:
				bad_input = False
			else:
				print "\n--- Answer not appropriate, outside range - please try again!\n"
				
		except ValueError:
			print "\n--- Answer not appropriate, please try again!\n"
	
	return input_float

def read_good_int(message,min_value,max_value):
		
	bad_input = True
	input_value = 0.
	
	while bad_input:
		try:
			input_value = raw_input(message)
			input_float = int(input_value)
			if input_float >= min_value and input_float <=max_value:
				bad_input = False
			else:
				print "\n--- Answer not appropriate, outside range - please try again!\n"
				
		except ValueError:
			print "\n--- Answer not appropriate, please try again!\n"
	
	return input_float


			
