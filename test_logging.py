import logging


def fract(x,y):
	logging.info('Calculating fraction of %4.2f and %4.2f', x, y)
	return float(x)/float(y)


logging.basicConfig(format = '%(asctime)s %(levelname)s:%(message)s', datefmt='%H:%M:%S', level = logging.INFO)

logging.debug('This is a debug message.')
#logging.warning('This is a warning.')
logging.info('START')



logging.info('FINISH')


logging.captureWarnings(True)
formatter = logging.Formatter('%(asctime)s %(levelname)s:%(message)s', datefmt='%H:%M:%S')
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(formatter)
warnings_logger = logging.getLogger()
warnings_logger.addHandler(console_handler)

print 'fraction = ',fract(3,0)

