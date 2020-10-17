from matplotlib.pylab import *
import matplotlib.gridspec as gridspec

def fancylayout(help=False, istex=True, tickright=True, ticktop=True, xtick=13, ytick=13, axeslw=1., \
				tickdirection='in', majorticklw=1., minorticklw=1, majortickll=6, minortickll=2):
	if help:
		print ('		istex\t\t\t\t = \tTrue/False')
		print ('		tickright\t\t\t = \tTrue/False')
		print ('		ticktop\t\t\t = \tTrue/False')
		print ('		xtick\t\t\t\t = \tlabel size on x/y ticks')
		print ('		axeslw\t\t\t\t = \taxes linewidth')
		print ('		tickdirection\t\t\t = \tin/out')
		print ('		majorticklw/minorticklw\t = \tmajor/minor ticks linewidth')
		print ('		majortickll/minortickll\t = \tmajor/minor ticks size')
		return 

	matplotlib.rc('xtick', labelsize=xtick) 
	matplotlib.rc('ytick', labelsize=ytick)
	 
	plt.rc('text', usetex=istex); plt.rc('font', family='serif')
	plt.rcParams['mathtext.default']='regular'

	for string in ['xtick.direction', 'ytick.direction']:
		plt.rcParams[string]	= tickdirection
	for string in ['ytick.right']:
		plt.rcParams[string]	= tickright
	for string in ['xtick.top']:
		plt.rcParams[string]	= ticktop
	for string in ['axes.linewidth']:
		plt.rcParams[string]	= axeslw
	for string in ['xtick.major.width', 'ytick.major.width']:
		plt.rcParams[string] 	= majorticklw
	for string in ['xtick.minor.width', 'ytick.minor.width']:
		plt.rcParams[string] 	= minorticklw
	for string in ['xtick.major.size', 'ytick.major.size']:
		plt.rcParams[string]	= majortickll
	for string in ['xtick.minor.size', 'ytick.minor.size']:
		plt.rcParams[string]	= minortickll		


def fmt(x, pos):                                                                                  
	a, b = '{:.2e}'.format(x).split('e')                                                                                          
	b = int(b)                                                                                                          
	return r'${} \times 10^{{{}}}$'.format(a, b) 

