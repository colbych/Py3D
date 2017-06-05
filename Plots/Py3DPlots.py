import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie



X = raw_input('Load Movie? \n ') 
if ((X == 'Y') or (X == 'y')):
	d = load_movie()
	print("Plotting...")
	
	Page1 = plt.figure(1)
	Page1.subplots_adjust(hspace = .5)
	Page1.set_size_inches(8.5,11)
	
	sp11 = Page1.add_subplot(321)
	plt.pcolormesh(d['bx'])
	plt.xticks([2000,8000])
	sp11.set_title("Bx Comp")
	sp11.set_aspect(1)
	
	sp12 = Page1.add_subplot(322)
	plt.pcolormesh(d['by'])
	plt.xticks([2000,8000])
	sp12.set_title("By Comp")
	sp12.set_aspect(1)
	
	sp13 = Page1.add_subplot(323)
	plt.pcolormesh(d['bz'])
	plt.xticks([2000,8000])
	sp13.set_title("Bz Comp")
	sp13.set_aspect(1)
	
	Mb  = np.sqrt(d['bx']**2 + d['by']**2 + d['bz']**2)
	sp14 = Page1.add_subplot(324)
	plt.pcolormesh(Mb)
	plt.xticks([2000,8000])
	sp14.set_title("|B|")
	sp14.set_aspect(1)
	
	sp15 = Page1.add_subplot(325)
	plt.pcolormesh(d['ni'])
	plt.xticks([2000,8000])
	sp15.set_title("NumDen ions")
	sp15.set_aspect(1)
	
	sp16 = Page1.add_subplot(326)
	plt.pcolormesh(d['ne'])
	plt.xticks([2000,8000])
	sp16.set_title("NumDen Elec")
	sp16.set_aspect(1)
	
	print("Saving...")
	
	Page1.show()
	
	# savefig takes FOREVER for a PDF
	# PNG is a little quicker but still slow
	# is there a better way?
	#Page1.savefig('Py3D_Page_1.png')
	
else:
	print("Load Data Cancelled")

