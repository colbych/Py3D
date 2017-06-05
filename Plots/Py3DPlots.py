import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie


def P_1_6_Plots():
	d = load_movie()

#-----------------------------------------------

	print("Plotting Page 1...")
	
	Page1 = plt.figure(1)
	Page1.subplots_adjust(hspace = .5)
	# resize not working???
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
	
	print("Saving Page 1...")
	
	Page1.show()
	
	# savefig takes FOREVER for a PDF
	# PNG is a little quicker but still slow
	# is there a better way?
	#Page1.savefig('Py3D_Page_1.png')

#------------------------------------------------

        print("Plotting Page 2...")

        Page2 = plt.figure(2)
        Page2.subplots_adjust(hspace = .5)
        Page2.set_size_inches(8.5,11)

        sp21 = Page2.add_subplot(321)
        plt.pcolormesh(d['ex'])
        plt.xticks([2000,8000])
        sp21.set_title("Ex Comp")
        sp21.set_aspect(1)

        sp22 = Page2.add_subplot(322)
        plt.pcolormesh(d['ey'])
        plt.xticks([2000,8000])
        sp22.set_title("Ey Comp")
        sp22.set_aspect(1)

        sp23 = Page2.add_subplot(323)
        plt.pcolormesh(d['ez'])
        plt.xticks([2000,8000])
        sp23.set_title("Ez Comp")
        sp23.set_aspect(1)

        Me  = np.sqrt(d['ex']**2 + d['ey']**2 + d['ez']**2)
        sp24 = Page2.add_subplot(324)
        plt.pcolormesh(Me)
        plt.xticks([2000,8000])
        sp24.set_title("|E|")
        sp24.set_aspect(1)

        EdB = d['ex']*d['bx'] + d['ey']*d['by'] + d['ez']*d['bz']
        sp25 = Page2.add_subplot(325)
        plt.pcolormesh(EdB)
        plt.xticks([2000,8000])
        sp25.set_title("E dot B")
        sp25.set_aspect(1)

        sp26 = Page2.add_subplot(326)
        plt.pcolormesh(d['rho'])
        plt.xticks([2000,8000])
        sp26.set_title("Rho")
        sp26.set_aspect(1)

        print("Saving Page 2...")

        Page2.show()

        #Page2.savefig('Py3D_Page_2.png')

#----------------------------------------------

        print("Plotting Page 3...")

        Page3 = plt.figure(3)
        Page3.subplots_adjust(hspace = .5)
        Page3.set_size_inches(8.5,11)

        sp31 = Page3.add_subplot(321)
        plt.pcolormesh(d['jx'])
        plt.xticks([2000,8000])
        sp31.set_title("Jx Comp")
        sp31.set_aspect(1)

        sp32 = Page3.add_subplot(322)
        plt.pcolormesh(d['jy'])
        plt.xticks([2000,8000])
        sp32.set_title("Jy Comp")
        sp32.set_aspect(1)

        sp33 = Page3.add_subplot(323)
        plt.pcolormesh(d['jz'])
        plt.xticks([2000,8000])
        sp33.set_title("Jz Comp")
        sp33.set_aspect(1)

        sp34 = Page3.add_subplot(324)
        plt.pcolormesh(d['jix'])
        plt.xticks([2000,8000])
        sp34.set_title("Jix Comp")
        sp34.set_aspect(1)

        sp35 = Page3.add_subplot(325)
        plt.pcolormesh(d['jiy'])
        plt.xticks([2000,8000])
        sp35.set_title("Jiy Comp")
        sp35.set_aspect(1)

        sp36 = Page3.add_subplot(326)
        plt.pcolormesh(d['jiz'])
        plt.xticks([2000,8000])
        sp36.set_title("Jiz Comp")
        sp36.set_aspect(1)

        print("Saving Page 3...")

        Page3.show()

        #Page3.savefig('Py3D_Page_3.png')

#---------------------------------------------

        print("Plotting Page 4...")

        Page4 = plt.figure(4)
        Page4.subplots_adjust(hspace = .5)
        Page4.set_size_inches(8.5,11)

        sp41 = Page4.add_subplot(321)
        plt.pcolormesh(d['jex'])
        plt.xticks([2000,8000])
        sp41.set_title("Jex Comp")
        sp41.set_aspect(1)

        sp42 = Page4.add_subplot(322)
        plt.pcolormesh(d['jey'])
        plt.xticks([2000,8000])
        sp42.set_title("Jey Comp")
        sp42.set_aspect(1)

        sp43 = Page4.add_subplot(323)
        plt.pcolormesh(d['jez'])
        plt.xticks([2000,8000])
        sp43.set_title("Jez Comp")
        sp43.set_aspect(1)

        Vex = -d['jex']/d['ne']
        sp44 = Page4.add_subplot(324)
        plt.pcolormesh(Vex)
        plt.xticks([2000,8000])
        sp44.set_title("Vex Comp")
        sp44.set_aspect(1)

        Vey = -d['jey']/d['ne']
        sp45 = Page4.add_subplot(325)
        plt.pcolormesh(Vey)
        plt.xticks([2000,8000])
        sp45.set_title("Vey Comp")
        sp45.set_aspect(1)

        Vez = -d['jez']/d['ne']
        sp46 = Page4.add_subplot(326)
        plt.pcolormesh(Vez)
        plt.xticks([2000,8000])
        sp46.set_title("Vez Comp")
        sp46.set_aspect(1)

        print("Saving Page 4...")

        Page4.show()

        #Page4.savefig('Py3D_Page_4.png')

#-------------------------------------------

        print("Plotting Page 5...")

        Page5 = plt.figure(5)
        Page5.subplots_adjust(hspace = .5)
        Page5.set_size_inches(8.5,11)

        Vix = d['jix']/d['ni']
        sp51 = Page5.add_subplot(321)
        plt.pcolormesh(Vix)
        plt.xticks([2000,8000])
        sp51.set_title("Vix Comp")
        sp51.set_aspect(1)

        Viy = d['jiy']/d['ni']
        sp52 = Page5.add_subplot(322)
        plt.pcolormesh(Viy)
        plt.xticks([2000,8000])
        sp52.set_title("Viy Comp")
        sp52.set_aspect(1)

        Viz = d['jiz']/d['ni']
        sp53 = Page5.add_subplot(323)
        plt.pcolormesh(Viz)
        plt.xticks([2000,8000])
        sp53.set_title("Viz Comp")
        sp53.set_aspect(1)

        Tixx = d['pixx']/d['ni']
        sp54 = Page5.add_subplot(324)
        plt.pcolormesh(Tixx)
        plt.xticks([2000,8000])
        sp54.set_title("Tixx Comp")
        sp54.set_aspect(1)

        Tiyy = d['piyy']/d['ni']
        sp55 = Page5.add_subplot(325)
        plt.pcolormesh(Tiyy)
        plt.xticks([2000,8000])
        sp55.set_title("Tiyy Comp")
        sp55.set_aspect(1)

        Tizz = d['pizz']/d['ni']
        sp56 = Page5.add_subplot(326)
        plt.pcolormesh(Tizz)
        plt.xticks([2000,8000])
        sp56.set_title("Tizz Comp")
        sp56.set_aspect(1)

        print("Saving Page 5...")

        Page5.show()

        #Page5.savefig('Py3D_Page_5.png')

#----------------------------------------------

        print("Plotting Page 6...")

        Page6 = plt.figure(6)
        Page6.subplots_adjust(hspace = .5)
        Page6.set_size_inches(8.5,11)

        Texx = d['pexx']/d['ne']
        sp61 = Page6.add_subplot(321)
        plt.pcolormesh(Texx)
        plt.xticks([2000,8000])
        sp61.set_title("Texx Comp")
        sp61.set_aspect(1)

        Teyy = d['peyy']/d['ne']
        sp62 = Page6.add_subplot(322)
        plt.pcolormesh(Teyy)
        plt.xticks([2000,8000])
        sp62.set_title("Teyy Comp")
        sp62.set_aspect(1)

        Tezz = d['pezz']/d['ne']
        sp63 = Page6.add_subplot(323)
        plt.pcolormesh(Tezz)
        plt.xticks([2000,8000])
        sp63.set_title("Tezz Comp")
        sp63.set_aspect(1)

        Texy = d['pexy']/d['ne']
        sp64 = Page6.add_subplot(324)
        plt.pcolormesh(Texy)
        plt.xticks([2000,8000])
        sp64.set_title("Texy Comp")
        sp64.set_aspect(1)

        Texz = d['pexz']/d['ne']
        sp65 = Page6.add_subplot(325)
        plt.pcolormesh(Texz)
        plt.xticks([2000,8000])
        sp65.set_title("Texz Comp")
        sp65.set_aspect(1)

        Teyz = d['peyz']/d['ne']
        sp66 = Page6.add_subplot(326)
        plt.pcolormesh(Teyz)
        plt.xticks([2000,8000])
        sp66.set_title("Teyz Comp")
        sp66.set_aspect(1)

        print("Saving Page 6...")

        Page6.show()

        #Page6.savefig('Py3D_Page_6.png')

#--------------------------------------------

i = 0
while(i == 0):
    X = raw_input('Load Movie? Y or N \n ')
    if ((X == 'Y') or (X == 'y')):
        i = 1
        P_1_6_Plots()
    elif((X == 'N') or (X == 'n')):
        i = 1
        print("Load Data Cancelled")
    else:
        print("invalid input \n")

