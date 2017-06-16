import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie

# sub-sampling of 1D or 2D data, rate being number of points to skip
def SubSample(data, rate):
     if len(data.shape) > 1:
         return data[::rate, ::rate]
     else:
         return data[::rate]
    
# sub-sampling of entire dictionary
def SS_Dict(dictionary):
    return {k: SubSample(dictionary[k], 6) for k in dictionary}

# primary plotting function
def Py3D_Plots():
    # loads initial dictionary
    Q = load_movie()
    print("Sub-sampling....")
    # sub-sample dictionary
    d = SS_Dict(Q)
    # set X and Y for clarity in plotting functions
    X = d['xx']
    Y = d['yy']
#-----------------------------------------------

    print("Plotting Page 1...")
    
    # 6 Plots per page
    Page1 = plt.figure(1)
    # formatting page
    Page1.set_size_inches(8.5,11, forward = True)
    Page1.subplots_adjust(hspace = .5)
    
    # first plot
    sp11 = Page1.add_subplot(321)
    sp11.pcolormesh(X,Y,d['bx'])
    # formatting subplot
    sp11.locator_params(nbins = 6)
    sp11.set_title("$B_x$", fontsize=20)
    sp11.set_aspect(1)
    sp11.set_xlim([0,102.5])
    sp11.set_ylim([0,102.5])    
    
    sp12 = Page1.add_subplot(322)
    sp12.pcolormesh(X,Y,d['by'])
    sp12.locator_params(nbins = 6)
    sp12.set_title("$B_y$", fontsize=20)
    sp12.set_aspect(1)
    sp12.set_xlim([0,102.5])
    sp12.set_ylim([0,102.5])

    sp13 = Page1.add_subplot(323)
    sp13.pcolormesh(X,Y,d['bz'])
    sp13.locator_params(nbins = 6)
    sp13.set_title("$B_z$", fontsize=20)
    sp13.set_aspect(1)
    sp13.set_xlim([0,102.5])
    sp13.set_ylim([0,102.5])
    
    # calculating magnitude of Mag field
    Mb  = np.sqrt(d['bx']**2 + d['by']**2 + d['bz']**2)
    sp14 = Page1.add_subplot(324)
    sp14.pcolormesh(X,Y,Mb)
    sp14.locator_params(nbins = 6)
    sp14.set_title("$\mid B\mid$", fontsize=20)
    sp14.set_aspect(1)
    sp14.set_xlim([0,102.5])
    sp14.set_ylim([0,102.5])

    sp15 = Page1.add_subplot(325)
    sp15.pcolormesh(X,Y,d['ni'])
    sp15.locator_params(nbins = 6)
    sp15.set_title("$N_i$", fontsize=20)
    sp15.set_aspect(1)
    sp15.set_xlim([0,102.5])
    sp15.set_ylim([0,102.5])

    sp16 = Page1.add_subplot(326)
    sp16.pcolormesh(X,Y,d['ne'])
    sp16.locator_params(nbins = 6)
    sp16.set_title("$N_e$", fontsize=20)
    sp16.set_aspect(1)
    sp16.set_xlim([0,102.5])
    sp16.set_ylim([0,102.5])
    
    Page1.show()
    
    print("Saving Page 1...")

    # savefig takes FOREVER for a PDF
    # PNG is a little quicker but still slow
    # is there a better way?
    #Page1.savefig('Py3D_Page_1.png', dpi = 49)

#------------------------------------------------

    print("Plotting Page 2...")

    Page2 = plt.figure(2)
    Page2.subplots_adjust(hspace = .5)
    Page2.set_size_inches(8.5,11, forward = True)

    sp21 = Page2.add_subplot(321)
    sp21.pcolormesh(X,Y,d['ex'])
    sp21.locator_params(nbins = 6)
    sp21.set_title("$E_x$", fontsize=20)
    sp21.set_aspect(1)
    sp21.set_xlim([0,102.5])
    sp21.set_ylim([0,102.5])

    sp22 = Page2.add_subplot(322)
    sp22.pcolormesh(X,Y,d['ey'])
    sp22.locator_params(nbins = 6)
    sp22.set_title("$E_y$", fontsize=20)
    sp22.set_aspect(1)
    sp22.set_xlim([0,102.5])
    sp22.set_ylim([0,102.5])

    sp23 = Page2.add_subplot(323)
    sp23.pcolormesh(X,Y,d['ez'])
    sp23.locator_params(nbins = 6)
    sp23.set_title("$E_z$", fontsize=20)
    sp23.set_aspect(1)
    sp23.set_xlim([0,102.5])
    sp23.set_ylim([0,102.5])

    # calculating magnitude of Elec Field
    Me  = np.sqrt(d['ex']**2 + d['ey']**2 + d['ez']**2)
    sp24 = Page2.add_subplot(324)
    sp24.pcolormesh(X,Y,Me)
    sp24.locator_params(nbins = 6)
    sp24.set_title("$\mid E\mid$", fontsize=20)
    sp24.set_aspect(1)
    sp24.set_xlim([0,102.5])
    sp24.set_ylim([0,102.5])
    
    # calculating magnitude of E dot B
    EdB = d['ex']*d['bx'] + d['ey']*d['by'] + d['ez']*d['bz']
    sp25 = Page2.add_subplot(325)
    sp25.pcolormesh(X,Y,EdB)
    sp25.locator_params(nbins = 6)
    sp25.set_title("$E\cdot B$", fontsize=20)
    sp25.set_aspect(1)
    sp25.set_xlim([0,102.5])
    sp25.set_ylim([0,102.5])

    sp26 = Page2.add_subplot(326)
    sp26.pcolormesh(X,Y,d['rho'])
    sp26.locator_params(nbins = 6)
    sp26.set_title("$P$", fontsize=20)
    sp26.set_aspect(1)
    sp26.set_xlim([0,102.5])
    sp26.set_ylim([0,102.5])

    Page2.show()
    
    print("Saving Page 2...")

    #Page2.savefig('Py3D_Page_2.png', dpi = 49)

#----------------------------------------------

    print("Plotting Page 3...")

    Page3 = plt.figure(3)
    Page3.subplots_adjust(hspace = .5)
    Page3.set_size_inches(8.5,11, forward = True)

    sp31 = Page3.add_subplot(321)
    sp31.pcolormesh(X,Y,d['jx'])
    sp31.locator_params(nbins = 6)
    sp31.set_title("$J_x$", fontsize=20)
    sp31.set_aspect(1)
    sp31.set_xlim([0,102.5])
    sp31.set_ylim([0,102.5])

    sp32 = Page3.add_subplot(322)
    sp32.pcolormesh(X,Y,d['jy'])
    sp32.locator_params(nbins = 6)
    sp32.set_title("$J_y$", fontsize=20)
    sp32.set_aspect(1)
    sp32.set_xlim([0,102.5])
    sp32.set_ylim([0,102.5])

    sp33 = Page3.add_subplot(323)
    sp33.pcolormesh(X,Y,d['jz'])
    sp33.locator_params(nbins = 6)
    sp33.set_title("$J_z$", fontsize=20)
    sp33.set_aspect(1)
    sp33.set_xlim([0,102.5])
    sp33.set_ylim([0,102.5])

    sp34 = Page3.add_subplot(324)
    sp34.pcolormesh(X,Y,d['jix'])
    sp34.locator_params(nbins = 6)
    sp34.set_title("$J_{ix}$", fontsize=20)
    sp34.set_aspect(1)
    sp34.set_xlim([0,102.5])
    sp34.set_ylim([0,102.5])

    sp35 = Page3.add_subplot(325)
    sp35.pcolormesh(X,Y,d['jiy'])
    sp35.locator_params(nbins = 6)
    sp35.set_title("$J_{iy}$", fontsize=20)
    sp35.set_aspect(1)
    sp35.set_xlim([0,102.5])
    sp35.set_ylim([0,102.5])

    sp36 = Page3.add_subplot(326)
    sp36.pcolormesh(X,Y,d['jiz'])
    sp36.locator_params(nbins = 6)
    sp36.set_title("$J_{iz}$", fontsize=20)
    sp36.set_aspect(1)
    sp36.set_xlim([0,102.5])
    sp36.set_ylim([0,102.5])

    Page3.show()
    
    print("Saving Page 3...")

    #Page3.savefig('Py3D_Page_3.png', dpi = 49)

#---------------------------------------------

    print("Plotting Page 4...")

    Page4 = plt.figure(4)
    Page4.subplots_adjust(hspace = .5)
    Page4.set_size_inches(8.5,11, forward = True)

    sp41 = Page4.add_subplot(321)
    sp41.pcolormesh(X,Y,d['jex'])
    sp41.locator_params(nbins = 6)
    sp41.set_title("$J_{ex}$", fontsize=20)
    sp41.set_aspect(1)
    sp41.set_xlim([0,102.5])
    sp41.set_ylim([0,102.5])

    sp42 = Page4.add_subplot(322)
    sp42.pcolormesh(X,Y,d['jey'])
    sp42.locator_params(nbins = 6)
    sp42.set_title("$J_{ey}$", fontsize=20)
    sp42.set_aspect(1)
    sp42.set_xlim([0,102.5])
    sp42.set_ylim([0,102.5])

    sp43 = Page4.add_subplot(323)
    sp43.pcolormesh(X,Y,d['jez'])
    sp43.locator_params(nbins = 6)
    sp43.set_title("$J_{ez}$", fontsize=20)
    sp43.set_aspect(1)
    sp43.set_xlim([0,102.5])
    sp43.set_ylim([0,102.5])

    # calculating velocities
    Vex = -d['jex']/d['ne']
    sp44 = Page4.add_subplot(324)
    sp44.pcolormesh(X,Y,Vex)
    sp44.locator_params(nbins = 6)
    sp44.set_title("$V_{ex}$", fontsize=20)
    sp44.set_aspect(1)
    sp44.set_xlim([0,102.5])
    sp44.set_ylim([0,102.5])

    Vey = -d['jey']/d['ne']
    sp45 = Page4.add_subplot(325)
    sp45.pcolormesh(X,Y,Vey)
    sp45.locator_params(nbins = 6)
    sp45.set_title("$V_{ey}$", fontsize=20)
    sp45.set_aspect(1)
    sp45.set_xlim([0,102.5])
    sp45.set_ylim([0,102.5])

    Vez = -d['jez']/d['ne']
    sp46 = Page4.add_subplot(326)
    sp46.pcolormesh(X,Y,Vez)
    sp46.locator_params(nbins = 6)
    sp46.set_title("$V_{ez}$", fontsize=20)
    sp46.set_aspect(1)
    sp46.set_xlim([0,102.5])
    sp46.set_ylim([0,102.5])

    Page4.show()
    
    print("Saving Page 4...")

    #Page4.savefig('Py3D_Page_4.png', dpi = 49)

#-------------------------------------------

    print("Plotting Page 5...")

    Page5 = plt.figure(5)
    Page5.subplots_adjust(hspace = .5)
    Page5.set_size_inches(8.5,11, forward = True)

    Vix = d['jix']/d['ni']
    sp51 = Page5.add_subplot(321)
    sp51.pcolormesh(X,Y,Vix)
    sp51.locator_params(nbins = 6)
    sp51.set_title("$V_{ix}$", fontsize=20)
    sp51.set_aspect(1)
    sp51.set_xlim([0,102.5])
    sp51.set_ylim([0,102.5])
    
    Viy = d['jiy']/d['ni']
    sp52 = Page5.add_subplot(322)
    sp52.pcolormesh(X,Y,Viy)
    sp52.locator_params(nbins = 6)
    sp52.set_title("$V_{iy}$", fontsize=20)
    sp52.set_aspect(1)
    sp52.set_xlim([0,102.5])
    sp52.set_ylim([0,102.5])

    Viz = d['jiz']/d['ni']
    sp53 = Page5.add_subplot(323)
    sp53.pcolormesh(X,Y,Viz)
    sp53.locator_params(nbins = 6)
    sp53.set_title("$V_{iz}$", fontsize=20)
    sp53.set_aspect(1)
    sp53.set_xlim([0,102.5])
    sp53.set_ylim([0,102.5])

    # calculating temperatures
    Tixx = d['pixx']/d['ni']
    sp54 = Page5.add_subplot(324)
    sp54.pcolormesh(X,Y,Tixx)
    sp54.locator_params(nbins = 6)
    sp54.set_title("$T_{ixx}$", fontsize=20)
    sp54.set_aspect(1)
    sp54.set_xlim([0,102.5])
    sp54.set_ylim([0,102.5])

    Tiyy = d['piyy']/d['ni']
    sp55 = Page5.add_subplot(325)
    sp55.pcolormesh(X,Y,Tiyy)
    sp55.locator_params(nbins = 6)
    sp55.set_title("$T_{iyy}$", fontsize=20)
    sp55.set_aspect(1)
    sp55.set_xlim([0,102.5])
    sp55.set_ylim([0,102.5])

    Tizz = d['pizz']/d['ni']
    sp56 = Page5.add_subplot(326)
    sp56.pcolormesh(X,Y,Tizz)
    sp56.locator_params(nbins = 6)
    sp56.set_title("$T_{izz}$", fontsize=20)
    sp56.set_aspect(1)
    sp56.set_xlim([0,102.5])
    sp56.set_ylim([0,102.5])

    Page5.show()
    
    print("Saving Page 5...")

    #Page5.savefig('Py3D_Page_5.png', dpi = 49)

#----------------------------------------------

    print("Plotting Page 6...")

    Page6 = plt.figure(6)
    Page6.subplots_adjust(hspace = .5)
    Page6.set_size_inches(8.5,11, forward = True)

    Texx = d['pexx']/d['ne']
    sp61 = Page6.add_subplot(321)
    sp61.pcolormesh(X,Y,Texx)
    sp61.locator_params(nbins = 6)
    sp61.set_title("$T_{exx}$", fontsize=20)
    sp61.set_aspect(1)
    sp61.set_xlim([0,102.5])
    sp61.set_ylim([0,102.5])

    Teyy = d['peyy']/d['ne']
    sp62 = Page6.add_subplot(322)
    sp62.pcolormesh(X,Y,Teyy)
    sp62.locator_params(nbins = 6)
    sp62.set_title("$T_{eyy}$", fontsize=20)
    sp62.set_aspect(1)
    sp62.set_xlim([0,102.5])
    sp62.set_ylim([0,102.5])

    Tezz = d['pezz']/d['ne']
    sp63 = Page6.add_subplot(323)
    sp63.pcolormesh(X,Y,Tezz)
    sp63.locator_params(nbins = 6)
    sp63.set_title("$T_{ezz}$", fontsize=20)
    sp63.set_aspect(1)
    sp63.set_xlim([0,102.5])
    sp63.set_ylim([0,102.5])

    Texy = d['pexy']/d['ne']
    sp64 = Page6.add_subplot(324)
    sp64.pcolormesh(X,Y,Texy)
    sp64.locator_params(nbins = 6)
    sp64.set_title("$T_{exy}$", fontsize=20)
    sp64.set_aspect(1)
    sp64.set_xlim([0,102.5])
    sp64.set_ylim([0,102.5])

    Texz = d['pexz']/d['ne']
    sp65 = Page6.add_subplot(325)
    sp65.pcolormesh(X,Y,Texz)
    sp65.locator_params(nbins = 6)
    sp65.set_title("$T_{exz}$", fontsize=20)
    sp65.set_aspect(1)
    sp65.set_xlim([0,102.5])
    sp65.set_ylim([0,102.5])

    Teyz = d['peyz']/d['ne']
    sp66 = Page6.add_subplot(326)
    sp66.pcolormesh(X,Y,Teyz)
    sp66.locator_params(nbins = 6)
    sp66.set_title("$T_{eyz}$", fontsize=20)
    sp66.set_aspect(1)
    sp66.set_xlim([0,102.5])
    sp66.set_ylim([0,102.5])

    Page6.show()
    
    print("Saving Page 6...")
#
#    #Page6.savefig('Py3D_Page_6.png', dpi = 49)

#--------------------------------------------

# user interface
i = 0
while(i == 0):
    # catch for non-string entry
    while True:
        try:
            X = str(raw_input('Load Movie? Y or N \n '))
            break
        except ValueError:
            print("invalid input try again...\n")
    # handling string entries
    if ((X == 'Y') or (X == 'y')):
        i = 1        
        Py3D_Plots()
    elif((X == 'N') or (X == 'n')):
        i = 1
        print("Load Data Cancelled")
    else:
        print("invalid input try again...\n")
        
        
        
