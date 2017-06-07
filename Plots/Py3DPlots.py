import numpy as np
import matplotlib.pyplot as sp
from Py3D.sub import load_movie

def TwoD_SubSample64(Par_Orig):
    print("Subsampling....")
    Par_Sub = np.empty([1024, 1024])
    for i in range (1,1024):
        for n in range(1,1024):
            Par_Sub[i,n] = Par_Orig[((i*8)-1),((n*8)-1)]
    return Par_Sub
    
def OneD_SubSample64(Par_Orig):
    print("Subsampling....")
    Par_Sub = np.empty([1024])
    for i in range (1,1024):
        Par_Sub[i] = Par_Orig[((i*8)-1)]
    return Par_Sub

def SubSampleDictionary(dict):
    New_dict = {'bx': TwoD_SubSample64(dict['bx']),'by': TwoD_SubSample64(dict['by']),'bz': TwoD_SubSample64(dict['bz']),
    'ex': TwoD_SubSample64(dict['ex']),'ey': TwoD_SubSample64(dict['ey']),'ez': TwoD_SubSample64(dict['ez']),
    'jx': TwoD_SubSample64(dict['jx']),'jy': TwoD_SubSample64(dict['jy']),'jz': TwoD_SubSample64(dict['jz']),
    'ni': TwoD_SubSample64(dict['ni']),'ne': TwoD_SubSample64(dict['ne']),'rho': TwoD_SubSample64(dict['rho']),
    'jex': TwoD_SubSample64(dict['jex']),'jey': TwoD_SubSample64(dict['jey']),'jez': TwoD_SubSample64(dict['jez']),
    'jix': TwoD_SubSample64(dict['jix']),'jiy': TwoD_SubSample64(dict['jiy']),'jiz': TwoD_SubSample64(dict['jiz']),
    'pexx': TwoD_SubSample64(dict['pexx']),'peyy': TwoD_SubSample64(dict['peyy']),'pezz': TwoD_SubSample64(dict['pezz']),
    'pixx': TwoD_SubSample64(dict['pixx']),'piyy': TwoD_SubSample64(dict['piyy']),'pizz': TwoD_SubSample64(dict['pizz']),
    'pexy': TwoD_SubSample64(dict['pexy']),'peyz': TwoD_SubSample64(dict['peyz']),'pexz': TwoD_SubSample64(dict['pexz']),
    'pixy': TwoD_SubSample64(dict['pixy']),'piyz': TwoD_SubSample64(dict['piyz']),'pixz': TwoD_SubSample64(dict['pixz']),
    'xx': OneD_SubSample64(dict['xx']),'yy': OneD_SubSample64(dict['yy'])}
    return New_dict


def P_1_6_Plots():
    Q = load_movie()
    d = SubSampleDictionary(Q)
    X = d['xx']
    Y = d['yy']
#-----------------------------------------------

    print("Plotting Page 1...")

    Page1 = sp.figure(1)
    Page1.subplots_adjust(hspace = .5)
    # resize not working???
    Page1.set_size_inches(8.5,11)

    sp11 = Page1.add_subplot(321)
    sp11.pcolormesh(X,Y,d['bx'])
    sp11.locator_params(nbins = 6)
    sp11.set_title("$B_x$")
    sp11.set_aspect(1)
    sp11.set_xlim([0,102.5])
    sp11.set_ylim([0,102.5])    
    
    sp12 = Page1.add_subplot(322)
    sp12.pcolormesh(X,Y,d['by'])
    sp12.locator_params(nbins = 6)
    sp12.set_title("$B_y$")
    sp12.set_aspect(1)
    sp12.set_xlim([0,102.5])
    sp12.set_ylim([0,102.5])

    sp13 = Page1.add_subplot(323)
    sp13.pcolormesh(X,Y,d['bz'])
    sp13.locator_params(nbins = 6)
    sp13.set_title("$B_z$")
    sp13.set_aspect(1)
    sp13.set_xlim([0,102.5])
    sp13.set_ylim([0,102.5])

    Mb  = np.sqrt(d['bx']**2 + d['by']**2 + d['bz']**2)
    sp14 = Page1.add_subplot(324)
    sp14.pcolormesh(X,Y,Mb)
    sp14.locator_params(nbins = 6)
    sp14.set_title("$\mid B\mid$")
    sp14.set_aspect(1)
    sp14.set_xlim([0,102.5])
    sp14.set_ylim([0,102.5])

    sp15 = Page1.add_subplot(325)
    sp15.pcolormesh(X,Y,d['ni'])
    sp15.locator_params(nbins = 6)
    sp15.set_title("$N_i$")
    sp15.set_aspect(1)
    sp15.set_xlim([0,102.5])
    sp15.set_ylim([0,102.5])

    sp16 = Page1.add_subplot(326)
    sp16.pcolormesh(X,Y,d['ne'])
    sp16.locator_params(nbins = 6)
    sp16.set_title("$N_e$")
    sp16.set_aspect(1)
    sp16.set_xlim([0,102.5])
    sp16.set_ylim([0,102.5])

    print("Saving Page 1...")
    
    Page1.show()

    # savefig takes FOREVER for a PDF
    # PNG is a little quicker but still slow
    # is there a better way?
    #Page1.savefig('Py3D_Page_1.png', dpi = 49)

#------------------------------------------------

    print("Plotting Page 2...")

    Page2 = sp.figure(2)
    Page2.subplots_adjust(hspace = .5)
    Page2.set_size_inches(8.5,11)

    sp21 = Page2.add_subplot(321)
    sp21.pcolormesh(X,Y,d['ex'])
    sp21.locator_params(nbins = 6)
    sp21.set_title("$E_x$")
    sp21.set_aspect(1)
    sp21.set_xlim([0,102.5])
    sp21.set_ylim([0,102.5])

    sp22 = Page2.add_subplot(322)
    sp22.pcolormesh(X,Y,d['ey'])
    sp22.locator_params(nbins = 6)
    sp22.set_title("$E_y$")
    sp22.set_aspect(1)
    sp22.set_xlim([0,102.5])
    sp22.set_ylim([0,102.5])

    sp23 = Page2.add_subplot(323)
    sp23.pcolormesh(X,Y,d['ez'])
    sp23.locator_params(nbins = 6)
    sp23.set_title("$E_z$")
    sp23.set_aspect(1)
    sp23.set_xlim([0,102.5])
    sp23.set_ylim([0,102.5])

    Me  = np.sqrt(d['ex']**2 + d['ey']**2 + d['ez']**2)
    sp24 = Page2.add_subplot(324)
    sp24.pcolormesh(X,Y,Me)
    sp24.locator_params(nbins = 6)
    sp24.set_title("$\mid E\mid$")
    sp24.set_aspect(1)
    sp24.set_xlim([0,102.5])
    sp24.set_ylim([0,102.5])

    EdB = d['ex']*d['bx'] + d['ey']*d['by'] + d['ez']*d['bz']
    sp25 = Page2.add_subplot(325)
    sp25.pcolormesh(X,Y,EdB)
    sp25.locator_params(nbins = 6)
    sp25.set_title("$E\cdot B$")
    sp25.set_aspect(1)
    sp25.set_xlim([0,102.5])
    sp25.set_ylim([0,102.5])

    sp26 = Page2.add_subplot(326)
    sp26.pcolormesh(X,Y,d['rho'])
    sp26.locator_params(nbins = 6)
    sp26.set_title("$P$")
    sp26.set_aspect(1)
    sp26.set_xlim([0,102.5])
    sp26.set_ylim([0,102.5])

    print("Saving Page 2...")

    Page2.show()

    #Page2.savefig('Py3D_Page_2.png', dpi = 49)

#----------------------------------------------

    print("Plotting Page 3...")

    Page3 = sp.figure(3)
    Page3.subplots_adjust(hspace = .5)
    Page3.set_size_inches(8.5,11)

    sp31 = Page3.add_subplot(321)
    sp31.pcolormesh(X,Y,d['jx'])
    sp31.locator_params(nbins = 6)
    sp31.set_title("$J_x$")
    sp31.set_aspect(1)
    sp31.set_xlim([0,102.5])
    sp31.set_ylim([0,102.5])

    sp32 = Page3.add_subplot(322)
    sp32.pcolormesh(X,Y,d['jy'])
    sp32.locator_params(nbins = 6)
    sp32.set_title("$J_y$")
    sp32.set_aspect(1)
    sp32.set_xlim([0,102.5])
    sp32.set_ylim([0,102.5])

    sp33 = Page3.add_subplot(323)
    sp33.pcolormesh(X,Y,d['jz'])
    sp33.locator_params(nbins = 6)
    sp33.set_title("$J_z$")
    sp33.set_aspect(1)
    sp33.set_xlim([0,102.5])
    sp33.set_ylim([0,102.5])

    sp34 = Page3.add_subplot(324)
    sp34.pcolormesh(X,Y,d['jix'])
    sp34.locator_params(nbins = 6)
    sp34.set_title("$J_{ix}$")
    sp34.set_aspect(1)
    sp34.set_xlim([0,102.5])
    sp34.set_ylim([0,102.5])

    sp35 = Page3.add_subplot(325)
    sp35.pcolormesh(X,Y,d['jiy'])
    sp35.locator_params(nbins = 6)
    sp35.set_title("$J_{iy}$")
    sp35.set_aspect(1)
    sp35.set_xlim([0,102.5])
    sp35.set_ylim([0,102.5])

    sp36 = Page3.add_subplot(326)
    sp36.pcolormesh(X,Y,d['jiz'])
    sp36.locator_params(nbins = 6)
    sp36.set_title("$J_{iz}$")
    sp36.set_aspect(1)
    sp36.set_xlim([0,102.5])
    sp36.set_ylim([0,102.5])

    print("Saving Page 3...")

    Page3.show()

    #Page3.savefig('Py3D_Page_3.png', dpi = 49)

#---------------------------------------------

    print("Plotting Page 4...")

    Page4 = sp.figure(4)
    Page4.subplots_adjust(hspace = .5)
    Page4.set_size_inches(8.5,11)

    sp41 = Page4.add_subplot(321)
    sp41.pcolormesh(X,Y,d['jex'])
    sp41.locator_params(nbins = 6)
    sp41.set_title("$J_{ex}$")
    sp41.set_aspect(1)
    sp41.set_xlim([0,102.5])
    sp41.set_ylim([0,102.5])

    sp42 = Page4.add_subplot(322)
    sp42.pcolormesh(X,Y,d['jey'])
    sp42.locator_params(nbins = 6)
    sp42.set_title("$J_{ey}$")
    sp42.set_aspect(1)
    sp42.set_xlim([0,102.5])
    sp42.set_ylim([0,102.5])

    sp43 = Page4.add_subplot(323)
    sp43.pcolormesh(X,Y,d['jez'])
    sp43.locator_params(nbins = 6)
    sp43.set_title("$J_{ez}$")
    sp43.set_aspect(1)
    sp43.set_xlim([0,102.5])
    sp43.set_ylim([0,102.5])

    Vex = -d['jex']/d['ne']
    sp44 = Page4.add_subplot(324)
    sp44.pcolormesh(X,Y,Vex)
    sp44.locator_params(nbins = 6)
    sp44.set_title("$V_{ex}$")
    sp44.set_aspect(1)
    sp44.set_xlim([0,102.5])
    sp44.set_ylim([0,102.5])

    Vey = -d['jey']/d['ne']
    sp45 = Page4.add_subplot(325)
    sp45.pcolormesh(X,Y,Vey)
    sp45.locator_params(nbins = 6)
    sp45.set_title("$V_{ey}$")
    sp45.set_aspect(1)
    sp45.set_xlim([0,102.5])
    sp45.set_ylim([0,102.5])

    Vez = -d['jez']/d['ne']
    sp46 = Page4.add_subplot(326)
    sp46.pcolormesh(X,Y,Vez)
    sp46.locator_params(nbins = 6)
    sp46.set_title("$V_{ez}$")
    sp46.set_aspect(1)
    sp46.set_xlim([0,102.5])
    sp46.set_ylim([0,102.5])

    print("Saving Page 4...")

    Page4.show()

    #Page4.savefig('Py3D_Page_4.png', dpi = 49)

#-------------------------------------------

    print("Plotting Page 5...")

    Page5 = sp.figure(5)
    Page5.subplots_adjust(hspace = .5)
    Page5.set_size_inches(8.5,11)

    Vix = d['jix']/d['ni']
    sp51 = Page5.add_subplot(321)
    sp51.pcolormesh(X,Y,Vix)
    sp51.locator_params(nbins = 6)
    sp51.set_title("$V_{ix}$")
    sp51.set_aspect(1)
    sp51.set_xlim([0,102.5])
    sp51.set_ylim([0,102.5])
    
    Viy = d['jiy']/d['ni']
    sp52 = Page5.add_subplot(322)
    sp52.pcolormesh(X,Y,Viy)
    sp52.locator_params(nbins = 6)
    sp52.set_title("$V_{iy}$")
    sp52.set_aspect(1)
    sp52.set_xlim([0,102.5])
    sp52.set_ylim([0,102.5])

    Viz = d['jiz']/d['ni']
    sp53 = Page5.add_subplot(323)
    sp53.pcolormesh(X,Y,Viz)
    sp53.locator_params(nbins = 6)
    sp53.set_title("$V_{iz}$")
    sp53.set_aspect(1)
    sp53.set_xlim([0,102.5])
    sp53.set_ylim([0,102.5])

    Tixx = d['pixx']/d['ni']
    sp54 = Page5.add_subplot(324)
    sp54.pcolormesh(X,Y,Tixx)
    sp54.locator_params(nbins = 6)
    sp54.set_title("$T_{ixx}$")
    sp54.set_aspect(1)
    sp54.set_xlim([0,102.5])
    sp54.set_ylim([0,102.5])

    Tiyy = d['piyy']/d['ni']
    sp55 = Page5.add_subplot(325)
    sp55.pcolormesh(X,Y,Tiyy)
    sp55.locator_params(nbins = 6)
    sp55.set_title("$T_{iyy}$")
    sp55.set_aspect(1)
    sp55.set_xlim([0,102.5])
    sp55.set_ylim([0,102.5])

    Tizz = d['pizz']/d['ni']
    sp56 = Page5.add_subplot(326)
    sp56.pcolormesh(X,Y,Tizz)
    sp56.locator_params(nbins = 6)
    sp56.set_title("$T_{izz}$")
    sp56.set_aspect(1)
    sp56.set_xlim([0,102.5])
    sp56.set_ylim([0,102.5])

    print("Saving Page 5...")

    Page5.show()

    #Page5.savefig('Py3D_Page_5.png', dpi = 49)

#----------------------------------------------

    print("Plotting Page 6...")

    Page6 = sp.figure(6)
    Page6.subplots_adjust(hspace = .5)
    Page6.set_size_inches(8.5,11)

    Texx = d['pexx']/d['ne']
    sp61 = Page6.add_subplot(321)
    sp61.pcolormesh(X,Y,Texx)
    sp61.locator_params(nbins = 6)
    sp61.set_title("$T_{exx}$")
    sp61.set_aspect(1)
    sp61.set_xlim([0,102.5])
    sp61.set_ylim([0,102.5])

    Teyy = d['peyy']/d['ne']
    sp62 = Page6.add_subplot(322)
    sp62.pcolormesh(X,Y,Teyy)
    sp62.locator_params(nbins = 6)
    sp62.set_title("$T_{eyy}$")
    sp62.set_aspect(1)
    sp62.set_xlim([0,102.5])
    sp62.set_ylim([0,102.5])

    Tezz = d['pezz']/d['ne']
    sp63 = Page6.add_subplot(323)
    sp63.pcolormesh(X,Y,Tezz)
    sp63.locator_params(nbins = 6)
    sp63.set_title("$T_{ezz}$")
    sp63.set_aspect(1)
    sp63.set_xlim([0,102.5])
    sp63.set_ylim([0,102.5])

    Texy = d['pexy']/d['ne']
    sp64 = Page6.add_subplot(324)
    sp64.pcolormesh(X,Y,Texy)
    sp64.locator_params(nbins = 6)
    sp64.set_title("$T_{exy}$")
    sp64.set_aspect(1)
    sp64.set_xlim([0,102.5])
    sp64.set_ylim([0,102.5])

    Texz = d['pexz']/d['ne']
    sp65 = Page6.add_subplot(325)
    sp65.pcolormesh(X,Y,Texz)
    sp65.locator_params(nbins = 6)
    sp65.set_title("$T_{exz}$")
    sp65.set_aspect(1)
    sp65.set_xlim([0,102.5])
    sp65.set_ylim([0,102.5])

    Teyz = d['peyz']/d['ne']
    sp66 = Page6.add_subplot(326)
    sp66.pcolormesh(X,Y,Teyz)
    sp66.locator_params(nbins = 6)
    sp66.set_title("$T_{eyz}$")
    sp66.set_aspect(1)
    sp66.set_xlim([0,102.5])
    sp66.set_ylim([0,102.5])

    print("Saving Page 6...")

    Page6.show()
#
#    #Page6.savefig('Py3D_Page_6.png', dpi = 49)

#--------------------------------------------

i = 0
while(i == 0):
    while True:
        try:
            X = str(raw_input('Load Movie? Y or N \n '))
            break
        except ValueError:
            print("invalid input try again...\n")
    if ((X == 'Y') or (X == 'y')):
        i = 1        
        P_1_6_Plots()
    elif((X == 'N') or (X == 'n')):
        i = 1
        print("Load Data Cancelled")
    else:
        print("invalid input try again...\n")
