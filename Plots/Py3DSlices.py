import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie

# primary slicing and plotting function for constant X
def sliceX(index):
    # loading data
    d = load_movie()
    # this is a slice of constant X so plot against Y
    Y = d['yy']
    print("slicing...")
    
    # slicing
    Bx_Xslice = d['bx'][:,index]
    By_Xslice = d['by'][:,index]
    Bz_Xslice = d['bz'][:,index]
    # calculating magnitude of Mag field
    Mb_Xslice = np.sqrt((d['bx'][:,index])**2 + (d['by'][:,index])**2 + (d['bz'][:,index])**2)
    
    Ex_Xslice = d['ex'][:,index]
    Ey_Xslice = d['ey'][:,index]
    Ez_Xslice = d['ez'][:,index]
    # calculating magnitude of Elec field
    Me_Xslice = np.sqrt((d['ex'][:,index])**2 + (d['ey'][:,index])**2 + (d['ez'][:,index])**2)    
    
    Ni_Xslice = d['ni'][:,index]
    Ne_Xslice = d['ne'][:,index]

    Jx_Xslice = d['jx'][:,index]
    Jy_Xslice = d['jy'][:,index]
    Jz_Xslice = d['jz'][:,index]
    
    # calculating velocities
    Vix_Xslice = (d['jix'][:,index])/(d['ni'][:,index])
    Viy_Xslice = (d['jiy'][:,index])/(d['ni'][:,index])
    Viz_Xslice = (d['jiz'][:,index])/(d['ni'][:,index])

    Vex_Xslice = (-d['jex'][:,index])/(d['ne'][:,index])
    Vey_Xslice = (-d['jey'][:,index])/(d['ne'][:,index])   
    Vez_Xslice = (-d['jez'][:,index])/(d['ne'][:,index])
    
    # calculating E X B
    EcrBI = (d['ey']*d['bz']) - (d['ez']*d['by'])
    EcrBJ = (d['ez']*d['bx']) - (d['ex']*d['bz'])
    EcrBK = (d['ex']*d['by']) - (d['ey']*d['bx'])
    MEcrB = np.sqrt(EcrBI**2 + EcrBJ**2 + EcrBK**2)
    EBX_Xslice = EcrBI[:,index]
    EBY_Xslice = EcrBJ[:,index]
    EBZ_Xslice = EcrBK[:,index]
    Meb_Xslice = MEcrB[:,index]
    
    # calculating total pressure
    Pb_Xslice = (Mb_Xslice**2)/2
    Piyy_Xslice = d['piyy'][:,index]
    Peyy_Xslice = d['peyy'][:,index]
    Ptot_Xslice = Pb_Xslice + Piyy_Xslice + d['pezz'][:,index]
    
    # calculating temperatures
    Tixx_Xslice = (d['pixx'][:,index])/(d['ni'][:,index])
    Tiyy_Xslice = (d['piyy'][:,index])/(d['ni'][:,index])
    Tizz_Xslice = (d['pizz'][:,index])/(d['ni'][:,index])
    
    Texx_Xslice = (d['pexx'][:,index])/(d['ne'][:,index])
    Teyy_Xslice = (d['peyy'][:,index])/(d['ne'][:,index])
    Tezz_Xslice = (d['pezz'][:,index])/(d['ne'][:,index])
    
    print("plotting...")
    
    Page7 = plt.figure(7)
    Page7.set_size_inches(8.5,11, forward = True)
    Page7.subplots_adjust(hspace = 1.5)
    Page7.suptitle('Slice With X = ' + str(index), fontsize=20)
        
    sp71 = Page7.add_subplot(521)
    sp71.plot(Y,Bx_Xslice)
    sp71.plot(Y,By_Xslice)
    sp71.plot(Y,Bz_Xslice)
    sp71.plot(Y,Mb_Xslice)
    #  Latex, formatting, positioning of legend
    sp71.legend(['$B_x$', '$B_y$', '$B_z$', '$\mid B\mid$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    # axes formatting 
    sp71.locator_params(nbins = 5, axis = 'y')
    sp71.locator_params(nbins = 6, axis = 'x')
    sp71.set_xlim([0,102.5])

    sp72 = Page7.add_subplot(522)
    sp72.plot(Y,Ex_Xslice)
    sp72.plot(Y,Ey_Xslice)
    sp72.plot(Y,Ez_Xslice)
    sp72.plot(Y,Me_Xslice)
    sp72.legend(['$E_x$', '$E_y$', '$E_z$', '$\mid E\mid$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp72.locator_params(nbins = 5, axis = 'y')
    sp72.locator_params(nbins = 6, axis = 'x')
    sp72.set_xlim([0,102.5])
    
    sp73 = Page7.add_subplot(523)
    sp73.plot(Y,Ni_Xslice)
    sp73.plot(Y,Ne_Xslice)
    sp73.legend(['$N_i$', '$N_e$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp73.locator_params(nbins = 5, axis = 'y')
    sp73.locator_params(nbins = 6, axis = 'x')
    sp73.set_xlim([0,102.5])
    
    sp74 = Page7.add_subplot(524)
    sp74.plot(Y,Jx_Xslice)
    sp74.plot(Y,Jy_Xslice)
    sp74.plot(Y,Jz_Xslice)
    sp74.legend(['$J_x$', '$J_y$', '$J_z$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp74.locator_params(nbins = 5, axis = 'y')
    sp74.locator_params(nbins = 6, axis = 'x')
    sp74.set_xlim([0,102.5])
    
    sp75 = Page7.add_subplot(525)
    sp75.plot(Y,Vix_Xslice)
    sp75.plot(Y,Viy_Xslice)
    sp75.plot(Y,Viz_Xslice)
    sp75.legend(['$V_{ix}$', '$V_{iy}$', '$V_{iz}$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp75.locator_params(nbins = 5, axis = 'y')
    sp75.locator_params(nbins = 6, axis = 'x')
    sp75.set_xlim([0,102.5])
    
    sp76 = Page7.add_subplot(526)
    sp76.plot(Y,Vex_Xslice)
    sp76.plot(Y,Vey_Xslice)
    sp76.plot(Y,Vez_Xslice)
    sp76.legend(['$V_{ex}$', '$V_{ey}$', '$V_{ez}$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp76.locator_params(nbins = 5, axis = 'y')
    sp76.locator_params(nbins = 6, axis = 'x')
    sp76.set_xlim([0,102.5])
    
    sp77 = Page7.add_subplot(527)
    sp77.plot(Y,EBX_Xslice)
    sp77.plot(Y,EBY_Xslice)
    sp77.plot(Y,EBZ_Xslice)
    sp77.plot(Y,Meb_Xslice)
    sp77.legend([r'$(E\times B)_x$', r'$(E\times B)_y$', r'$(E\times B)_z$', r'$\mid E\times B\mid$'], ncol = 2, loc = 'upper center', 
                 bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp77.locator_params(nbins = 5, axis = 'y')
    sp77.locator_params(nbins = 6, axis = 'x')
    sp77.set_xlim([0,102.5])
    
    sp78 = Page7.add_subplot(528)
    sp78.plot(Y,Pb_Xslice)
    sp78.plot(Y,Piyy_Xslice)
    sp78.plot(Y,Peyy_Xslice)
    sp78.plot(Y,Ptot_Xslice)
    sp78.legend(['$P_b$', '$P_{iyy}$', '$P_{eyy}$', '$P_{tot}$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp78.locator_params(nbins = 5, axis = 'y')
    sp78.locator_params(nbins = 6, axis = 'x')
    sp78.set_xlim([0,102.5])
    
    sp79 = Page7.add_subplot(529)
    sp79.plot(Y,Tixx_Xslice)
    sp79.plot(Y,Tiyy_Xslice)
    sp79.plot(Y,Tizz_Xslice)
    sp79.legend(['$T_{ixx}$', '$T_{iyy}$', '$T_{izz}$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp79.locator_params(nbins = 5, axis = 'y')
    sp79.locator_params(nbins = 6, axis = 'x')
    sp79.set_xlim([0,102.5])
    
    sp710 = Page7.add_subplot(5,2,10)
    sp710.plot(Y,Texx_Xslice)
    sp710.plot(Y,Teyy_Xslice)
    sp710.plot(Y,Tezz_Xslice)
    sp710.legend(['$T_{exx}$', '$T_{eyy}$', '$T_{ezz}$'], ncol = 3, loc = 'upper center', 
                 bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp710.locator_params(nbins = 5, axis = 'y')
    sp710.locator_params(nbins = 6, axis = 'x')
    sp710.set_xlim([0,102.5])
    
    
    Page7.show()


# primary slicing and plotting function for constant Y
def sliceY(index):
    d = load_movie()
    X = d['xx']
    
    print("slicing...")    
    
    Bx_Yslice = d['bx'][index,:]
    By_Yslice = d['by'][index,:]
    Bz_Yslice = d['bz'][index,:]
    Mb_Yslice = np.sqrt((d['bx'][index,:])**2 + (d['by'][index,:])**2 + (d['bz'][index,:])**2)
    
    Ex_Yslice = d['ex'][index,:]
    Ey_Yslice = d['ey'][index,:]
    Ez_Yslice = d['ez'][index,:]
    Me_Yslice = np.sqrt((d['ex'][index,:])**2 + (d['ey'][index,:])**2 + (d['ez'][index,:])**2)    
    
    Ni_Yslice = d['ni'][index,:]
    Ne_Yslice = d['ne'][index,:]

    Jx_Yslice = d['jx'][index,:]
    Jy_Yslice = d['jy'][index,:]
    Jz_Yslice = d['jz'][index,:]
    
    Vix_Yslice = (d['jix'][index,:])/(d['ni'][index,:])
    Viy_Yslice = (d['jiy'][index,:])/(d['ni'][index,:])
    Viz_Yslice = (d['jiz'][index,:])/(d['ni'][index,:])

    Vex_Yslice = (-d['jex'][index,:])/(d['ne'][index,:])
    Vey_Yslice = (-d['jey'][index,:])/(d['ne'][index,:])   
    Vez_Yslice = (-d['jez'][index,:])/(d['ne'][index,:])
    
    EcrBI = (d['ey']*d['bz']) - (d['ez']*d['by'])
    EcrBJ = (d['ez']*d['bx']) - (d['ex']*d['bz'])
    EcrBK = (d['ex']*d['by']) - (d['ey']*d['bx'])
    MEcrB = np.sqrt(EcrBI**2 + EcrBJ**2 + EcrBK**2)
    EBX_Yslice = EcrBI[index,:]
    EBY_Yslice = EcrBJ[index,:]
    EBZ_Yslice = EcrBK[index,:]
    Meb_Yslice = MEcrB[index,:]
    
    Pb_Yslice = (Mb_Yslice**2)/2
    Piyy_Yslice = d['piyy'][index,:]
    Peyy_Yslice = d['peyy'][index,:]
    Ptot_Yslice = Pb_Yslice + Piyy_Yslice + d['pezz'][index,:]    
    
    Tixx_Yslice = (d['pixx'][index,:])/(d['ni'][index,:])
    Tiyy_Yslice = (d['piyy'][index,:])/(d['ni'][index,:])
    Tizz_Yslice = (d['pizz'][index,:])/(d['ni'][index,:])
    
    Texx_Yslice = (d['pexx'][index,:])/(d['ne'][index,:])
    Teyy_Yslice = (d['peyy'][index,:])/(d['ne'][index,:])
    Tezz_Yslice = (d['pezz'][index,:])/(d['ne'][index,:])
    
    print("plotting...")    
    
    Page7 = plt.figure(7)
    Page7.set_size_inches(8.5,11, forward = True)
    Page7.subplots_adjust(hspace = 1.5)
    Page7.suptitle('Slice  With Y = ' + str(index), fontsize=20)
        
    sp71 = Page7.add_subplot(521)
    sp71.plot(X,Bx_Yslice)
    sp71.plot(X,By_Yslice)
    sp71.plot(X,Bz_Yslice)
    sp71.plot(X,Mb_Yslice)
    sp71.legend(['$B_x$', '$B_y$', '$B_z$', '$\mid B\mid$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp71.locator_params(nbins = 5, axis = 'y')
    sp71.locator_params(nbins = 6, axis = 'x')
    sp71.set_xlim([0,102.5])

    sp72 = Page7.add_subplot(522)
    sp72.plot(X,Ex_Yslice)
    sp72.plot(X,Ey_Yslice)
    sp72.plot(X,Ez_Yslice)
    sp72.plot(X,Me_Yslice)
    sp72.legend(['$E_x$', '$E_y$', '$E_z$', '$\mid E\mid$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp72.locator_params(nbins = 5, axis = 'y')
    sp72.locator_params(nbins = 6, axis = 'x')
    sp72.set_xlim([0,102.5])
    
    sp73 = Page7.add_subplot(523)
    sp73.plot(X,Ni_Yslice)
    sp73.plot(X,Ne_Yslice)
    sp73.legend(['$N_i$', '$N_e$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp73.locator_params(nbins = 5, axis = 'y')
    sp73.locator_params(nbins = 6, axis = 'x')
    sp73.set_xlim([0,102.5])
    
    sp74 = Page7.add_subplot(524)
    sp74.plot(X,Jx_Yslice)
    sp74.plot(X,Jy_Yslice)
    sp74.plot(X,Jz_Yslice)
    sp74.legend(['$J_x$', '$J_y$', '$J_z$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp74.locator_params(nbins = 5, axis = 'y')
    sp74.locator_params(nbins = 6, axis = 'x')
    sp74.set_xlim([0,102.5])
    
    sp75 = Page7.add_subplot(525)
    sp75.plot(X,Vix_Yslice)
    sp75.plot(X,Viy_Yslice)
    sp75.plot(X,Viz_Yslice)
    sp75.legend(['$V_{ix}$', '$V_{iy}$', '$V_{iz}$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp75.locator_params(nbins = 5, axis = 'y')
    sp75.locator_params(nbins = 6, axis = 'x')
    sp75.set_xlim([0,102.5])
    
    sp76 = Page7.add_subplot(526)
    sp76.plot(X,Vex_Yslice)
    sp76.plot(X,Vey_Yslice)
    sp76.plot(X,Vez_Yslice)
    sp76.legend(['$V_{ex}$', '$V_{ey}$', '$V_{ez}$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp76.locator_params(nbins = 5, axis = 'y')
    sp76.locator_params(nbins = 6, axis = 'x')
    sp76.set_xlim([0,102.5])
    
    sp77 = Page7.add_subplot(527)
    sp77.plot(X,EBX_Yslice)
    sp77.plot(X,EBY_Yslice)
    sp77.plot(X,EBZ_Yslice)
    sp77.plot(X,Meb_Yslice)
    sp77.legend([r'$(E\times B)_x$', r'$(E\times B)_y$', r'$(E\times B)_z$', r'$\mid E\times B\mid$'], ncol = 2, loc = 'upper center', 
                 bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp77.locator_params(nbins = 5, axis = 'y')
    sp77.locator_params(nbins = 6, axis = 'x')
    sp77.set_xlim([0,102.5])
    
    sp78 = Page7.add_subplot(528)
    sp78.plot(X,Pb_Yslice)
    sp78.plot(X,Piyy_Yslice)
    sp78.plot(X,Peyy_Yslice)
    sp78.plot(X,Ptot_Yslice)
    sp78.legend(['$P_b$', '$P_{iyy}$', '$P_{eyy}$', '$P_{tot}$'], ncol = 2, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp78.locator_params(nbins = 5, axis = 'y')
    sp78.locator_params(nbins = 6, axis = 'x')
    sp78.set_xlim([0,102.5])
    
    sp79 = Page7.add_subplot(529)
    sp79.plot(X,Tixx_Yslice)
    sp79.plot(X,Tiyy_Yslice)
    sp79.plot(X,Tizz_Yslice)
    sp79.legend(['$T_{ixx}$', '$T_{iyy}$', '$T_{izz}$'], ncol = 3, loc = 'upper center', 
                bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp79.locator_params(nbins = 5, axis = 'y')
    sp79.locator_params(nbins = 6, axis = 'x')
    sp79.set_xlim([0,102.5])
    
    sp710 = Page7.add_subplot(5,2,10)
    sp710.plot(X,Texx_Yslice)
    sp710.plot(X,Teyy_Yslice)
    sp710.plot(X,Tezz_Yslice)
    sp710.legend(['$T_{exx}$', '$T_{eyy}$', '$T_{ezz}$'], ncol = 3, loc = 'upper center', 
                 bbox_to_anchor = (0.5, -0.2), fontsize = 12)
    sp710.locator_params(nbins = 5, axis = 'y')
    sp710.locator_params(nbins = 6, axis = 'x')
    sp710.set_xlim([0,102.5])
    
    
    Page7.show()

  
# user interface  
i = 0
while(i == 0):
    # catch for non string input
    while True:
        try:
            X = str(raw_input('Take Cuts? Y or N\n '))
            break
        except ValueError:
            print("invalid input try again... \n")
    # handling string entries to first question
    if ((X == 'Y') or (X == 'y')):
        i = 1
        p = 0
        while(p == 0):
            while True:
                try:
                    Axis = str(raw_input('Slice X or Y axis? \n'))
                    break
                except ValueError:
                    print("invalid input try again... \n")
            # handling string entries for axes question
            if ((Axis == 'X') or (Axis == 'x')):
                p = 1
                q = 0
                while(q == 0):
                    # catch for non int entries
                    while True:
                        try:
                            R = int(raw_input('Enter # from 0 to 8191 \n'))
                            break
                        except ValueError:
                            print("invalid input try again... \n")
                    # handling slice position
                    if (0 <= R <= 8191):
                        q = 1
                        sliceX(R)
                    else:
                        print("invalid input try again...\n")
            elif ((Axis == 'Y') or (Axis == 'y')):
                p = 1
                q = 0
                while(q == 0):
                    while True:
                        try:
                            R = int(raw_input('Enter # from 0 to 8191 \n'))
                            break
                        except ValueError:
                            print("invalid input try again... \n")
                    if (0 <= R <= 8191):
                        q = 1
                        sliceY(R)
                    else:
                        print("invalid input try again...\n")
            else:
                print("invalid input try again...\n")
    elif ((X == 'N') or (X == 'n')):
        i = 1
        print("Cuts Cancellled \n")
    else:
        print("invalid input try again...\n")
        

