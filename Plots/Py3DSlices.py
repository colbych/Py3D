import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie


def sliceX(index):
    d = load_movie()
    
    Bx_Xslice = d['bx'][:,index]
    By_Xslice = d['by'][:,index]
    Bz_Xslice = d['bz'][:,index]
    Mb_Xslice = np.sqrt((d['bx'][:,index])**2 + (d['by'][:,index])**2 + (d['bz'][:,index])**2)
    
    Ex_Xslice = d['ex'][:,index]
    Ey_Xslice = d['ey'][:,index]
    Ez_Xslice = d['ez'][:,index]
    Me_Xslice = np.sqrt((d['ex'][:,index])**2 + (d['ey'][:,index])**2 + (d['ez'][:,index])**2)    
    
    Ni_Xslice = d['ni'][:,index]
    Ne_Xslice = d['ne'][:,index]

    Jx_Xslice = d['jx'][:,index]
    Jy_Xslice = d['jy'][:,index]
    Jz_Xslice = d['jz'][:,index]
    
    Vix_Xslice = (d['jix'][:,index])/(d['ni'][:,index])
    Viy_Xslice = (d['jiy'][:,index])/(d['ni'][:,index])
    Viz_Xslice = (d['jiz'][:,index])/(d['ni'][:,index])

    Vex_Xslice = (-d['jex'][:,index])/(d['ne'][:,index])
    Vey_Xslice = (-d['jey'][:,index])/(d['ne'][:,index])   
    Vez_Xslice = (-d['jez'][:,index])/(d['ne'][:,index])
    
    EcrBI = (d['ey']*d['bz']) - (d['ez']*d['by'])
    EcrBJ = (d['ez']*d['bx']) - (d['ex']*d['bz'])
    EcrBK = (d['ex']*d['by']) - (d['ey']*d['bx'])
    MEcrB = np.sqrt(EcrBI**2 + EcrBJ**2 + EcrBK**2)
    EBX_Xslice = EcrBI[:,index]
    EBY_Xslice = EcrBJ[:,index]
    EBZ_Xslice = EcrBK[:,index]
    Meb_Xslice = MEcrB[:,index]
    
    Tixx_Xslice = (d['pixx'][:,index])/(d['ni'][:,index])
    Tiyy_Xslice = (d['piyy'][:,index])/(d['ni'][:,index])
    Tizz_Xslice = (d['pizz'][:,index])/(d['ni'][:,index])
    
    Texx_Xslice = (d['pexx'][:,index])/(d['ne'][:,index])
    Teyy_Xslice = (d['peyy'][:,index])/(d['ne'][:,index])
    Tezz_Xslice = (d['pezz'][:,index])/(d['ne'][:,index])
    
    
    Page7 = plt.figure(7)
    Page7.subplots_adjust(hspace = .3)
    Page7.subplots_adjust(wspace = .2)   
        
    sp71 = Page7.add_subplot(581)
    plt.plot(Bx_Xslice)
    plt.xticks([2000,8000])
    sp71.set_title("Bx X Slice")
    
    sp72 = Page7.add_subplot(582)
    plt.plot(By_Xslice)
    plt.xticks([2000,8000])
    sp72.set_title("By X Slice")
    
    sp73 = Page7.add_subplot(583)
    plt.plot(Bz_Xslice)
    plt.xticks([2000,8000])
    sp73.set_title("Bz X Slice")
    
    sp74 = Page7.add_subplot(584)
    plt.plot(Mb_Xslice)
    plt.xticks([2000,8000])
    sp74.set_title("|B| X Slice")

    sp75 = Page7.add_subplot(585)
    plt.plot(Ex_Xslice)
    plt.xticks([2000,8000])
    sp75.set_title("Ex X Slice")
    
    sp76 = Page7.add_subplot(586)
    plt.plot(Ey_Xslice)
    plt.xticks([2000,8000])
    sp76.set_title("Ey X Slice")
    
    sp77 = Page7.add_subplot(587)
    plt.plot(Ez_Xslice)
    plt.xticks([2000,8000])
    sp77.set_title("Ez X Slice")
    
    sp78 = Page7.add_subplot(588)
    plt.plot(Me_Xslice)
    plt.xticks([2000,8000])
    sp78.set_title("|E| X Slice")
    
    sp79 = Page7.add_subplot(589)
    plt.plot(Ni_Xslice)
    plt.xticks([2000,8000])
    sp79.set_title("Ni X Slice")

    sp710 = Page7.add_subplot(5,8,10)
    plt.plot(Ne_Xslice)
    plt.xticks([2000,8000])
    sp710.set_title("Ne X Slice")
    
    sp711 = Page7.add_subplot(5,8,11)
    plt.plot(Jx_Xslice)
    plt.xticks([2000,8000])
    sp711.set_title("Jx X Slice")
    
    sp712 = Page7.add_subplot(5,8,12)
    plt.plot(Jy_Xslice)
    plt.xticks([2000,8000])
    sp712.set_title("Jy X Slice")
    
    sp713 = Page7.add_subplot(5,8,13)
    plt.plot(Jz_Xslice)
    plt.xticks([2000,8000])
    sp713.set_title("Jz X Slice")
    
    sp717 = Page7.add_subplot(5,8,17)
    plt.plot(Vix_Xslice)
    plt.xticks([2000,8000])
    sp717.set_title("Vix X Slice")
    
    sp718 = Page7.add_subplot(5,8,18)
    plt.plot(Viy_Xslice)
    plt.xticks([2000,8000])
    sp718.set_title("Viy X Slice")

    sp719 = Page7.add_subplot(5,8,19)
    plt.plot(Viz_Xslice)
    plt.xticks([2000,8000])
    sp719.set_title("Viz X Slice")
    
    sp720 = Page7.add_subplot(5,8,20)
    plt.plot(Vex_Xslice)
    plt.xticks([2000,8000])
    sp720.set_title("Vex X Slice")
    
    sp721 = Page7.add_subplot(5,8,21)
    plt.plot(Vey_Xslice)
    plt.xticks([2000,8000])
    sp721.set_title("Vey X Slice")
    
    sp722 = Page7.add_subplot(5,8,22)
    plt.plot(Vez_Xslice)
    plt.xticks([2000,8000])
    sp722.set_title("Vez X Slice")
    
    sp725 = Page7.add_subplot(5,8,25)
    plt.plot(EBX_Xslice)
    plt.xticks([2000,8000])
    sp725.set_title("EcrBx X Slice")
    
    sp726 = Page7.add_subplot(5,8,26)
    plt.plot(EBY_Xslice)
    plt.xticks([2000,8000])
    sp726.set_title("EcrBy X Slice")
    
    sp727 = Page7.add_subplot(5,8,27)
    plt.plot(EBZ_Xslice)
    plt.xticks([2000,8000])
    sp727.set_title("EcrBz X Slice")
    
    sp728 = Page7.add_subplot(5,8,28)
    plt.plot(Meb_Xslice)
    plt.xticks([2000,8000])
    sp728.set_title("|EcrB| X Slice")
    
    sp733 = Page7.add_subplot(5,8,33)
    plt.plot(Tixx_Xslice)
    plt.xticks([2000,8000])
    sp733.set_title("Tixx X Slice")
    
    sp734 = Page7.add_subplot(5,8,34)
    plt.plot(Tiyy_Xslice)
    plt.xticks([2000,8000])
    sp734.set_title("Tiyy X Slice")
    
    sp735 = Page7.add_subplot(5,8,35)
    plt.plot(Tizz_Xslice)
    plt.xticks([2000,8000])
    sp735.set_title("Tizz X Slice")
    
    sp736 = Page7.add_subplot(5,8,36)
    plt.plot(Texx_Xslice)
    plt.xticks([2000,8000])
    sp736.set_title("Texx X Slice")
    
    sp737 = Page7.add_subplot(5,8,37)
    plt.plot(Teyy_Xslice)
    plt.xticks([2000,8000])
    sp737.set_title("Teyy X Slice")
    
    sp738 = Page7.add_subplot(5,8,38)
    plt.plot(Tezz_Xslice)
    plt.xticks([2000,8000])
    sp738.set_title("Tezz X Slice")
    
    
    Page7.show()

    
def sliceY(index):
    d = load_movie()
    
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
    
    Tixx_Yslice = (d['pixx'][index,:])/(d['ni'][index,:])
    Tiyy_Yslice = (d['piyy'][index,:])/(d['ni'][index,:])
    Tizz_Yslice = (d['pizz'][index,:])/(d['ni'][index,:])
    
    Texx_Yslice = (d['pexx'][index,:])/(d['ne'][index,:])
    Teyy_Yslice = (d['peyy'][index,:])/(d['ne'][index,:])
    Tezz_Yslice = (d['pezz'][index,:])/(d['ne'][index,:])
    
    
    Page7 = plt.figure(7)
    Page7.subplots_adjust(hspace = .3)
    Page7.subplots_adjust(wspace = .2) 
        
        
    sp71 = Page7.add_subplot(581)
    plt.plot(Bx_Yslice)
    plt.xticks([2000,8000])
    sp71.set_title("Bx Y Slice")
    
    sp72 = Page7.add_subplot(582)
    plt.plot(By_Yslice)
    plt.xticks([2000,8000])
    sp72.set_title("By Y Slice")
    
    sp73 = Page7.add_subplot(583)
    plt.plot(Bz_Yslice)
    plt.xticks([2000,8000])
    sp73.set_title("Bz Y Slice")
    
    sp74 = Page7.add_subplot(584)
    plt.plot(Mb_Yslice)
    plt.xticks([2000,8000])
    sp74.set_title("|B| Y Slice")

    sp75 = Page7.add_subplot(585)
    plt.plot(Ex_Yslice)
    plt.xticks([2000,8000])
    sp75.set_title("Ex Y Slice")
    
    sp76 = Page7.add_subplot(586)
    plt.plot(Ey_Yslice)
    plt.xticks([2000,8000])
    sp76.set_title("Ey Y Slice")
    
    sp77 = Page7.add_subplot(587)
    plt.plot(Ez_Yslice)
    plt.xticks([2000,8000])
    sp77.set_title("Ez Y Slice")
    
    sp78 = Page7.add_subplot(588)
    plt.plot(Me_Yslice)
    plt.xticks([2000,8000])
    sp78.set_title("|E| Y Slice")
    
    sp79 = Page7.add_subplot(589)
    plt.plot(Ni_Yslice)
    plt.xticks([2000,8000])
    sp79.set_title("Ni Y Slice")

    sp710 = Page7.add_subplot(5,8,10)
    plt.plot(Ne_Yslice)
    plt.xticks([2000,8000])
    sp710.set_title("Ne Y Slice")
    
    sp711 = Page7.add_subplot(5,8,11)
    plt.plot(Jx_Yslice)
    plt.xticks([2000,8000])
    sp711.set_title("Jx Y Slice")
    
    sp712 = Page7.add_subplot(5,8,12)
    plt.plot(Jy_Yslice)
    plt.xticks([2000,8000])
    sp712.set_title("Jy Y Slice")
    
    sp713 = Page7.add_subplot(5,8,13)
    plt.plot(Jz_Yslice)
    plt.xticks([2000,8000])
    sp713.set_title("Jz Y Slice")
    
    sp717 = Page7.add_subplot(5,8,17)
    plt.plot(Vix_Yslice)
    plt.xticks([2000,8000])
    sp717.set_title("Vix Y Slice")
    
    sp718 = Page7.add_subplot(5,8,18)
    plt.plot(Viy_Yslice)
    plt.xticks([2000,8000])
    sp718.set_title("Viy Y Slice")

    sp719 = Page7.add_subplot(5,8,19)
    plt.plot(Viz_Yslice)
    plt.xticks([2000,8000])
    sp719.set_title("Viz Y Slice")
    
    sp720 = Page7.add_subplot(5,8,20)
    plt.plot(Vex_Yslice)
    plt.xticks([2000,8000])
    sp720.set_title("Vex Y Slice")
    
    sp721 = Page7.add_subplot(5,8,21)
    plt.plot(Vey_Yslice)
    plt.xticks([2000,8000])
    sp721.set_title("Vey Y Slice")
    
    sp722 = Page7.add_subplot(5,8,22)
    plt.plot(Vez_Yslice)
    plt.xticks([2000,8000])
    sp722.set_title("Vez Y Slice")
    
    sp725 = Page7.add_subplot(5,8,25)
    plt.plot(EBX_Yslice)
    plt.xticks([2000,8000])
    sp725.set_title("EcrBx Y Slice")
    
    sp726 = Page7.add_subplot(5,8,26)
    plt.plot(EBY_Yslice)
    plt.xticks([2000,8000])
    sp726.set_title("EcrBy Y Slice")
    
    sp727 = Page7.add_subplot(5,8,27)
    plt.plot(EBZ_Yslice)
    plt.xticks([2000,8000])
    sp727.set_title("EcrBz Y Slice")
    
    sp728 = Page7.add_subplot(5,8,28)
    plt.plot(Meb_Yslice)
    plt.xticks([2000,8000])
    sp728.set_title("|EcrB| Y Slice")
    
    sp733 = Page7.add_subplot(5,8,33)
    plt.plot(Tixx_Yslice)
    plt.xticks([2000,8000])
    sp733.set_title("Tixx Y Slice")
    
    sp734 = Page7.add_subplot(5,8,34)
    plt.plot(Tiyy_Yslice)
    plt.xticks([2000,8000])
    sp734.set_title("Tiyy Y Slice")
    
    sp735 = Page7.add_subplot(5,8,35)
    plt.plot(Tizz_Yslice)
    plt.xticks([2000,8000])
    sp735.set_title("Tizz Y Slice")
    
    sp736 = Page7.add_subplot(5,8,36)
    plt.plot(Texx_Yslice)
    plt.xticks([2000,8000])
    sp736.set_title("Texx Y Slice")
    
    sp737 = Page7.add_subplot(5,8,37)
    plt.plot(Teyy_Yslice)
    plt.xticks([2000,8000])
    sp737.set_title("Teyy Y Slice")
    
    sp738 = Page7.add_subplot(5,8,38)
    plt.plot(Tezz_Yslice)
    plt.xticks([2000,8000])
    sp738.set_title("Tezz Y Slice")
    
    
    Page7.show()

    
    
i = 0
while(i == 0):
    while True:
        try:
            X = str(raw_input('Take Cuts? Y or N\n '))
            break
        except ValueError:
            print("invalid input try again... \n")
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
            if ((Axis == 'X') or (Axis == 'x')):
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
        

