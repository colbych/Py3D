import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie


def sliceX(index):
    d = load_movie()
    
    d['bx'] = np.array(d['bx'])
    Bx_Xslice = d['bx'][:,index]
    
    plt.figure(1)
    plt.plot(Bx_Xslice)
    plt.show()
    
def sliceY(index):
    d = load_movie()
    
    d['bx'] = np.array(d['bx'])
    Bx_Yslice = d['bx'][index,:]
    
    plt.figure(2)
    plt.plot(Bx_Yslice)
    plt.show()

i = 0
while(i == 0):
    while True:
        try:
            X = str(raw_input('Take Cuts? Y or N\n '))
            break
        except ValueError:
            print("invalid raw_input try again... \n")
    if ((X == 'Y') or (X == 'y')):
        i = 1
        p = 0
        while(p == 0):
            while True:
                try:
                    Axis = str(raw_input('Slice X or Y axis? \n'))
                    break
                except ValueError:
                    print("invalid raw_input try again... \n")
            if ((Axis == 'X') or (Axis == 'x')):
                p = 1
                q = 0
                while(q == 0):
                    while True:
                        try:
                            R = int(raw_input('Enter # from 0 to 8191 \n'))
                            break
                        except ValueError:
                            print("invalid raw_input try again... \n")
                    if (0 <= R <= 8191):
                        q = 1
                        sliceX(R)
                    else:
                        print("invalid raw_input try again...\n")
            elif ((Axis == 'Y') or (Axis == 'y')):
                p = 1
                q = 0
                while(q == 0):
                    while True:
                        try:
                            R = int(raw_input('Enter # from 0 to 8191 \n'))
                            break
                        except ValueError:
                            print("invalid raw_input try again... \n")
                    if (0 <= R <= 8191):
                        q = 1
                        sliceX(R)
                    else:
                        print("invalid raw_input try again...\n")
            else:
                print("invalid raw_input try again...\n")
    elif ((X == 'N') or (X == 'n')):
        i = 1
        print("Cuts Cancellled \n")
    else:
        print("invalid raw_input try again...\n")
        
        

