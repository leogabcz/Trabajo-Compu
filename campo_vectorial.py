import numpy as np
import matplotlib.pyplot as plt

glocx = np.arange(0, 1, 0.01)
glocy = np.arange(0, 1, 0.01)
k = 5.
l = 4.
ndiv = 5

filtro = np.zeros((len(glocx),len(glocy)), dtype=bool)
for i in range(len(filtro)):
        for j in range(len(filtro)):
            if (i+j) % ndiv == 0:
                filtro[i,j] = True

            else:
                filtro[i,j] = False
                      
posx = np.zeros((len(glocx),len(glocy)))
posy = np.zeros((len(glocx),len(glocy)))
vx = np.zeros((len(glocx),len(glocy)))
vy = np.zeros((len(glocx),len(glocy)))
#vz = np.zeros((len(glocx),len(glocy)))
#dirv = np.zeros((len(glocx),len(glocy)))

for i in range(len(glocx)):
    posx[i,:] = glocx[i]
    
for j in range(len(glocx)):
    posy[:,j] = glocy[j]

            

tRange = np.linspace(0, 5, 101)
for t in tRange:
    for i in range(len(glocx)):
        for j in range(len(glocy)):
            vx[i,j] = np.sin(l*glocx[i]*np.pi - k*t)
            vy[i,j] = np.sin(l*glocy[j]*np.pi - k*t)

    
    plt.quiver(posx[filtro],posy[filtro],vx[filtro],vy[filtro],scale=40.0)
    plt.grid()
    plt.ylim(-0.1, 1.1)
    plt.xlim(glocx[0], glocx[-1])
    plt.pause(0.1) # pausa 0.1s
    plt.cla() # cierra los plots anteriores