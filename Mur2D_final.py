import numpy as np
import matplotlib.pyplot as plt
from simulation2Dfdtd import fdtd2D


sm = fdtd2D(201, 201, 0.95, 'mur', 'gaussrad')

glocx = np.linspace(0, sm.Lx*sm.dx, sm.Lx-1)    #Genero la dimensión X del espacio de nuestro sistema
glocy = np.linspace(0, sm.Ly*sm.dy, sm.Ly-1)    #Genero la dimensión Y del espacio de nuestro sistema
#1 unidad menos de posicion en ambas direcciones porque se representa la célula de hz
k = 5.
l = 4.
ndiv = 5

filtro = np.zeros((len(glocx),len(glocy)), dtype=bool) #Filtro se encarga de quitar las felchas del plot para que no esté
for i in range(len(filtro)):                           #sobrecargado. Es una matriz de verdadero o falso.
        for j in range(len(filtro)):
            if (i % ndiv == 0 and j % ndiv == 0):
                filtro[i,j] = True

            else:
                filtro[i,j] = False
posx = np.zeros((len(glocx),len(glocy)))   #Genera una matriz en la que se muestra la posición de cada celda en el eje X
posy = np.zeros((len(glocx),len(glocy)))   #Genera una matriz en la que se muestra la posición de cada celda en el eje Y
vx = np.zeros((len(glocx),len(glocy)))
vy = np.zeros((len(glocx),len(glocy)))     #Son las componentes X, Y, Z de los vectores, en nuestro caso un vector tendrá 
vz = np.zeros((len(glocx),len(glocy)))     #sólo X e Y y el otro Z, por lo que no hace falta crear más.

vx[:,:] = 1./2*( sm.ex[:,1:] + sm.ex[:,:-1] )
vy[:,:] = 1./2*( sm.ey[1:,:] + sm.ey[:-1,:] )
vz[:,:] = sm.hz[:,:]

#vz = np.zeros((len(glocx),len(glocy)))
#dirv = np.zeros((len(glocx),len(glocy)))

for i in range(len(glocx)):     #Asigno valores de las posiciones en X
    posx[i,:] = glocx[i]
    
for j in range(len(glocy)):     #Asigno valores de las posiciones en Y
    posy[:,j] = glocy[j]

            
figure = plt.figure(figsize=(12,7))  #Aquí genero la figura
figure.tight_layout()                #Esto deja todo compactado (evita que se pisen mutuamente las gráficas)

ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)             #Esto declara los plot que se pintarán por serparado, así se pueden quitar o añadir 
                                     #fácilmente
        
tRange = np.linspace(0, 10, 401)
for t in tRange:
    
    sm.sim(1)
    vx[:,:] = 1/2*( sm.ex[:,1:] + sm.ex[:,:-1] )
    vy[:,:] = 1/2*( sm.ey[1:,:] + sm.ey[:-1,:] )
    vz[:,:] = sm.hz[:,:]
    minim = np.nanmin(vz)
    maxim = np.nanmax(vz)*1.2
    #Esto es el cálculo de componentes, no forma 
    #No forma parte de la representación en sí
    
    #ax1 es el plot vectorial con el color de fondo, ax2 sólo son las flechas y ax3 los colores con líneal de nivel
    
    e_mod = np.sqrt(vx**2 + vy**2)
    ax1.quiver(posx[filtro],posy[filtro],vx[filtro],vy[filtro],scale=18,angles='xy')
    im=ax1.imshow(e_mod.T,extent=(0.0,1.0,0.0,1.0),origin='lower',cmap='viridis',vmax=maxim,vmin=-0.3) 
#    cb1 = plt.colorbar(im, ax=ax1,label=r'|$\vec{E}$|')
#    ax1.grid()
    
    im=ax2.imshow(vz.T,extent=(0.0,1.0,0.0,1.0),cmap='turbo',origin='lower')
    ax2.contour(posx,posy,vz,4,colors='white',alpha=0.7)
#    cb2 = plt.colorbar(im, ax=ax2,label=r'H$_z$')
#    ax2.grid()
#    plt.grid()
#    plt.ylim(-0.1, 1.1)
#    plt.xlim(glocx[0], glocx[-1])
    plt.pause(0.001) # pausa 0.1s
    ax1.cla() # cierra los plots anteriores
    ax2.cla()
#    ax3.cla()
#    cb1.remove()
#    cb2.remove()
print("END")
