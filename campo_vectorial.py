import numpy as np
import matplotlib.pyplot as plt

glocx = np.arange(0, 1, 0.01)    #Genero la dimensión X del espacio de nuestro sistema
glocy = np.arange(0, 1, 0.01)    #Genero la dimensión Y del espacio de nuestro sistema
k = 5.
l = 4.
ndiv = 5

filtro = np.zeros((len(glocx),len(glocy)), dtype=bool) #Filtro se encarga de quitar las felchas del plot para que no esté
for i in range(len(filtro)):                           #sobrecargado. Es una matriz de verdadero o falso.
        for j in range(len(filtro)):
            if (i+j) % ndiv == 0:
                filtro[i,j] = True

            else:
                filtro[i,j] = False
                      
posx = np.zeros((len(glocx),len(glocy)))   #Genera una matriz en la que se muestra la posición de cada celda en el eje X
posy = np.zeros((len(glocx),len(glocy)))   #Genera una matriz en la que se muestra la posición de cada celda en el eje Y
vx = np.zeros((len(glocx),len(glocy)))
vy = np.zeros((len(glocx),len(glocy)))     #Son las componentes X, Y, Z de los vectores, en nuestro caso un vector tendrá 
vz = np.zeros((len(glocx),len(glocy)))     #sólo X e Y y el otro Z, por lo que no hace falta crear más.
#vz = np.zeros((len(glocx),len(glocy)))
#dirv = np.zeros((len(glocx),len(glocy)))

for i in range(len(glocx)):     #Asigno valores de las posiciones en X
    posx[i,:] = glocx[i]
    
for j in range(len(glocx)):     #Asigno valores de las posiciones en Y
    posy[:,j] = glocy[j]

            
figure = plt.figure(figsize=(18,6))  #Aquí genero la figura
figure.tight_layout()                #Esto deja todo compactado (evita que se pisen mutuamente las gráficas)

ax1 = plt.subplot(1,3,1)
ax2 = plt.subplot(1,3,2)             #Esto declara los plot que se pintarán por serparado, así se pueden quitar o añadir 
ax3 = plt.subplot(1,3,3)             #fácilmente
        
tRange = np.linspace(0, 10, 101)
for t in tRange:
    for i in range(len(glocx)):
        for j in range(len(glocy)):
            vx[i,j] = np.sin(l*glocx[i]*np.pi - 2*t)
            vy[i,j] = np.sin(l*glocy[j]*np.pi - 2*t)
            vz[i,j] = np.sin(l*np.sqrt(glocx[i]**2+glocy[j]**2)*np.pi - t)     #Esto es el cálculo de componentes, no forma 
                                                                              #No forma parte de la representación en sí
    
    #ax1 es el plot vectorial con el color de fondo, ax2 sólo son las flechas y ax3 los colores con líneal de nivel
    
    
    ax1.quiver(posx[filtro],posy[filtro],vx[filtro],vy[filtro],scale=40.0)
    im=ax1.imshow(vz,vmax=2,vmin=-2,extent=(0.0,1.0,0.0,1.0),cmap='turbo',origin='lower') 
    cb1 = plt.colorbar(im,ax=ax1)
    ax1.grid()
    
    ax2.quiver(posx[filtro],posy[filtro],vx[filtro],vy[filtro],scale=40.0)
    ax2.grid()
    
    im=ax3.imshow(vz,vmax=2,vmin=-2,extent=(0.0,1.0,0.0,1.0),cmap='turbo',origin='lower')
    ax3.contour(posx,posy,vz,5,colors='white',alpha=0.5)
    cb3 = plt.colorbar(im,ax=ax3)
    ax3.grid()
    
#    plt.grid()
#    plt.ylim(-0.1, 1.1)
#    plt.xlim(glocx[0], glocx[-1])
    plt.pause(0.0001) # pausa 0.1s
    ax1.cla() # cierra los plots anteriores
    ax2.cla()
    ax3.cla()
    cb1.remove()
    cb3.remove()