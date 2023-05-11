import numpy as np
import matplotlib.pyplot as plt

##--- Declaracion  de las variables iniciales-----##
mu = 1.0
epsilon = 1.0
c0 = 1.0/np.sqrt(epsilon*mu)

CFL = 1.0
tfinal = 20.0
L = 10.0

dx = 0.1
dy = 0.1

#--Creacion de la maya y de la dual
x = np.arange(0, 10+dy, dy)
y = np.arange(0, 10+dy, dy)
xv, yv = np.meshgrid(x, y, sparse=True) #El yv es porque meshgrid devuelve dos grids, pero no lo uso para nada

xDual = (x[1:] + x[:-1])/2
yDual = (y[1:] + y[:-1])/2
xvDual, yvDual = np.array(np.meshgrid(xDual, yDual))

#--Inicializo a gaussiana
x0, y0 = (5, 5)
s0 = 1.

ez =  np.exp( -(xv-x0)**2 / (2*s0**2) )*np.exp( -(yv-y0)**2 / (2*s0**2) )

ez[:,0] = 0.0 
ez[:,-1] = 0.0
ez[0,:] = 0.0 
ez[-1,:] = 0.0

epsilon = np.ones((len(x), len(y)))

hx = np.zeros((len(xDual),len(yDual)))
hy = np.zeros((len(xDual),len(yDual)))

dt = CFL * dx / c0 

tRange = np.arange(0,10,dt)

for t in tRange: 

    ez[1:-1,1:-1] = ez[1:-1,1:-1] + (dt/dx/epsilon[1:-1,1:-1])*(hx[1:,:-1] - hx[1:,1:] + hy[1:,1:] - hy[:-1,1:])
        
    #Actualizamos condiciones de forntera, ahora mismo las tengo puesto a 0
    ez[:,0] = 0.0 
    ez[:,-1] = 0.0
    ez[0,:] = 0.0 
    ez[-1,:] = 0.0

    #Actualizamos las h
    hx = hx + (-dt/dx/mu)*(ez[1:,:-1] - ez[1:,1:])
    hy = hy + (-dt/dx/mu)*(ez[1:,1:] - ez[:-1,1:])

    plt.contourf(x, y, ez)
    plt.show()
    plt.pause(0.1)
    plt.cla()

plt.show()