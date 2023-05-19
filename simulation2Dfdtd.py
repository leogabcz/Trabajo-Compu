import numpy as np
import matplotlib.pyplot as plt

class fdtd2D:
    Lx = 201 #numero de celdas horizontales
    Ly = 201 #numero de celdas verticales
    ex = np.zeros(Lx-1, Ly)
    ey = np.zeros(Lx, Ly-1)   #Se declaran todos del mismo tamaño
    hz = np.zeros(Lx-1, Ly-1)   #luego simplemente se omitiran los ultimos
    bcl = ''    
    CFL = 0.9

    def __init__(self, Lx, Ly, CFL, bcl = ['mur', 'pbc', 'pec']) -> None:
        
        self.Lx = Lx
        self.Ly = Ly
        self.CFL = CFL

        self.x = np.linspace(0, 10, Lx)
        self.y = np.linspace(0, 10, Ly)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.xDual = (self.x[1:] + self.x[:-1])/2
        self.yDual = (self.y[1:] + self.y[:-1])/2
        self.initGauss(self, Lx/20, Lx/20, Lx/2, Ly/2)
        self.dt = CFL*np.sqrt(self.dx**2 + self.dy**2)
        self.mu = 1
        self.eps = 1
        self.sigma = 0
        self.bcl = 'mur'

    def sim(self, T):
        
        exaux=np.zeros(self.x.shape**2)
        eyaux=np.zeros(self.x.shape**2)
        hzaux=np.zeros(self.xDual.shape)
        alpha = 1 - 1/2*self.sigma*self.dt/self.eps
        beta = 1 + 1/2*self.sigma*self.dt/self.eps
        #Condicones periodicas por ahora
        self.ex[:,0]= np.zeros(self.lx-1)
        self.ey[0,:]= np.zeros(self.ly-1)
        self.ex[:,-1]= np.zeros(self.lx-1)
        self.ey[-1,:]= np.zeros(self.ly-1)
        for t in T:
            #Hay un rotacional, las Ex y dy van alternadas
            self.hz[:,:] = self.hz[:,:] - self.dt/self.mu/self.dx*( self.ey[1:,:] - self.ey[:-1,:]) - self.dt/self.mu/self.dy*( self.ex[:,1:] - self.ex[:,:-1]) 
            self.ex[:,1:-1] = alpha/beta*self.ex[:,1:-1] + self.dt/beta/self.eps/self.dy*(self.hz[:,1:] - self.hz[:,:-1])
            self.ey[1:-1,:] = alpha/beta*self.ey[1:-1,:] + self.dt/beta/self.eps/self.dx*(self.hz[1:,:] - self.hz[:-1,:])
        
    def initGauss(self, sx, sy, x0, y0):
        
        lx = self.Lx
        ly = self.Ly
        i0 = x0/self.dx
        j0 = y0/self.dy
        #Modo tangente eléctrico
        for j in range( ly ):
            #El campo x adelantado 1/2 en x (1 posicion en x menos)
            for i in range( lx-1 ):
                self.ex[i,j] = np.exp( -1/2*self.dx*( i + 1/2 - i0 )**2/sx**2  + -1/2*self.dy*( j - j0 )**2/sy**2 )
        for j in range( ly-1 ):
            #El campo y adelantado 1/2 en y (1 posicion menos en y)
            for i in range( lx ):
                self.ey[i,j] = 0
        for j in range( ly-1 ):
            #El campo magnetico en z adelantado diagonalmente
            for i in range( lx-1 ):
                self.hz[i,j] = np.exp( -1/2*self.dx*( i + 1/2 - i0 )**2/sx**2  + -1/2*self.dy*( j + 1/2 - j0 )**2/sy**2 )
        
        