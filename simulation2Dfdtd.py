import numpy as np

class fdtd2D:

    def initGauss(self, sx, sy, x0, y0):
        
        lx = self.Lx
        ly = self.Ly
        i0 = x0/self.dx
        j0 = y0/self.dy
        #Modo tangente eléctrico
        for j in range( ly ):
            #El campo x adelantado 1/2 en x (1 posicion en x menos)
            for i in range( lx-1 ):
                self.ex[i,j] = 2*np.exp( -1/2*self.dx*( i + 1/2 - i0 )**2/sx**2  + -1/2*self.dy*( j - j0 )**2/sy**2 )
        for j in range( ly-1 ):
            #El campo y adelantado 1/2 en y (1 posicion menos en y)
            for i in range( lx ):
                self.ey[i,j] = 0
        for j in range( ly-1 ):
            #El campo magnetico en z adelantado diagonalmente
            for i in range( lx-1 ):
                self.hz[i,j] = 2.5*np.exp( -1/2*self.dx*( i + 1/2 - i0 )**2/sx**2  + -1/2*self.dy*( j + 1/2 - j0 )**2/sy**2 )
        
        
    def __init__(self, Lx, Ly, CFL, bcl = ['mur', 'pbc', 'pec']):
        
        self.Lx = Lx
        self.Ly = Ly
        self.CFL = CFL
        self.mu = 1
        self.eps = 1
        self.sigma = 0
        self.c = 1/np.sqrt(self.eps*self.mu)
        self.bcl = 'mur'

        self.x = np.linspace(0, 1, Lx)
        self.y = np.linspace(0, 1, Ly)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dt = CFL/np.sqrt(1/self.dx**2 + 1/self.dy**2)/self.c
        self.ex = np.zeros((Lx-1,Ly))
        self.ey = np.zeros((Lx,Ly-1))
        self.hz = np.zeros((Lx-1,Ly-1))
        self.initGauss(0.4, 0.4, 0.5, 0.5)
        
    def sim(self, T):
        
        mu = self.mu
        eps = self.eps
        dt = self.dt
        dx = self.dx
        dy = self.dy
        assert(dt != 0)
        assert(dx != 0)
        assert(dy != 0)
        assert(eps != 0)
        assert(mu != 0)
        c = 1./np.sqrt(eps*mu)
        sigma = self.sigma
        alpha = 1 - 1/2*sigma*dt/eps
        beta = 1 + 1/2*sigma*dt/eps
        #Metemos las fronteras
        frxmin = int(self.Lx/6)
        frxmax = self.Lx- int(self.Lx/6)-1
        frymin = int(self.Ly/6)
        frymax = self.Ly- int(self.Ly/6)-1
        assert(frxmin >= 0)
        assert(frxmax < self.Lx and frxmax > frxmin)
        assert(frymin >= 0)
        assert(frymax < self.Ly and frymax > frymin)

        if self.bcl == 'pec':
            self.ex[:,0]= np.zeros(self.Lx-1)
            self.ey[0,:]= np.zeros(self.Ly-1)
            self.ex[:,-1]= np.zeros(self.Lx-1)
            self.ey[-1,:]= np.zeros(self.Ly-1)
        for t in range(T):
            if self.bcl == 'mur':
                aux1x = self.ex[:,frymin] #Frontera OY inferior
                aux2x = self.ex[:,frymax] #Frontera OY superior
                aux3x = self.ex[:,frymin + 1] #Antes de frontera OY inf 
                aux4x = self.ex[:,frymax-1] #Antes de frontera OY sup
                aux1y = self.ey[frxmin,:] #Frontera OX izda
                aux2y = self.ey[frxmax,:] # Frontera OX dhca
                aux3y = self.ey[frxmin + 1,:] #Antes de frontera OX izda 
                aux4y = self.ey[frxmax-1,:] #Antes de frontera OX dcha
            elif self.bcl == 'pec':
                self.ex[:,0]= np.zeros(self.Lx-1)
                self.ey[0,:]= np.zeros(self.Ly-1)
                self.ex[:,-1]= np.zeros(self.Lx-1)
                self.ey[-1,:]= np.zeros(self.Ly-1)
            #Hay un rotacional, las Ex y dy van alternadas
            self.hz[:,:] = self.hz[:,:] - dt/mu/dx*( self.ey[1:,:] - self.ey[:-1,:]) + dt/mu/dy*( self.ex[:,1:] - self.ex[:,:-1]) 
            self.ex[:,1:-1] = alpha/beta*self.ex[:,1:-1] + dt/beta/eps/dy*(self.hz[:,1:] - self.hz[:,:-1])
            self.ey[1:-1,:] = alpha/beta*self.ey[1:-1,:] - dt/beta/eps/dx*(self.hz[1:,:] - self.hz[:-1,:])
            if self.bcl == 'mur':
                self.ex[:,frymin] = aux3x[:] + (c * dt - dy )/(c * dt + dy)*(self.ex[:,frymin+1] - aux1x[:])
                self.ex[:,frymax] = aux4x[:] + (c * dt - dy )/(c * dt + dy)*(self.ex[:,frymax-1] - aux2x[:])
                self.ey[frxmin,:] = aux3y[:] + (c * dt - dx )/(c * dt + dx)*(self.ey[frxmin+1,:] - aux1y[:])
                self.ey[frxmax,:] = aux4y[:] + (c * dt - dx )/(c * dt + dx)*(self.ey[frxmax-1,:] - aux2y[:])
    
    def energia(self):
        
        eps = self.eps
        mu = self.mu
        
        e2x = self.ex*self.ex
        e2y = self.ey*self.ey
        h2z = self.hz*self.hz

        self.energy = 1/2 * eps *(np.sum(e2x)) + 1/2 * eps *(np.sum(e2y)) + (1 /2 /nu) * (np.sum(h2z))
                            