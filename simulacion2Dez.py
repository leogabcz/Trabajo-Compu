import numpy as np

class fdtd2Dez:

    def initGauss(self, sx, sy, x0, y0):
        
        lx = self.Lx
        ly = self.Ly
        i0 = x0/self.dx
        j0 = y0/self.dy
        #Modo tangente eléctrico
        for j in range( ly ):
            #El campo x adelantado 1/2 en x (1 posicion en x menos)
            for i in range( lx ):
                self.ez[i,j] = 2*np.exp( -1/2*self.dx*( i -i0 )**2/sx**2  + -1/2*self.dy*( j - j0 )**2/sy**2 )
        for j in range( ly -1):
            #El campo magnetico en z adelantado diagonalmente
            for i in range( lx ):
                self.hx[i,j] = 0
        for j in range( ly ):
            #El campo magnetico en z adelantado diagonalmente
            for i in range( lx -1 ):
                self.hy[i,j] = 0



        
        
    def __init__(self, Lx, Ly, CFL, bcl = ['mur', 'pbc', 'pec', 'mur2'], condini = ['gausstang', 'gaussrad']):
        
        self.Lx = Lx
        self.Ly = Ly
        self.CFL = CFL
        self.mu = 1
        self.eps = 1
        self.sigma = 0
        self.c = 1/np.sqrt(self.eps*self.mu)
        self.bcl = bcl
        self.condini = condini 
        
        self.x = np.linspace(0, 1, Lx)
        self.y = np.linspace(0, 1, Ly)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dt = CFL/np.sqrt(1/self.dx**2 + 1/self.dy**2)/self.c
        self.hx = np.zeros((Lx,Ly-1))
        self.hy = np.zeros((Lx-1,Ly))
        self.ez = np.zeros((Lx,Ly))
        
        if self.condini == 'gausstang':
            self.initGauss(0.4, 0.4, 0.5, 0.5)
        elif self.condini == 'gaussrad':
            self.initGaussRadial(0.4, 0.4, 0.5, 0.5)
        
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
        frxmin = int(self.Lx/1000)
        frxmax = self.Lx - 1
        frymin = int(self.Ly/1000) 
        frymax = self.Ly - 1
        assert(frxmin >= 0)
        assert(frxmax < self.Lx and frxmax > frxmin)
        assert(frymin >= 0)
        assert(frymax < self.Ly and frymax > frymin)

        if self.bcl == 'pec':
            self.ez[:,0]= np.zeros(self.Ly )
            self.ez[-1,:]= np.zeros(self.Lx )
            self.ez[:,-1]= np.zeros(self.Ly )
            self.ez[0,:]= np.zeros(self.Lx )
        for t in range(T):
            if self.bcl == 'mur' or self.bcl == 'mur2':
                aux1 = self.ez[:,frymin] #Frontera OY inferior
                aux2 = self.ez[:,frymax] #Frontera OY superior
                aux3 = self.ez[:,frymin + 1] #Antes de frontera OY inf 
                aux4 = self.ez[:,frymax - 1] #Antes de frontera OY sup
                aux5 = self.ez[frxmin,:] #Frontera OX izda
                aux6 = self.ez[frxmax,:] # Frontera OX dhca
                aux7 = self.ez[frxmin + 1,:] #Antes de frontera OX izda 
                aux8 = self.ez[frxmax - 1,:] #Antes de frontera OX dcha
            elif self.bcl == 'pec':
                 self.ez[:,0]= np.zeros(self.Lx)
                 self.ez[-1,:]= np.zeros(self.Ly)
                 self.ez[:,-1]= np.zeros(self.Lx)
                 self.ez[0,:]= np.zeros(self.Ly)
            
            #Hay un rotacional, las Ex y dy van alternadas
            self.hx[:,:] =  self.hx[:,:] - dt/mu/dx*( self.ez[:,1:] - self.ez[:,:-1]) 
            self.hy[:,:] =  self.hy[:,:] + dt/mu/dx*( self.ez[1:,:] - self.ez[:-1,:]) 
            self.ez[1:-1,1:-1] = alpha/beta*self.ez[1:-1,1:-1] + dt/beta/eps/dy*(self.hy[1:,1:-1] - self.hy[:-1,1:-1]) - dt/beta/eps/dy*(self.hx[1:-1,1:] - self.hx[1:-1,:-1])
            
            
            if self.bcl == 'mur':
                self.ez[:,frymin] = aux3[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[:,frymin+1] - aux1[:])
                self.ez[:,frymax] = aux4[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[:,frymax-1] - aux2[:])
                self.ez[frxmin,:] = aux7[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[frxmin,:] - aux5[:])
                self.ez[frxmax,:] = aux8[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[frxmax,:] - aux6[:])

            if self.bcl == 'mur2':
                self.ez[:,frymin] = aux3[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[:,frymin+1] - aux1[:]) + mu*c/ 2*(dy + c*dt) * (self.hy[1:,frymin] - self.hy[:-1,frymin] + self.hy[1:,frymin+1] - self.hy[:-1,frymin+1])
                self.ez[:,frymax] = aux4[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[:,frymax-1] - aux2[:]) + mu*c/ 2*(dy + c*dt) * (self.hy[1:,frymax] - self.hy[:-1,frymax] + self.hy[1:,frymax-1] - self.hy[:-1,frymax-1])
                self.ez[frxmin,:] = aux7[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[frxmin,:] - aux5[:]) - mu*c/ 2*(dy + c*dt) * (self.hx[frxmin,1:] - self.hx[frxmin,:-1] + self.hx[frxmin+1,1:] - self.hx[frxmin+1,:-1])
                self.ez[frxmax,:] = aux8[:] + (c * dt - dy )/(c * dt + dy)*(self.ez[frxmax,:] - aux6[:]) - mu*c/ 2*(dy + c*dt) * (self.hx[frxmax,1:] - self.hx[frxmax,:-1] + self.hx[frxmax-1,1:] - self.hx[frxmax-1,:-1])

     
        self.energy = self.energia(frxmin, frxmax, frymin, frymax)        
    
    
    def energia(self, frxmin, frxmax, frymin, frymax):
        
        eps = self.eps
        mu = self.mu
        
        h2x = self.hx*self.hx
        h2y = self.hy*self.hy
        e2z = self.ez*self.ez
        energy = 1/2 * eps *(np.sum(e2z[frxmin:frxmax,frymin:frymax])) + (1/ 2 / mu) *(np.sum(h2y)) + (1 /2 /mu) * (np.sum(h2x))
       
        return energy               
        
        
    
