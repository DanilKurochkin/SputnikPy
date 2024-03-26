import numpy as np
import pandas
from scipy import linalg
from itertools import product

sigma = np.float64(5.67*10**(-8))
sigma = np.float64(5.67*10**(-8))
earth_radius = np.float64(6371000)
earth_mass = np.float64(5.976 * 10**24)
G = np.float64(6.67 * 10**(-11))
albedo_coef = 0.3
earth_radiation = 240
qs = np.float64(1367)

class Solver:
    def __init__(self, width, Lx, Ly, Lz, c, p, R, orbit_height, generation_heat) -> None:
        self.generationHeat = generation_heat

        self.width = 0.005
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.R = R * np.ones(shape=(6,6))
        for i in range(6):
            self.R[i][i] = 0
        for i in range(3):
            self.R[2*i+1][2*i] = 0
            self.R[2*i][2*i+1] = 0

        self.areas = np.array(
            [
            Ly*Lz, Ly*Lz,
            Lx*Lz, Lx*Lz,
            Lx*Ly, Lx*Ly
            ]
        )
        self.heat_capacity = c*p*width*self.areas[:]
        
        self.radius = orbit_height*1000 + earth_radius # в метры
        velocity = np.sqrt(G * earth_mass / self.radius) # скорость 
        length = 2*np.pi*self.radius #длина орбиты
        self.period = length/velocity #период полного оборота
        self.angleSunSet = np.arccos(earth_radius/(self.radius))
        self.tDay = self.period * (np.pi/2+self.angleSunSet)/np.pi
        
        self.labels = ["Зенит", "Надир", "Перпенд 1", "Перпенд 2", "Теневая", "Солнечная"]
    
    def temperature_diapozone(self, As, e):
        Qav = np.zeros(6)
        Qav += self.generationHeat
        Qav

        # От солнца
        Qav[0]+=self.areas[0]*As[0]*qs/np.pi #зенит
        Qav[1]+=self.areas[1]*As[1]*(1-np.cos(self.angleSunSet))*qs/np.pi #надир после 6 часов
        Qav[5]+=self.areas[5]*As[5]*(1+np.sin(self.angleSunSet))*qs/np.pi #после до 6 часов

        #От солнечного альбедо
        Qav[1]+=As[1]*self.areas[1]*albedo_coef*qs/np.pi

        #От земли
        Qav[1]+=e[1]*self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        
        Tav = 300*np.ones(6)
        eps = 0.001
        eps_curr = 100
        iter = 0
        alpha = 0.2
        while eps < eps_curr:
            eps_curr = -1
            iter += 1
            for i in range(6):
                index = self.R[i] > 0

                newT = (1.0-alpha)*Tav[i] + alpha*(Qav[i] + np.sum(Tav[index]/self.R[i][index])) / (sigma*Tav[i]**3*self.areas[i]*e[i]+np.sum(1/self.R[i][index]))
                eps_curr = max(np.abs(newT - Tav[i]), eps_curr)
                Tav[i] = newT

        Qex = np.zeros(6)
        Qex[0] += qs * As[0] * self.areas[0] * np.sqrt(2)/2
        Qex[5] += qs * As[1] * self.areas[1] * np.sqrt(2.5)/2
        Qex[1] += e[1]*self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        Qex[1] += 2*albedo_coef*qs*As[1]*self.areas[1]/np.pi
        Qex += self.generationHeat
        tDayCalculation = self.tDay/2
        
        m = np.zeros(self.R.shape)
        for i in range(6):
            for j in range(6):
                if i == j:
                    index = self.R[i] > 0
                    m[i][i] = np.abs(np.sum(1/self.R[i][index])) + self.heat_capacity[i]/tDayCalculation
                else:
                    if self.R[i][j] != 0:
                        m[i][j] = -1/self.R[i][j]
        d = Qex[:] - Qav[:] + Tav[:]*self.heat_capacity[:]/tDayCalculation
        Tmax = linalg.solve(m, d)
        
        tNight = (self.period - self.tDay)/2

        Qex = self.generationHeat
        Qex[1]+=e[1]*self.areas[1]*earth_radiation*earth_radius**2/self.radius**2  
        for i in range(6):
            for j in range(6):
                if i == j:
                    index = self.R[i] > 0
                    m[i][i] = np.abs(np.sum(1/self.R[i][index]/2)) + self.heat_capacity[i]/tNight
                else:
                    if self.R[i][j] != 0:
                        m[i][j] = -1/self.R[i][j]/2

        d = Qex[:] - Qav[:] + Tav[:]*self.heat_capacity[:]/tNight
        Tmin = linalg.solve(m, d)
        
        return Tmin, Tav, Tmax
    
s = Solver(
    width=0.005,
    Lx=1,
    Ly=1,
    Lz=1,
    c=800,
    p=2700,
    R=2,
    orbit_height=500,
    generation_heat=50*np.ones(6)
)