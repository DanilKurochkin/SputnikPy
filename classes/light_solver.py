import numpy as np
import pandas
from scipy import linalg
from scipy import optimize

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
            self.R[i][i] = 10**46
        for i in range(3):
            self.R[2*i+1][2*i] = 10**46
            self.R[2*i][2*i+1] = 10**46

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
        self.angle_sunset = np.arccos(earth_radius/(self.radius))
        self.time_day = self.period * (np.pi/2+self.angle_sunset)/np.pi
        
        self.labels = ["Зенит", "Надир", "Перпенд 1", "Перпенд 2", "Теневая", "Солнечная"]
        
        self.Qav = np.zeros(6)
        self.Qex = np.zeros(6)
        self.m = np.zeros(self.R.shape)
        self.result = np.zeros(6*3)
        
        self.x0 = 300*np.ones(6)
        self.lambdas = [None]*6

    def make_lambda(self, i, e):
        Qav_curr = self.Qav[i]
        R_reverse_sum = np.sum(1/self.R[i])
        product_rad = sigma*self.areas[i]*e
        R_reverse = self.R[i]
        result = lambda x: (Qav_curr + np.sum(x/R_reverse)) / (product_rad*x[i]**3+R_reverse_sum) - x[i]
        return result
    
    def func(self, x):
        return list(map(lambda f: f(x), self.lambdas))
        
    def temperature_diapozone(self, As, e):
        self.Qav *= 0
        self.Qav += self.generationHeat

        # От солнца
        self.Qav[0]+=self.areas[0]*As[0]*qs/np.pi #зенит
        self.Qav[1]+=self.areas[1]*As[1]*(1-np.cos(self.angle_sunset))*qs/np.pi #надир после 6 часов
        self.Qav[5]+=self.areas[5]*As[5]*(1+np.sin(self.angle_sunset))*qs/np.pi #после до 6 часов

        #От солнечного альбедо
        self.Qav[1]+=As[1]*self.areas[1]*albedo_coef*qs/np.pi

        #От земли
        self.Qav[1]+=e[1]*self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        
        for i in range(6):
            self.lambdas[i] = self.make_lambda(i, e[i])
        
        Tav = optimize.root(self.func, x0=self.x0).x
        
        #расчёт максимума
        self.Qex *= 0
        self.Qex[0] += qs * As[0] * self.areas[0] * np.sqrt(2)/2
        self.Qex[5] += qs * As[1] * self.areas[1] * np.sqrt(2.5)/2
        self.Qex[1] += e[1]*self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        self.Qex[1] += 2*albedo_coef*qs*As[1]*self.areas[1]/np.pi
        self.Qex += self.generationHeat
        time_day_calc = self.time_day/2
        
        self.m = -1/self.R
        self.m[np.diag_indices_from(self.m)] = np.sum(1/self.R, axis=1) + self.heat_capacity[:]/time_day_calc        
        d = self.Qex[:] - self.Qav[:] + Tav[:]*self.heat_capacity[:]/time_day_calc
        Tmax = linalg.solve(self.m, d, overwrite_a=True, overwrite_b=True)
        
        #расчёт минимума
        self.Qex *= 0
        self.Qex += self.generationHeat
        self.Qex[1]+=e[1]*self.areas[1]*earth_radiation*earth_radius**2/self.radius**2
        time_night_calc = (self.period - self.time_day)/2

        self.m = -1/self.R
        self.m[np.diag_indices_from(self.m)] = np.sum(1/self.R, axis=1) + self.heat_capacity[:]/time_night_calc
        d = self.Qex[:] - self.Qav[:] + Tav[:]*self.heat_capacity[:]/time_night_calc
        Tmin = linalg.solve(self.m, d, overwrite_a=True, overwrite_b=True)
        
        np.concatenate([Tmin,Tav,Tmax], out=self.result)
        return self.result