import numpy as np
from classes.models import Sputnik

earth_mass = np.float64(5.976 * 10**24) # в килограммах
G = np.float64(6.67 * 10**(-11)) #гравитационная постоянная
earth_radius = np.float64(6371000) #в метрах

class ClassicOrbit(): #орбита
    def __init__(self, radiusAboveEarth):
        self.radius = radiusAboveEarth*1000 + earth_radius # в метры
        self.velocity = np.sqrt(G * earth_mass / self.radius) # скорость 
        self.length = 2*np.pi*self.radius #длина орбиты
        self.period = self.length/self.velocity #период полного оборота
        self.position = np.float64(0)
        self.sunset = np.array([np.arccos(earth_radius/self.radius) + np.pi/2, 3*np.pi/2 - np.arccos(earth_radius/self.radius)])
        
    def getAlpha(self):
        return self.position/self.radius
    
    def InShadow(self):
        a = self.getAlpha()
        if a>self.sunset[0] and a<self.sunset[1]:
            return True
        return False
    
    def Move(self, ht, sputnik : Sputnik):
        self.position += ht*self.velocity
        if(self.position >= self.length):
            self.position -= self.length
         
class SunLookingOrbit(ClassicOrbit):
    def __init__(self, radiusAboveEarth):
        super().__init__(radiusAboveEarth)
        self.ZenitRotate = False
        self.NadirRotate = True
    
    def Move(self, ht, sputnik : Sputnik):
        self.position += ht*self.velocity
        if(self.position >= self.length):
            self.position -= self.length
        
        if self.NadirRotate and self.position>=self.length/2:
            self.LookAtSun(sputnik)
            self.NadirRotate = False
            self.ZenitRotate = True
            return
        
        if self.position < self.length/2 and self.ZenitRotate and self.position>=np.float64(0):
            self.LookAtSun(sputnik)
            self.NadirRotate = True
            self.ZenitRotate = False
            return
        
    def LookAtSun(self, sputnik : Sputnik):
        matrix = np.array([[1, 0, 0],
                           [0, np.cos(np.pi), -np.sin(np.pi)],
                           [0, np.sin(np.pi), np.cos(np.pi)]])
        for box in sputnik.boxes:
            new_orientation = matrix.dot(box.orientation)
            box.orientation = new_orientation  