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
        self.sunset = np.array([np.arccos(earth_radius/self.radius) + np.pi/2, 3*np.pi/2 - np.arccos(earth_radius/self.radius)]) # угол захода в тень
        
    def getAlpha(self): # угол между линией Земля-Солнце и спутником
        return self.position/self.radius
    
    def InShadow(self): #проверка на нахожеднии в тени
        a = self.getAlpha()
        if a>self.sunset[0] and a<self.sunset[1]:
            return True
        return False
    
    def Move(self, ht, sputnik : Sputnik): #измение позиции спутника
        self.position += ht*self.velocity
        if(self.position >= self.length):
            self.position -= self.length
         
class SunLookingOrbit(ClassicOrbit): #разница между этой орбитой и классической в том, что спутник разворачивается по оси надир-зенит в полдень и полночь
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
        
    def LookAtSun(self, sputnik : Sputnik): #разворот производится с помощью умножения матрицы вращения и нормали поверхностей пластин
        matrix = np.array([[1, 0, 0],
                           [0, np.cos(np.pi), -np.sin(np.pi)],
                           [0, np.sin(np.pi), np.cos(np.pi)]])
        for box in sputnik.boxes:
            new_orientation = matrix.dot(box.orientation)
            box.orientation = new_orientation  