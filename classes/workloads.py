from abc import ABC, abstractmethod
import numpy as np

sigma = np.float64(5.67*10**(-8))

class Load(ABC):
    def __init__(self) -> None:
        pass

    def rotate(self, alpha):
        pass

    def heatFlux(self, box, T):
        pass
    
    def heat(self, box, T):
        Q = self.heatFlux(box, T) * box.area
        
        return Q
    
    def derivative(self, box, T):
        return 0

class Sun(Load): #солнышко и его вращение**
    
    def __init__(self, q):
        self.q = np.float64(q)
        self.default_orientation = np.array([1, 0, 0])
        self.orientation = []
        self.rotate(0)
    
    def rotate(self, alpha):
        matrix = np.array([[np.cos(alpha), 0, np.sin(alpha)],
                  [0, 1, 0],
                  [-np.sin(alpha), 0, np.cos(alpha)]])
    
        self.orientation = matrix.dot(self.default_orientation)

    def heatFlux(self, box, T): #метод радиатион нужен чтобы вычислить, чтобы вычислить плотность теплового потока
        if box.parent.orbit.InShadow():
            return np.float64(0)
        
        dot = np.dot(self.orientation, box.orientation)
        if(dot > 0):
            return dot* self.q * box.parent.coat.As
        
        return np.float64(0)

class Radiaton(Load): #Простое излучение с пластины
    def _init_(self):
        pass
    
    def heatFlux(self, box, T):
        q = np.float64(-sigma*box.parent.coat.epsilon*T**4)
        return q
    
    def derivative(self, box, T):
        return -4*sigma*box.parent.coat.epsilon*T**3
    
class Isolated(Load): #изолированная сторона(вообще по большей части это placeholder)
    def __init__(self):
        pass
    
    def heatFlux(self, box, T):
        return np.float64(0)
    
class Conditions(): #класс для нагрузок действующих на спутник
    def __init__(self):
        self.external = []
        self.ethernal = []
        
    def addEx(self, *objects):
        for obj in objects:
            self.external.append(obj)
    
    def addEt(self, *objects):
        for obj in objects:
            self.ethernal.append(obj)

class Connect():
    def byNum(sputnik, num1, num2, R):
        sputnik.boxes[num1].connections.append(Connection(R, sputnik.boxes[num1],sputnik.boxes[num2]))
        sputnik.boxes[num2].connections.append(Connection(R, sputnik.boxes[num2],sputnik.boxes[num1]))
    
    def neighbours(sputnik, R):
        for box in sputnik.boxes:
            for neighbour in box.neighbours:
                box.connections.append(Connection(R ,box, neighbour))
    
    def byHeatPipe(sputnik, num1, num2, R, maxHeatFlux):
        sputnik.boxes[num1].connections.append(HeatPipe(R, sputnik.boxes[num1], sputnik.boxes[num2], maxHeatFlux))
        sputnik.boxes[num2].connections.append(HeatPipe(R, sputnik.boxes[num2], sputnik.boxes[num1], maxHeatFlux))

class Connection(Load):
    
    def __init__(self, R, box, connectedBox) -> None:
        self.R = R
        self.box = box
        self.connectedBox = connectedBox
    
    def heatFlux(self, T):
        n = self.connectedBox.prevIterT.size - 1
        q = (self.connectedBox.prevIterT[n]- T)/self.R
        
        return q
    
    def heat(self, T):
        Q = self.heatFlux(T) * self.box.area
        
        return Q
    
    def derivative(self):
        return -1/self.R

class EarthRadiation(Load):
    def __init__(self, q, radius) -> None:
       self.q = q
       self.radius = radius
       self.orientation = [-1, 0, 0]
    
    def heatFlux(self, box, T):
        dot = np.dot(self.orientation, box.orientation)
        
        if(dot > 0):
            q = box.coat.epsilon * dot * self.q * self.radius**2/box.parent.orbit.radius**2
            return q
        
        return 0

class VoidRadiation(Load):
    def __init__(self, T) -> None:
        self.T = T
    
    def heatFlux(self, box, T):
        q = box.coat.epsilon * sigma * self.T ** 4

        return q

class ConstantHeatFlux(Load):
    def __init__(self, q) -> None:
        self.q = q
    
    def heatFlux(self, box, T):
        return self.q

class EarthAlbedo(Load):
    def __init__(self, k, q) -> None:
        self.k = k
        self.q = q
        self.orientation = [-1, 0, 0]
    
    def heatFlux(self, box, T):
        if box.parent.orbit.InShadow():
            return np.float64(0)
        
        dot = np.dot(self.orientation, box.orientation)
        
        if dot > 0 and np.cos(box.parent.orbit.getAlpha()) > 0:
            q = dot * np.cos(box.parent.orbit.getAlpha()) * self.q * box.coat.As * self.k
            
            return q

        return 0

class HeatPipe(Connection):
    def __init__(self, R, box, connectedBox, maxHeatFlux) -> None:
        super().__init__(R, box, connectedBox)
        self.maxHeatFlux = maxHeatFlux
    
    def heatFlux(self, T):
        q = super().heatFlux(T)
        
        if q > self.maxHeatFlux:
            return self.maxHeatFlux
        
        if q < -self.maxHeatFlux:
            return -self.maxHeatFlux
        
        return q
        
        

#albedo 0.3, тень, солнце, солнечные панели