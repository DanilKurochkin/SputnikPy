from abc import ABC
import numpy as np
from classes.models import Box, Sputnik

sigma = np.float64(5.67*10**(-8))

class Load(ABC): #Базовый функционал для любой нагрузки
    def __init__(self) -> None:
        pass

    def rotate(self, alpha): #возможность вращения, определяется наследником
        pass

    def heatFlux(self, box : Box, T): #тепловой поток, определяется наследником
        pass
    
    def heat(self, box : Box, T): #общий для всех наследников метод вычисления количества теплоты попадающего на пластину за секунду
        if not self.timeDependent():
            Q = self.heatFlux(box, T) * box.area
        else:
            Q = self.heatFlux(box, box.parent.time) * box.area
        
        return Q
    
    def timeDependent(self):
        return False
    
    def derivative(self, box : Box, T): #Производная теплового потока потока по температуре, определяется наследником
        return 0

class Sun(Load): #Солнце
    
    def __init__(self, q):
        self.q = np.float64(q)
        self.default_orientation = np.array([1, 0, 0]) #По умолчанию мы считаем, что спутник в 0 момент времени находится между Землёй и Солнцем
        self.orientation = []
        self.rotate(0)
    
    def rotate(self, alpha): #Вращаем вектор солнечного излучения с помощью матрицы вращения на соответствующий угол
        matrix = np.array([[np.cos(alpha), 0, np.sin(alpha)],
                  [0, 1, 0],
                  [-np.sin(alpha), 0, np.cos(alpha)]])
    
        self.orientation = matrix.dot(self.default_orientation)

    def heatFlux(self, box, T): #вычислям тепловой поток с помощью скалярного произведения векторов нормали поверхности и вектора солнечного излучения
        if box.parent.orbit.InShadow(): #Обращаемся к орбите спутника, находится ли он в тени. Если в тени, неожиданно, тепловой поток от солнца равен нулю
            return np.float64(0)
        
        dot = np.dot(self.orientation, box.orientation)
        if(dot > 0):
            return dot* self.q * box.parent.coat.As
        
        return np.float64(0)

class Radiaton(Load): #Тепловое излучения с пластины
    def _init_(self):
        pass
    
    def heatFlux(self, box, T): #Закон Стефана-Больцмана
        q = np.float64(-sigma*box.parent.coat.epsilon*T**4)
        return q
    
    def derivative(self, box, T):
        return -4*sigma*box.parent.coat.epsilon*T**3
    
class Isolated(Load): #изолированная сторона
    def __init__(self):
        pass
    
    def heatFlux(self, box, T):
        return np.float64(0)

class Connect(): #Служебный класс для соединения пластин
    def byNum(sputnik : 'Sputnik', num1, num2, R):
        sputnik.boxes[num1].conditions.addEt(Connection(R, sputnik.boxes[num1],sputnik.boxes[num2]))
        sputnik.boxes[num2].conditions.addEt(Connection(R, sputnik.boxes[num2],sputnik.boxes[num1]))
    
    def byBox(box1 : 'Box', box2 : 'Box', R):
        if(box1.parent is not box2.parent):
            raise ValueError
        box1.conditions.addEt(Connect(R, box1, box2))
        box1.conditions.addEt(Connect(R, box2, box1))
    
    def neighbours(sputnik : Sputnik, R):
        for box in sputnik.boxes:
            for neighbour in box.neighbours:
                box.conditions.addEt(Connection(R ,box, neighbour))
    
    def byHeatPipeByNum(sputnik : Sputnik, num1, num2, R, maxHeatFlux):
        sputnik.boxes[num1].conditions.addEt(HeatPipe(R, sputnik.boxes[num1], sputnik.boxes[num2], maxHeatFlux))
        sputnik.boxes[num2].conditions.addEt(HeatPipe(R, sputnik.boxes[num2], sputnik.boxes[num1], maxHeatFlux))
    
    def byHeatPipeByBox(box1 : 'Box', box2 : 'Box', R):
        if(box1.parent is not box2.parent):
            raise ValueError
        box1.conditions.addEt(HeatPipe(R, box1, box2))
        box1.conditions.addEt(HeatPipe(R, box2, box1))

class Connection(Load): #Класс термического контакта
    
    def __init__(self, R, box, connectedBox) -> None:
        self.R = R
        self.box = box
        self.connectedBox = connectedBox
    
    def heatFlux(self, box, T): #Температура берётся с предыдущий итерации, так как если брать с текущей, то закон сохранения нарушится
        n = self.connectedBox.prevIterT.size - 1
        q = (self.connectedBox.prevIterT[n]- T)/self.R * (self.connectedBox.area / self.box.area)
        
        return q
    
    def heat(self, box, T):
        Q = self.heatFlux(box, T) * self.box.area
        
        return Q
    
    def derivative(self, box, t):
        return -1/self.R

class EarthRadiation(Load): #Земное излучение
    def __init__(self, q, radius) -> None:
       self.q = q
       self.radius = radius
       self.orientation = [-1, 0, 0] #Поток по умолчанию направлен в надир
    
    def heatFlux(self, box, T):
        dot = np.dot(self.orientation, box.orientation)
        
        if(dot > 0): #аналогично Солнечному произведению, вычисление производится с помощью скалярного произведения
            q = box.coat.epsilon * dot * self.q * self.radius**2/box.parent.orbit.radius**2
            return q
        
        return 0

class VoidRadiation(Load): #Излучение пустоты
    def __init__(self, T) -> None:
        self.T = T
    
    def heatFlux(self, box, T): #Стефан-Больцман
        q = box.coat.epsilon * sigma * self.T ** 4

        return q

class ConstantHeatFlux(Load): #Постоянный тепловой поток
    def __init__(self, q) -> None:
        self.q = q
    
    def heatFlux(self, box, T):
        return self.q

class EarthAlbedo(Load): #Солнечное излучение отраженной от земной поверхности
    def __init__(self, k, q) -> None:
        self.k = k
        self.q = q
        self.orientation = [-1, 0, 0]
    
    def heatFlux(self, box, T): #Аналогично тому что происходит для солнечного излучения, но добавляется ограничение, так как земная поверхность должна быть не в тени
        if box.parent.orbit.InShadow():
            return np.float64(0)
        
        dot = np.dot(self.orientation, box.orientation)
        
        if dot > 0 and np.cos(box.parent.orbit.getAlpha()) > 0:
            q = dot * np.cos(box.parent.orbit.getAlpha()) * self.q * box.coat.As * self.k
            
            return q

        return 0

class HeatPipe(Connection): #Тепловая трубка
    def __init__(self, R, box, connectedBox, maxHeatFlux) -> None:
        super().__init__(R, box, connectedBox)
        self.maxHeatFlux = maxHeatFlux
    
    def heatFlux(self, T): #ограничение по максимально тепловому потоку появляется из-за физических свойств тепловой трубки
        q = super().heatFlux(T)
        
        if q > self.maxHeatFlux:
            return self.maxHeatFlux
        
        if q < -self.maxHeatFlux:
            return -self.maxHeatFlux
        
        return q
    
class TableFunctionLoadFlux(Load): #Тепловой поток задаваемый табличной функцией
    def __init__(self, heatPoints, timePoints): 
        self.heatPoints = np.array(heatPoints)
        self.timePoints = np.array(timePoints)
    
    def timeDependent(self):
        return True
    
    def heatFlux(self, box, t):
        result = np.interp(t, self.timePoints, self.heatPoints)
        return result

class TableFunctionPeriodicalFlux(TableFunctionLoadFlux): # таблица закцикленная периодечки сама на себя
    def __init__(self, heatPoints, timePoints, period):
        super().__init__(heatPoints, timePoints)
        self.period = period
    
    def heatFlux(self, box, t):
        t = t % self.period
        result = super().heatFlux(box, t)
        return result        

class TableFunctionLoadHeat(TableFunctionLoadFlux): # тоже самое но для теплоты в секунду
    def __init__(self, heatPoints, timePoints):
        super().__init__(heatPoints, timePoints)
    
    def heatFlux(self, box, t):
        return super().heatFlux(box, t)/box.area

class TableFunctionPeriodicalHeat(TableFunctionPeriodicalFlux): # тоже самое но для теплоты в секунду
    def __init__(self, heatPoints, timePoints, period):
        super().__init__(heatPoints, timePoints, period)
    
    def heatFlux(self, box, t):
        return super().heatFlux(box, t)/box.area