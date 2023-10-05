import numpy as np
from classes.material import Material, Coating
import classes.elemath as SMath
import numpy.typing as npt

class PlateCollector():
        def __init__(self ,sputnik : 'Sputnik'):
            self.zenit = sputnik.boxes[0]
            self.zenit : 'Box'
            self.nadir = sputnik.boxes[1]
            self.nadir : 'Box'
            self.ortogonal1 = sputnik.boxes[2]
            self.ortogonal1 : 'Box'
            self.ortogonal2 = sputnik.boxes[3]
            self.ortogonal2 : 'Box'
            self.shade = sputnik.boxes[4]
            self.shade : 'Box'
            self.solar = sputnik.boxes[5]
            self.solar : 'Box'

class Sputnik(): # спутник
    
    default_orientation = np.array([np.array([1, 0 , 0]),
                                    np.array([-1, 0, 0]),
                                    np.array([0, 1, 0]),
                                    np.array([0, -1, 0]),
                                    np.array([0, 0, 1]),
                                    np.array([0, 0, -1])]) #дефолтная ориентация пластин для куба
    
    def __init__(self, Lx : np.float64, Ly : np.float64, Lz : np.float64, width : np.float64, material, coat, orbit):
        self.width = width #толщина спутника
        self.size = [Ly*Lz, Ly*Lz,
                     Lx*Lz, Lx*Lz,
                     Lx*Ly, Lx*Ly] #площадь пластинок
        self.boxes = np.empty(6, dtype=Box) 
        self.boxes : dict[Box, Box]
        self.coat = coat
        self.coat : 'Coating'
        self.orbit = orbit
        self.time = np.float64(0)
        self.externalConditions = []
        self.dtime = None
        #инициалируем стенки спутника
        for i in np.arange(self.boxes.size):
            self.boxes[i] = Box(0, width, self.size[i], material, i, self.default_orientation[i], coat, self)
        self.plate = PlateCollector(self)
        self.plate : PlateCollector
        
    
    def knitPlates(self): #связываем пластины, чтобы знать какая с какой соприкасается
        for box in self.boxes:
            box : Box
            box.neighbours = np.empty(4, dtype=Box)
        
        for i in np.arange(self.boxes.size):
            z = 0
            for j in np.arange(self.boxes.size):
                if i != j:
                    if np.dot(self.boxes[i].orientation, self.boxes[j].orientation) == 0:
                        self.boxes : dict[Box, Box]
                        self.boxes[i].neighbours[z] = self.boxes[j]
                        z += 1

    def createVolumes(self, h : np.float64): #нарезаем пластины спутника на конечные объёмы
        for box in self.boxes:
            box : Box
            box.createVolumes(h)
    
    def writeResult(self, file, ht, i, j): #запись результата в отдельный файл, результат это распределение температур внутри пластины, угол на котором находится объект
        format1 = '{0} {1} {2} {3}\n'
        format2 = '{0:13.3f} '
        file.write(format1.format(i, j, ht*j+i*self.orbit.period, self.orbit.getAlpha()))
        for box in self.boxes:
            file.write('\t')
            for temper in box.T[1 : box.T.size - 1]:
                file.write(format2.format(temper))
            file.write('\n')
    
    def writeHeat(self, file, ht): #запись полученного/потерянного тепла/секунду в отдельный файл, учитывая что шаг по времени постойнный, с помощью этого можно осуществлять проверку по сохранению энергии
        format2 = '{0:16.3f} '
        for box in self.boxes:
            for external in box.conditions.external:
                file.write(format2.format(external.heat(box, box.T[0]) * ht))
            
            for ethernal in box.conditions.ethernal:
                file.write(format2.format(ethernal.heat(box, box.T[box.T.size-1]) * ht))
            
            for connect in box.connections:
                file.write(format2.format(connect.heat(box.T[box.T.size-1]) * ht))
            file.write('\n')

        file.write('\n')
    
    def writeInnerEnergy(self, file, startT): #записывает итоговое изменение внутренней энергии
        format2 = '{0:13.3f} '
        for box in self.boxes:
            for i in np.arange(1, box.T.size - 1):
                V = box.volumes[i].area * box.volumes[i].length
                dU = box.volumes[i].material.p * box.volumes[i].material.c * V * (startT-box.T[i])
                file.write(format2.format(dU))
            
            file.write('\n')
    
    def boxesNextT(self, ht, a0, b, c, d, a, P, Q): # для улучшения сходимости итерируемся по пластинам последовательно
        for i in np.arange(self.boxes.size):
            self.boxes[i].iterT = self.boxes[i].T
        disperancy = 1000
        
        new_disp = np.empty(self.boxes.size, dtype=np.float64)

        while disperancy > 10**(-3): #вычисляем с заданной точностью
            for i in np.arange(self.boxes.size):
                self.boxes[i].prevIterT = self.boxes[i].iterT
            
            for i in np.arange(self.boxes.size):
                self.boxes[i].iter(ht, a0, b, c, d, a, P, Q)

            for i in np.arange(self.boxes.size):
                new_disp[i] = SMath.discrepancy(self.boxes[i].iterT, self.boxes[i].prevIterT)
            disperancy = np.max(new_disp)
        
        for i in np.arange(self.boxes.size):
            self.boxes[i].T = self.boxes[i].iterT 
    
    def solve(self, amountOfRounds : int, pointsInRounds : int, save_every : int, startT = 300,filePath = 'output.txt', radiation_check = False, HeatCheckPath = 'outputheat.txt'): #решаем численно всё для всех пластинок в спутнике
        n = self.boxes[0].T.size

        a0 = np.empty(n, dtype=np.float64) #чтобы лишний раз много памяти не выделять
        b = np.empty(n, dtype=np.float64)
        c = np.empty(n, dtype=np.float64)
        d = np.empty(n, dtype=np.float64)
        a = np.empty(n, dtype=np.float64)
        P = np.empty(n, dtype=np.float64)
        Q = np.empty(n, dtype=np.float64)
        
        self.dtime = self.orbit.period/pointsInRounds
        
        self.SetStartT(startT)
        
        file = open(filePath, 'w+')
        if radiation_check:
            file2 = open(HeatCheckPath, 'w+')
        
        for i in np.arange(amountOfRounds):
            for j in np.arange(pointsInRounds):
                self.boxesNextT(self.dtime, a0, b, c, d, a, P, Q)
                    
                if radiation_check:
                    self.writeHeat(file2, self.dtime) 

                if j % save_every == 0:
                    self.writeResult(file, self.dtime, i, j)
                
                self.orbit.Move(self.dtime, self)
                for externalConditon in self.externalConditions:
                    externalConditon.rotate(self.orbit.getAlpha())
                
                self.time += self.dtime
        
        if radiation_check:
            self.writeInnerEnergy(file2, startT)
                    
    def newContactT(self): # мы усредняем температуру по каждой пластине, чтобы понять как пойдёт теплообмен в соединёненных пластинах        
        for box in self.boxes:
            box.contactT = np.mean(box.T)
        
    def SetStartT(self, startT): # для численных методов необходима стартовая температура, здесь мы её устанавливаем
        for box in self.boxes:
            box : Box
            for i in np.arange(box.T.size):
                box.T[i] = startT
        self.newContactT()
    
    def addCondition(self, condition): #добавление граничного условия для всех пластин, то есть того, которое будет действовать на все пластины
        for box in self.boxes:
            for ethernal in condition.ethernal:
                box.conditions.addEt(ethernal)
                
            for external in condition.external:
                box.conditions.addEx(external)
        
        for external in condition.external:
            self.externalConditions.append(external)
    
    def addConditionByNum(self, condition, *nums): # аналогично выше описанному методу, но добавляет только для конкретных пластин
        for num in nums:
            for ethernal in condition.ethernal:
                self.boxes[num].conditions.addEt(ethernal)
            
            for external in condition.external:
                self.boxes[num].conditions.addEx(external)
                
        for external in condition.external:
            self.externalConditions.append(external)
    
class Box(): #родная коробочка
    
    def __init__(self, x: np.float64, length: np.float64, area : np.float64, material: Material, number, orientation, coat : 'Coating', parent : Sputnik ):
        self.length = length
        self.material = material
        self.x = x
        self.area = area
        self.number = number
        self.orientation = orientation
        self.coat = coat
        self.parent = parent

        self.conditions = Conditions()
        self.T = []
        self.iterT = []
        self.prevIterT = []
        self.contactT = []
        self.h = []
        self.neighbours = []
        self.connections = []
        self.neighbours : dict[Box, Box]

    def createVolumes(self, amount : int): #нарезам всё на конечные объёмы
        self.h = self.length / amount
        n = amount + 2
        self.volumes = np.empty(n, dtype=FiniteVolume)
        self.volumes : dict[FiniteVolume, FiniteVolume]
        self.T = np.empty(n, dtype=np.float64)
        
        x = self.x - self.h/2
        self.volumes[0] = FiniteVolume(x, self.h/2, self.area, self.material , self)
        for i in range(1, n-1):
            x += self.h
            self.volumes[i] = FiniteVolume(x, self.h, self.area, self.material , self)
        self.volumes[n-1] = FiniteVolume(x + self.h, self.h/2, self.area, self.material , self)
        Box.knitVolumes(self)

    def knitVolumes(self): #связываем конечные объёмы
        n = self.volumes.size
        self.volumes[0].knit([], self.volumes[1], True, False)
        for i in np.arange(1, n-1):
            self.volumes[i].knit(self.volumes[i-1], self.volumes[i+1], False, False)
        self.volumes[n-1].knit(self.volumes[n-2], [], False, True)

    def iter(self, ht, a0, b, c, d, a, P, Q): #итерация
        def FindKoef(self : Box ,a0 : npt.NDArray[np.float64], b : npt.NDArray[np.float64], c : npt.NDArray[np.float64], d : npt.NDArray[np.float64] , a : npt.NDArray[np.float64]):
            for  i in np.arange(1, a0.size-1):
                a0[i] = self.volumes[i].material.p*self.volumes[i].material.c*self.volumes[i].length/ht
                b[i] = self.volumes[i].rightNeighbour.material.k/self.h
                c[i] = self.volumes[i].leftNeighbour.material.k/self.h
                d[i] = a0[i]*self.T[i]
                a[i] = b[i] + c[i] + a0[i]
        
        #вычисляем коэициенты для вектора температур self.T
        c[0] = 0 #Так всегда, из-за того что матрица трёхдиагональная
        b[0] = self.volumes[0].rightNeighbour.material.k/self.volumes[0].parent.h
        a[0], d[0] = BoundaryCondition.FindCoefsEt(self, self.conditions.external, self.iterT[0], b[0])
        FindKoef(self, a0, b, c, d, a)
        size = a0.size
        b[size-1] = 0 #Так всегда, аналогично тому что выше
        c[size-1] = self.volumes[size-1].leftNeighbour.material.k/self.volumes[size-1].parent.h
        a[size-1], d[size-1] = BoundaryCondition.FindCoefsEx(self, self.conditions.ethernal, self.iterT[size-1], c[size-1])
        #вычислили коэфициенты 
        
        self.iterT = SMath.TDMA(a, b, c, d, P, Q) # решили методом простой прогонки и обновили вектор температур
    
class FiniteVolume(): #конечный объём и его характеристики

    def __init__(self, x: np.float64, length : np.float64, area : np.float64, material: Material, parent: Box): 
        self.x = x
        self.length = length
        self.material = material
        self.parent = parent
        self.area = area

    def knit(self, leftNeightbour : 'FiniteVolume', rightNeighbour : 'FiniteVolume', onLeftEdge: bool, onRightEdge: bool):
        #связываем соседние объёмы чтобы проще было к ним обраться в решателе
        self.leftNeighbour = leftNeightbour
        self.rightNeighbour = rightNeighbour
        self.onLeftEdge = onLeftEdge
        self.onRightEdge = onRightEdge

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
            
class BoundaryCondition(): #граничное условие

    def ConnectionCoefs(box, T):
        fT = 0
        fp = 0
        for connection in box.connections: #линеаризация
            fT += connection.heatFlux(T)
            fp += connection.derivative()
        
        return fT, fp
    
    def FindCoefsEx(box, conditions, T, a1): #вычисление коэфициентов для граничных условий
        fT = 0
        fp = 0
        for cond in conditions: #линеаризация
            fT += cond.heatFlux(box, T)    
            fp += cond.derivative(box, T)
        fc = fT - fp * T
        
        a = a1 - fp
        d = fc
        return a, d
    
    def FindCoefsEt(box, conditions, T,a1): #вычисление коэфициентов для граничных условий
        fT = 0
        fp = 0
        for cond in conditions: #линеаризация
            fT += cond.heatFlux(box, T)
            fp += cond.derivative(box, T)
        res = BoundaryCondition.ConnectionCoefs(box, T)
        fT += res[0]
        fp += res[1]
        
        fc = fT - fp * T
        
        a = a1 -fp
        d = fc
        return a, d