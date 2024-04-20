import numpy as np
from classes.material import Coating, Material
import classes.elemath as SMath
import numpy.typing as npt
from typing import List

class Sputnik(): # спутник

    default_orientation = np.array([np.array([1, 0 , 0]),
                                    np.array([-1, 0, 0]),
                                    np.array([0, 1, 0]),
                                    np.array([0, -1, 0]),
                                    np.array([0, 0, 1]),
                                    np.array([0, 0, -1])]) #дефолтная ориентация пластин для куба
    
    def __init__(self, Lx, Ly, Lz, width : np.float64, material : Material, coat : 'Coating', orbit):
        self.width = width #толщина спутника
        self.size = np.array([Ly*Lz, Ly*Lz,
                                Lx*Lz, Lx*Lz,
                                Lx*Ly, Lx*Ly]) #площадь пластинок
        self.boxes : List[Box] = []
        self.coat = None
        self.orbit = orbit
        self.boxes
        self.ht = None
        self.time = None
        
        self.externalConditions = []
        #инициалируем стенки спутника
        for i in np.arange(6):
            self.boxes.append(Box(0, width, self.size[i], material, i, self.default_orientation[i], coat, self))
        self._knitPlates()
        
    def _knitPlates(self): #связываем пластины, чтобы знать какая с какой соприкасается
        for box in self.boxes:
            box.neighbours = []
        
        for i in np.arange(len(self.boxes)):
            z = 0
            for j in np.arange(len(self.boxes)):
                if i != j:
                    if np.dot(self.boxes[i].orientation, self.boxes[j].orientation) == 0:
                        self.boxes : list[Box]
                        self.boxes[i].neighbours.append(self.boxes[j])
                        z += 1

    def createVolumes(self, n : int): #нарезаем пластины спутника на конечные объёмы
        for box in self.boxes:
            box.createVolumes(n)
    
    def writeResult(self, file, ht, i, j): #запись результата в отдельный файл, результат это распределение температур внутри пластины, угол на котором находится объект
        format1 = '{0} {1} {2} {3}\n'
        format2 = '{0:13.3f} '
        j += 1
        file.write(format1.format(i, j, ht*j+i*self.orbit.period, self.orbit.getAlpha()))
        for box in self.boxes:
            file.write('\t')
            for temper in box.T:
                file.write(format2.format(temper))
            file.write('\n')
    
    def writeHeat(self, file, ht): #запись полученного/потерянного тепла/секунду в отдельный файл, учитывая что шаг по времени постойнный, с помощью этого можно осуществлять проверку по сохранению энергии
        format2 = '{0:16.3f} '
        for box in self.boxes:
            for external in box.conditions.external:
                value = external.heat(box, box.T[0])
                file.write(format2.format(value * ht))
                if value > 0:
                    box.integr.append(value)
            
            for ethernal in box.conditions.ethernal:
                value = ethernal.heat(box, box.T[0])
                file.write(format2.format(value * ht))
                if (not ethernal.is_connection()) and value > 0:
                    box.integr.append(value)
            file.write('\n')    
                
        file.write('\n\n')
    
    def writeInnerEnergy(self, file, startT):
        format2 = '{0:13.3f} '
        for box in self.boxes:
            for i in np.arange(box.T.size):
                V = box.volumes[i].area * box.volumes[i].length
                dU = box.volumes[i].material.p * box.volumes[i].material.c * V * (startT-box.T[i])
                file.write(format2.format(dU))
            
            file.write('\n')
    
    def boxesNextT(self, ht, a0, b, c, d, a, P, Q):
        for i in np.arange(len(self.boxes)):
            self.boxes[i].iterT = self.boxes[i].T
        disperancy = 1000
        
        new_disp = np.empty(len(self.boxes), dtype=np.float64)
        
        while disperancy > 10**(-2):
            for i in np.arange(len(self.boxes)):
                self.boxes[i].prevIterT = self.boxes[i].iterT
            
            for i in np.arange(len(self.boxes)):
                self.boxes[i].iter(ht, a0, b, c, d, a, P, Q)

            for i in np.arange(len(self.boxes)):
                new_disp[i] = SMath.discrepancy(self.boxes[i].iterT, self.boxes[i].prevIterT)
            disperancy = np.max(new_disp)
        
        for i in np.arange(len(self.boxes)):
            self.boxes[i].T = self.boxes[i].iterT
    
    def solve(
        self,
        amountOfRounds : int,
        pointsInRounds : int,
        save_every,
        startT = 300,
        filePath = 'result/output.txt',
        radiation_check = False,
        HeatCheckPath = 'result/outputheat.txt'
        ): #решаем численно всё для всех пластинок в спутнике
        
        n = self.boxes[0].T.size

        a0 = np.empty(n, dtype=np.float64) #чтобы лишний раз много памяти не выделять
        b = np.empty(n, dtype=np.float64)
        c = np.empty(n, dtype=np.float64)
        d = np.empty(n, dtype=np.float64)
        a = np.empty(n, dtype=np.float64)
        P = np.empty(n, dtype=np.float64)
        Q = np.empty(n, dtype=np.float64)
        
        self.ht = self.orbit.period/pointsInRounds
        self.time = 0
        self.SetStartT(startT)
        
        file = open(filePath, 'w+')
        if radiation_check:
            file2 = open(HeatCheckPath, 'w+')
            
        file.write(f'{pointsInRounds}\n')
        for i in np.arange(amountOfRounds):
            for j in np.arange(pointsInRounds):
                print(f"Solving {i} round {j} point")
                self.boxesNextT(self.ht, a0, b, c, d, a, P, Q)
                    
                if radiation_check:
                    self.writeHeat(file2, self.ht) 
                if j % save_every == 0:
                    self.writeResult(file, self.ht, i, j)
                
                self.orbit.Move(self.ht, self)
                for externalConditon in self.externalConditions:
                    externalConditon.rotate(self.orbit.getAlpha())
                self.time += self.ht
        
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
    
    def addCondition(self, condition : 'Conditions'): #добавление граничного условия для всех пластин, то есть того, которое будет действовать на все пластины
        for box in self.boxes:
            for ethernal in condition.ethernal:
                box.conditions.addEt(ethernal)
                
            for external in condition.external:
                box.conditions.addEx(external)
        
        for external in condition.external:
            self.externalConditions.append(external)
    
    def addConditionByNum(self, condition : 'Conditions', *nums): # аналогично выше описанному методу, но добавляет только для конкретных пластин
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
        self.coat : Coating = coat
        self.parent : Sputnik = parent
        self.conditions : Conditions = Conditions()
        self.T = []
        self.iterT = []
        self.prevIterT = []
        self.contactT = []
        self.neighbours : list[Box] = []
        self.integr = []

    def createVolumes(self, n : int): #нарезам всё на конечные объёмы
        h = self.length/(n-1)
        self.volumes : list[FiniteVolume] = []
        self.T = np.empty(n, dtype=np.float64)

        x = self.x
        self.volumes.append(FiniteVolume(x, h/2, self.area, self.material, self))
        for i in range(1, n-1):
            self.volumes.append(FiniteVolume(x + h*i, h, self.area, self.material , self))
        self.volumes.append(FiniteVolume(self.x+self.length, h/2, self.area, self.material , self))
        self._knitVolumes()

    def _knitVolumes(self): #связываем конечные объёмы
        n = len(self.volumes)
        self.volumes[0].knit([], self.volumes[1], True, False)
        for i in np.arange(1, n-1):
            self.volumes[i].knit(self.volumes[i-1], self.volumes[i+1], False, False)
        self.volumes[n-1].knit(self.volumes[n-2], [], False, True)

    def iter(self, ht, a0, b, c, d, a, P, Q): #итерация
        def FindKoef(self : Box ,a0 : npt.NDArray[np.float64], b : npt.NDArray[np.float64], c : npt.NDArray[np.float64], d : npt.NDArray[np.float64] , a : npt.NDArray[np.float64]):
            for  i in np.arange(1, a0.size-1):
                volume : FiniteVolume = self.volumes[i]
                a0[i] = volume.material.p*volume.material.c*volume.length/ht
                b[i] = volume.rightNeighbour.material.k/volume.distantW
                c[i] = volume.leftNeighbour.material.k/volume.distantE
                d[i] = a0[i]*self.T[i]
                a[i] = b[i] + c[i] + a0[i]
        
        #вычисляем коэициенты для вектора температур self.T
        c[0] = 0 #Так всегда, из-за того что матрица трёхдиагональная
        volume = self.volumes[0]
        a0[0] = volume.material.p*volume.material.c*volume.length/ht
        b[0] = volume.rightNeighbour.material.k/volume.distantE
        fp, fc = BoundaryCondition.linearize(self, self.conditions.external, self.iterT[0])
        a[0] = a0[0] + b[0] - fp
        d[0] = a0[0]*self.T[0] + fc
        
        FindKoef(self, a0, b, c, d, a)
        
        b[-1] = 0 #Так всегда, из-за того что матрица трёхдиагональная
        volume = self.volumes[-1]
        a0[-1] = volume.material.p*volume.material.c*volume.length/ht
        c[-1] = volume.leftNeighbour.material.k/volume.distantW
        fp, fc = BoundaryCondition.linearize(self, self.conditions.ethernal, self.iterT[-1])
        a[-1] = a0[-1] + c[-1] - fp
        d[-1] = a0[-1]*self.T[-1] + fc
        #вычислили коэфициенты 
        
        self.iterT = SMath.TDMA(a, b, c, d, P, Q) # решили методом простой прогонки и обновили вектор температур
    
class FiniteVolume(): #конечный объём и его характеристики

    def __init__(self, x: np.float64, length : np.float64, area : np.float64, material: Material, parent: Box): 
        self.x = x
        self.length = length
        self.material = material
        self.parent = parent
        self.distantE = None
        self.distantW = None
        self.area = area

    def knit(self, leftNeightbour : 'FiniteVolume', rightNeighbour : 'FiniteVolume', onLeftEdge: bool, onRightEdge: bool):
        #связываем соседние объёмы чтобы проще было к ним обраться в решателе
        self.leftNeighbour = leftNeightbour
        self.rightNeighbour = rightNeighbour
        self.onLeftEdge = onLeftEdge
        self.onRightEdge = onRightEdge
        
        if self.onLeftEdge:
            self.distantW = None
        else:
            self.distantW = np.abs(self.x - self.leftNeighbour.x)
        
        if self.onRightEdge:
            self.distantE = None
        else:
            self.distantE = np.abs(self.x - self.rightNeighbour.x)
        

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
    def integrateTime(condtion, box, time, dtime):
        result = condtion.heatFlux(box, time + dtime) + condtion.heatFlux(box, time)
        result /= 2
        
        return result
    
    def linearize(box : Box, conditions, T):
        fT = 0
        fp = 0
        for cond in conditions:
            if not cond.timeDependent():
                fT += cond.heatFlux(box, T)    
                fp += cond.derivative(box, T)
            else:
                fT += BoundaryCondition.integrateTime(cond, box, box.parent.time, box.parent.ht)
        
        fc = fT - fp * T
        
        return fp, fc