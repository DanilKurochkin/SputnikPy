import numpy as np
import numpy.typing as npt
import elemath as SMath
import workloads as wl
import time

#средневитковая температура

earth_mass = np.float64(5.976 * 10**24) # в килограммах
G = np.float64(6.67 * 10**(-11)) #гравитационная постоянная
earth_radius = np.float64(6371000) #в метрах
sigma = np.float64(5.67 * 10**(-8))

class Material(): #материал пластины

    def __init__(self, c: np.float64, p: np.float64, k: np.float64):
        self.c = c #теплоёмкость
        self.p = p #плотность
        self.k = k #температуропроводность
 
class Sputnik(): # спутник

    default_orientation = np.array([np.array([1, 0 , 0]),
                                    np.array([-1, 0, 0]),
                                    np.array([0, 1, 0]),
                                    np.array([0, -1, 0]),
                                    np.array([0, 0, 1]),
                                    np.array([0, 0, -1])]) #дефолтная ориентация пластин для куба
    
    def __init__(self, size : np.float64, width : np.float64, material : Material, coat : 'Coating', orbit : 'ClassicOrbit'):
        self.width = width #толщина спутника
        self.size = size #площадь пластинок
        self.boxes = np.empty(6, dtype=Box) 
        self.boxes : dict[Box, Box]
        self.coat = coat
        self.orbit = orbit
        
        self.externalConditions = []
        #инициалируем стенки спутника
        for i in np.arange(self.boxes.size):
            self.boxes[i] = Box(0, width, size, material, i, self.default_orientation[i], coat, self)
    
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
            for temper in box.T:
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
    
    def writeInnerEnergy(self, file, startT):
        format2 = '{0:13.3f} '
        for box in self.boxes:
            for i in np.arange(box.T.size):
                V = box.volumes[i].area * box.volumes[i].length
                dU = box.volumes[i].material.p * box.volumes[i].material.c * V * (startT-box.T[i])
                file.write(format2.format(dU))
            
            file.write('\n')
    
    def boxesNextT(self, ht, a0, b, c, d, a, P, Q):
        for i in np.arange(self.boxes.size):
            self.boxes[i].iterT = self.boxes[i].T
        disperancy = 1000
        
        new_disp = np.empty(self.boxes.size, dtype=np.float64)
        
        while disperancy > 10**(-3):
            for i in np.arange(self.boxes.size):
                self.boxes[i].prevIterT = self.boxes[i].iterT
            
            for i in np.arange(self.boxes.size):
                self.boxes[i].iter(ht, a0, b, c, d, a, P, Q)

            for i in np.arange(self.boxes.size):
                new_disp[i] = SMath.discrepancy(self.boxes[i].iterT, self.boxes[i].prevIterT)
            disperancy = np.max(new_disp)
        
        for i in np.arange(self.boxes.size):
            self.boxes[i].T = self.boxes[i].iterT
    
    def solve(self, amountOfRounds : int, pointsInRounds : int, save_every, startT = 300,filePath = 'output.txt', radiation_check = False, HeatCheckPath = 'outputheat.txt'): #решаем численно всё для всех пластинок в спутнике
        n = self.boxes[0].T.size

        a0 = np.empty(n, dtype=np.float64) #чтобы лишний раз много памяти не выделять
        b = np.empty(n, dtype=np.float64)
        c = np.empty(n, dtype=np.float64)
        d = np.empty(n, dtype=np.float64)
        a = np.empty(n, dtype=np.float64)
        P = np.empty(n, dtype=np.float64)
        Q = np.empty(n, dtype=np.float64)
        
        ht = self.orbit.period/pointsInRounds
        
        self.SetStartT(startT)
        
        file = open(filePath, 'w+')
        if radiation_check:
            file2 = open(HeatCheckPath, 'w+')
        
        for i in np.arange(amountOfRounds):
            for j in np.arange(pointsInRounds):
                self.boxesNextT(ht, a0, b, c, d, a, P, Q)
                    
                if radiation_check:
                    self.writeHeat(file2, ht) 

                if j % save_every == 0:
                    self.writeResult(file, ht, i, j)
                
                self.orbit.Move(ht, self)
                for externalConditon in self.externalConditions:
                    externalConditon.rotate(self.orbit.getAlpha())
        
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
    
    def addCondition(self, condition : wl.Conditions): #добавление граничного условия для всех пластин, то есть того, которое будет действовать на все пластины
        for box in self.boxes:
            for ethernal in condition.ethernal:
                box.conditions.addEt(ethernal)
                
            for external in condition.external:
                box.conditions.addEx(external)
        
        for external in condition.external:
            self.externalConditions.append(external)
    
    def addConditionByNum(self, condition : wl.Conditions, *nums): # аналогично выше описанному методу, но добавляет только для конкретных пластин
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

        self.conditions = wl.Conditions()
        self.conditions : wl.Conditions
        self.T = []
        self.iterT = []
        self.prevIterT = []
        self.contactT = []
        self.h = []
        self.neighbours = []
        self.connections = []
        self.connections : list[wl.Connection]
        self.neighbours : dict[Box, Box]

    def createVolumes(self, h : np.float64): #нарезам всё на конечные объёмы
        self.h = h
        n = (np.arange(0, self.length+h/2, h)).size
        self.volumes = np.empty(n, dtype=FiniteVolume)
        self.volumes : dict[FiniteVolume, FiniteVolume]
        self.T = np.empty(n, dtype=np.float64)
        
        x = self.x
        self.volumes[0] = FiniteVolume(x, h/2, self.area, self.material , self)
        for i in range(1, n-1):
            x += self.h
            self.volumes[i] = FiniteVolume(x, h, self.area, self.material , self)
        self.volumes[n-1] = FiniteVolume(x + self.h, h/2, self.area, self.material , self)
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
        a[0], d[0] = wl.BoundaryCondition.FindCoefsEt(self, self.conditions.external, self.iterT[0], b[0])
        FindKoef(self, a0, b, c, d, a)
        size = a0.size
        b[size-1] = 0 #Так всегда, аналогично тому что выше
        c[size-1] = self.volumes[size-1].leftNeighbour.material.k/self.volumes[size-1].parent.h
        a[size-1], d[size-1] = wl.BoundaryCondition.FindCoefsEx(self, self.conditions.ethernal, self.iterT[size-1], c[size-1])
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
    
class Coating(): #покрытие спутника
    def __init__(self, As, epsilon):
        self.As = As #коэфициент поглащения солнечного излучения, то есть высокочастотного
        self.epsilon = epsilon #степень черноты, а также коэфициент поглащения в низкочастотном спектре

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

class HelioStationarOrbit(ClassicOrbit):
    def __init__(self, radiusAboveEarth, localTime):
        super().__init__(radiusAboveEarth)
        res = -np.pi
        
    def Move(self, ht, sputnik: Sputnik):
        pass
        
def main():
    
    A6061 = Material(800, 2700, 202) #алюминий 6061
    coat = Coating(0.9, 0.91) #ЭКОМ-1П a05 e04
    orbit = SunLookingOrbit(500)
    cond = wl.Conditions()
    
    sp = Sputnik(4, 0.002, A6061, coat, orbit)
    sp.createVolumes(0.0002)
    sp.knitPlates()
    
    cond.addEx(wl.Radiaton(), wl.EarthRadiation(240, earth_radius), wl.Sun(1500), wl.EarthAlbedo(0.3, 1500))
    cond.addEx()
    cond.addEt(wl.ConstantHeatFlux(200))
    sp.addCondition(cond)

    wl.Connect.neighbours(sp, 0.5)

    start = time.time()
    sp.solve(30, 100, 1, 290, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()