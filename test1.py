from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
mat = Material(c=800, p=2700, k=100) #материал спутника
coat = Coating(As=0.6, epsilon=0.4) #покрытие
orbit = SunLookingOrbit(radiusAboveEarth=500) #орбита

#создание спутника
sp = Sputnik(Lx=1,Ly=1, Lz=1, width=0.005, material=mat, coat=coat, orbit=orbit)  
#задание КО в пластинах
sp.createVolumes(10)

#Создание и подключение граничных условий
cond = Conditions()
cond.addEx(wl.Sun(1367), wl.Radiaton())
sp.addCondition(cond)

#связывание пластинок
wl.Connect.neighbours(sputnik=sp, R=0.01)
#моделирование

start = time.time()
sp.solve(30, 100, 1, 300, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.radiationCheck() * 100)
print(Checker.averageTCheck(sp.size, 100000))
print(f'Теоретическое:290.89')

import graph