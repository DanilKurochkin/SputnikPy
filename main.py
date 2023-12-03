from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
mat = Material(c=800, p=2700, k=100)
coat = Coating(As=0.1, epsilon=0.1)
orbit = SunLookingOrbit(radiusAboveEarth=500)

sp = Sputnik(Lx=1,Ly=1, Lz=1, width=0.005, material=mat, coat=coat, orbit=orbit)  
sp.createVolumes(10)

sp.boxes[0].coat = Coating(As=0.1, epsilon=0.9) # ЭКОМ-ж2
sp.boxes[5].coat = Coating(As=0.1, epsilon=0.9)

cond = Conditions()
cond.addEx(wl.Sun(1367), wl.Radiaton())
cond.addEt(wl.TableFunctionLoadHeat([50, 50], [0, 10]))

cond2 = Conditions()
cond2.addEt(wl.TableFunctionPeriodicalHeat([0,0,300,300,0,0], 
                                           [0, orbit.period/4, orbit.period/4 + 1, orbit.period/4 + 1000, orbit.period/4 + 1001, orbit.period],
                                           orbit.period))

sp.addCondition(cond)
sp.addConditionByNum(cond2, 0)
wl.Connect.neighbours(sp, 0.001)

start = time.time()
sp.solve(30, 100, 1, 314, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.averageTCheck(sp.size, 80000))
print(Checker.radiationCheck() * 100)