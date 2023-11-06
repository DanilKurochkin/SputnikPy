from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
mat = Material(8000, 2700, 100)
coat = Coating(0.9, 0.5)
orbit = SunLookingOrbit(500)

sp = Sputnik(Lx=1,Ly=2, Lz=3, width=0.001, material=mat, coat=coat, orbit=orbit)  
sp.createVolumes(10)

cond = Conditions()
cond.addEx(wl.Sun(1500), wl.Radiaton())
cond2 = Conditions()
cond2.addEt(wl.TableFunctionPeriodicalHeat([0, 0, 1000, 1000, 0, 0],
                                           [0, orbit.period/4, orbit.period/4+1,
                                            3*orbit.period/4, 3*orbit.period/4+1, orbit.period],
                                           orbit.period))

sp.addCondition(cond)
sp.addConditionByNum(cond2, 0)
wl.Connect.neighbours(sp, 0.1)

start = time.time()
sp.solve(60, 350, 1, 300, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.averageTCheck(sp.size, 20000))
print(Checker.radiationCheck() * 100)