from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
A6061 = Material(1000, 1000, 202)
coat = Coating(0.3, 0.3)
orbit = SunLookingOrbit(500)
cond = Conditions()
sp = Sputnik(1,2,3, 0.001, A6061, coat, orbit)

sp.createVolumes(3)
sp.knitPlates()

#cond.addEx(wl.Sun(1500), wl.Radiaton())
cond.addEt(wl.TableFunctionPeriodicalHeat([1000, 1000], [0, 1000], orbit.period))
sp.addCondition(cond)
#wl.Connect.neighbours(sp, 0.1)

start = time.time()
sp.solve(1, 100, 1, 300, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.averageTCheck(sp.size, 60000))
print(Checker.radiationCheck() * 100)