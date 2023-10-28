from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
A6061 = Material(8000, 2700, 100)
coat = Coating(0.9, 0.5)
orbit = SunLookingOrbit(500)

sp = Sputnik(1,2,3, 0.001, A6061, coat, orbit)  
sp.createVolumes(7)

cond = Conditions()
cond.addEx(wl.Sun(1500), wl.Radiaton())
cond2 = Conditions()
cond2.addEt(wl.TableFunctionPeriodicalFlux([0,0,1000,1000,0,0],
                                     [0,orbit.period/4, orbit.period/4+1,
                                      3*orbit.period/4,3*orbit.period/4+1, orbit.period],
                                     orbit.period))

sp.addCondition(cond)
sp.addCondition(cond2)
wl.Connect.neighbours(sp, 0.01)

start = time.time()
sp.solve(50, 333, 1, 300, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.averageTCheck(sp.size, 20000))
print(Checker.radiationCheck() * 100)