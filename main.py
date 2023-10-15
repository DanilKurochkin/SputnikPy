from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
A6061 = Material(800, 2700, 202) #алюминий 6061
coat = Coating(0.3, 0.3) #ЭКОМ-1П a05 e04
orbit = SunLookingOrbit(500)
cond = Conditions()
sp = Sputnik(1,2,3, 0.001, A6061, coat, orbit)

sp.createVolumes(12)
sp.knitPlates()

cond.addEx(wl.Radiaton(), wl.Sun(1500))
cond.addEt(wl.TableFunctionPeriodicalHeat([0, 1000, 1000, 0],[0,orbit.period/4,3*orbit.period/4, orbit.period],orbit.period))
sp.addCondition(cond)
wl.Connect.neighbours(sp, 0.1)

start = time.time()
sp.solve(40, 800, 1, 305, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.averageTCheck(sp.size, 300000))
print(Checker.radiationCheck())