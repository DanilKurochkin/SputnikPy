from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Cheker
import time
        
A6061 = Material(10000, 3000, 202) #алюминий 6061
coat = Coating(0.3, 0.3) #ЭКОМ-1П a05 e04
orbit = SunLookingOrbit(500)
cond = Conditions()
sp = Sputnik(1,2,3, 0.001, A6061, coat, orbit)

sp.createVolumes(12)
sp.knitPlates()

cond.addEx(wl.Radiaton(), wl.Sun(1500))
cond.addEt(wl.TableFunctionLoadHeat([0, 1000, 0], [0,orbit.period/2,orbit.period], orbit.period))
sp.addCondition(cond)
wl.Connect.neighbours(sp, 0.1)

start = time.time()
sp.solve(300, 100, 1, 305, radiation_check=True)
end = time.time()
print(end-start)

print(Cheker.averageTCheck(sp.size, 40000))
print(Cheker.radiationCheck())