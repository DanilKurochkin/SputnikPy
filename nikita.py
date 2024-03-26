from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time
        
mat = Material(c=800, p=2700, k=100)
coat = Coating(As=0.92, epsilon=0.95)
orbit = SunLookingOrbit(radiusAboveEarth=766)

sp = Sputnik(Lx=1,Ly=1, Lz=1, width=0.0005, material=mat, coat=coat, orbit=orbit)  
sp.createVolumes(10)

cond = Conditions()
cond.addEx(wl.Sun(1367), wl.Radiaton())

sp.addCondition(cond)

print(orbit.period)
start = time.time()
sp.solve(10, 500, 1, 273, radiation_check=True)
end = time.time()
print(end-start)

print(Checker.averageTCheck(sp.size, 80000))
print(Checker.radiationCheck() * 100)