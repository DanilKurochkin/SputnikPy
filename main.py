from classes.orbits import SunLookingOrbit, earth_radius
from classes.models import *
import classes.workloads as wl
from classes.material import Material, Coating
from classes.results import Checker
import time

points_in_orbit = 56
orbit_sim_amount = 40
    
mat = Material(c=800, p=2700, k=100)
coat = Coating(As=0.6, epsilon=0.4)
orbit = SunLookingOrbit(radiusAboveEarth=500)

sp = Sputnik(Lx=1,Ly=1, Lz=1, width=0.005, material=mat, coat=coat, orbit=orbit)  
sp.createVolumes(10)

cond = Conditions()
cond.addEx(wl.Sun(1367), wl.Radiaton(), wl.EarthAlbedo(0.3, 1367), wl.EarthRadiation(240, earth_radius))
cond.addEt(wl.TableFunctionLoadHeat([50, 50], [0, 10]))

sp.addCondition(cond)
wl.Connect.neighbours(sp, 2)

print(orbit.period)
start = time.time()
sp.solve(orbit_sim_amount, points_in_orbit, 1, 273, radiation_check=True)
end = time.time()
print(end-start)

for box in sp.boxes:
    print(np.sum(box.integr)/(orbit_sim_amount*points_in_orbit))

print(Checker.averageTCheck(sp.size, orbit.period*20))
print(Checker.radiationCheck() * 100)