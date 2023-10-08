import numpy as np
import classes.workloads as wl
import time
from classes.material import Coating, Material
from classes.models import *
from classes.orbits import SunLookingOrbit, earth_radius

def main():
    A6061 = Material(15000, 2700, 209) #алюминий 6061
    coat = Coating(0.4, 0.2)
    orbit = SunLookingOrbit(500)
    cond = Conditions()

    sp = Sputnik(1, 1, 1, 0.001, A6061, coat, orbit)
    sp.createVolumes(20)
    sp.knitPlates()
    
    cond.addEx(wl.Sun(1500), wl.Radiaton())
    
    sp.addCondition(cond)
    wl.Connect.neighbours(sp, 0.01)
    start = time.time()
    sp.solve(40, 100, 1, 300, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()