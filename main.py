import numpy as np
import classes.workloads as wl
import time
from classes.material import Coating, Material
from classes.models import *
from classes.orbits import SunLookingOrbit

def main():
    A6061 = Material(9000, 2700, 209) #алюминий 6061
    coat = Coating(0.5, 0.6)
    orbit = SunLookingOrbit(500)
    cond = Conditions()

    sp = Sputnik(1, 2, 3, 0.001, A6061, coat, orbit)
    sp.createVolumes(20)
    sp.knitPlates()
    
    cond.addEx(wl.Sun(1500), wl.Radiaton())
    sp.addCondition(cond)
    wl.Connect.neighbours(sp, 0.1)
    start = time.time()
    sp.solve(50, 1000, 1, 250, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()