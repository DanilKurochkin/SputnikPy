import numpy as np
import classes.workloads as wl
import time
from classes.material import Coating, Material
from classes.models import *
from classes.orbits import SunLookingOrbit, earth_radius

def main():
    A6061 = Material(1000, 1000, 202) #алюминий 6061
    coat = Coating(1, 0.0)
    orbit = SunLookingOrbit(500)
    cond = Conditions()

    sp = Sputnik(1, 2, 3, 0.001, A6061, coat, orbit)
    sp.createVolumes(5)
    sp.knitPlates()
    
    cond.addEx(wl.Isolated())
    cond.addEt(wl.ConstantHeatFlux(1000))
    sp.addCondition(cond)

    start = time.time()
    sp.solve(4, 100, 1, 300, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()