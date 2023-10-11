import numpy as np
from classes.orbits import SunLookingOrbit
from classes.models import Sputnik
import classes.workloads as wl
from classes.materials import Material, Coating
import time
        
def main():
    
    A6061 = Material(8000, 2700, 202) #алюминий 6061
    coat = Coating(0.4, 0.2) #ЭКОМ-1П a05 e04
    orbit = SunLookingOrbit(500)
    cond = wl.Conditions()
    
    sp = Sputnik(4, 0.002, A6061, coat, orbit)
    sp.createVolumes(0.0002)
    sp.knitPlates()
    
    cond.addEx(wl.Radiaton(), wl.Sun(1500))
    sp.addCondition(cond)
    wl.Connect.neighbours(sp, 0.01)

    start = time.time()
    sp.solve(60, 1000, 1, 240, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()