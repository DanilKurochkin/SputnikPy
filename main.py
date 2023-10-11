import numpy as np
from classes.orbits import SunLookingOrbit
from classes.models import *
import classes.workloads as wl
from classes.materials import Material, Coating
import time
        
def main():
    A6061 = Material(1000, 3000, 202) #алюминий 6061
    coat = Coating(0.4, 0.2) #ЭКОМ-1П a05 e04
    orbit = SunLookingOrbit(500)
    cond = wl.Conditions()
    
    sp = Sputnik(1,2,3, 0.005, A6061, coat, orbit)
    sp.createVolumes(7)
    sp.knitPlates()
    
    cond.addEx(wl.Radiaton(), wl.Sun(1500))
    cond.addEt(wl.ConstantHeatFlux(500))
    sp.addCondition(cond)
    wl.Connect.neighbours(sp, 0.1)
    start = time.time()
    sp.solve(60, 1000, 1, 450, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()