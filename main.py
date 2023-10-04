import numpy as np
import classes.workloads as wl
import time
from classes.material import Coating, Material
from classes.models import *
from classes.orbits import SunLookingOrbit, earth_radius

def main():
    A6061 = Material(800, 2700, 202) #алюминий 6061
    coat = Coating(0.5, 0.4)
    orbit = SunLookingOrbit(500)
    cond = Conditions()

    sp = Sputnik(1, 2, 1, 0.002, A6061, coat, orbit)
    sp.createVolumes(0.0002)
    sp.knitPlates()
    
    cond.addEx(wl.Radiaton(), wl.EarthRadiation(240, earth_radius), wl.Sun(1500), wl.EarthAlbedo(0.3, 1500))
    cond.addEt(wl.TableFunctionLoadPeriodical(np.array([0, 0, 1000, 1000], dtype=np.float64),
                                            np.array([0, orbit.period/2, orbit.period/2+1, orbit.period], dtype=np.float64),
                                            orbit.period))
    sp.addCondition(cond)
    wl.Connect.neighbours(sp, 0.5)

    start = time.time()
    sp.solve(4, 40, 1, 290, radiation_check=True)
    end = time.time()
    
    print(end-start)

main()