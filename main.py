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
    
    cond.addEx()
    cond.addEt(wl.ConstantHeatFlux(1000))
    sp.addCondition(cond)

    start = time.time()
    sp.solve(4, 100, 1, 300, radiation_check=True)
    end = time.time()
    
    print(end-start)
    # orbit = SunLookingOrbit(500)
    # timeDots = np.array([0, orbit.period/2 - 0.5, orbit.period/2 + 0.5, orbit.period], dtype=np.float64)
    # heatDots = np.array([0, 0, 1000, 1000], dtype=np.float64)
    # tabl = wl.TableFunctionLoadPeriodical(heatDots,
    #                                         timeDots,
    #                                         orbit.period)
    # n = timeDots[(0<timeDots) & (timeDots<orbit.period)]
    # print(str(tabl.calculator(orbit.period/2-0.5, 0.5)))

main()