import classes.workloads as wl
import classes.models as model
import matplotlib.pyplot as plt
import numpy as np

box1 = model.Box(0, 1, 1, 2, 2, 2, 2, 2)
table = wl.TableFunctionPeriodicalFlux([0, 0, 1000, 1000, 0, 0], [0, 99.9, 100, 200, 200.1, 300], 300)

dt = 10.0
listQ = []
listTime = []

for t in np.arange(0, 600, 10):
    listQ.append(table.heatFlux(box1, t))
    listTime.append(t)

plt.plot(listTime, listQ)
plt.show()