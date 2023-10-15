import classes.workloads as wl
import classes.models as model
import matplotlib.pyplot as plt
import numpy as np

box1 = model.Box(0, 1, 4, 2, 2, 2, 2, 2)
box2 = model.Box(0,1,2, 3,3,3,3,3)
table = wl.TableFunctionPeriodicalHeat([0,1000, 1000, 0], [0, 100, 200, 300], 300)

t = np.linspace(0, 300*300, 300*300)

q1 = table.heatFlux(box1, t)
q2 = table.heatFlux(box2, t)
plt.plot(t, q1)
plt.plot(t, q2)
plt.show()