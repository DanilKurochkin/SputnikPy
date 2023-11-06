import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import numpy as np
t = []
angle = []
labels = ['Зенит', 'Надир', 'Перпендикулярная 1', 'Перпендикулярная 2', 'Теневая', 'Солнечная']
colors = ['red', 'blue', 'green', 'grey', 'purple', 'yellow']
T = [[],[],[],[],[],[]]

with open('output.txt', 'r') as f:
    string = f.readline()
    period = int(string.rstrip().rsplit()[0])
    while True:
        string = f.readline()
        if string == '':
            break
        string = string.rstrip().split()
        t.append(float(string[2]))
        angle.append(float(string[3]))
        for i in range(6):
            T[i].append([float(num) for num in f.readline().rstrip().split()])

def extractcol(l):
    n_l = []
    for items in l:
        n_l.append(np.average(items))
    return np.array(n_l)

cut = period
for i in range(len(T)):
    T[i] = extractcol(T[i])
    newT = np.empty(len(T[i]) - period)
    for j in range(len(T[i]) - period - 1):
        newT[j] = T[i][j+period+1] - T[i][j]
    T[i] = newT
    # T1 = T[i][cut:]
    # T2 = T[i][:-cut]
    # T[i] = T1 - T2

for i in range(len(labels)):
    plt.plot(t[cut:], T[i], color = colors[i], label = labels[i])

plt.xlabel('Время [c]')
plt.ylabel('Температура [K]')

plt.legend()

plt.show()

plt.show()