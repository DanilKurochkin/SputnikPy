import matplotlib.pyplot as plt
import numpy as np
t = []
angle = []
labels = ['Зенит', 'Надир', 'Перпендикулярная 1', 'Перпендикулярная 2', 'Теневая', 'Солнечная']
colors = ['red', 'blue', 'green', 'grey', 'purple', 'yellow']
T = [[],[],[],[],[],[]]

with open('output.txt', 'r') as f:
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
    return n_l

for i in range(len(labels)):
    plt.plot(t, extractcol(T[i]), color = colors[i], label = labels[i])

plt.xlabel('Время [c]')
plt.ylabel('Температура [K]')

plt.legend()

plt.show()