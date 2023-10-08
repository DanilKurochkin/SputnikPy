import matplotlib.pyplot as plt

# file.write(str(self.conditions.external[1].radiation(self.boxes[5], 0)) + ' ' + str(ht*j) + '\n')
# file.write(str(j) + '\n')
# for box in self.boxes:
#     file.write(str(box.number))
#     for T in box.T:
#         file.write(" {0:16.4f}".format(float(T)))
#     file.write('\n')

t = []
angle = []
T1 = []
T2 = []
T3 = []
T4 = []
T5 = []
T6 = []

areas = [6,6,3,3,2,2]
area = sum(areas)
with open('output.txt', 'r') as f:
    while True:
        string = f.readline()
        if string == '':
            break
        string = string.rstrip().split()
        t.append(float(string[2]))
        angle.append(float(string[3]))
        T1.append([float(num) for num in f.readline().rstrip().split()])
        T2.append([float(num) for num in f.readline().rstrip().split()])
        T3.append([float(num) for num in f.readline().rstrip().split()])
        T4.append([float(num) for num in f.readline().rstrip().split()])
        T5.append([float(num) for num in f.readline().rstrip().split()])
        T6.append([float(num) for num in f.readline().rstrip().split()])

def extractcol(l):
    n_l = []
    for items in l:
        n_l.append(items[0])
    return n_l

T1 = extractcol(T1)
T2 = extractcol(T2)
T3 = extractcol(T3)
T4 = extractcol(T4)
T5 = extractcol(T5)
T6 = extractcol(T6)

plt.plot(t, T1, color = 'red', label = 'Зенит')
plt.plot(t, T2, color = 'blue', label = 'Надир')
plt.plot(t, T3, color = 'green', label = 'Перпендикулярная 1', )
plt.plot(t, T4, color = 'grey', label = 'Перпендикулярная 2')
plt.plot(t, T5, color = 'purple', label = 'Теневая')
plt.plot(t, T6, color = 'yellow', label = 'Солнечная')
plt.xlabel('Время [c]')
plt.ylabel('Температура [K]')

plt.legend()

sumT = 0
sumTrad = 0
for i in range(len(T1)//2, len(T1)):
    sumT += areas[0]*T1[i] + areas[1]*T2[i] + areas[2]*T3[i] + areas[3]*T4[i] + areas[4]*T5[i] + areas[5]*T6[i]
    sumTrad += areas[0]*T1[i]**4 + areas[1]*T2[i]**4 + areas[2]*T3[i]**4 + areas[3]*T4[i]**4 + areas[4]*T5[i]**4 + areas[5]*T6[i]**4
sumT /=(len(T1)//2)
sumT /= area
sumTrad /= (area*len(T1)//2)
sumTrad = sumTrad**0.25

print(str(sumT),str(sumTrad))

plt.show()