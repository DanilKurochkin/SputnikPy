import numpy as np
import scipy as sp
import classes.elemath as SMath

l = 0.0012
size = 7
p = 1000.0
cp = 1000.0
k = 202.0

dh = l / (size - 1)
dt = 56.68153

a = [0] * size
b = [0] * size
c = [0] * size
d = [0] * size

b[0] = k/dh
a[0] = b[0]
d[0] = 0

for i in range(1, size -1):
    b[i] = k/dh
    c[i] = k/dh
    a0 = (p*cp*dh)/dt
    d[i] = a0 * 300.0
    a[i] = a0 + b[i] + c[i]

c[-1] = k/dh
d[-1] = 1000
a[-1] = c[-1]
a = np.array(a)
b = np.array(b)
c = np.array(c)
d = np.array(d)

print(SMath.TDMA(a, b, c, d, np.empty(size), np.empty(size)))