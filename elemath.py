import numpy as np
import numpy.typing as npt

def TDMA(a : npt.NDArray[np.float64], b : npt.NDArray[np.float64], c : npt.NDArray[np.float64], d : npt.NDArray[np.float64], P  : npt.NDArray[np.float64] , Q  : npt.NDArray[np.float64]):
    n = a.size
    x = np.empty(n, dtype=np.float64)
    
    P[0] = b[0]/a[0]
    Q[0] = d[0]/a[0]
    for i in np.arange(1, n):
        P[i] = b[i]/(a[i]-c[i]*P[i-1])
        Q[i] = (d[i]+c[i]*Q[i-1])/(a[i]-c[i]*P[i-1])

    x[n-1] = Q[n - 1]
    for i in np.arange(n-1, -1, -1):
        x[i - 1] = P[i - 1]*x[i] + Q[i - 1]
    return x

def discrepancy(array1 : npt.NDArray[np.float64], array2 : npt.NDArray[np.float64]):
    n = array1.size
    
    epsilon = 200
    
    for i in np.arange(n):
        new_epsilon = abs(array1[i] - array2[i])
        if new_epsilon < epsilon:
            epsilon = new_epsilon
    
    return epsilon