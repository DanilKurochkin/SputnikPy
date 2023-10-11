import numpy as np

class Material(): #материал пластины

    def __init__(self, c: np.float64, p: np.float64, k: np.float64):
        self.c = c #теплоёмкость
        self.p = p #плотность
        self.k = k #температуропроводность
        
class Coating(): #покрытие спутника
    def __init__(self, As, epsilon):
        self.As = As #коэфициент поглащения солнечного излучения, то есть высокочастотного
        self.epsilon = epsilon #степень черноты, а также коэфициент поглащения в низкочастотном спектре