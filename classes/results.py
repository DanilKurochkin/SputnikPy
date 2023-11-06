import numpy as np

class Checker():
    def radiationCheck():
        with open('outputheat.txt', 'r') as f:
            r = f.read().rsplit()

        res = 0
        ans = 0
        for item in r:
            res += float(item)
            ans += abs(float(item))

        result = res/ans
        
        return result
    
    def averageTCheck(areas, ignoreVal = 0):
        def byPass(list, val):
            n = len(list)
            for i in range(n):
                if list[i] > val:
                    return i
            
            return n
        
        t = []
        T = [[],[],[],[],[],[]]

        with open('output.txt', 'r') as f:
            f.readline()
            while True:
                string = f.readline()
                if string == '':
                    break
                string = string.rstrip().split()
                t.append(float(string[2]))
                for i in range(6):
                    T[i].append([float(num) for num in f.readline().rstrip().split()])
        
        ignoreIndex = byPass(t, ignoreVal)
        result = 0
        for i in range(len(areas)):
            result += np.average(T[i][ignoreIndex:]) * areas[i]
        result = result/np.sum(areas)
        
        return result