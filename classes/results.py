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
                    T[i].append(np.average([float(num) for num in f.readline().rstrip().split()]))
        
        ignoreIndex = byPass(t, ignoreVal)
        labels = ["Зенит", "Надир", "Перпенд 1", "Перпенд 2", "Теневая", "Солнечная"]
        for i in range(len(areas)):
            T[i] = np.array(T[i])
            result = np.average(T[i][ignoreIndex:]**4)
            print(labels[i])
            print(f"avg:{result**0.25:5.2f}, min: {np.min(T[i][ignoreIndex:]):5.2f}, max: {np.max(T[i][ignoreIndex:]):5.2f}")
        
        return result