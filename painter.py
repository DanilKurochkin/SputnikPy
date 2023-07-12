min = 0
max = 100

def floatToRGB(num):
    delta = max - min
    num -= min
    
    if num <= delta * 0.25:
        num /= (delta * 0.25)
        return (0, num, 1)
    
    elif num <= delta * 0.5:
        num -= delta * 0.25
        num /= (delta * 0.25)
        return (0, 1, 1-num)
    
    elif num <= delta * 0.75:
        num -= (delta * 0.5)
        num /= (delta * 0.25)
        return (num, 1, 0)
    
    else:
        num -= (delta * 0.75)
        num /= (delta * 0.25)
        return (1, 1-num, 0)
