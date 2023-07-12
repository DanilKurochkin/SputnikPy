with open('outputheat.txt', 'r') as f:
    r = f.read().rsplit()

res = 0
ans = 0
for item in r:
    res += float(item)
    ans += abs(float(item))


print(res/ans*100, end = ' %')