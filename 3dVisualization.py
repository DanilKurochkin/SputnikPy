import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import painter as pter

cube_size = 3
colors = [[],
          [],
          [],
          [],
          [],
          []]
orbit_radius = 13
earth_radius = 8
rotate = True
cube_size /= 2

t = []
angle = []

with open('output.txt', 'r') as f:
    while True:
        string = f.readline()
        if string == '':
            break
        string = string.rstrip().split()
        t.append(float(string[2]))
        angle.append(float(string[3]))
        for T in colors:
            T.append([float(num) for num in f.readline().rstrip().split()])

def extractcol(l):
    n_l = []
    for items in l:
        n_l.append(items[0])
    return n_l

for i in np.arange(len(colors)):
    colors[i] = extractcol(colors[i])

pter.min = np.min(colors) - 1
pter.max = np.max(colors) + 1

def init():
    Z1 = np.linspace(-cube_size, cube_size, 2) #зенит
    Y1 = np.linspace(-cube_size, cube_size, 2)
    Z1, Y1 = np.meshgrid(Z1, Y1)
    X1 = Z1*0 + Y1*0 + cube_size

    Z2 = np.linspace(-cube_size, cube_size, 2) #надир
    Y2 = np.linspace(-cube_size, cube_size, 2)
    Z2, Y2 = np.meshgrid(Z2, Y2)
    X2 = Z2*0 + Y2*0 - cube_size

    X3 = np.linspace(-cube_size, cube_size, 2) #солнечная
    Y3 = np.linspace(-cube_size, cube_size, 2)
    X3, Y3 = np.meshgrid(X3, Y3)
    Z3 = X3*0 + Y3*0 - cube_size

    X4 = np.linspace(-cube_size, cube_size, 2) #теневая
    Y4 = np.linspace(-cube_size, cube_size, 2)
    X4, Y4 = np.meshgrid(X4, Y4)
    Z4 = X4*0 + Y4*0 + cube_size

    X5 = np.linspace(-cube_size, cube_size, 2) #перпендикулярые
    Z5 = np.linspace(-cube_size, cube_size, 2)
    X5, Z5 = np.meshgrid(X5, Z5)
    Y5 = X5*0 + Z5*0 - cube_size

    X6 = np.linspace(-cube_size, cube_size, 2)
    Z6 = np.linspace(-cube_size, cube_size, 2)
    X6, Z6 = np.meshgrid(X6, Z6)
    Y6 = X6*0 + Z6*0 + cube_size
    
    X = [X1, X2, X5, X6, X4, X3]
    Y = [Y1, Y2, Y5, Y6, Y4, Y3]
    Z = [Z1, Z2, Z5, Z6, Z4, Z3]
    
    return X, Y, Z

def rotateOX(x, y, z, X, Y, Z, angle):
    matrix = np.array([[1, 0, 0],
                  [0, np.cos(angle), -np.sin(angle)],
                  [0, np.sin(angle), np.cos(angle)]])
    
    
    for i in range(len(X)):
        for j in range(len(X[i])):
            vec = [X[i,j], Y[i,j], Z[i,j]]
            vec = matrix.dot(vec)
            
            x[i,j] = vec[0]
            y[i,j] = vec[1]
            z[i,j] = vec[2]

def rotateOY(x, y, z, X, Y, Z, angle):
    matrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                  [0, 1, 0],
                  [-np.sin(angle), 0, np.cos(angle)]])
    
    
    for i in range(len(X)):
        for j in range(len(X[i])):
            vec = [X[i,j], Y[i,j], Z[i,j]]
            vec = matrix.dot(vec)
            
            x[i,j] = vec[0]
            y[i,j] = vec[1]
            z[i,j] = vec[2]
            
def rotateOZ(x, y, z, X, Y, Z, angle):
    matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                  [np.sin(angle), np.cos(angle), 0],
                  [0, 0, 1]])
    
    
    for i in range(len(X)):
        for j in range(len(X[i])):
            vec = [X[i,j], Y[i,j], Z[i,j]]
            vec = matrix.dot(vec)
            
            x[i,j] = vec[0]
            y[i,j] = vec[1]
            z[i,j] = vec[2]

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

X, Y, Z = init()

x = np.empty(X[0].shape, dtype=np.float32)
y = np.empty(X[0].shape, dtype=np.float32)
z = np.empty(X[0].shape, dtype=np.float32)

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x1 = earth_radius * np.cos(u)*np.sin(v)
y1 = earth_radius * np.sin(u)*np.sin(v)
z1 = earth_radius * np.cos(v)

firstPhase = True
secondPhase = False

def animate(i):
    ax.cla()
    global x, y, z
    global firstPhase, secondPhase
    if rotate:
        if angle[i]%(2*np.pi) > np.pi and firstPhase:
            for j in range(len(X)):
                rotateOX(X[j], Y[j], Z[j], X[j], Y[j], Z[j], np.pi)
            firstPhase = False
            secondPhase = True
        elif angle[i]%(2*np.pi) > 0 and angle[i]%(2*np.pi) < np.pi and secondPhase:
            for j in range(len(X)):
                rotateOX(X[j], Y[j], Z[j], X[j], Y[j], Z[j], np.pi)
            firstPhase = True
            secondPhase = False
      
    for j in range(len(X)):
        rotateOY(x, y, z, X[j], Y[j], Z[j], -angle[i])
        x = x + orbit_radius*np.cos(angle[i])
        z = z + orbit_radius*np.sin(angle[i])
        new_color = pter.floatToRGB(colors[j][i])
        ax.plot_surface(x, y, z, color=new_color, shade=False)
    
    ax.quiver(20, 0, 0, -5, 0, 0, color = '#000')
    ax.quiver(20, -5, 0, -5, 0, 0, color = '#000')
    ax.quiver(20, 5, 0, -5, 0, 0, color = '#000')

    ax.plot_surface(x1, y1, z1, color = '#ff00f7', shade=False)
    ax.axes.set_xlim3d(left=-20, right=20) 
    ax.axes.set_ylim3d(bottom=-20, top=20) 
    ax.axes.set_zlim3d(bottom=-20, top=20)
    ax.set_xlabel('Ось X')
    ax.set_ylabel('Ось Y')
    ax.set_zlabel('Ось Z')
    return 0,

ani = animation.FuncAnimation(
    fig, animate, interval=1, frames = len(angle))

plt.show()