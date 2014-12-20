import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

f = open('data.log', 'r')
lines = f.readlines()
f.close()
arrs = []
idx = 0
while(idx < len(lines)):
    nToRead = int(lines[idx])
    idx += 1
    arr = np.zeros((nToRead, 2))
    for i in range(nToRead):
        arr[i, :] = map(float, filter(None, lines[idx+i].split(' ')))
    idx += nToRead
    arrs.append(arr)

fig = plt.figure()
ax = plt.axes(xlim=(0,1), ylim=(0,1))
line, = ax.plot([], [], 'ko')

def init():
    line.set_data([],[])
    return line,

def animate(i):
    line.set_data(arrs[i][:,0], arrs[i][:,1])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=100, interval=200, blit=True)
plt.show()
