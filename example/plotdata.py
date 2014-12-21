import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm

colors = [
    (1., 0., 0., 1.),
    (0., 1., 0., 1.),
    (0., 0., 1., 1.),
    (.5, .5, 0., 1.),
    (.5, 0., .5, 1.),
    (0., .5, .5, 1.),
    (0., 1., .5, 1.),
    (0., .5, 1., 1.),
    (1., 0., .5, 1.),
    (.5, 0., 1., 1.),
    (1., .5, 0., 1.),
    (.5, 1., 0., 1.)]


f = open('data.log', 'r')
dlines = f.readlines()
f.close()
f = open('lbls.log', 'r')
llines = f.readlines()
f.close()
arrs = []
lbls = []
idx = 0
while(idx < len(dlines)):
    nToRead = int(dlines[idx])
    idx += 1
    arr = np.zeros((nToRead, 2))
    lbl = np.zeros(nToRead)
    for i in range(nToRead):
        arr[i, :] = map(float, filter(None, dlines[idx+i].split(' ')))
        lbl[i] = float(llines[idx+i])
    idx += nToRead
    arrs.append(arr)
    lbls.append(lbl)

#fig = plt.figure()
#ax = plt.axes(xlim=(0,1), ylim=(0,1))
#scat = plt.scatter([], [], 'ko')
#
#def animate(i, fig, scat):
#    scat.set_offsets( tuple(map(list, arrs[i])) )
#    return scat
#
#anim = animation.FuncAnimation(fig, animate, fargs=(fig, scat), frames=len(arrs), interval=200, blit=True)
#plt.show()

#plot using lines
#fig = plt.figure()
#ax = plt.axes(xlim=(0,1), ylim=(0,1))
#line, =ax.plot([], [], 'ko')
#
#def init():
#    line.set_data([], [])
#    return line,
#
#def animate(i):
#    line.set_data(arrs[i][:, 0], arrs[i][:, 1])
#    return line,
#
#anim = animation.FuncAnimation(fig, animate, frames=len(arrs), init_func=init, interval=200, blit=True)
#plt.show()


#plot using lines
fig = plt.figure()
ax = plt.axes(xlim=(0,1), ylim=(0,1))
lns = []
for i in range(arrs[0].shape[0]):
    line, =ax.plot([], [], 'ko')
    lns.append(line)

def init():
    for line in lns:
        line.set_data([], [])
    return lns

def animate(i):
    idx = 0
    for line in lns:
        line.set_data(arrs[i][idx, 0], arrs[i][idx, 1])
        line.set_color(colors[int(lbls[i][idx]) % len(colors)])
        idx += 1
    return lns

anim = animation.FuncAnimation(fig, animate, frames=len(arrs), init_func=init, interval=200, blit=True)
plt.show()
