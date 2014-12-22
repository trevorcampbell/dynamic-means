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
maxToPlot = 0
while(idx < len(dlines)):
    nToRead = int(dlines[idx])
    idx += 1
    arr = np.zeros((nToRead, 2))
    lbl = np.zeros(nToRead)
    for i in range(nToRead):
        arr[i, :] = map(float, filter(None, dlines[idx+i].split(' ')))
        lbl[i] = float(llines[idx+i])
    idx += nToRead
    if (nToRead > maxToPlot):
        maxToPlot = nToRead
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





if (True):
    fig = plt.figure()
    for i in range(len(arrs)):
        ax = plt.axes(xlim=(0,1), ylim=(0,1))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        lns = []
        for j in range(maxToPlot):
            line, =ax.plot([], [], 'wo', lw=0, markeredgecolor='w', markeredgewidth=0, markersize=10)
            lns.append(line)
        for j in range(arrs[i].shape[0]):
            lns[j].set_data(arrs[i][j, 0], arrs[i][j, 1])
            lns[j].set_color(colors[int(lbls[i][j]) % len(colors)])
        for j in range(arrs[i].shape[0], len(lns)):
            lns[j].set_data([], [])
            lns[j].set_color( (1.0, 1.0, 1.0, 0.0) )
        plt.savefig('frames/fr-'+str(i)+'.pdf')
        plt.cla()
    plt.close()


#plot using lines
fig = plt.figure()
ax = plt.axes(xlim=(0,1), ylim=(0,1))
lns = []
for i in range(maxToPlot):
    line, =ax.plot([], [], 'wo', lw=0, markeredgecolor='w', markeredgewidth=0, markersize=10)
    lns.append(line)

def init():
    for line in lns:
        line.set_data([], [])
    return lns

def animate(i):
    for j in range(arrs[i].shape[0]):
        lns[j].set_data(arrs[i][j, 0], arrs[i][j, 1])
        lns[j].set_color(colors[int(lbls[i][j]) % len(colors)])
    for j in range(arrs[i].shape[0], len(lns)):
        lns[j].set_data([], [])
        lns[j].set_color( (1.0, 1.0, 1.0, 0.0) )
    return lns

anim = animation.FuncAnimation(fig, animate, frames=len(arrs), init_func=init, interval=200, blit=True)
plt.show()



