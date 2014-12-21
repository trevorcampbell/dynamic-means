import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm


f = open('tree.log', 'r')
dlines = f.readlines()
f.close()
nToRead = int(dlines[0])
arr = np.zeros((nToRead, 2))
pa = np.zeros(nToRead)
for i in range(nToRead):
    rw = map(float, filter(None, dlines[1+i].split(' ')))
    pa[i] = int(rw[0])
    arr[i, :] = rw[1:]


fig = plt.figure()
for i in range(1, pa.shape[0]):
    plt.plot( [arr[i,0], arr[pa[i], 0]], [arr[i,1], arr[pa[i], 1]], c='b')

plt.show()


