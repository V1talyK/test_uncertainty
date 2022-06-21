import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance

Lx, Ly = 1000, 1000
k, h, mu = 0.1*8.31, 10, 1
Pk = 10
x = np.linspace(0, 1000, 20)
y = np.linspace(0, 1000, 20)

grid = np.meshgrid(x, y)
W = 2*math.pi*k*h/mu

wxy = [[600,400],[400,600]]
dist = np.zeros([len(wxy),len(x)*len(y)])
dstw = np.zeros([len(wxy),len(wxy)])
tt = range(1,24)
P = np.zeros([len(tt),len(x)*len(y)])
qw = np.zeros((len(tt), len(wxy)))
pw = np.zeros([len(tt), len(wxy)])
ppl = np.zeros([len(tt), len(wxy)])

for t in tt:
    qw[t-1][0] = 1/(1 + math.exp(-t/16))*32
    qw[t-1][1] = 1/(1 + math.exp(-t/2))*16 - 1/(1 + math.exp(-t/8))*48

for i in range(len(wxy)):
    for j in range(len(wxy)):
        dstw[i][j] = distance.euclidean(wxy[i], wxy[j])
        if dstw[i][j]<0.05:
            dstw[i][j] = 0.05
    for k0, gxy0 in enumerate(zip(grid[0], grid[1])):
         for k, gxy in enumerate(zip(gxy0[0], gxy0[1])):
            dist[i][k0*20+k] = distance.euclidean(gxy, wxy[i])
            if dist[i][k0*len(x)+k]<0.05:
                dist[i][k0 * len(y) + k] = 0.05

for t in range(len(tt)):
    for i in range(len(x)):
        for j in range(len(y)):
            # print(i,j,dist[0][i*20+j],dist[1][i*20+j])
            P[t][i*20+j] = Pk + qw[t][0]*math.log(0.05/dist[0][i*20+j])/W+qw[t][1]*math.log(0.05/dist[1][i*20+j])/W

    for i in range(len(wxy)):
        pw[t][i] = Pk
        ppl[t][i] = np.mean(P[t][dist[i]<100])
        for j in range(len(wxy)):
            pw[t][i] = pw[t][i]+qw[t][j]*math.log(dstw[i][j]/1000)/W

print(qw)
# P[0] = grid[0].reshape((1,400))*grid[1].reshape((1,400))

fig, axs = plt.subplots(2,2)
axs[0,0].contourf(x, y, P[0].reshape(20,20))
axs[0,0].set_title('Карта давления')
axs[0,1].plot(tt, qw)
axs[0,1].set_title('Дебит')
axs[1,0].plot(tt, pw)
axs[1,0].set_title('Забойное')
axs[1,1].plot(tt, ppl)
axs[1,1].set_title('Пластовое')

plt.show()
