import numpy as np
import math
import matplotlib.pyplot as plt
from libs import calc_sim, calc_dist, make_fact

Lx, Ly = 1000, 1000
k, h, mu = 0.1*8.31, 10, 1
Pk = 10
x = np.linspace(0, Lx, 20)
y = np.linspace(0, Ly, 20)

grid = np.meshgrid(x, y)
W = 2*math.pi*k*h/mu

wxy = [[700, 300], [300, 700]]

tt = range(1, 24)
P = np.zeros([len(tt),len(x)*len(y)])
qw = np.zeros((len(tt), len(wxy)))
pw = np.zeros([len(tt), len(wxy)])
ppl = np.zeros([len(tt), len(wxy)])
qw0 = [1.1,1]
for t in tt:
    qw[t-1][0] = (1/(1 + math.exp(-t/16))*32)*qw0[1]
    qw[t-1][1] = (1/(1 + math.exp(-t/2))*16 - 1/(1 + math.exp(-t/8))*48)*qw0[1]

dist, dstw = calc_dist(wxy, grid, x, y)
for t in range(len(tt)):
    P[t], pw[t], ppl[t] = calc_sim(Pk, qw[t], dist, dstw, W, x, y)

e1r = 0.05
pplf = make_fact(tt, Pk, qw, dist, dstw, W, x, y, e1r)

ttf = [[],[]]
pplfp = [[],[]]
JMAPE = 0;
for v in pplf:
    ttf[int(v[0])].append(v[1])
    pplfp[int(v[0])].append(v[2])
    JMAPE = JMAPE + abs(v[2]-ppl[v[1]-1][int(v[0])])/ppl[v[1]-1][int(v[0])]

print(JMAPE/16)

fig, axs = plt.subplots(2,2)
axs[0,0].contourf(x, y, P[0].reshape(20,20))
axs[0,0].set_title('Карта давления')
axs[0,1].plot(tt, qw)
axs[0,1].set_title('Дебит')
axs[1,0].plot(tt, pw)
axs[1,0].set_title('Забойное')
axs[1,1].plot(tt, ppl)
axs[1,1].set_title('Пластовое')
axs[1,1].plot(ttf[0], pplfp[0], 'o')
axs[1,1].plot(ttf[1], pplfp[1], 'o')


plt.show()
