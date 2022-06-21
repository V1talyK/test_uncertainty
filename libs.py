import math
import numpy as np
from scipy.spatial import distance
import random

def calc_sim(Pk,qw,dist, dstw, W, x, y):
    P = np.zeros(len(dist[0]))
    pw = np.zeros(len(qw))
    ppl = np.zeros(len(qw))

    for i in range(len(x)):
        for j in range(len(y)):
            P[i * 20 + j] = Pk
            for k in range(len(qw)):
                P[i*20+j] = P[i * 20 + j] + qw[k]*math.log(0.05/dist[k][i*20+j])/W

    for i in range(len(qw)):
        pw[i] = Pk
        ppl[i] = np.mean(P[dist[i]<50])
        for j in range(len(qw)):
            pw[i] = pw[i]+qw[j]*math.log(dstw[i][j]/1000)/W

    return P, pw, ppl

def calc_dist(wxy, grid, x, y):
    dist = np.zeros([len(wxy), len(x) * len(y)])
    dstw = np.zeros([len(wxy), len(wxy)])
    for i in range(len(wxy)):
        for j in range(len(wxy)):
            dstw[i][j] = distance.euclidean(wxy[i], wxy[j])
            if dstw[i][j] < 0.05:
                dstw[i][j] = 0.05
        for k0, gxy0 in enumerate(zip(grid[0], grid[1])):
            for k, gxy in enumerate(zip(gxy0[0], gxy0[1])):
                dist[i][k0 * 20 + k] = distance.euclidean(gxy, wxy[i])
                if dist[i][k0 * len(x) + k] < 0.05:
                    dist[i][k0 * len(y) + k] = 0.05
    return dist, dstw

def make_fact(tt, Pk, qw, dist, dstw, W, x, y, e1r):
    pw = np.zeros([len(tt), len(qw[0])])
    ppl = np.zeros([len(tt), len(qw[0])])

    for t in range(len(tt)):
        P, pw[t], ppl[t] = calc_sim(Pk, qw[t], dist, dstw, W, x, y)

    pplf = []
    for i in range(8):
        for j in range(len(qw[0])):
            ia = random.randint(0, len(qw)-1)
            shft = (1+e1r*random.randint(-5, 5)/5)
            pplf.append([j, ia+1, ppl[ia][j]*shft])
    return pplf