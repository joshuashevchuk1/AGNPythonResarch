import pencil_old as pc
import numpy as np
import matplotlib.pyplot as plt

q = 0
ecc_int = 0
ts = None
t = None
name = 'test'


def initPars():
    global q
    global ecc_int
    global ts
    global t
    print('entering getCutOff')
    ts = pc.read_ts()
    t = ts.t

    radius = ts.xq2
    LinearVelocity = ts.vxq2
    AngularVelocity = ts.vyq2

    par = pc.read_param()

    if (par.iprimary == 1):
        q = par.pmass[1]
    else:
        q = par.pmass[0]

    v2 = LinearVelocity ** 2 + AngularVelocity ** 2
    semi_major = 1. / (2 / radius - v2)
    DArclength = radius ** 2 * (AngularVelocity / radius)
    ep1 = (DArclength ** 2) / semi_major
    eccentricity = (1 - ep1) ** 0.5
    ecc_int = par.eccentricity
    ecc = eccentricity

    indexTimeCutOff = 0

    for i in range(len(ecc)):
        if ecc_int != 0:
            if ecc[i] <= 0.01:
                indexTimeCutOff = i
                break

    timeCutOff = t[indexTimeCutOff] / (np.pi * 2)
    print('timeCutOff is ', np.round(timeCutOff))
    print('leaving get cutoff')


def plotTemp():
    global q
    global ecc_int
    global ts
    global t
    global name
    temp = ts.TTm
    temp_max = ts.TTmax
    temp_min = ts.TTmin

    t = t / (2 * np.pi)

    A = 1

    ticks = np.arange(0, 200, 1)
    line = A * ticks

    plt.plot(t, temp)
    plt.xlabel(r'$t/T_0$')
    plt.ylabel('Temperature')
    plt.tight_layout()
    plt.grid(True)
    plt.xlim([0, t.max()])
    plt.title("T(q=" + str(q) +" , " +r'$\varepsilon$' + "="+ str(ecc_int) + ")")
    # save figure
    plt.savefig('PAPER_TTm' + str(name) + '.png',bbox_inches="tight")
    plt.close()

initPars()
plotTemp()
