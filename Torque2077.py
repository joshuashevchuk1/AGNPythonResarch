import numpy as np
import pencil_old as pc
import math
import sys

def initData():
    # get only torque data and eccentricity data
    ts = pc.read_ts()
    t = ts.t / 2 * math.pi
    N = 850

    kernel = np.ones((N,)) / N
    time = np.convolve(t, kernel, mode='valid')
    torqint = np.convolve(ts.torqint_2, kernel, mode='valid')
    torqext = np.convolve(ts.torqext_2, kernel, mode='valid')
    torqtotal = torqext[:] + torqint[:]

    indexTimeCutOff = 0

    radius = ts.xq2
    Omega = ts.yq2
    LinearVelocity = ts.vxq2
    AngularVelocity = ts.vyq2

    v2 = LinearVelocity ** 2 + AngularVelocity ** 2
    semi_major = 1. / (2 / radius - v2)
    DArclength = radius ** 2 * (AngularVelocity / radius)
    ep1 = (DArclength ** 2) / semi_major
    eccentricity = (1 - ep1) ** 0.5
    ecc_int = par.eccentricity

    ecc = eccentricity

    for i in range(len(ecc)):
        if ecc_int != 0:
            if ecc[i] <= 0.01:
                indexTimeCutOff = i
                break

    timeCutOff = time[indexTimeCutOff] / (np.pi * 2)

    vars_dict = {'ts':ts,'time':time,'torqint':torqint,'torqext':torqext,'torqtotal':torqtotal,'timeCutOff':timeCutOff}
    return vars_dict

def plotData():
    data_frame = initData()
    time = data_frame[n]['time']

    # init for cut_off
    cut_off = data_frame[n]['indexTimeCutOff']
    time = time[:cut_off]
    torqint = data_frame[n]['torqint']
    torqext = data_frame[n]['torqext']
    torqtotal = data_frame[n]['torqtotal']
    #

    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    fig.subplots_adjust(bottom=0.07, top=0.95)
    ax1.plot(time, torqint[:len(time)], '--', label='Inner')
    ax1.plot(time, torqext[:len(time)], '--', label='Outer')
    ax1.plot(time, torqtotal[:len(time)], label='Total')
    ax1.set_aspect('auto')
    plt.xlim([0, time.max()])
    plt.tight_layout()
    plt.grid(True)

    # =========================
    # legend handles
    # =========================

    # use this information when adding important parameters to the run

    eccentricity = data_frame[n]['eccentricity']
    sound_speed = data_frame[n]['cs']
    initial_pressure = data_frame[n]['initial_pressure']
    Sigma = data_frame[n]['Sigma']
    rsmooth = data_frame[n]['rsmooth']
    gamma = data_frame[n]['gamma']
    alpha = data_frame[n]['alpha']
    beta = data_frame[n]['beta']

    DirMass = data_frame[n]['par1']
    DirEcc = (round(eccentricity[0], 1))

    DirMass_patches = mpatches.Patch(
        color='white', label='q :' + str(DirMass))
    DirEcc_patches = mpatches.Patch(
        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
    Dirsound_speed_patches = mpatches.Patch(
        color='white', label='sound speed :' + str(sound_speed))
    Dirinitial_pressure_patches = mpatches.Patch(
        color='white', label='inital pressure :' + str(initial_pressure))
    DirSigma_patches = mpatches.Patch(
        color='white', label=r'$\Sigma$ :' + str(Sigma))
    Dir_rsmooth_patches = mpatches.Patch(
        color='white', label='potential smoothing :' + str(rsmooth))
    Dir_gamma_patches = mpatches.Patch(
        color='white', label=r'$\gamma$ :' + str(gamma))
    Dir_alpha_patches = mpatches.Patch(
        color='white', label=r'$\alpha$ :' + str(alpha))
    Dir_beta_patches = mpatches.Patch(
        color='white', label=r'$\beta$ :' + str(beta))
    Dir_cutoff_patches = mpatches.Patch(
        color='white', label=r'$t_{c}$ :' + str(time[cut_off - 1]))

    plt.legend(handles=[DirMass_patches,
                        DirEcc_patches,
                        Dirsound_speed_patches,
                        Dirinitial_pressure_patches,
                        DirSigma_patches,
                        Dir_rsmooth_patches,
                        Dir_gamma_patches,
                        Dir_alpha_patches,
                        Dir_beta_patches,
                        Dir_cutoff_patches], loc=2)

    # =========================

    plt.title(r'Torques')
    plt.xlabel(r'$t/T_0$')
    plt.ylabel(r'$\Gamma$')
    plt.xlim([0, time[cut_off - 1]])
    plt.subplots_adjust(bottom=0.05, top=0.95)
    plt.savefig('Standard_Torque_' + '2077' + '.png')
    plt.close(fig)

plotData()