from itertools import islice
from random import random
from time import perf_counter
from scipy import signal
from scipy.signal import *
import numpy as np
import scipy
import scipy.integrate
from scipy.interpolate import interp1d
import triggerfishvcv
from plotting import *

"""
w in radian / sec
"""
def vanderPol_numpy(t,xt, mu, w):
    x = interp1d(t,xt, fill_value="extrapolate")
    def func(y,t):
        input = x(t)
        return np.array([y[1], mu * (1.0 - y[0] ** 2) * y[1] * w - y[0] * w ** 2 + w**2 * input])

    return func


def getModelOutputs_numpy(time, x,  mus, w):
    outputs = []
    initVals = [0,1.0]
    for mu in mus:
        f = vanderPol_numPy(time, x, mu, w)
        y = scipy.integrate.odeint(f, initVals, time)[:,0]
        outputs = outputs + [y]

    return outputs

def getModelOutputs(time, x,  mu, ws, fs):
    outputs = []
    for m in mu:
        mus = np.full(len(time), m)
        y = triggerfishvcv.vdpO(x, mus, ws, fs)
        outputs = outputs + [y]

    return outputs

if __name__ == "__main__":
    #print(np.get_include())
    fs = 48000 * 1
    #N = 10 * fs
    N = 10 * fs
    time = np.arange(N) / float(fs)

    fmax = 4200

    x = 0.0*chirp(time, f0 = 0, t1 = time[-1], f1 = fmax)
    #ws = 2* np.pi * np.linspace(fmax/3, fmax/3, len(time))
    ws = 2* np.pi * np.linspace(0, fmax, len(time))
    #ws = 2* np.pi * np.append(np.linspace(0.0, fmax, len(time)/2), np.linspace(fmax, 0.0, len(time)/2))
    #x = np.where(time >= 1.0, x, 0)

    #muVals = np.linspace(1.0,3.0,4)
    muVals = [0.1, 0.5, 3.0,8.53]

    models = ['mu=' + str(m) for m in muVals]
    f = 1200

    ys = getModelOutputs(time, x, muVals, ws, fs)
    MultiPlotBokeh(spectrogramBokeh, fs, time, x, ys, models)
    #MultiPlotBokeh(powerSpectralDensity, fs, time, x, ys, models)

    N = fs * 1
    time = np.arange(N) / float(fs)
    #x = 10* np.where(time >0.0, 0.0, 1.0) # Impulse
    #x = 5* np.where(time >0.0, 1.0, 0.0) # Heaviside
    f = 100
    x = 5.0* np.sin(time * 2 * np.pi * f/10)
    ys = getModelOutputs(time, x, muVals, 2 * np.pi * f)
    MultiPlotBokeh(timeResponse, fs, time, x, ys, models)

