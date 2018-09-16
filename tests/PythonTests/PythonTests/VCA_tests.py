from itertools import islice
from random import random
from time import perf_counter
from scipy import signal
from scipy.signal import *
import triggerfishvcv
import numpy as np
from plotting import *


def getModelOutputs(time,x, models,fs):
    cv = np.full(len(x),1.)
    responses = []
    for m in models:
        vcaMethod = getattr(triggerfishvcv, m)
        y = vcaMethod(x,cv,fs)
        responses += [y]
    return responses



if __name__ == "__main__":
    #print(np.get_include())
    fs = 48000 * 1
    N = 10 * fs
    time = np.arange(N) / float(fs)

    x = 5*chirp(time, f0 = 0, t1 = time[-1], f1 = 24000)
    #x = np.where(time >= 1.0, x, 0)
    models = ["vca_OTA_noOversampling", "vca_OTA_butterworth5", "vca_OTA_cheby7", "vca_OTA_x4_cheby7"]
    models += ["vca_Transistor_noOversampling", "vca_Transistor_butterworth5", "vca_Transistor_cheby7", "vca_Transistor_x4_cheby7"]

    ys = getModelOutputs(time, x,models,fs)

    MultiPlotBokeh(spectrogramBokeh, fs, time, x, ys, models)
    #MultiPlotBokeh(powerSpectralDensity, fs, time, x, ys, models)

    N = 1000
    time = np.arange(N) / float(fs)
    #x = 10* np.where(time >0.0, 0.0, 1.0) # Impulse
    x = 3* np.where(time >0.0, 1.0, 0.0) # Heaviside
    ys = getModelOutputs(time, x,models)
    MultiPlotBokeh(timeResponse, fs, time, x, ys, models)
