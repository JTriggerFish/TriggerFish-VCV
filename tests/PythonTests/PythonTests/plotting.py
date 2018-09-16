from itertools import islice
from random import random
from time import perf_counter
from scipy import signal
from scipy.signal import *
import matplotlib.pyplot as plt
import numpy as np
from bokeh.plotting import figure, show, output_file, gridplot
import bokeh
import bokeh.models
import bokeh.models.mappers

def spectrogram(name, time, x, y, fs, i):
    f, t, Sxx = signal.spectrogram(y, fs ,  window = signal.get_window('blackman',1024))
    Sxx = 10*np.log(Sxx) / np.log(10)
    min_display_dB = -80
    #max_display_dB = 0
    plt.subplot(2,2,i)
    plt.pcolormesh(t, f, Sxx, vmin= min_display_dB)#, vmax = max_display_dB)
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.gca().set_title(name + " " + str(fs))

def spectrogramBokeh(name, time, x, y, fs):
    f, t, Sxx = signal.spectrogram(y, fs ,  window = signal.get_window('blackman',1024))
    Sxx = 10*np.log(Sxx) / np.log(10)
    min_display_dB = -80
    #max_display_dB = 0

    p = figure(x_range = (t[0], t[-1]), y_range = (f[0], f[-1]), plot_width = 600, plot_height = 400,
              title = name + " %d" % fs,
              x_axis_label = "Time [sec]",
              y_axis_label = "Freq [Hz]")
             
    p.image(image = [Sxx], x=0, y=0, dw = t[-1], dh = f[-1], color_mapper =bokeh.models.mappers.LinearColorMapper(palette = "Magma11", low = -80))

    return p

def timeResponse(name, time, x, y, fs):
    p = figure(plot_width = 600, plot_height = 400,
              title = name + " %d" % fs,
              #x_axis_label = "Time [sec]",
              x_axis_label = "Sample",
              y_axis_label = "VCA out")
             
    p.multi_line([time*fs,time*fs], [x,y], color=["firebrick", "navy"], alpha=[0.3, 0.9])
    #p.line(time,y)

    return p

def freqResponse(name, time, x, y, fs):
    n = len(y)
    freqs = np.arange(n) * float(fs /n)
    freqs = freqs[range(int(n/2))]
    Y = np.fft.fft(y) / n
    Y = Y[range(int(n/2))] # Keep positive freqs only
    amplitude = 20*np.log10(np.abs(Y))

    p = figure(plot_width = 600, plot_height = 400,
              title = name + " %d" % fs,
              #x_axis_label = "Time [sec]",
              y_range=(0,22000),
              x_axis_label = "Freq [Hz]",
              y_axis_label = "amplitude [dB]")
             
    p.line(freqs, power)

    return p

def powerSpectralDensity(name, time, x, y, fs):
    f, Pxx_den = signal.welch(y, fs, nperseg=1024)
    Pxx_den = 10*np.log10(Pxx_den)

    p = figure(plot_width = 600, plot_height = 400,
              title = name + " %d" % fs,
              #x_axis_label = "Time [sec]",
              x_axis_label = "Freq [Hz]",
              x_axis_type = "log",
              y_axis_label = "Psd [dB]")
             
    p.line(f, Pxx_den)

    return p


def MultiPlotBokeh(plotFunction, fs, time, x, yArray, namesArray):
    output_file(plotFunction.__name__ + ".html")

    numCols = 2
    rows = []
    i=0
    row = []
    for y,name in zip(yArray, namesArray):
        fig = plotFunction(name, time, x, y, fs)
        row  = row + [fig]
        i +=1 
        if i >= numCols:
            i = 0
            rows = rows + [row]
            row = []
    if len(row) >0:
        rows = rows + [row]

    p = gridplot(rows)
    show(p)
