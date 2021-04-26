# fft_sine_wave.py
# This example file demonstrates FFT analysis of a signal that is generated using Python functions

import matplotlib.pyplot as plt
import scipy.fftpack
import numpy as np
import math, cmath
import time

# Generate test signal
N = 500 # set the number of point for the signal to be analyzed
freq = 0.02 # this is the frequency of a sine wave to be generated
phase = math.pi/2 # phase of the sine wave

dt = 1  # time increment for the signal
t = [dt*i for i in range(N)] # time data
x = [math.sin(2*math.pi*t[i]*freq + phase) for i in range(N)] # Build the sine wave data set

# FFT analysis
# use the scipy.fftpack.fft
Xk = scipy.fftpack.fft(x) # Note that Xk is a complex number
Nk = len(Xk) # the number of total points; note only Nk/2 are useful
print(Nk) # echo the Nk

# You want to now determine the spacing between the points in the frequency domain 'delta f'
# This is an important relation to remember because it defines the resolution of the
# frequency spectrum in terms of the two KEY parameters in the data set, the time interval and N
df = 1/(dt*N)
# generate the discrete frequency
fk = [k*df for k in range(N)]
print(df) # echo the df (my resolution in the frequency domain)

# Now, find the real and imaginary parts of the FFT result, Xk
ReXk = [Xk[i].real for i in range(N)]
ImXk = [Xk[i].imag for i in range(N)]

# Now, here's how you can calculate a PSD = power spectral density
# we call this one a "one-sided PSD"
# NOTE: sometimes the term autospectral density is also used for a signal PSD
# here 'auto' basically refers to how values at each frequency for a given signal
# relate to values at other frequencies
# this distinguishes the 'auto' from the 'cross-spectral density'
# which is one where you examine the relation between two signals

# estimate the PSD (or autospectral density)
Gxx = [(2.0/N) * abs(Xk[i]) for i in range(N)]

# Plot the results
# plot the sine wave
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
ax1.plot(t,x)
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel("X")
# ax1.set_ylim([-ymax,ymax])

# plot the real part
ax2.plot(fk, ReXk, label='Real Xk', alpha=0.5)
ax2.legend()
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel("ReXk")

# plot the impaginary part
ax3.plot(fk, ImXk, label='Imag Xk', alpha=0.5)
ax3.legend()
ax3.set_xlabel('Frequency (Hz)')
ax3.set_ylabel("ImXk")

# plot the PSD
ax4.plot(fk, Gxx, label='Gxx', alpha=0.5)
ax4.legend()
ax4.set_xlabel('Frequency (Hz)')
ax4.set_ylabel("Gxx")
# ax2.set_ylim([-1,psdmax])
# ax2.set_xlim([-1,fmaxp])

fig.tight_layout()
plt.show()