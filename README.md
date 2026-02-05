# MySNR Package

MySNR is a Python package designed to calculate Signal-to-Noise Ratio (SNR) and Signal-to-Noise and Distortion Ratio (SNDR) from Power Spectral Density (PSD). This package provides functions for analyzing frequency data, estimating signal power, harmonic power, and noise power, and computing SNR and SNDR metrics.

## Installation

You can install MySNR using pip:

```bash
pip install mysnr
```

## Usage

### Importing the Package

To use the MySNR package, you first need to import it into your Python script or Jupyter notebook.

```python
import mysnr
```

### Example Code

Here's an example of how to use the MySNR package:

```python
import numpy as np
from scipy.signal import periodogram
import matplotlib.pyplot as plt

# Generate a sample signal with noise and harmonics
N = int(1E6)    # total number of samples
T_samp = 100E-9 # sampling interval
f = 3.0E3       # fundamental signal frequency

signal_power = 0.5
harmonic_power = 0.045
noise_power = 0.1

t = np.linspace(0.0, N*T_samp, N, endpoint=False)                           # time vector
y = np.sqrt(signal_power * 2) * np.sin(f * 2.0 * np.pi*t)                  # base signal, power = 0.5 * 1**2 = 0.5
y += np.sqrt(harmonic_power * 2) * np.sin(2*f * 3.0 * np.pi *t)             # third harmonic, power = 0.5 * 0.3**2 = 0.045
y += np.random.normal(0, np.sqrt(noise_power), len(t))                      # add white noise over the full band

# Compute the PSD using periodogram
py_xf, py_yf = periodogram(x=y, fs=1/T_samp, window='hann', scaling='density')

# Find harmonics and calculate SNR and SNDR
harmonic_bins = mysnr.findBinnedHarmonics(py_yf, py_xf, f, 0.5/T_samp, 3)
signal_peak_bin = mysnr.freqToBin(py_xf, f)
signal_bins = mysnr.getIndizesAroundPeak(py_yf, signal_peak_bin)

print(f"Expected Signal Power: {signal_power}")
print(f"Calculated Signal power: {mysnr.bandpower(py_yf[signal_bins], tau=T_samp*N)}")

print(f"Expected Harm Power: {harmonic_power}")
print(f"Calculated Harmonic Power: {mysnr.bandpower(py_yf[harmonic_bins], tau=T_samp*N)}")

print(f"Expected SNR: {10*np.log10(signal_power/noise_power)}")
print(f"Calculated SNR: {10*np.log10(mysnr.snr(py_yf, py_xf, T_samp=T_samp, N=N, bw=None))}")

print(f"Expected SNDR: {10*np.log10(signal_power/(noise_power + harmonic_power))}")
print(f"Calculated SNDR: {10*np.log10(mysnr.sndr(py_yf, py_xf, T_samp=T_samp, N=N, bw=None))}")

# Plot the PSD
fig, fftax = plt.subplots()
py_yfdb = 10*np.log10(np.abs(py_yf))

fftax.semilogx(py_xf, py_yfdb)
fftax.semilogx(py_xf[harmonic_bins], py_yfdb[harmonic_bins])
fftax.semilogx(py_xf[signal_bins], py_yfdb[signal_bins])

plt.grid(which='minor', linestyle=':', linewidth=0.5, color='lightgray')  # Minor grid lines
plt.show()
```

## Functions

- **bandpower(psd, mode='psd', tau=1)**: Estimate the bandpower from a power spectral density.
- **getIndizesAroundPeak(psd, peakIndex, searchWidth=1000)**: Find indices around a peak in the PSD.
- **freqToBin(fAxis, f)**: Convert a frequency to its corresponding bin index.
- **findBinnedHarmonics(psd, faxis, fund, bw, nHarmonics=6, aliased=False)**: Find binned harmonics of a fundamental frequency.
- **snr(psd, faxis, T_samp, N, bw=None)**: Calculate the Signal-to-Noise Ratio (SNR).
- **sndr(psd, faxis, T_samp, N, bw=None)**: Calculate the Signal-to-Noise and Distortion Ratio (SNDR).

## Contributing

Contributions are welcome! If you have any suggestions or bug reports, please open an issue on the [GitHub repository](https://github.com/kl-david/mysnr).