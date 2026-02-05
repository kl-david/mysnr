import numpy as np
import scipy

def bandpower(ps, mode='psd', tau=1):
    """
    estimate bandpower, see https://de.mathworks.com/help/signal/ref/bandpower.html
    """
    if mode=='time':
        x = ps
        l2norm = linalg.norm(x)**2./len(x)
        return l2norm
    elif mode == 'psd':
        return sum(ps)/tau # we need to normalize to signal time, i.e. nr_of_samples * sampling_time

def getIndizesAroundPeak(psd, peakIndex,searchWidth=1000):
    '''
    takes in a psd and finds all monotonicly falling bins around it
    i.e. all bins that belong to a certain peak, including bleading
    works entirely bin-based
    '''
    peakBins = []
    magMax = psd[peakIndex]
    curVal = magMax
    for i in range(searchWidth):
        newBin = peakIndex+i
        newVal = psd[newBin]
        if newVal>curVal:
            break
        else:
            peakBins.append(int(newBin))
            curVal=newVal
    curVal = magMax
    for i in range(searchWidth):
        newBin = peakIndex-i
        newVal = psd[newBin]
        if newVal>curVal:
            break
        else:
            peakBins.append(int(newBin))
            curVal=newVal
    return np.array(list(set(peakBins)))

def freqToBin(fAxis, f):
    '''
    finds the closest bin in a frequency axis
    '''
    return np.argmin(abs(fAxis-f))

def getPeakInArea(psd, faxis, estimation, searchWidthHz = 10):
    """
    returns bin and frequency of the maximum in an area
    """
    binLow = freqToBin(faxis, estimation-searchWidthHz)
    binHi = freqToBin(faxis, estimation+searchWidthHz)
    peakbin = binLow+np.argmax(psd[binLow:binHi])
    return peakbin, faxis[peakbin]

def getHarmonics(fund,bw,nHarmonics=6,aliased=False):
    '''
    return the n Harmonics of a fundamental frequency
    returns actual frequencies, not bins
    '''
    harmonicMultipliers = np.arange(2,nHarmonics+2)
    harmonicFs = fund*harmonicMultipliers
    if not aliased:
        harmonicFs[harmonicFs>bw] = -1
        harmonicFs = np.delete(harmonicFs,harmonicFs==-1)
    else:
        nyqZone = floor(harmonicFs/(bw))
        oddEvenNyq = nyqZone%2  
        harmonicFs = mod(harmonicFs,bw)
        harmonicFs[oddEvenNyq==1] = (bw)-harmonicFs[oddEvenNyq==1]
    return harmonicFs 

def findBinnedHarmonics(psd, faxis, fund, bw, nHarmonics=6, aliased=False):
    '''
    finds all bins belonging to a harmonic
    returns bins, not frequencies
    '''
    harmonicFs = getHarmonics(fund, bw, nHarmonics=nHarmonics, aliased=aliased)

    harmonicBorders = np.zeros([2, nHarmonics], dtype=np.int16).T
    fullHarmonicBins = np.array([], dtype=np.int16)
    fullHarmonicBinList = []
    harmPeakFreqs=[]
    harmPeakBins=[]

    for i, harmonic in enumerate(harmonicFs):
        searchWidth = 0.1*fund
        estimation = harmonic

        harmBin, harmFreq = getPeakInArea(psd, faxis, estimation, searchWidth)
        harmPeakFreqs.append(harmFreq)
        harmPeakBins.append(harmBin)

        allHarmBins = getIndizesAroundPeak(psd, harmBin)
        fullHarmonicBins = np.append(fullHarmonicBins, allHarmBins)
        fullHarmonicBinList.append(allHarmBins)

        harmonicBorders[i,:] = [allHarmBins[0], allHarmBins[-1]]
    
    return fullHarmonicBins

def findFundamental(psd, faxis, bw=None):
    '''
    finds the highest peak in psd (up to bw)
    returns the bin-aligned frequency and the bin
    '''
    if bw == None:
        fbin = np.argmax(psd)
    else:
        bwbin = freqToBin(faxis, bw)
        fbin = np.argmax(psd[:bwbin])
    return faxis[fbin], fbin

def snr(psd, faxis, T_samp, N, bw=None) -> float:
    '''
    Returns the snr of a psd.

    A power spectral density (psd) and its frequency axis (faxis) are
    separated into signal, harmonic and noise components.
    
    Signal and harmonic peaks are identified via a monotonicity-criterion.

    Bins identified as part of a harmonic are replaced by the mean of the (non-harmonic) noise.

    Returns signal power divided by noise power
    '''
    if bw == None:
        bw = faxis[-1]
    tau = T_samp * N
    fund, fund_bin = findFundamental(psd, faxis, bw)

    bw_bin = freqToBin(faxis, bw)

    # get bins of signal and harmonics
    signal_bins = getIndizesAroundPeak(psd, fund_bin)
    harm_bins = findBinnedHarmonics(psd, faxis, fund, bw)

    # estimate signal power
    signal_power = bandpower(psd[signal_bins], tau=tau)
    
    # estimate noise power:
    psd_noise = np.copy(psd)
    psd_noise[signal_bins] = 0
    psd_noise[harm_bins] = 0
    noise_mean = np.median(psd_noise[psd_noise != 0])
    psd_noise[signal_bins] = noise_mean
    psd_noise[harm_bins] = noise_mean
    noise_power = bandpower(psd_noise[:bw_bin], tau=tau)

    return signal_power/noise_power

def sndr(psd, faxis, T_samp, N, bw=None):
    '''
    Returns the sndr of a psd.

    A power spectral density (psd) and its frequency axis (faxis) are
    separated into signal, harmonic and noise components.
    
    Signal and harmonic peaks are identified via a monotonicity-criterion.

    Returns signal power divided by the sum of noise power and harmonic power.
    '''
    if bw == None:
        bw = faxis[-1]
    tau = T_samp * N
    fund, fund_bin = findFundamental(psd, faxis, bw)

    bw_bin = freqToBin(faxis, bw)

    # get bins of signal and harmonics
    signal_bins = getIndizesAroundPeak(psd, fund_bin)
    harm_bins = findBinnedHarmonics(psd, faxis, fund, bw)

    # estimate signal power
    signal_power = bandpower(psd[signal_bins], tau=tau)

    # estimate harm power
    harm_power = bandpower(psd[harm_bins], tau=tau)
    
    # estimate noise power:
    psd_noise = np.copy(psd)
    psd_noise[signal_bins] = 0
    psd_noise[harm_bins] = 0
    noise_mean = np.median(psd_noise[psd_noise != 0])
    psd_noise[signal_bins] = noise_mean
    psd_noise[harm_bins] = noise_mean
    noise_power = bandpower(psd_noise[:bw_bin], tau=tau)

    return signal_power/(noise_power + harm_power)