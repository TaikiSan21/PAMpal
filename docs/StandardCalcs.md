Items that are **bolded** are parameters output by the function, items in *italics* are input parameters to the function.

## standardClickCalcs
The process for calculating the included click parameters is as follows:

First a Butterworth filter is applied, either a highpass filter starting at *filterfrom_khz*, or bandpass if *filterto_khz*  is also specified (bandpass was added in more recent versions of PAMpal, previously only a highpass was available, uses additional parameter *filterto_khz*)

**peakTime** is the time (seconds) after the start of the click clip where the maximum value of the waveform is located.

After **peakTime** is calculated, the signal is shortened to be of size equal to that specified by *winLen_sec*. The windowed clip is created by centering the window around the maximum value of the waveform, so if the original clip was 1000 samples with a peak at 250 samples, a 300 sample clip would be created from samples 100 to 399. This shortened clip is used for all further calculations.
Next the Teager Kaiser energy of the signal is calculated. The **noiseLevel** is the median of the TK energy (reported in dB, 10*Log10(TKEnergy). duration is defined by counting the number of samples above 100 times the 40th Percentile of the TKenergy level, duration is reported in microseconds. TK energy calculations follow methods in Soldevilla et al 2008

**dBPP** is calculated as 20 * log10 of the difference between the maximum and minimum value of the waveform. If a calibration function has been supplied, then the calibration value at the peak frequency will also be added to this value.

Next we calculate the spectrum of the shortened signal clip from earlier using FFT length equal to *winLen_sec*, which at this point is the same as the length of the clip. We use the `spec` function from the `seewave` package. Signal spectrum values are converted to dB (20*log10), and then if a calibration function has been supplied it is applied by fitting a GAM to the sensitivity curve of the calibration file then adding these values to the spectrum (Methods C in Griffiths et al 2020).

This calibrated spectrum is used for the next set of calculations. First peak and “trough” frequencies of the spectrum are calculated using the the function `peakTrough` in the package `PAMmisc`. This attempts to find multiple peaks in the spectrum, and if multiple peaks are found it will find the “trough” or valley point between them. To do this, first the spectrum is smoothed slightly using a local rolling average of 5 points. Then the frequency value with the highest dB level is defined as the **peak**. From here, potential candidates for a second peak are identified. A candidate point must be greater than both its neighbors, it must be more then 10kHz from the first peak, but not greater than 30kHz, and its dB value must be no more than 15dB below the value of the peak frequency. The point with the highest dB value that meets these criteria is labeled **peak2**, and then this process is repeated for **peak3**. If No points meet these criteria at either step, then the values will be set to 0. The value of **trough** is the frequency value between peak and peak2 with the lowest dB value, and **trough2** is similarly calculated between **peak2** and **peak3**. If **peak2** or **peak3** are 0, then **trough** and/or **trough2** will also be 0. **peakToPeak2** is the difference between the frequency values of **peak** and **peak2**, if **peak2** is 0 then this will also be 0. **peakToPeak3** and **peak2ToPeak3** are calculated similarly.

The calibrated spectrum is then used to calculate six measures at the -3 and -10 dB threshold values. The minimum frequency (**fmin_3dB**), maximum frequency (**fmax_3dB**), frequency bandwidth (**BW_3dB**),  resonant quality factor (**Q_3dB**), and center frequency (**centerHz_3dB**). From Griffiths et al 2020: Q estimates “the frequency pureness of a time wave at a specific dB level. Q is calculated by dividing the center frequency by the bandwidth, such that a higher Q indicates a lower rate of energy loss relative to the stored energy of the resonator”.

Clicks may also have parameters **angle** and **angleError**, but these are read in directly from the Pamguard binary files and are not calculated by this function.


Griffiths, E. T. et al (2020) “Detection and classification of narrow-band high frequency echolocation clicks from drifting recorders” J. Acoust. Soc. Am. 147 3511-3522

Soldevilla, M. et al (2008) “Classification of Risso’s and Pacific white-sided dolphins using spectral properties of echolocation clicks” J. Acoust. Soc. Am. 124 609-624

## roccaWhistleCalcs

Whistle calculations are a reimplementation of the ROCCA calculations currently present in Pamguard with
permission from Julie and Michael Oswald. More details can be found in:

Oswald et al (2007) "A tool for real-time acoustic species identification of delphinid whistles", J. Acoust. Soc. Am. 122 587

## standardCepstrumCalcs
**ici** is the median ICI value across the detected contour.
**duration** is the length in seconds of the contour
**iciSlope** is the slope of the ICI contour, given by fitting a line to the contour data using the `lm` function