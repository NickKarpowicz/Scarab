# <p style="text-align: center;">Scarab</p>
<p style="text-align: center;">Spectrometer Control and Recording for Attosecond Beamlines</p>
<p style="text-align: center;">Nick Karpowicz</p>
<p style="text-align: center;">Max Planck Institute of Quantum optics</p>

<p style="text-align: center;"><img src="Documentation/screenshot.png"></p>

This is a simple c++ application for controlling spectrometers from Ocean Insight (nee Ocean Optics).

GUI with GTK4, should work on Linux, Mac, and Windows, but only tested on Windows so far.


# Quick start guide
To launch the program, run "Scarab.exe" in the "bin" folder of the unzipped archive.

The controls are quite basic at the moment, so let's just look at what you see in the screenshot above:

## Spectrometer selector
This is the pulldown menu that reads 1:USB4000 in the screenshot. If you have multiple spectrometers connected to your PC, it will let you switch between them

## Mode selector

This is the pulldown menu that reads "Spectrometer live (nm)" in the screenshot. It holds the different view an operation modes you can put the program in, including:
- _Spectrometer live (nm)_: view the live spectrum coming from the spectrometer, with a wavelength x-axis.
- _Spectrometer live (THz)_: view the live spectrum with a THz x-axis.
- _All spectrometers_: simultaneously show the live output of all connected spectrometers (THz-axis only). Note that this mode won't set exposure time, you have to do that in the individual spectrometer's live view.
- _Interferometry: Spectrum_: View the live spectrum and the reference spectra for performing white-light interferometry.
- _Interferometry: Time_: View the Fourier-transformed interferometry data, allowing you to pick an appropriate filter for processing the data
- _Interferometry: Phase_: View the interferometer's retrieved spectral phase as it accumulates.
- _Interferometry: Group delay_: View the group delay curve determined by the spectral phase as it accumulates.

## Run button
This toggles the spectrometer live view on and off.

## Line  color controls
Under the Mode selector, there are 12 text boxes. These let you put in red, green, and blue values (from 0 to 1) to set the colors of the main data line (top row) and the overlays (next 3 rows).

## Acquire dark spectrum (🕯️)
Press this button to acquire a dark spectrum, which will automatically be subtracted from the spectra you acquire

## Delete dark spectrum (🗑️)
Press the trash can next to the candle to delete your dark spectrum.

## Acquire overlay (📈)
Press this button to acquire an overlay spectrum to plot together with your live spectrum. There are three slots for overlays, with their own acquire 📈 and delete 🗑️ buttons.

## Output file path
This reads "DefaultOutput.txt" in the screenshot. This is where the data will be saved when you press one of the file saving buttons (Acquire and Save, more on them in a bit).

## File path search button (...)
Press this to open a dialog box to set your save path.

## Exposture (ms)
This reads 40 in the screenshot. This is the exposure time of the spectrometer - you can set it in the live views.

## Acquire button
Press this and the program will take a series of spectra and save them to the file you indicate in the output file path. The number of acquisitions is set by the text box directly below it (reading 10 in the screenshot). You can optionally set length of time to pause between captures with in the box to its left (reading 0 in the screenshot).
The file will be in plain text with the following columns:
1. Wavelength (nm)
2. First spectrum
3. Second spectrum

... and so on

## Timestamp checkbox (⌚)
If you check this, a timestamp will automatically be appended to the filename when you save. This way you can just keep mashing the Acquire button when you want to save a bunch of measurements, and your files will have unique and sequentially-ordered names. The timestamp is milliseconds since the Unix epoch.

## Frequency grid controls
The spectral interferometry mode is done on a uniformly spaced frequency grid, which you set up here. The first box, which reads 2048 in the screenshot is the number of points in the grid (note: a power of two may give slightly better performance because there are a lot of FFTs of this size). The next two boxes are the minimum and maximum frequency of the grid. If they are left empty, the grid will automatically span the full range of the spectrometer data. If you want to set them manually (better performance when the interference only occupies a small range of the spectrum), the values are in THz.

## Interferometry filter (Filter (t0, sig., ord))
This filter applies a supergaussian filter to the data in the pseudo-time domain (FFT of the spectral interferogram). The filter should be placed such that it encapsulates the interferogram data but excludes artifacts in the spectrum -- go to the "Inferferometry: Time" view to see what you're doing. The three boxes are for the central time, the width, and the order of the supergaussian, respectively.

## Interferometry references - Ref. A, Ref. B buttons
Press these buttons to acquire a reference spectrum from the two arms of your interferometer. The references will be automatically averaged over the number of acquisitions set for the acquisition mode - since these will be used in every phase you calculate, it worth taking a few to reduce the noise.
Procedure:
1. Block arm B.
2. Press "Ref. A"
3. Block arm A.
4. Press "Ref. B" Now you should have the required references to do interferometry. You'll have to re-take them if you change your exposure time. Next steps:
5. Go to the "Interferometry: Time" view to setup your filter.
6. Go to the "Interferometry: Phase" or "Interferometry: Group delay" view to start acquiring the data.

## Reset button
Once you're acquiring phase, it will take a running average of all the values. If you want to reset the running average (e.g. if you change something and want to see what it did) press this.

## Save button
This will save the interferometry data.
The columns are:
1. Frequency (THz)
2. Phase (rad)
3. Standard deviation of phase (rad)

## Plot controls
Beneath the plot, there are a few more buttons and controls. These are:
- _SVG button_: save a .svg file of the current plot
- _xlim button, and next two text boxes_: set the x-axis limits manually. If empty, data will be autoscaled.
- _ylim button, and next two text boxes_: set the y-axis limits manually. If empty, data will be autoscaled.
- _Log checkbox_: Plot on a log scale. Note that if you check this, you will probably need to set a the y-limits manually (the lower limit should be greater than zero!).

