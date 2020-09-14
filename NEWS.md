## PAMr 0.9.2

* Lots of prep for CRAN submission

* Changes to `standardClickCalcs` so it doesnt crash if peak frequency is 0, happens
if lower filter bound is 0

* Added `updateFiles` function to update file locations of a study object if you have
changed computers or moved folders around

* Added `checkStudy` function to do some sanity checking for possible issues after
data is processed. Currently only checks for peak freq of 0 which would mean you
probably want a different filter value or add a calibration function

* Added test datasets to inst/extdata to go along with added testthat support

## PAMr 0.9.1

* NEW FUN FUNCTION `calculateAverageSpectra`! Calcultes and plots average spectra from
an event! Woohoo!

* Fixed issue with `export_banter` not finding ancillary measures properly

* `matchEnvData` compatible with update in `PAMmisc` package, works with custom
functions now

* `setSpecies` handles things more gracefully without provided method argument

* Changed `setSpecies` to choose insted of stop when multiple methods provided (default)

## PAMr 0.9.0

* Added `filter` method that works like dplyr's filter for detections in an object. Has
a special case when specifying `species` or `Species` that will filter an
`AcousticStudy` object by the species listed in the `$id` spot of each species slot

* Also added first version of `addAnnotation` function, it asks a lot of questions for now

* `AcousticStudy` `[` accessor changed to return another `AcousticStudy` after subsetting
the events in that study

## PAMr 0.8.3

* Minor bug fix for `writeEventClips`

## PAMr 0.8.2

* Fixed bug in `addGps` that could create duplicate rows after adding gps to your data.
Would occur if two GPS entries occurred within the same second, now adds milliseconds 
before matching so that this is avoided.

* Added support for SoundTrap files for `writeEventClips` with `format='soundtrap'`

## PAMr 0.8.1

* Added better progress bar for `mode='db'` that is based on binary files not databases

* Added `writeEventClips` function for creating wav clips of events. Currently needs
some adjustments but should work for Pamguard data if your selection of wav files
definitely covers all of your events

* Fixed bug where click detectors were being improperly named if a matched classifier was
present in the binary file

## PAMr 0.8.0

* Finally not just a bug fix! New function `matchEnvData` extends a function of the
same name in the "PAMmisc" package (as of PAMmisc v1.4.1). This lets you download
environmental data and match it to your `AcousticStudy` and `AcousticEvent` objects

## PAMr 0.7.18

* Bug in `export_banter` if an event had exactly 1 detection fixed

## PAMr 0.7.17

* Fixed bug in `export_banter` for list of `AcousticStudy` objects not propagated event 
names properly

## PAMr 0.7.16

* `export_banter` can take list of `AcousticStudy` objects now, and fixed a bug where 
NA values were being introduced

## PAMr 0.7.15

* Skip over binaries with no UIDs in `mode='db'`

* `mode='db'` was not properly passing along an ID if you supply it

* `export_banter` now lists which measures were `NA` in the `$na` dataframe output

## PAMr 0.7.14

* First pass fixing `egClicks`

## PAMr 0.7.13

* Error message in `processPgDetections` for `mode='db'` now reports db and binary file
failed on better

* Changed folder choosing funciton for `addBinaries`, previous version was windows-only
and would not properly start in current working directory

* `addDatabase` now prints added databases for reference

## PAMr 0.7.12

* Fixed calibration loading file if first column is just row names

## PAMr 0.7.11

* `setSpecies` now has `method = 'reassign'` that can reassign species id
to new ones provided in a dataframe

## PAMr 0.7.10

* `processPgDetctions` bug fixed when encountering a binary with 0 detections

## PAMr 0.7.9

* `processPgDetections` shouldnt crash on empty binary file that also didnt finish
writing the header and footer for `mode='time'`

## PAMr 0.7.8

* `processPgDetections` handles binary files that never finished writing without
a mysterious crash now. Should only affect `mode='time'`.

## PAMr 0.7.6 and 0.7.7 because oops

* Just adde a debugger for JKs problem

## PAMr 0.7.5

* `standardClickCalcs` changed to have a `filterfrom_khz` and `filterto_khz` option
so that a bandpass filter can be applied. If only a highpass filter is desired simply
leave `filterto_khz` as the default `NULL`

## PAMr 0.7.4

* `addCalibration` reworked slightly again, created calibration function now takes
frequency in Hz as input and outputs dB correction to add instead of taking in
result of seewave::spec and outputting corrected spec. This makes it more flexible
and means we aren't calculating spectrum twice for no reason

* `standardClickCalcs` adjusted to use new calibration, and temporarily has a new
parameter calculated - `dBPP` - which is the peak to peak value as traditionally 
defined in acoustics. This is mostly for testing comparability to other Matlab
based code of new calibration stuff. Also the spectrum is never rescaled to a
max of 0 like it was previously, hopefully should represent accurate SPL levels
with new calibration code

## PAMr 0.7.3

* `addCalibration` now asks for units of your calibration, and should properly add
or subtract this value from calculated spectrum depending on units selected

## PAMr 0.7.2

* `export_banter` now has an option to split your data into training and test sets by 
specifying a proportion for the `training` argument, and now puts out a pretty table
showing a summary of your events and species

## PAMr 0.7.1

* `processPgDetections` has some minor updates to improve clarity when binary files or 
databases aren't found using `mode='time'`

## PAMr 0.7.0

### Functions Renamed

* `getPgDetections` has now been re-named to `processPgDetections` to more accurately
reflect what it is doing. Functions starting with `get___` will be used just for 
accessing data

* `showWaveform` and related functions have been renamed to `plotWaveform`

* General naming consistency overhaul - `sr` and `db` should now be used in place of
`sampleRate` and `database` wherever these were previously used, and this should be
the naming convention going forward. This means included functions will need to be 
re-added (they previously relied on it being named `sampleRate`), and any custom 
functions should expect data to contain the item `sr` and not `sampleRate`

### Major Changes

* `AcousticStudy` class has been added. This will now be the class of object
returned by `processPgDetections`, it stores your list of `AcousticEvent` objects with
other important data. 

* Any data processed in previous versions will need to be re-run, and all functions in
older `PAMrSettings` files will need to be removed and re-added 

* New Function: `calculateICI` calculates the inter-click interval of click data, and 
stores it in the `ancillary` slot of an `AcousticEvent`. 

* New Function: `plotDataExplorer` creates an interactive plot to allow users to 
explore the data in an `AcousticStudy`, `AcousticEvent`, or `data.frame`. Plot
can be colored or facetted by any columns that are characters or factors, and any
numeric columns can be plotted with `ggplot`'s `geom_density`

* New Function: `getDetectorData` gathers all detector data into single data frames
for each detector type (`'click'`, `'whistle'`, and `'cepstrum'`) for easier manipulation

### Minor Changes

* Detector dataframes in `AcousticEvents` now have a `'calltype'` associated with them,
currently one of '`whistle'`, `'click'`, or `'cepstrum'`. This allows for future fun stuff
to happen

* `getBinaryData` will now attempt to get the appropriate sample rate for each
data point, either from the settings or matching by time using the database file
if more than one sample rate was in your data.

* Many speed improvements - click calculations should take about half the time, and
`processPgDetections` with `mode='time'` will now skip over binaries that are outside
of the time range of specified events

* `getPgDetections` no longer has an option `mode='all'`, and it is no longer necessary
to supply a sample rate for `mode='time'`. It will now attempt to match events to a
correct database used on time stamps, the match the appropriate sample rate the same
way that `mode='db'` does. If a manual sample rate needs to be applied it should now be
entered into the `grouping` table

* Functions `addBinaries` and `addDatabase` can now add files from another `PAMrSettings`
object, and will report how many files have been added when finishing

* `export_banter` will now attempt to find event level measures in the `'measures'` item in
the `ancillary` slot of each `AcousticEvent`. 

* `export_banter` no longer has a `reportNA` option, any `NA` values that are removed are now
always saved to a separate `'na'` item in the list output. The function also reports how many
detections were processed.

* `addGps` now also stores all of the gps data loaded into the `gps` slot of the `AcousticStudy`
instead of only adding coordinates to detections

## PAMr 0.6.8

* `getPgDetections` will name events with database appended for `method = 'db'` instead
of just event ID number to ensure uniqueness across multiple databases

## PAMr 0.6.7

* I don't remember what happened here. I bet it was important

## PAMr 0.6.6

* `setSpecies` can now use a dataframe for `method = 'manual'`, and has
a top secret option for SR

## PAMr 0.6.5

* Updated `export_banter` with options to exclude certain species and to
export data without species codes to use for prediction only instead
of training a banter model

## PAMr 0.6.4

* Better error tracking when functions cause `getPgDetections` to crash

* Now reads in angles and angleErrors from click data

* `export_banter` allows you to specify certain columns to not export

* `'time'` mode for `getPgDetections` will now report a sample event time
so you can see if times are being converted properly from your csv before
proceeding with calculations

## PAMr 0.6.3

* Fixed an issue where `getPgDetections` would not work if both Detection Group Localizer
and Offline Click events were present in a database.

* Renamed `eventType` and `Text_Annotation` columns from event databases to `eventLabel`
within detection dataframes so there is consistency between the two

* Removed some unnecessary columns from detection dataframes, including `detectorName`,
`sampleRate`, `Id`, `parentUID`, and `comment`

* Fixed a bug in how `getBinaryData` was checking for multiple matches on the same UID

## PAMr 0.6.2

* Fixed an issue with repeated entries in click detections for modes other than `'db'`

## PAMr 0.6.1

* Sometimes whistles would not get proper decimated sample rate, fixed

## PAMr 0.6.0

* Added an `id` slot to `AcousticEvent` objects. Note that this will cause existing `AcousticEvent`
objects to behave poorly until they have their `id` slot created / updated using the new
`setIdSlot` function.

* Added `setIdSlot` function to update older `AcousticEvent` objects to the new format

## PAMr 0.5.9

* Fixed bug in SR / FFT parameter calculation in whistles if there was a gap in
the whistle contour

## PAMr 0.5.8

* Changed event naming for `mode='time'` in `getPgDetections`. Will now only append
numbers if event is not unique, and will also insert an underscore before the number

* `export_banter` now properly checks for cases when there are no detectors in an event
instead of crashing confusingly

## PAMr 0.5.7

* Rocca whistle calcs minor change - boundary settings for sweeps set to "flat",
should reduce number of inflection points

## PAMr 0.5.6

* minor change to `export_banter` no longer using list names to index and create unique event names

## PAMr 0.5.5

* `export_banter` now removes NA rows, and has reportNA option to see which ones are NA

* Detector names have spaces replaced with _ to avoid weird issues later

## PAMr 0.5.4

* Dealing with zero detection events better for `export_banter`

* Temporary fix for click calculations with lots of zeroes in wave form - no more NA

## PAMr 0.5.3

* `standardClickCalcs` now supports manual input of sample rate. Default argument is `'auto'`,
which will read from the database. User can supply a numeric value for sample rate in hertz instead.

## PAMr 0.5.2

* Added a check in `addDatabase` to see if all files are actually .sqlite3 databases

* Added a `tryCatch` in `getPgDetections` for mode `db` so that it shouldn't 
stop completely when encountering an error and lose all previously analysed DBs

## PAMr 0.5.1

* Minor bug fix when using `seewave::spec`, it can produce NA values for the frequency
if the input wave is long. Adjusted parts `standardClickCalcs` and `addCalibration` to
work around this.

## PAMr 0.5.0

* `getPgDetections` changed to work by specifying a `mode` as an argument instead of
calling separate functions. Can now create events using a csv or dataframe with 
start and end times specified.

* Changed `export_banter` to export a list with named item `detectors` instead of 
`detections`, no other functional change.

## PAMr 0.4.0

* Added `setSpecies` functions for assigning species classifications to 
AcousticEvents.

* Changed `getDbData` to default to looking for both OfflineEvents tables 
and DetectionGroupLocaliser tables, will still only load one if a specific
type is provided for `grouping`.

## PAMr 0.3.0

* Added `addCalibration`, `applyCalibration`, and `findCalibration` functions,
as well as a `plot` method that will show the calibration function used. See
the  *Calibration* section above for more details.

* Fixed a bug in `removeFunction` that would cause the incorrect number of
functions to show in certain cases, and that would cause all functions to
be removed when only one was selected.

* `standardClickCalcs` has been adjusted to work with the new calibration
methods. See *Calibration* section above for more details.

## PAMr 0.2.2

* Added `addGps` function for adding matching GPS data to your detections. This 
allows you to supply a dataframe of Lat/Long locations with timestamps to 
match to your detections.

* Added `showWaveform`, `showSpectrogram`, and `showWigner` functions that
allow you to easily plot the waveform, spectrogram, or wigner plot of a detection
in an `AcousticEvent` object by selecting the UID(s) you want to investigate further.
`getBinaryData` is also added as a helper function for these, lets you easily
get the binary file data for a single detection.

* `standardClickCalcs` now supports supplying a `Wave` class object as input

## PAMr 0.2.0

* Rocca (`roccaWhistleCalcs`) and cepstrum (`standardCepstrumCalcs`) 
functions added. These are also added by default to a new PRS.

* changed `AcousticEvent` class slot name from `specClass` to `species`

