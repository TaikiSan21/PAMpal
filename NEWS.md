## PAMpal 0.20.1

- Exporting `timeJoin` from `matchTimeData` so other packages can use

## PAMpal 0.20.0

- Adding support for Click Train Detector with `grouping='clickTrain'` when processing
with `mode='db'`. Events have prefix CT instead of OE or DGL. By default `mode='db'` will
process all three of `grouping=c('event', 'detGroup', 'clickTrain')`

## PAMpal 0.19.5

- Adding `sampleDetector` function to subsample by detector

## PAMpal 0.19.4

- Adding `addFPOD` and `getFPODData` functions to deal with adding/getting FPOD 
data from AcousticStudies

## PAMpal 0.19.3

- Updating `getClipData` and `writeEventClips` to have a `fixLength` option that 
allows users to create clips of a fixed length using just the `buffer` argument

- Allowing `processPgDetections` to take multiple `label` arguments to specify
different table columns if using both DG and OE events in same db

## PAMpal 0.19.2

- Updating `markAnnotated` to give warnings on NA values instead of failing

- Updating `roccaWhistleCalcs` to avoid divide by zero NA/Inf values

## PAMpal 0.19.1

- Adding `markAnnotated` function to flag whether or not detections are within annotation
boxes

## PAMpal 0.19.0

- Adding `matchTimeData` function to add arbitrary data to AcousticStudies

- Bug fixes for some random stuff

## PAMpal 0.18.1

- Updated to warnings in `addGps` to now be less overwhelming

- Fixed bug for displaying some notes after combining with `bindStudies`

- `addGps` now adds `gpsUncertainty` which says how many meters since to
the closest GPS point. Useful to catch points where we have interpolated a lot.

## PAMpal 0.18.0

- Adding `addAnnotation` and related functions for annomate exporting

- Explicitly marking joins that are expected to have multiple matches re:
`dplyr::join`'s new behavior

## PAMpal 0.17.2

- Smol change to not grab unnecessary columns in `addGps`

- Big speed improvements for processing with `mode='time'` and `mode='recording'`

## PAMpal 0.17.1

- Converting all paths to use '/' instead of '\\'

- Found out `updateFiles` was unuseably slow, so I fixed it. I dunno, maybe 1000x faster
depending on how many files you have? 

## PAMpal 0.17.0

-   Removing the `$source` item in the `settings` slot for each `AcousticEvent`. These were not used for anything, and when the source was wav files it would end up with a huge repetitive list for every event. Ended up taking up a huge amount of memory, removing this results in `AcousticStudy`s that are half the size.

-   Updated `filter` function for a 2x speed increase. Turns out quos(...) is expensive when you need to do it thousands of times, so now passing along the character version of `dots` when we can instead of recomputing it.

-   JK pdate again `filter` is actually 10x faster. Way better to pull out all the data into a single df and then just reassign everything to the proper event

-   Updated `localizations` slot to use more standardized names for imported PGTargetMotion locs. Currently using latitude, longitude, perpDist, perpDistErr, depth, depthErr. These seem like good standard localization columns going forward. Changed the `ppVars()$tarMoCols` to pull in only columns corresponding to these, and added new `ppVars()$locCols` to reflect new standard column names

-   Added `addNote`, `getNotes`, `removeNote` functions for adding and viewing notes. Can be added to study or event, and will print when they are printed.

## PAMpal 0.16.8

-   `addDatabase` not working properly for directory of dbs

## PAMpal 0.16.7

-   Bug fixes for `addGps` when using a dataframe as the source, `thresh` argument was not doing anything

## PAMpal 0.16.6

-   Changing `myGram` function used for power spectrum calculations to normalize by average window *energy* instead of average window *amplitude* because we are almost always using power spectra

## PAMpal 0.16.5

-   Updating `addHydrophoneDepth` to be able to take numeric input for constant depth

## PAMpal 0.16.4

-   Changing added column fro `depth` to `hpDepth` in `addHydrophoneDepth`

-   `writeEventClips` now has a `filter` option to apply low/high/bandpass filters

-   `getSr` and `matchSR` functions were not working properly in some edge cases

-   changed `getDetectorData` to return `db` for more usefulness

-   `getBinaryData` returns UIDs in the order given

## PAMpal 0.16.3

-   Minor adjustment to `mode='db'` to not repeat UIDs within a single event

## PAMpal 0.16.2

-   Missed a bug in new `addRecordings`

-   `addBinaries` now works with vector of folders

## PAMpal 0.16.1

-   Changing `addRecordings` to return NA values when reads are not possible instead of omitting

-   Some fixes to `plotGram` - previously whistle contours could plot at wrong frequnecy and time scale could be thrown off

## PAMpal 0.16.0

-   `addDatabase` can now accept a folder containing databases

-   cepstrum calculations handled better internally

-   Detection Group Localiser had changed in PAMGuard - reworked code to work with new naming scheme and fail less often if some tables are missing

-   `filter` has more regular behavior - no longer returns `NULL` if no detections are present, was previously causing odd behavior in loops

-   Event order is now preserved when doing `mode='time'`

-   Added functions `nDetections`, `nClicks`, `nWhistles`, `nCepstrum`, `nGPL` that get total number of detections in an event or study

-   `getICI` for `type='data'` returns one big dataframe instead of lists of lists of dataframes

-   Added `showBreaks` option to `calculateAverageSpectra` to toggle showing line breaks for multiple events in concatenated spectrogram

-   Added `bindStudies` function to combine multiple `AcousticStudy` objects

## PAMpal 0.15.3

-   Trying to stabilize processing against minor errors - should no longer crash entire processing if one piece of binary data fails, will try again instead

-   Trying to fix some crashes with missing Sound Acquisition data

## PAMpal 0.15.2

-   Bug with `mode='time'` where sample rate could be read in as a character, causing crashes

-   `calculateAverageSpectra` would have problems when doing multiple events that had detections with repeated UIDs

-   `matchEnvData` now has a `depth` parameter to specify depth or range of depths to match

-   Store and print total processing time for `processPgDetections`

## PAMpal 0.15.1

-   Added ability to add functions during `PAMpalSettings` call

## PAMpal 0.15.0

-   Added support for GPL Detector data!

-   Added `processDate` to ancillary info of each study

-   Removing some unnecessary package dependencies

## PAMpal 0.14.6

-   Changed included `testClick` from 192kHz to 500kHz to avoid unnecessary errors when trying to add filters \> 96kHz

## PAMpal 0.14.5

-   Some bug fixes for weird SoundTrap click detector stuff

## PAMpal 0.14.4

-   Added `brightness` and `contrast` parameters to `calculateAverageSpectra`. These only affect the concatenated spectra plot. Also added a `ylim` parameter for the mean spectra plot.

## PAMpal 0.14.3

-   Bug fixed in `calculateAverageSpectra`, previously was not working with event names as an argument for `evNum`

## PAMpal 0.14.2

-   `calculateAverageSpectra` was always computing noise spectrum cuz I'm bad at stuff

-   Bug fix for `addRecordings` with certain SoundAcquisition table setups

-   `getBinaryData` was improperly complaining about multiple matches

## PAMpal 0.14.1

-   `checkStudy` complains about too many things, shush already

## PAMpal 0.14.0

-   Normalizing some paths internally

-   `checkStudy` checks for NA values and reports

-   Small bug in `addRecordings` that wasn't assigning `startSample` properly for separate wav files

-   Trying to re-work how we can consistently access SR after processing. Functions should now be able to automatically tell if a decimated sample rate is appropriate (`calculateAverageSpectra` and `writeEventClips` mainly) for clicks in binaries as long as XML settings have been added with `addSettings` before processing.

-   `PAMpalSettings` has new parameter `settings` to add XML settings without an additional call to `addSettings`

-   `filter` works better with detector names

-   Grouping in `mode='recording'` wasn't working well when recordings were processed with "Merge contiguous files" checked, fixed now

-   `writeEventClips` works with decimated data now with `useSample = TRUE`

## PAMpal 0.13.0

-   Added more info to `loadPamguardXML`. Now loads in click sweep classifier type data as well as pre-filter and trigger filter data

-   `standardClickCalcs` no longer fails if upper end of filter \>= to Nyquist

-   `processPgDetections` with `mode='time'` can run successfully if no Sound_Acquisition table is present by providing columns `db` and `sr` the grouping file (`db` not needed if processing a single database)

-   `processPgDetections` with `mode='db'` can run successfully if no Sound_Acquisition table is present, will ask for sample rate during processing

-   Missed some function warnings with new warning system

## PAMpal 0.12.8

-   `calculateAverageSpectra` now draws lines breaking up events in concatenated spectrogram if plotting multiple events at once

-   `plotGram` now uses start and end times based on the start of the event, not the recording

## PAMpal 0.12.7

-   `standardClickCalcs` incorrectly had outputs labelled `centerHz_3dB` and `centerHz_10dB` when units were actually in kilohertz. Names have been changed to `centerkHz_3dB` and `centerkHz_10dB`

-   `calculateAverageSpectra` has a few new features. `plot` can be a vector of lenght two specifying which of the two plots to create, useful when only one is desired. New parameters `mode` and `decimate` added. `mode` allows users to specify `'spec'`(default) or `'ceps'` to calculate either spectrum or cepstrum of signal. `decimate` reduces the samplerate of the signal by an integer factor before calculations, this can be especially useful for cepstrum data.

-   Fixed bug in `processPgDetetions` for detection group localiser data crashing with an error related to `tarMoCols`, and also updated code to automatically check for event labels. Will look for (in order) `Text_Annotation` column, column with `species` in the name, column with `label` in the name, or column with `id` in the name. If none of these exist, defaults to the first non-standard column present in the database.

## PAMpal 0.12.6

-   `filter` has new special case to filter by detector names

-   `processPgDetections` rare bug fix with processing bad event grouping files for `mode='time'` and `mode='recording'`

## PAMpal 0.12.5

-   `getBinaryData` can now select for detection type, other functions updated to use this

-   `filter` works better with species, in particular can do `is.na(species)` now

## PAMpal 0.12.4

-   `addHydrophoneDepth` function added that will try to read depth data from a database

-   `processPgDetections` for `mode='db'` now adds target motion localizations if they are present and stores in the `localizations` slot of events

## PAMpal 0.12.3

-   `addRecordings` stores some sample number and samplerate data

-   `writeEventClips` has new option `useSample` to try and use the startSample values stored in binary files instead of just UTC times

## PAMpal 0.12.2

-   `addRecordings` tries to use to get start sample data from database

-   Fixed bug where `processPgDetections` would crash if no comment field was present

## PAMpal 0.12.1

-   `calculateICI` error fixed when exactly two detections were present for a detector

-   `export_banter` reports names of event level measures it finds

-   Switching from using eventUID to event ID for all labeling and `id()` purposes

## PAMpal 0.12.0

-   Added ability to work with Pamguard's XML settings files with functions `addSettings`, `removeSettings`, and `loadPamguardXML`. Currently this lets `PAMpal` tell if a decimator was in use and adjust the sample rate accordingly, no longer need to manually set the `sr_hz` parameter for `standardClickCalcs` if XML settings have been added to the `pps`.

-   `PAMpalSettings` objects now have a `@settings` slot to accommodate this addition, older objects will need to have this slot manually added or they will give warnings of `Not a validObject(), no slot of name "settings" for this object of class "PAMpalsettings"`

-   Added `updatePamObject` function to deal with updating S4 classes created in older versions of `PAMpal`

-   Added message on load to warn about updating objects

-   Added new warning system that captures warnings during processing and function calls and stores them in the `ancillary` slot of the `AcousticStudy` object. These can be accessed with the new `getWarnings` function

## PAMpal 0.11.2

-   `calculateICI` keeps UID info with the data

## PAMpal 0.11.1

-   `calculateAverageSpectra` has a `sort` option to sort clicks by peak frequency in the concatenated spectrogram

-   `processPgDetections` with `mode='time'` does a better job of parsing sample rates from database

## PAMpal 0.11.0

-   Added `mode='recording'` to `processPgDetections` that will group events by the start and end times of wav files used

-   `processPgDetections` will try to auotmatically determine an appropriate `mode` if none specified

## PAMpal 0.10.6

-   `processPgDetections` with `mode='db'` now reads in the comment column of the database and stores in the `ancillary()` slot of each event

-   `addRecordings` reads in an extra format of dates [0-9]{14}\_[0-9]{3}

-   Channels for clicks are read in separately from `standardClickCalcs`, now these will be saved even if no calc functions are present along with UID and UTC

-   `checkStudy` behaves better if run on a study with no events

## PAMpal 0.10.5

-   `calculateAverageSpectra` doesn't fail on high SNR values anymore

-   `plotWaveform` and friends have a length argument and will now attempt to get sample rate values from the `AcousticStudy` object if possible

-   `getBinaryData` tries to get sample rate better

-   `plotGram` added in preliminary form. Plots spectrograms and cepstrograms of events, can also overlay WMD detections onto plot.

## PAMpal 0.10.4

-   `calculateAverageSpectra` - more improvements. Noise shouldnt fail anymore, added SNR threshold argument, returns all noise data and UIDs. Output names changed from `all` and `average` to `allSpec` and `avgSpec`.

## PAMpal 0.10.3

-   `calculateAverageSpectra` has a channel argument, better color scaling for concatenated spectrogram, and a beta version of average noise floor

-   `processPgDetections` will now automatically use `mode='time'` if a `grouping` argument is provided and is either a dataframe or a valid filepath

-   `processPgDetections` now properly updates progress bar for `mode='time'` when some binary files have no data or are not within any events, previously would not get to 100% or would jump

## PAMpal 0.10.2

-   `calculateAverageSpectra` now averages in exponential space rather than in dB space, and plots are created in separate plotting windows. The `norm` option now just scales to a maximum dB value of 0

## PAMpal 0.10.1

-   Multiple functions that tried to match files against a list would fail if filepaths contained "\\" as file separators, these have been changed to use `fixed=TRUE` so that this no longer happens

-   Can now set function parameters for functions being adding during the initial call to `PAMpalSettings` and `addFunction` instead of being forced to go through interactive menus

## PAMpal 0.9.14

-   Cepstrum module can now be named "BurstPulse" or "Burst Pulse" in addition to "Cepstrum" to be recognized by `PAMpal`

-   Function checker was set up poorly and would unintentially write to global env on some failed checks

-   `calculateAverageSpectra` was not passing along manual `sr` input properly

-   `processPgDetections` with `mode='time'` handles date converting better, will no longer fail with Excel rounding 12:10:00 to 12:10. Also removed some unnecessary interactive steps, will now stop and ask user to provide instead of asking user to provide during call

## PAMpal 0.9.13

-   Fixed bug in `calculateICI` if large values were present rounding would lead to bad answers. Now filters results by z-score \< 2 before taking ICI to get rid of outlier large values

-   Fixed bug in `writeEventClips` that was undoing buffer in the wrong direction if the buffered time was not present in any of your wav files, would previously get a false warning of "cant find wav file"

-   Updated how warnings work for `writeEventClips`, now always reported even on failure

-   Adding filter options to `calculateAverageSpectra` similar to `standardClickCalcs`

-   Fixed bug where processing with `mode='time'` was leaving in full binary path names instead of just file basename, messing up other stuff like `getBinaryData`. Fixed and temporary adjustment to `getBinaryData` made

## PAMpal 0.9.12

-   `getBinaryData` was not doing a good job if binaries had the same name. It do better now

-   `processPgDetections` will only store file names that it actually used instead of all possible matching file names, related to other binary fix.

-   Added `getICI`, `getClickData`, `getWhistleData`, and `getCepstrumData` as easier accessors

-   Updating webpage with `calculateICI` and `getDetectorData` guides

## PAMpal 0.9.11

-   `writeEventClips` now attaches clip times with \_YYYYMMDDHHMMSS_mmm to the end of file names

## PAMpal 0.9.10

-   `writeWignerData` warns about bad directories instead of crashing

-   `filter` fixed for bug in species and database filtering if these keywords were present on the right hand side of the call

## PAMpal 0.9.9

-   Added `writeWignerData` function

## PAMpal 0.9.8

-   newUIDs not getting used properly for some users because they were integer64

-   fixed bug with `writeEventClips` not reporting all file names if multiple databases present

## PAMpal 0.9.7

-   `addRecordings` prints some failure warnings

-   `filter` for databases also removes those databasess from the `files(x)$db`

## PAMpal 0.9.6

-   `writeEventClips` updated to use channel argument

-   R and package version now stored in studies when created

-   `standardClickCalcs` now uses "channelMap" from binaries to assign correct channel if it is present instead of defaulting to Channel 1 & 2

-   `export_banter` warning messagaes cleaned up

-   Removed extra "Peak 3dB" and "Peak 10dB" values from `standardClickCalcs`

-   `filter` now works with a special case for `database` and also works for any environmental variables added

-   Updating some docs for GitHub pages website

-   Adding `calculateAverageSpectra` back

## PAMpal 0.9.5

-   Modified `processPgDetections` for `mode='db'` to work with a possible "newUID" column from PAMmisc's updateUID function. Also does a better job of matching if there are multiple possible binary files.

-   `addGps` now asks if you want to add gps data to your database if it has none

-   Better warning messages for processing data

-   `writeEventClips` modified to work with events or individual detections

## PAMpal 0.9.4

-   Big update to documentation and testing for CRAN prep

## PAMpal 0.9.3

-   `addRecordings` function added to deal with keeping track of wav recordings in a study

-   Name changed to PAMpal because `pamr` was taken on CRAN

## PAMr 0.9.2

-   Lots of prep for CRAN submission

-   Changes to `standardClickCalcs` so it doesnt crash if peak frequency is 0, happens if lower filter bound is 0

-   Added `updateFiles` function to update file locations of a study object if you have changed computers or moved folders around

-   Added `checkStudy` function to do some sanity checking for possible issues after data is processed. Currently only checks for peak freq of 0 which would mean you probably want a different filter value or add a calibration function

-   Added test datasets to inst/extdata to go along with added testthat support

## PAMr 0.9.1

-   NEW FUN FUNCTION `calculateAverageSpectra`! Calcultes and plots average spectra from an event! Woohoo!

-   Fixed issue with `export_banter` not finding ancillary measures properly

-   `matchEnvData` compatible with update in `PAMmisc` package, works with custom functions now

-   `setSpecies` handles things more gracefully without provided method argument

-   Changed `setSpecies` to choose insted of stop when multiple methods provided (default)

## PAMr 0.9.0

-   Added `filter` method that works like dplyr's filter for detections in an object. Has a special case when specifying `species` or `Species` that will filter an `AcousticStudy` object by the species listed in the `$id` spot of each species slot

-   Also added first version of `addAnnotation` function, it asks a lot of questions for now

-   `AcousticStudy` `[` accessor changed to return another `AcousticStudy` after subsetting the events in that study

## PAMr 0.8.3

-   Minor bug fix for `writeEventClips`

## PAMr 0.8.2

-   Fixed bug in `addGps` that could create duplicate rows after adding gps to your data. Would occur if two GPS entries occurred within the same second, now adds milliseconds before matching so that this is avoided.

-   Added support for SoundTrap files for `writeEventClips` with `format='soundtrap'`

## PAMr 0.8.1

-   Added better progress bar for `mode='db'` that is based on binary files not databases

-   Added `writeEventClips` function for creating wav clips of events. Currently needs some adjustments but should work for Pamguard data if your selection of wav files definitely covers all of your events

-   Fixed bug where click detectors were being improperly named if a matched classifier was present in the binary file

## PAMr 0.8.0

-   Finally not just a bug fix! New function `matchEnvData` extends a function of the same name in the "PAMmisc" package (as of PAMmisc v1.4.1). This lets you download environmental data and match it to your `AcousticStudy` and `AcousticEvent` objects

## PAMr 0.7.18

-   Bug in `export_banter` if an event had exactly 1 detection fixed

## PAMr 0.7.17

-   Fixed bug in `export_banter` for list of `AcousticStudy` objects not propagated event names properly

## PAMr 0.7.16

-   `export_banter` can take list of `AcousticStudy` objects now, and fixed a bug where NA values were being introduced

## PAMr 0.7.15

-   Skip over binaries with no UIDs in `mode='db'`

-   `mode='db'` was not properly passing along an ID if you supply it

-   `export_banter` now lists which measures were `NA` in the `$na` dataframe output

## PAMr 0.7.14

-   First pass fixing `egClicks`

## PAMr 0.7.13

-   Error message in `processPgDetections` for `mode='db'` now reports db and binary file failed on better

-   Changed folder choosing funciton for `addBinaries`, previous version was windows-only and would not properly start in current working directory

-   `addDatabase` now prints added databases for reference

## PAMr 0.7.12

-   Fixed calibration loading file if first column is just row names

## PAMr 0.7.11

-   `setSpecies` now has `method = 'reassign'` that can reassign species id to new ones provided in a dataframe

## PAMr 0.7.10

-   `processPgDetctions` bug fixed when encountering a binary with 0 detections

## PAMr 0.7.9

-   `processPgDetections` shouldnt crash on empty binary file that also didnt finish writing the header and footer for `mode='time'`

## PAMr 0.7.8

-   `processPgDetections` handles binary files that never finished writing without a mysterious crash now. Should only affect `mode='time'`.

## PAMr 0.7.6 and 0.7.7 because oops

-   Just adde a debugger for JKs problem

## PAMr 0.7.5

-   `standardClickCalcs` changed to have a `filterfrom_khz` and `filterto_khz` option so that a bandpass filter can be applied. If only a highpass filter is desired simply leave `filterto_khz` as the default `NULL`

## PAMr 0.7.4

-   `addCalibration` reworked slightly again, created calibration function now takes frequency in Hz as input and outputs dB correction to add instead of taking in result of seewave::spec and outputting corrected spec. This makes it more flexible and means we aren't calculating spectrum twice for no reason

-   `standardClickCalcs` adjusted to use new calibration, and temporarily has a new parameter calculated - `dBPP` - which is the peak to peak value as traditionally defined in acoustics. This is mostly for testing comparability to other Matlab based code of new calibration stuff. Also the spectrum is never rescaled to a max of 0 like it was previously, hopefully should represent accurate SPL levels with new calibration code

## PAMr 0.7.3

-   `addCalibration` now asks for units of your calibration, and should properly add or subtract this value from calculated spectrum depending on units selected

## PAMr 0.7.2

-   `export_banter` now has an option to split your data into training and test sets by specifying a proportion for the `training` argument, and now puts out a pretty table showing a summary of your events and species

## PAMr 0.7.1

-   `processPgDetections` has some minor updates to improve clarity when binary files or databases aren't found using `mode='time'`

## PAMr 0.7.0

### Functions Renamed

-   `getPgDetections` has now been re-named to `processPgDetections` to more accurately reflect what it is doing. Functions starting with `get___` will be used just for accessing data

-   `showWaveform` and related functions have been renamed to `plotWaveform`

-   General naming consistency overhaul - `sr` and `db` should now be used in place of `sampleRate` and `database` wherever these were previously used, and this should be the naming convention going forward. This means included functions will need to be re-added (they previously relied on it being named `sampleRate`), and any custom functions should expect data to contain the item `sr` and not `sampleRate`

### Major Changes

-   `AcousticStudy` class has been added. This will now be the class of object returned by `processPgDetections`, it stores your list of `AcousticEvent` objects with other important data.

-   Any data processed in previous versions will need to be re-run, and all functions in older `PAMrSettings` files will need to be removed and re-added

-   New Function: `calculateICI` calculates the inter-click interval of click data, and stores it in the `ancillary` slot of an `AcousticEvent`.

-   New Function: `plotDataExplorer` creates an interactive plot to allow users to explore the data in an `AcousticStudy`, `AcousticEvent`, or `data.frame`. Plot can be colored or facetted by any columns that are characters or factors, and any numeric columns can be plotted with `ggplot`'s `geom_density`

-   New Function: `getDetectorData` gathers all detector data into single data frames for each detector type (`'click'`, `'whistle'`, and `'cepstrum'`) for easier manipulation

### Minor Changes

-   Detector dataframes in `AcousticEvents` now have a `'calltype'` associated with them, currently one of '`whistle'`, `'click'`, or `'cepstrum'`. This allows for future fun stuff to happen

-   `getBinaryData` will now attempt to get the appropriate sample rate for each data point, either from the settings or matching by time using the database file if more than one sample rate was in your data.

-   Many speed improvements - click calculations should take about half the time, and `processPgDetections` with `mode='time'` will now skip over binaries that are outside of the time range of specified events

-   `getPgDetections` no longer has an option `mode='all'`, and it is no longer necessary to supply a sample rate for `mode='time'`. It will now attempt to match events to a correct database used on time stamps, the match the appropriate sample rate the same way that `mode='db'` does. If a manual sample rate needs to be applied it should now be entered into the `grouping` table

-   Functions `addBinaries` and `addDatabase` can now add files from another `PAMrSettings` object, and will report how many files have been added when finishing

-   `export_banter` will now attempt to find event level measures in the `'measures'` item in the `ancillary` slot of each `AcousticEvent`.

-   `export_banter` no longer has a `reportNA` option, any `NA` values that are removed are now always saved to a separate `'na'` item in the list output. The function also reports how many detections were processed.

-   `addGps` now also stores all of the gps data loaded into the `gps` slot of the `AcousticStudy` instead of only adding coordinates to detections

## PAMr 0.6.8

-   `getPgDetections` will name events with database appended for `method = 'db'` instead of just event ID number to ensure uniqueness across multiple databases

## PAMr 0.6.7

-   I don't remember what happened here. I bet it was important

## PAMr 0.6.6

-   `setSpecies` can now use a dataframe for `method = 'manual'`, and has a top secret option for SR

## PAMr 0.6.5

-   Updated `export_banter` with options to exclude certain species and to export data without species codes to use for prediction only instead of training a banter model

## PAMr 0.6.4

-   Better error tracking when functions cause `getPgDetections` to crash

-   Now reads in angles and angleErrors from click data

-   `export_banter` allows you to specify certain columns to not export

-   `'time'` mode for `getPgDetections` will now report a sample event time so you can see if times are being converted properly from your csv before proceeding with calculations

## PAMr 0.6.3

-   Fixed an issue where `getPgDetections` would not work if both Detection Group Localizer and Offline Click events were present in a database.

-   Renamed `eventType` and `Text_Annotation` columns from event databases to `eventLabel` within detection dataframes so there is consistency between the two

-   Removed some unnecessary columns from detection dataframes, including `detectorName`, `sampleRate`, `Id`, `parentUID`, and `comment`

-   Fixed a bug in how `getBinaryData` was checking for multiple matches on the same UID

## PAMr 0.6.2

-   Fixed an issue with repeated entries in click detections for modes other than `'db'`

## PAMr 0.6.1

-   Sometimes whistles would not get proper decimated sample rate, fixed

## PAMr 0.6.0

-   Added an `id` slot to `AcousticEvent` objects. Note that this will cause existing `AcousticEvent` objects to behave poorly until they have their `id` slot created / updated using the new `setIdSlot` function.

-   Added `setIdSlot` function to update older `AcousticEvent` objects to the new format

## PAMr 0.5.9

-   Fixed bug in SR / FFT parameter calculation in whistles if there was a gap in the whistle contour

## PAMr 0.5.8

-   Changed event naming for `mode='time'` in `getPgDetections`. Will now only append numbers if event is not unique, and will also insert an underscore before the number

-   `export_banter` now properly checks for cases when there are no detectors in an event instead of crashing confusingly

## PAMr 0.5.7

-   Rocca whistle calcs minor change - boundary settings for sweeps set to "flat", should reduce number of inflection points

## PAMr 0.5.6

-   minor change to `export_banter` no longer using list names to index and create unique event names

## PAMr 0.5.5

-   `export_banter` now removes NA rows, and has reportNA option to see which ones are NA

-   Detector names have spaces replaced with \_ to avoid weird issues later

## PAMr 0.5.4

-   Dealing with zero detection events better for `export_banter`

-   Temporary fix for click calculations with lots of zeroes in wave form - no more NA

## PAMr 0.5.3

-   `standardClickCalcs` now supports manual input of sample rate. Default argument is `'auto'`, which will read from the database. User can supply a numeric value for sample rate in hertz instead.

## PAMr 0.5.2

-   Added a check in `addDatabase` to see if all files are actually .sqlite3 databases

-   Added a `tryCatch` in `getPgDetections` for mode `db` so that it shouldn't stop completely when encountering an error and lose all previously analysed DBs

## PAMr 0.5.1

-   Minor bug fix when using `seewave::spec`, it can produce NA values for the frequency if the input wave is long. Adjusted parts `standardClickCalcs` and `addCalibration` to work around this.

## PAMr 0.5.0

-   `getPgDetections` changed to work by specifying a `mode` as an argument instead of calling separate functions. Can now create events using a csv or dataframe with start and end times specified.

-   Changed `export_banter` to export a list with named item `detectors` instead of `detections`, no other functional change.

## PAMr 0.4.0

-   Added `setSpecies` functions for assigning species classifications to AcousticEvents.

-   Changed `getDbData` to default to looking for both OfflineEvents tables and DetectionGroupLocaliser tables, will still only load one if a specific type is provided for `grouping`.

## PAMr 0.3.0

-   Added `addCalibration`, `applyCalibration`, and `findCalibration` functions, as well as a `plot` method that will show the calibration function used. See the *Calibration* section above for more details.

-   Fixed a bug in `removeFunction` that would cause the incorrect number of functions to show in certain cases, and that would cause all functions to be removed when only one was selected.

-   `standardClickCalcs` has been adjusted to work with the new calibration methods. See *Calibration* section above for more details.

## PAMr 0.2.2

-   Added `addGps` function for adding matching GPS data to your detections. This allows you to supply a dataframe of Lat/Long locations with timestamps to match to your detections.

-   Added `showWaveform`, `showSpectrogram`, and `showWigner` functions that allow you to easily plot the waveform, spectrogram, or wigner plot of a detection in an `AcousticEvent` object by selecting the UID(s) you want to investigate further. `getBinaryData` is also added as a helper function for these, lets you easily get the binary file data for a single detection.

-   `standardClickCalcs` now supports supplying a `Wave` class object as input

## PAMr 0.2.0

-   Rocca (`roccaWhistleCalcs`) and cepstrum (`standardCepstrumCalcs`) functions added. These are also added by default to a new PRS.

-   changed `AcousticEvent` class slot name from `specClass` to `species`
