## Next Processing Steps

## Assigning Species IDs

## Adding GPS Data

## Gathering Detection Data in a Dataframe

The `AcousticStudy` and `AcousticEvent` classes that `PAMpal` creates can be awkward
to work with if you need to do something that doesn't have a built-in function. 
In order to get your data into a format that is easier to work with, `PAMpal` has functions
that will gather your data into dataframes. 

The function `getDetectorData` takes as input either an entire `AcousticStudy` or a single
`AcousticEvent`, and gathers all the detector data contained within into separate dataframes
for clicks, whistles, and cepstrum data. It returns a list of dataframes, named by these
detector types. Each dataframe within that list will contain all the parameters calculated
by the processing functions, as well as the event ID, detector name, and the species ID (
species will be NA if it has not been set using `setSpecies`). In addition to 
`getDetectorData`, there are also three functions that do the exact same thing to get data
for only specific detectors, `getClickData`, `getWhistleData`, and `getCepstrumData`. These
have the exact same functionality, and are just convenient to directly output a dataframe
instead of needing to access it from a list.

```r
# Get data for your entire study
allDets <- getDetectorData(myStudy)
# this will contain $click, $whistle, and $cepstrum (if those are present in your data)
names(allDets)
# To get the actual dataframe, get it out of the list first
str(allDets$click)
str(allDets$whistle)
str(allDets$cepstrum)

# The functions for accessing just one type of detector directly
justClicks <- getClickData(myStudy)
str(justClicks)
identical(justClicks, allDets$click)
justWhistles <- getWhistleData(myStudy)
justCepstrums <- getCepstrumData(myStudy)

# These also works for a single event
oneDets <- getDetectorData(myStudy$`Example.OE1`)
oneDets <- getDetectorData(myStudy[[1]])
oneClick <- getClickData(myStudy[[1]])
```

## Calculating Inter-Click Interval (ICI)

`PAMpal` has a built in function for calculating the inter-click interval of your
data since this is a common step for a lot of analyses. The calculation is done by
simply sorting all the detections by time, and then for each detection taking the 
difference in seconds between it and the previous detection. Then from these values
the most common number is selected as the ICI value (it is slightly more complicated 
than this because the individual time differences are likely to be all slightly different values,
but this is the idea).

The function is called `calculateICI`, and is very straight forward.
There is only really one option to set, this controls what number to use
as the time of each detection. `time='time'` simply uses UTC time in Pamguard

```r
myStudy <- calculateICI(myStudy, time='time')
```

`time='peakTime'` adjusts this slightly by using the time of the highest value
in the waveform clip. So if the peak value of the waveform for a given
detection is 500 samples into a clip, then 'peakTime' will use the UTC
time plus 500 / SampleRate as the time of that click

```r
myStudy <- calculateICI(myStudy, time='peakTime')
```

This calculation is done for every event, and is done separately for each click detector
in the event (note that `PAMpal` splits click detections up by click classification number,
so you have Click_Detector_0, Click_Detector_1, etc.), and also calculated combining all the
detectors in an event. These data are stored within the `ancillary` slot of each event which can
be accessed using the `ancillary()` function, but the easiest way to get the data back out is using
the `getICI` function. This has one parameter, selecting the type of data you want to get.
`type='value'` will return the single ICI value calculated for each
detector as a list named by detector name.

```r

iciValues <- getICI(myStudy, type='value')
```

This returns a list of results for every single event, so to see the
ICI values for your first event:

```r
iciValues[[1]]
```

`type='data'` will return all the individual time differences used to calculate
the number returned by 'value', this can be useful for making plots or 
if you have your own way of doing things

```r
iciData <- getICI(myStudy, type='data')
```
These are similarly returned as a list for each event, and the result
is a list of dataframes for each click detector that just contain
the name of the detector and the time difference values used

```r
str(iciData[[1]])
```

Looking at the actual numbers for the ICI data that combines all the detectors,
the first value will always be 0 since there is no time between the first detection
and the previous detection. It can also appear that the ICI values are repeated,
especially for `time = 'time'`, but this is because the time difference calculations are
done separately for each channel. In fact for `time = 'time'` the values across channels
will be exactly the same since Pamguard does not store a separate detection time for
each channel, but the ICI values should be close but slightly different for `time ='peakTime'`

```r
iciData[[1]]$All
```

## Calculate and Plot Average Click Spectra

## Adding Environmental Data

## Filtering Data

AcousticStudy objects can be filtered with syntax similar to the `dplyr` package
using the function `filter`. There are currently four ways data can be filtered:
by database, species, environmental data values, and function parameter values.

Filtering by database leaves only events with databases in `files(event)$db` matching
the criteria provided, but database names must exactly match the full file path to the
database. Criteria are specified using either `database` or `Database`, the best way to provide
the full names is typically by indexing from `files(myStudy)$db`.

```r
# This won't work because it needs the full file name and path
oneDb <- filter(myStudy, database == 'FirstDb.sqlite3')
# This is best, events from only first database
oneDb <- filter(myStudy, database == files(myStudy)$db[1])
# Events from all databases other than first
notFirstDb <- filter(myStudy, database != files(myStudy)$db[1])
# To specify multiple, use %in%
# Events with first two dbs only
twoDb <- filter(myStudy, database %in% files(myStudy)$db[1:2])
# Events from all dbs other than first two
notTwoDb <- filter(myStudy, !(database %in% files(myStudy)$db[1:2]))
```

Filtering by species leaves only events matching the species criteria provided, and thus
should only be done after species are assigned using `setSpecies`. Criteria are specified 
using either `species` or `Species`.

```r
# Only 'ZC' events
zcOnly <- filter(myStudy, species == 'ZC')
# Not ZC events
notZc <- filter(myStudy, species != 'ZC')
# To specify multiple species, use %in%
ZCGG <- filter(myStudy, species %in% c('ZC', 'GG'))
notZCGG <- filter(myStudy, !(species %in% c('ZC', 'GG')))
```

Filtering by environmental data leaves only events matching the criteria provided, and
the names of the criteria must exactly match the names of variables listed in
`ancillary(myStudy[[1]])$environmental`, so it is usually best to double check these names
before filtering.

```r
# This probably wont work because name is not exact match
shallowOnly <- filter(myStudy, sea_floor_depth > -500)
# Environmental parameters usually have mean or median added to the name, this works
shallowOnly <- filter(myStudy, sea_floor_depth_mean > -500)
```

Filtering by function parameters works slightly differently than the above three methods.
The interface is the same, but instead of removing entire events it will remove all detections
that do not match the criteria supplied. Parameter names must exactly match the names of parameters
calculated by some of the processing functions. If a name provided does not match the name of a parameter
in a detector, then all of that data will remain unfiltered. For example, `trough` is a value measured by
`standardClickCalcs`, so filtering by `trough > 10` will affect all click detections in your data, but will
leave all Whistle or Cepstrum detectors untouched. Any events that are left with 0 detections after filtering
are removed.

```r
less10Peak <- filter(myStudy, peak < 10)
peak10to20 <- filter(myStudy, peak > 10, peak < 20)
```

Multiple types of filters can also be combined into a single filter statement

```r
filterStudy <- filter(myStudy,
                        database != files(myStudy)$db[3],
                        species == 'OO',
                        sea_floor_depth_mean < -1000,
                        peak > 15)
```

If you want to filter out specific events by event id / event name, that is actually easiest
to accomplish without using the filter function at all, but rather just using `[]` to subset
your data. You can provide either indexes or full event names, event names are usually easiest
to provide by indexing into `names(events(myStudy))`

```r
firstOnly <- myStudy[1]
someOdds <- myStudy[c(1,3,5,7,9)]
byName <- myStudy[names(events(myStudy))[c(1,3,5)]]
```

**KNOWN ISSUES WITH FILTERING**

Currently there are two known issues with the `filter` function as implemented. First, if you
supply a long list of options for a single filter, it won't work and will likely give you an
error. As a workaround, the function works fine if you first assign these to a variable, then
filter.

```r
# This probably won't work
badFilt <- filter(myStudy, species %in% c("SPECIES1", "SPECIES2", "SPECIES3", "SPECIES4", "SPECIES5",
                             "SPECIES6", "SPECIES7", "SPECIES8", "SPECIES9", "SPECIES10"))
# This should be fine
mySpecies <- c("SPECIES1", "SPECIES2", "SPECIES3", "SPECIES4", "SPECIES5",
                "SPECIES6", "SPECIES7", "SPECIES8", "SPECIES9", "SPECIES10")
goodFilt <- filter(myStudy, species %in% mySpecies)
```

Second, the `filter` function currently does not behave well if you try to use it inside
other functions. Unfortunately there is not currently a work around for this, but I will
be looking into improving the function in the filter so that these issues do not occur.

```r
# This probably won't work
myFilter <- function(x, sp) {
    filter(x, species %in% sp)
}
# running this will give an error about "object 'sp' not found"
myFilter(myStudy, 'SPECIES1')
```

## Accessing Binary Data
