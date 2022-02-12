## Creating Events by Start and End Times

Rather than using Pamguard's built in tools, it is also possible to specify
events for your data using start and end times. To process your data this
way, run `processPgDetections` with `mode = 'time'` and provide the event time
data in the `grouping` argument. In this case *all* detections in the binary files
between the start and end times will be included in your data.

### Minimum Requirements

Event times can be provided either as a csv file or as a dataframe that is already
loaded into R. In either case there are three required columns for each event: 
start, end, and id. The start and end columns must contain these times in UTC
since this is the time that Pamguard works with, and the id column will be the
id name for that event, these should be unique. Times can overlap, and you can
specify as many events as you like.

**NOTE:** PAMpal needs the times to be in POSIXct, and will try to convert character
dates appropriately. This can be tricky, especially when reading from csv files 
(Excel has a fun habit of automatically changing your date formatting). You can 
either do this conversion yourself first in R (if supplying a dataframe), or use
the `format` argument to specify how the dates should be converted. If PAMpal tries
to convert the date using `format`, it will print out the start time of your first
event and ask you to confirm that it has been converted correctly. 
See `?strptime` for more details on how `format` should be specified.

**NOTE:** The id column should *not* be used for species id, it should be a unique
id to label each event. See below for species identification.

### Optional Columns

There are also three optional columns you can provide: `species`, `db`, and `sr`. The
`species` column is simply used to assign a species id to a specific event, if it
is known (this is entirely optional).

The `db` and `sr` columns are slightly more complicated, and usually do not need to
be provided manually. PAMpal uses the Pamguard database to match the appropriate 
sample rate to a detection based on the time of the detection. Since you can
add multiple databases to one study, PAMpal needs to sort out which databases
belong to which events.

It will attempt to do this based on the times of the 
events and the times in the Sound_Acquisition table of the databases (this is currently
the only table that has any relevant time information), and if it isn't able to sort it out
you will be prompted to select which event belongs to which database (this will
happen if you have two databases with overlapping time coverage, for example). If no
database matches an event (rare), then you will also be asked to supply the
sample rate (in hertz) to use for all detections in that event. 

PAMpal will save all of these database and samplerate selections, and the full
`grouping` dataframe used is stored in your resulting `AcousticStudy` object.
You can take a look, with the following code:

```r
myStudy <- processPgDetections(myPps, mode='time', grouping=myGrouping)
ancillary(myStudy)$grouping
```
And if you need to run this same data again in the future you can provide this
grouping file to avoid having to select databases again:

```r
newStudy <- processPgDetections(myPps, mode='time', grouping=ancillary(myStudy)$grouping)
```

### All Other Columns

Any other columns in your event grouping file will not be used for anything by PAMpal,
but they will be stored with the `AcousticStudy` object in case there is other
useful information there you wish to keep. The entire row of the grouping file for 
each event is stored in the `ancillary` slot of each `AcousticEvent`, similar to how
the full grouping dataframe is stored in the `AcousticStudy`.

```r
# This is just the relevant row for the first event
oneEventGrouping <- ancillary(myStudy[[1]])$grouping
```

One example of why you might want to do this is with the dataframe created by the
`readSpecAnno` function from the PAMmisc package. This function reads in the 
Spectrogram Annotation tables from a database, and these tables contain
frequency bounds `fmin` and `fmax` for your annotated boxes. When using this table
to `processPgDetections` with `mode='time'`, the frequency limits are not used to
filter out detections outside of that frequency range, but you can access them
later and use PAMpal's `filter` function to remove detections outside those
frequency ranges if desired.

```r
# This contains columns id, start, and end by default
sa <- readSpecAnno(db)
myStudy <- processPgDetections(myPps, mode='time', grouping=sa)
# Now we can use the stored grouping information to filter
filtStudy <- myStudy
for(i in seq_along(events(filtStudy))) {
    # Read out the grouping info for each event, which has fmin and fmax
    thisGrouping <- ancillary(filtStudy[[i]])$grouping
    # Filter using whistle detection parameters freqBeg and freqEnd (Hz)
    # and click detections using peak (kHz)
    filtStudy[[i]] <- filter(filtStudy[[i]],
                            freqBeg > thisGrouping$fmin,
                            freqBeg < thisGrouping$fmax,
                            freqEnd > thisGrouping$fmin,
                            freqBeg < thisGrouping$fmax,
                            peak > thisGrouping$fmin / 1e3,
                            peak < thisGrouping$fmax / 1e3,
                            UTC + duration < thisGrouping$end)
}
```
