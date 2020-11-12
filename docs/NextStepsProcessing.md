## Next Processing Steps

## Assigning Species IDs

## Adding GPS Data

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

## Accessing Binary Data

## Gathering Detection Data in a Dataframe
