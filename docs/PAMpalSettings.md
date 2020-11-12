## The PAMpalSettings object

The first step in using PAMpal is to create a PAMpalSettings object. This is an
S4 class created for this package, and will store all settings related to a
particular anaylsis. The goal of this object is to make it easy to share and
replicate results between users. All you need to do is send someone your
PAMpalSettings object and they can see exactly how you analysed your data, or 
even point the PPS at their own data and process it in the exact same way.

A PAMpalSettings object has four slots:

* db - This stores the full path to any SQLite databases that will be analysed
* binaries - This stores the folder names of any binaries, as well as the full
path to all individual Pamguard binary (.pgdf) files within those folders
* functions - This stores a set of functions that will be applied to the data 
when read in from the database/binary files. Stores functions for the 
'ClickDetector', 'WhistlesMoans', and 'Cepstrum' modules separately. When adding a new 
function, users are asked to supply values for any arguments of the function.
These values are saved with the function and cannot be changed without removing
the function entirely. 
* calibration - Stores a calibration function to correct for different hydrophone
characteristics. Optional, currently only affects certain 'ClickDetector' functions

Note that none of these slots should ever be manually edited, use the functions described
below if you want to add or remove anything.

A PAMpalSettings object can be created without supplying any arguments (as described in the
[quick start guide](README.md)), or it can be created by directly supplying the database and binary
file paths, although the user will still be asked to input parameters for the 
`standardClickCalcs` function.

```r
myDb <- './Data/TestDB.sqlite3'
myBinaryFolder <- './Data/Binaries'
myPps <- PAMpalSettings(db = myDb, binaries = myBinaryFolder)
```
#### Adding to Your PPS

After the initial set-up of your PPS, you may want to add to it. There are 
four functions that accomplish this: `addDatabase`, `addBinaries`,
`addFunction`, and `addCalibration`. The first two are simple, and can be called interactively
just like the initial PPS setup or by providing the paths.

```r
myPps <- addDatabase(myPps)
myPps <- addBinaries(myPps)
newDb <- './Data/NewDB.sqlite3'
newBinaries <- './Data/NewBinaries/'
myPps <- addDatabase(myPps, newDb)
myPps <- addBinaries(myPps, newBinaries)
```
##### Adding Functions

Adding a function is slightly more involved. First make sure the function (or
the package the function is in) is already sourced. Then add the function by
name, also specifying the module as either 'ClickDetector', 'WhistlesMoans',
or 'Cepstrum'.
If you do not specify, you will be asked to choose. `addFunction` will also
ask the user to set the value for any parameters that are arguments to the
function you provide, except for parameters named "data" or "calibration".
In the following example, the user would be ask to set a value for "a" and
would be told that the default for "a" is 1.

```r
meanAdd <- function(data, add=1) {
    result <- apply(data$wave, 2, mean) + add
    data.frame(Channel = 1:length(result),
               MeanAdd = result)
}

myPps <- addFunction(myPps, meanAdd, module = 'ClickDetector') 
```

PAMpal does some checking when a new function is added to a PPS to ensure
that the functions have output that is in the proper format. This is to 
hopefully catch any potentials errors in this step rather than at the 
processing step. If we try to add this same function to the "WhistlesMoans"
module, it will trigger a warning from these checks and the function
will not be added. This happens because the `data` for the "WhistlesMoans"
module does not have a $`wave` portion like the "ClickDetector" data, so the
function does not work properly. See [here][custom-functions] for more information
on the requirements for adding functions other than the built-in functions.

<a href="images/FnAddError.png" data-lightbox="fn-add-error" data-title="Added function successfully to ClickDetector but not WhistlesMoans">![](images/FnAddError.png)</a>

##### Adding Calibration

A calibration file can be added that will adjust caluclated dB values, this currently
mostly affects calculations within `standardClickCalcs`. The calibration can be supplied
either as a CSV file or as a dataframe, and will need to have two columns, with the first 
being Frequency in Hz and the second being sensitivity in dB. 

```r
myPps <- addCalibration(myPps, 'Calibration.csv')
```

This will bring up menu selections asking what units the calibration are in (see
`?addCalibration` for more details) and which functions this calibration should apply
to.

<a href="images/Calibration.png" data-lightbox="add-calibration" data-title="Adding calibration to a PPS">![](images/Calibration.png)</a>

### Adding from another PPS

A separate PAMpalSettings object can also be supplied as the source for any of the 
functions that add to your PPS. In this case everything from the corresponding
portion of the PPS will be added to the new PPS object.

```r
myPps <- addDatabase(myPps, otherPps)
myPps <- addBinaries(myPps, otherPps)
myPps <- addFunction(myPps, otherPps)
myPps <- addCalibration(myPps,otherPps)
```

#### Removing Things From Your PPS

There are four functions that remove items from your PPS, and it is recommended
that you use this instead of trying to alter the PPS manually. All of them can
be called interactively and provide menus for the user to select the item to 
remove. See the help pages for more info.

```r
myPps <- removeDatabase(myPps)
myPps <- removeBinaries(myPps)
myPps <- removeFunction(myPps)
myPps <- removeCalibration(myPps)
```

[custom-functions]: CustomFunctions.md
