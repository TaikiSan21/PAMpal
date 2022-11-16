db <- '../Data/BP_Orca_Example/BP_Orca.sqlite3'
bin <- '../Data/BP_Orca_Example/PAMBinary/'
library(PAMpal)
pps <- PAMpalSettings(db, bin)
data <- processPgDetections(pps, mode='recording')
data <- addRecordings(data, '../Data/BP_Orca_Example/')
loud <- filter(data, dBPP > -30)
hm <- getClipData(loud[2], mode='detection', FUN=makeOneSpec)
makeOneSpec <- function(wav, name, ...) {
    s <- PAMpal:::wavToGram(wav, wl=64, hop=.1, axes=T)
    png(filename = file.path('TestPngs', paste0(name, '.png')),
        width=4, height=3, units='in', res=300)
    image(s$mat, x=s$x, y=s$y)
    dev.off()
}
plot(hm[[1]], type='l')
writewav <- writeEventClips(loud[2], outDir='TestWav', mode='detection')
