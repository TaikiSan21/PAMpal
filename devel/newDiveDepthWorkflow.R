library(PAMpal)
data <- readRDS('../Data/3D2PAMr_test_files/AMStudy.rds')
data <- addGps(data)
data <- addRecordings(data, '../Data/3D2PAMr_test_files/Recordings/')
data <- addHydrophoneDepth(data)
data <- setSpecies(data, method='manual', value=c('cuviers', 'cuviers', 'sowerbys', 'cuviers'))

spFilt <- tribble(~species, ~from, ~to,
                  'cuviers', 18, 50,
                  'trues', 30, 70e3,
                  'gervais', 25, 70,
                  'mmme', 25, 70,
                  'sowerbys', 50, 80,
                  'blainvilles', 20, 50)

dd <- export_diveDepthNEFSC(data, filter= spFilt, wavFolder = './TestWav/', waveHeight=.2, overwrite = TRUE, clipLength = .03)


wavs <- list.files('./TestWav/', full.names=TRUE)

library(tuneR)

oneWav <- readWave(tail(wavs,1))
nextWav <- readWave(tail(wavs, 1))
par(mfrow=c(2,1))
plot(oneWav)
plot(nextWav)
par(mfrow=c(1,1))
seewave::spec(oneWav)
seewave::spec(nextWav)

