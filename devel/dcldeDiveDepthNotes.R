# Qs 5/2/2024 meeting

dt is diff between max and pair2/3, why?
dt is hmmmmm maybe dont need

wav height, array error, only used to estimate depth error?

min time is often in a non-productive part of wav,
where acf is just high bc close to signal. Maybe
this should be species related (length of click) or
just related to depth we care about (first couple hundred
meters of depth do not matter for perp dist error bc circles

max time does not incorp the errors currently? and sometimes
our echoes are out of bounds from previous memory?

One of big Qs moving forward is should we keep the folder of
wav files.
- No is bc its annoying and a lot of file sto move around
- Yes is bc they are useful for debugging/troubleshooting/
  making sure you trust what is happening
  and faster to rerun if we don't have to dig out the original
  wav files

grid plot is already in


# Answers

filtering by either duration or number of clicks

adjustable ACF threshold

adding hydrophone depths from different sources - currently
pg table but also spreadsheet (and bft/wav height)

ICI with the click metadata?

click milliseconds get dropped in the tables created

hold on to wav clips for now bc its not a big deal

min clip duration work on later - either based on minimum depth
where our method works

double check end of Annabels most recent R script it has shit added

# 5-15 Progress Notes

Next step is implement Shiny code as i want it - some kind of visual for 
"have I looked at this event already"

Filtering functions probably need a reset flag - we want to be able to keep the
filtering results from previous things but also 

Do more test to make sure all is saving properly

Show depth change if +-1/2 m on array - incorporate with general
error est methods

yaml/something input for species-based params

question - should the wav plots do first N, random N,
or equally distributed N through event?

Should wave height be only event level (not imporant)

Error estimation - what do we actually care about here. 
- Wave height is probably not actually mattering on an average scale. m
  may affect individual measurements, but should go back and forth? 
- Dt err I dont think matters at all bc we pick the betters through manual val
- slant measure err is real, % seems to propagate
- soundspeed is like 3%  min/max from spreadsheet (so 1.5% mean), does this matter?
  is this within measurement error?
- simulating possible depth errs is an option, but unclear at which level we need to
  understand errors. Detection? event? part of event track (since changing) ?

Hm..if depth different across dive then tarmo angles are incorrect in different
ways so single correction doesnt work, potentially recalc?

Fix areas marked hard coded in functions + test difference between decimated/not
for the sperm whales - currently I dont think I'm re-decimating so click detector
was running at 96e3 but clips are at 192e3. Does this need to be corrected at 
wav clip level or dd processing level?

Function to store new DD corrected localization in ASO?

# 5-21 Thoughts to discuss

Test between 192 wav clips and decimating down to 96 seemed to have same acorr results

Should wav plots do random N, sequence N, is first N okay?

Error estimation
- Does wav height actually matter? For an individual estimate, yes, but shouldn't
this go back and forth?
- Does DT error matter? We are selecting the best clicks through Shiny manual
- Does soundspeed matter? 3% diff from max to min, is this within measurement
error of whatever instrument?
- Definite sources are array error and TarMo error

Method should really only be accurate for a fixed(ish) depth, unclear how tarmo
loc works for a shifting depth (since this is shifting apparent surface point)
- Sim results show this on my fake dive profile
- I think only affects absolute depth/range values, but should be able to pick
  out track regardless, should we think about re-calculating?

# 5-22 practice run notes for me

### Note for presenting

I'm in charge of math lol

should we add obligatory cute animal photos to some slides

PAMpal bonus features to talk about - calculating stats, filtering, looking at events, connecting
wav files

Discuss clipLen depends on maximum potential echo time

Add example plots of "good" and "bad" wav clip summary plots and ACF plots

### Notes for todo list

Introduction last slide will have numbers for each section, so make sure RMD works

Break up gps/hydro depth steps to talk about different input sources

Move waveHeight or drop??

minDepth/maxRange params are confusing

Bottom mounted situation uses multiple reflection paths from bottom and everything
so we would have to output just the autocorrelation part and find the pattern
from that. Average and then multipeak?
