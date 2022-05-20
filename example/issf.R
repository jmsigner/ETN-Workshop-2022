# Getting started ----

# Load all packages
library(amt)
library(lubridate)
library(sf)
library(tidyverse)

# Set the working directory, it would be better to use an Rstudio project, bur
# since this session is rather short I decided to use `setwd()` here.
setwd("example/")

# Finally we should set a seed. We will generate random data later. 
set.seed(3912)

# Covariates ----

# We will use only one covariate: the distance to the lake shore. This is
# somewhat simplistic and likely does not reflect the ecology of pikes. But it
# serves as a good example to illustrate the method.

shore <- read_rds("shore.rds") 
ggplot(shore) + geom_sf()

# Lets first create a raster with a 5m resolution of the lake
r <- raster::raster(st_cast(shore, "POLYGON"), res = 5)

# Then rasterize the lake and set each pixel in the lake to `NA`. The
# `raster::distance()` function will calculate the distance from `NA` pixels to
# the next non-NA pixel.
shore.rast <- raster::rasterize(shore, r, background = 0, field = NA)
raster::plot(shore.rast)

# Now we can calculate the distances from all NA pixels to all non-NA pixels. 
shore.dist <- raster::distance(shore.rast)

# And mask this to the outline of the lake
shore.dist <- raster::mask(shore.dist, shore)
raster::plot(shore.dist)

# Tracking data ----

# Note: Henrik Baktoft kindly allowed me to use part of his data for this
# illustration. Please do not distribute or use the data set for other purposes,
# unless you have Henrik's permission. The data consists of tracking data of one
# pike in Denmark.

# Next we can read the tracking data of one pike. This is just a data.frame and
# we need to have at least 3 columns: x, y, and a timestamp.
p1 <- read_rds("pike1.rds")
nrow(p1)

# To speed the analysis up a bit, we will only use the first half of the data.
p1 <- p1[1:2.5e4, ]

head(p1)

# iSSF ----
# .. Preparing the data ----
# Now we can start with the actual analysis. 

# We start with selecting the relevant columns from the data.frame (this is not
# strictly necessary)
tr1 <- p1 %>% select(x, y, ts) %>% 
  # And pare the timestamp of the observation
  mutate(ts = ymd_hms(ts)) %>% 
  # And finally create a track. Tracks are the basic building blocks of tracking
  # data in `amt`. The column names for the coordinate and the timestamp are
  # mandatory, all other information optional. Tracks behave like other
  # data.frame, just that they have some additional information attached. With
  # the argument `crs` we provide the coordinate reference system of the data.
  make_track(x, y, ts, crs = 32632)

# Next, lets visually inspect the tracking data
ggplot() + geom_sf(data = shore) +
  geom_point(aes(x_, y_), data = tr1)  

# We an summarize the sampling rate, i.e. how often the position of the fish was
# recorded.
summarize_sampling_rate(tr1)

# Step selection are likely to not very meaningful at the interval of only a few
# minutes and also the data set become huge, we subsample to 30 minutes here,
# with a tolerance of 90 seconds.
tr2 <- track_resample(tr1, rate = minutes(30), tolerance = seconds(90))

ggplot() + geom_sf(data = shore) +
  geom_point(aes(x_, y_), data = tr2)  

# Now we can prepare the data set for the integrated step-selection analysis
tr3 <- tr2 %>% 
  # We only want to consider burts (= sequences of data with equal sampling
  # rates) with at least three relocations.
  filter_min_n_burst() %>% 
  # We then switch from a point to a step representation
  steps_by_burst() %>% 
  # Now we pair each observed step with 100 random steps. Some of these random
  # steps will be outside the lake, we will take of this later.
  random_steps(n_control = 100) %>% 
  # Now we can start with data annotation.
  
  # First we add to each step, if twas during day or the night. Using the
  # argument `where = "both"` calculate the time of day for the start and the
  # end of a step.
  time_of_day(where = "both") %>% 
  # Now we can extract the covariates from a raster layer or stack. Here we only
  # have a raster layer with the distance to shore, but we could easily add
  # additional layer and the would all be extracted at once. Again, with the
  # argument `where = "both"` we can get the value of the covariates at the
  # beginning and end of a step.
  extract_covariates(shore.dist, where = "both") %>% 
  # Finally, it easier for plotting, if we manually dummy-encode the categorical
  # covariate `tod_*_`
  mutate(
    night_end = ifelse(tod_end_ == "night", 1, 0), 
    night_start = ifelse(tod_start_ == "night", 1, 0), 
  )

# Now we have to care for random steps that are outside the lake. We oversampled
# random steps at the beginning. We can remove all steps that have an `NA` value
# for the covariate 
tr3.controll <- tr3 %>% filter(!is.na(layer_end), !case_) %>%
  nest(data = -step_id_) 
tr3.controll %>% mutate(n = map_int(data, nrow)) %>% pull(n) %>% table()
# and then sample for each strata 20 random steps
tr3.controll <- tr3.controll %>% mutate(data = map(data, slice_sample, n = 20)) %>% 
  unnest(cols = data)
# Finally, we have to bind the observed and controll steps together.
tr4 <- bind_rows(
  tr3 %>% filter(case_), 
  tr3.controll
)

# .. Model fitting ----

# We start by fitting three model: 

# The first model only models habitat selection and test if there is a selection
# for distance to shore.
m0 <- tr3 %>% 
  fit_issf(case_ ~ layer_end + strata(step_id_))

# For the second model, we also include `sl_` and interaction between `sl_` and
# the covariate distance to shore.
m1 <- tr3 %>% 
  fit_issf(case_ ~ layer_end + sl_ + sl_:layer_end + strata(step_id_))

# The third models includes a interaction with time of day. 
m2 <- tr3 %>% 
  fit_issf(case_ ~ 
             layer_end + layer_end:night_end + 
             sl_ +  sl_:night_start + sl_:layer_end + 
             strata(step_id_), model = TRUE)

# Lets compare different models with AIC
AIC(m0)
AIC(m1)
AIC(m2)

# .. Interpretation -----
# .. .. Habitat Selection ----

# The function `log_rss()` allows us to do an interpretation in terms of the
# relative selection strength. See also:

# - https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3122
# - https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13441 (especially the supplements)

# The ideas is to compare the position `x2` with one or more positions `x1` and
# see how likely the animal is to choose `x1` over `x2`.

# Lets see if the pike has a preference for for pixels closer to the shore
# during the day
x1 <- data.frame(
  layer_end = 5, # 5m from the shore
  night_end = 0, # 0 = day
  night_start = 0, 
  sl_ = 1
)

x2 <- data.frame(
  layer_end = 100,  # 100 m from the shore
  night_end = 0, # 0 = day
  night_start = 0, 
  sl_ = 1
)

exp(log_rss(m2, x1, x2)$df$log_rss)

# This means, if the pike is given the choice between two location (both during
# the day) , it is ~ 1.5 more likely to choose the one closer (5 m) from the
# shore than the one 100 m away from the shore.

# Lets see if the pike has a preference for day or night
x1 <- data.frame(
  layer_end = 0, 
  night_end = 1, # 1 = night
  night_start = 0, 
  sl_ = 1
)

x2 <- data.frame(
  layer_end = 0, 
  night_end = 0, # 0 = day
  night_start = 0, 
  sl_ = 1
)

exp(log_rss(m2, x1, x2)$df$log_rss)

# This means, if the pike is given the choice between the same location (100 m
# away from the shore), once during the day and once during the night, it ~ 2
# times as likely to find the pike there during the night

# We can also calculate the log-rss for a range of values of `x1`. 
x1 <- expand.grid(
  layer_end = 0:300, 
  night_end = 0:1, 
  sl_ = 1
)

x1$night_start <- x1$night_end

x2 <- data.frame(
  layer_end = 100, 
  night_end = 0,
  night_start = 0, 
  sl_ = 1
)

lrss1 <- log_rss(m2, x1, x2, ci = "se")

lrss1$df %>% 
  ggplot(aes(layer_end_x1, log_rss)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = factor(night_start_x1)), alpha = 0.1) +
  geom_line(aes(col = factor(night_start_x1))) +
  labs(x = "Distance to shore", y = "log RSS", col = "time", fill = "time") +
  theme_minimal()

# .. .. Movement ----
# Finally, since we also modeled the movement of animals, we can now compare if
# the expected displacement changes between day and night and if the distance to
# shore influences the expected displacement.

# You can find all the formulas here: 
# https://conservancy.umn.edu/bitstream/handle/11298/218272/AppC_iSSA_movement.html?sequence=10&isAllowed=y

x1$layer_start <- x1$layer_end

sc <- 1 / ((1 / sl_distr_params(m2)$scale) - coef(m2)["sl_"] + coef(m2)["sl_:night_start"] * x1$night_start +
             coef(m2)["layer_end:sl_"] * x1$layer_start)

movement <- x1 %>% mutate(scale = sc, shape = sl_distr_params(m2)$shape,
displacement = scale * shape)

movement %>% ggplot(aes(layer_start, displacement, col = factor(night_start))) +
  geom_line() +
  labs(x = "Distance to shore [m]", y = "Displacement [m/30 min]", col = "Time") +
  theme_minimal()
