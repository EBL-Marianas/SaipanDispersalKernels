# Note: this analysis requires JAGS to be installed.
# http://mcmc-jags.sourceforge.net/

###
### Functions
###

# a negative log-likelihood function for the rayleigh distribution
fit_rayleigh <- function(s, data) {
  -sum(log(drayleigh(data, scale = s)))
}

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


###
### Packages
###

# load (and install, if necessary) required packages
packages <- c("rjags", "bayesmeta")

ipak(packages)


###
### Working directories
###

setwd(top.wd) # ignore if error
top.wd <- getwd()
data.wd <- paste(top.wd, "data", sep = "/")
analysis.wd <- paste(top.wd, "analysis", sep = "/")


###
### User defined settings
###

# MCMC setings
# for testing:
# adapt   <- 1000
# samples <- 300
# thin    <- 10
# chains  <- 3

# for fitting
chains  <- 3
thin    <- 100
adapt   <- 10000
samples <- 1000 # number of samples stored


###
### Data wrangling
###

# load beacon test data
setwd(data.wd)
telem_test <- read.csv("final_beacon.csv")

# estimate scale parameter for Rayleigh distribution using beacon measurements
# with known locations
scale <- optim(par    = 0,
               fn     = fit_rayleigh,
               data   = telem_test$Diff.Z,
               method = "Brent",
               lower  = 0,
               upper  = 10E6)$par

# load all of the telemetry data, store as dataframe "D"
D <- read.csv("move_clean_NAs.csv")

# re-format columns in D
D$Date <- as.Date(D$Date,"%Y-%m-%d")  # change date to Dates
D$Time <- gsub(".*[ ]", "", D$Time)   # change times to hh:mm char
D$SID  <- tolower(D$SID)              # change to lowercase
D$SID  <- factor(D$SID, levels = c("mist", "mafd", "gowe", "brwe", "wtgd"))
D$sp   <- as.integer(D$SID)           # species as integer
D$iid  <- as.numeric(D$Name)          # unique numeric individual identifier

# re-format the date and time so that they're in "Y-m-d H:M" format
D$Date <- as.POSIXct(paste(D$Date, D$Time), format="%Y-%m-%d %H:%M")

# calculate net displacements so that they're always with respect to an
# individual's location at time t = 0
inds <- unique(D$ID2)                 # all unique individual/burst combos
D$net_dist <- NA                      # initialize a net_distance column
for (i in 1:length(inds)) {
  ind <- D[D$ID2 == inds[i],]         # grab individual/burst-level data

  # calculate net displacement and update D$net_dist appropriately
  tmp <- sqrt((ind$X[1] - ind$X)^2 + (ind$Y[1] - ind$Y)^2)
  D[D$ID2 == inds[i],]$net_dist <- tmp
}

# make a subset dataframe that only contains relevant columns
data <- D[,c("Date", "SID","sp", "Name","iid", "burst", "net_dist", "time.mov")]
colnames(data) <- c("obs_date", "SID", "sp", "name", "id", "burst", "disp", "t")


###
### Prep data for JAGS
###

# remove missing rows from data
data <- na.omit(data)

# For getting individual level data
nested.id.sp <- cbind(as.numeric(data$id), as.numeric(data$sp))[!duplicated(cbind(as.numeric(data$id), as.numeric(data$sp))),]

# prep JAGS input
data <- data[order(data$id),]
jags_data <- list(d_o   = as.numeric(data$disp),  # observed net displacements
                  t     = data$t,                 # time since burst start
                  scale = scale,                  # Rayleigh scale parameter
                  N_obs = nrow(data),             # number of observations
                  N_spp = 5,                      # number of species
                  z     = rep(0, nrow(data)),     # for "the zeros trick"
                  id    = as.numeric(data$id),
                  nested.id.sp = nested.id.sp,
                  N_id  = dim(nested.id.sp)[1])

# parameters to monitor
jags_params <- c("p.a", "p.b",
                 "a_spp_m", "a_spp_sd",
                 "b_spp_m", "b_spp_sd")


###
### Define model
###


sink("saturating_net_displacement.txt")
cat("
    model {



    ### likelihood

    for(i in 1:N_obs){

    # calculate expected distance at the current elapsed time (in seconds)
    d_e[i] <- a[id[i]] * t[i] / (1 + a[id[i]] * t[i] / b[id[i]])

    # find the magnitude of the difference between expected and observed values
    # (i.e. the residual)
    r[i] <- abs(d_e[i] - d_o[i])

    # find the likelihood of the residual given that error is expected to be
    # Rayleigh distributed.
    # Here we're using the zeros trick.
    L[i] <- (r[i] / scale^2) * exp(-r[i]^2 / (2 * scale^2)) + 1E-3
    z[i] ~ dpois(-log(L[i]))

    }


    ### priors

    for(i in 1:N_id){

    a[i] ~ dgamma(a_spp_sh[nested.id.sp[i,2]], a_spp_ra[nested.id.sp[i,2]])
    b[i] ~ dgamma(b_spp_sh[nested.id.sp[i,2]], b_spp_ra[nested.id.sp[i,2]])

    }


    for(s in 1:N_spp) {

    # reparameterizing the gamma by mode and sd
    a_spp_sh[s] <- 1 + a_spp_m[s] * a_spp_ra[s]
    a_spp_ra[s] <- (a_spp_m[s] + sqrt(a_spp_m[s]^2 + 4 * a_spp_sd[s]^2)) / (2 * a_spp_sd[s]^2)

    a_spp_m_unabs[s] ~ dt(a_m_overall_mu, a_m_overall_tau, 2) # Folded noncentral t following Gelman 2006
    a_spp_m[s] <- abs(a_spp_m_unabs[s])

    a_spp_sd_unabs[s] ~ dt(a_sd_overall_sdu, a_sd_overall_tau, 2) # Folded noncentral t following Gelman 2006
    a_spp_sd[s] <- abs(a_spp_sd_unabs[s])

    b_spp_sh[s] <- 1 + b_spp_m[s] * b_spp_ra[s]
    b_spp_ra[s] <- (b_spp_m[s] + sqrt(b_spp_m[s]^2 + 4 * b_spp_sd[s]^2)) / (2 * b_spp_sd[s]^2)

    b_spp_m_unabs[s] ~ dt(b_m_overall_mu, b_m_overall_tau, 2) # Folded noncentral t following Gelman 2006
    b_spp_m[s] <- abs(b_spp_m_unabs[s])

    b_spp_sd_unabs[s] ~ dt(b_sd_overall_sd, b_sd_overall_tau, 2) # Folded noncentral t following Gelman 2006
    b_spp_sd[s] <- abs(b_spp_sd_unabs[s])

    }


    a_m_overall_mu ~ dunif(0, 100)
    a_m_overall_tau <- pow(a_m_overall_sig, -2)
    a_m_overall_sig ~ dunif(0, 100)

    a_sd_overall_sdu ~ dunif(0, 100)
    a_sd_overall_tau <- pow(a_sd_overall_sig, -2)
    a_sd_overall_sig ~ dunif(0, 100)

    b_m_overall_mu ~ dunif(0, 2000)
    b_m_overall_tau <- pow(b_m_overall_sig, -2)
    b_m_overall_sig ~ dunif(0, 2000)

    b_sd_overall_sd ~ dunif(0, 2000)
    b_sd_overall_tau <- pow(b_sd_overall_sig, -2)
    b_sd_overall_sig ~ dunif(0, 2000)


    ### sample for dispersal distance predictions

    for(s in 1:N_spp){

    p.a[s] ~ dgamma(a_spp_sh[s], a_spp_ra[s])
    p.b[s] ~ dgamma(b_spp_sh[s], b_spp_ra[s])

    }


    }
    ",fill=TRUE)
sink()


###
### JAGS analysis
###

mod <- jags.model("saturating_net_displacement.txt", data = jags_data,
                  n.chains = chains, n.adapt = adapt)

update(mod, n.iter = adapt)

move.samp <- coda.samples(mod,
                      variable.names = jags_params,
                      n.iter = samples * thin, thin = thin)


###
### Save chains
###

setwd(analysis.wd)

move_bird_spp <- levels(data$SID)

rm(list = ls()[-which(ls() %in% c("move_bird_spp", "move.samp"))])

save.image(file = "movement_fit.RData")
