# Note: this analysis requires JAGS to be installed.
# http://mcmc-jags.sourceforge.net/

###
### Functions
###

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

factor.by.frequency <- function(x){
  tb <- table(x)
  factor(x, levels = names(tb[order(tb, decreasing = T)]))
}


###
### Packages
###

# Load libraries
packages <- c("plyr", "dplyr", "tidyr", "rjags")
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
thin    <- 200
adapt   <- 20000
samples <- 1000 # number of samples stored


###
### Load data
###

setwd(data.wd)
gpt.dat <- read.csv("GPT data.csv")

gpt.dat$bird.id <- factor(gsub(" ", "", gpt.dat$bird.id))

# Make it so each row represents a single fecal pile, not a single seed
gpt.dat <- gpt.dat[!duplicated(with(gpt.dat,(paste(date,bird.id,tree.sp,time.fruit.added,pile.id,time.pile)))), ]

# Little bit of data cleaning on n.seeds.pile
gpt.dat$n.seeds.pile <- as.numeric(gsub(">", "", gpt.dat$n.seeds.pile))
table(gpt.dat$n.seeds.pile)
gpt.dat <- subset(gpt.dat, n.seeds.pile != 0)

gpt.dat$n.seeds.minus.one <- gpt.dat$n.seeds.pile - 1


###
### Subset to 15 focal tree species
###

levels(gpt.dat$tree.sp)
gpt.dat <- gpt.dat[-which(gpt.dat$tree.sp %in% c("coccinia", "discocaylx", "elaeocarpus")),]
gpt.dat$tree.sp <- factor(gpt.dat$tree.sp)


###
### Data wrangling
###

head(gpt.dat)
str(gpt.dat)
gpt.dat$bird.sp <- factor.by.frequency(gpt.dat$bird.sp)
gpt.dat$tree.sp <- factor.by.frequency(gpt.dat$tree.sp)


###
### Set up for model
###

bt.combos <- ifelse(table(gpt.dat$bird.sp,
                          gpt.dat$tree.sp) > 0, 1, 0)
combo.mat <- t(ifelse(t(bt.combos) > 0, cumsum(t(bt.combos)), 0))

all.int <- which(bt.combos > 0, arr.ind = T)
all.int.but.first <- which(bt.combos > 0, arr.ind = T)[2:(nrow(which(bt.combos > 0, arr.ind = T))),]
all.int.but.last <- which(bt.combos > 0, arr.ind = T)[1:(nrow(which(bt.combos > 0, arr.ind = T))-1),]
last.int <- which(bt.combos > 0, arr.ind = T)[nrow(which(bt.combos > 0, arr.ind = T)),]
n.all.int <- nrow(all.int)

# Bird species of bird id
bird.sp.of.bird.id <- match(substr(gpt.dat$bird.id,1,4), levels(gpt.dat$bird.sp))
bird.sp.of.bird.id <- match(gsub('[0-9]+', '', levels(gpt.dat$bird.id)), levels(gpt.dat$bird.sp))


for(i in levels(gpt.dat$bird.sp)){
  gpt.dat[which(gpt.dat$bird.sp == i), "nested.bird.id"] <- as.numeric(factor(gpt.dat[which(gpt.dat$bird.sp == i), "bird.id"]))
}

bird.id.nested.list <- lapply(1:5, function(x) which(x == bird.sp.of.bird.id))
bird.id.nested.matrix <- as.matrix(rbind.fill(lapply(bird.id.nested.list, function(x)as.data.frame(t(x)))))
bird.id.nested.length <- sapply(bird.id.nested.list, length)



data <- list(

  log.time = log(gpt.dat$gpt),

  n.samples = nrow(gpt.dat),

  all.int = all.int,
  n.all.int = n.all.int,

  combo.mat = combo.mat,

  bird.sp = gpt.dat$bird.sp,
  tree.sp = gpt.dat$tree.sp,

  n.bird.sp = length(levels(gpt.dat$bird.sp)),
  n.tree.sp = length(levels(gpt.dat$tree.sp)),

  nested.bird.id = as.numeric(gpt.dat$nested.bird.id),
  bird.id.nested.length = bird.id.nested.length

)



###
### Define model
###

sink("gpt.txt")
cat("
    model {

    ### Time: Likelihood
    for(i in 1:n.samples){
    log.time[i] ~ dnorm(log.time.hat[i], time.tau[combo.mat[bird.sp[i], tree.sp[i]]])
    log.time.hat[i] <- ft.overall.est + ft.bird.sp.est[bird.sp[i]] + ft.tree.sp.est[tree.sp[i]] + ft.bird.tree.sp.est[combo.mat[bird.sp[i], tree.sp[i]]] + ft.bird.id[bird.sp[i], nested.bird.id[i]]
    }


    ### Time: Priors

    # Overall effects (intercepts)

    ft.overall.est ~ dnorm(0, 0.01)


    # Bird species effects

    ft.bird.sp.est[1] <- 0
    for(i in 2:n.bird.sp){
    ft.bird.sp.est[i] ~ dnorm(0, ft.bird.sp.tau)
    }

    ft.bird.sp.tau <- pow(ft.bird.sp.sigma, -2)
    ft.bird.sp.sigma ~ dunif(0, 10)


    # Tree species effects

    ft.tree.sp.est[1] <- 0
    for(i in 2:n.tree.sp){
    ft.tree.sp.est[i] ~ dnorm(0, ft.tree.sp.tau)
    }

    ft.tree.sp.tau <- pow(ft.tree.sp.sigma, -2)
    ft.tree.sp.sigma ~ dunif(0, 10)

    # Interaction effects

    ft.bird.tree.sp.est[1] <- 0
    for(i in 2:n.all.int){
    ft.bird.tree.sp.est[i] ~ dnorm(0, ft.bird.tree.sp.tau)
    }

    ft.bird.tree.sp.tau <- pow(ft.bird.tree.sp.sigma, -2)
    ft.bird.tree.sp.sigma ~ dunif(0, 10)

    # Time tau.
    for(i in 1:n.all.int){
    time.tau[i] <- pow(time.sigma[i], -2)
    time.sigma.unabs[i] ~ dt(mu.time.sigma, mu.time.tau, 2) # Folded noncentral t following Gelman 2006
    time.sigma[i] <- abs(time.sigma.unabs[i])
    }
    mu.time.sigma ~ dunif(0, 10)
    mu.time.tau <- pow(mu.time.sigma.overall, -2)
    mu.time.sigma.overall ~ dunif(0,10)

    # Bird individual effects

    for(i in 1:n.bird.sp){
    for(j in 1:bird.id.nested.length[i]){
    ft.bird.id[i,j] ~ dnorm(0, ft.bird.id.tau[i])
    }
    ft.bird.id.tau[i] <- pow(ft.bird.id.sigma[i], -2)
    ft.bird.id.sigma[i] ~ dunif(0, 10)
    }



    ### Predictions
    for(i in 1:n.all.int){
    log.time.pred[all.int[i, 1], all.int[i, 2]] ~ dnorm(ft.overall.est + ft.bird.sp.est[all.int[i, 1]] + ft.tree.sp.est[all.int[i, 2]] + ft.bird.tree.sp.est[combo.mat[all.int[i, 1], all.int[i, 2]]], time.tau[combo.mat[all.int[i, 1], all.int[i, 2]]])
    }


    } # End of model

    ",fill=TRUE)
sink()









###
### Run model
###


# Set to zero constraint inits
inits <- list(
  list(ft.bird.tree.sp.est = c(NA,rnorm(n.all.int-1, sd=0.1))),
  list(ft.bird.tree.sp.est = c(NA,rnorm(n.all.int-1, sd=0.1))),
  list(ft.bird.tree.sp.est = c(NA,rnorm(n.all.int-1, sd=0.1))))



mod <- jags.model("gpt.txt", data = data,
                  inits = inits,
                  n.chains = chains,
                  n.adapt = adapt)

update(mod, n.iter = adapt)

gpt.samp <- coda.samples(mod,
                     variable.names = c("log.time.pred"),
                     n.iter = samples * thin,
                     thin = thin)


###
### Save chains
###

setwd(analysis.wd)

# I think this will make it easier to pull these data into the kernel simulations
gut_tree_spp <- levels(gpt.dat$tree.sp)
gut_bird_spp <- levels(gpt.dat$bird.sp)

rm(list = ls()[-which(ls() %in% c("gut_tree_spp", "gut_bird_spp", "gpt.samp"))])

save.image(file = "gpt_fit.RData")
