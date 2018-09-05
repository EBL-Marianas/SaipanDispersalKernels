###
### Functions
###

rbindchain <- function(x){
  rbind(x[[1]], x[[2]], x[[3]])
}

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

###
### Working directories
###

setwd(top.wd) # ignore if error
top.wd <- getwd()
data.wd <- paste(top.wd, "data", sep = "/")
analysis.wd <- paste(top.wd, "analysis", sep = "/")
figures.wd <- paste(top.wd, "figures", sep = "/")


###
### Packages
###

# load (and install, if necessary) required packages
packages <- c("mapplots", "plotrix")

ipak(packages)


###
### Load posterior samples
###

setwd(analysis.wd)
load("gpt_fit.RData")
load("movement_fit.RData")


###
###
###

times <- rbindchain(gpt.samp)
moves <- rbindchain(move.samp)

dim(times)
dim(moves)


###
### Make these digestable for predictions
###

# times
time.cols <- grep("log.time.pred", colnames(times))
times <- exp(times[ , time.cols]) * 60 # to get time in seconds
head(times)

# a and b values of saturating function
a.cols <- grep("p.a", colnames(moves))
b.cols <- grep("p.b", colnames(moves))
a.vals <- moves[ , a.cols]
b.vals <- moves[ , b.cols]

# now plug in for dists

bird.sp.of.time.cols <- as.numeric(sapply(strsplit(gsub("log.time.pred[",
                                                        "", colnames(times),
                                                        fixed = T), ","), `[`, 1))

dists <- times
dists <- ifelse(dists>0, NA, NA)
for(i in 1:ncol(times)){
  dists[,i] <- a.vals[,bird.sp.of.time.cols[i]] * times[,i] / (1 + a.vals[,bird.sp.of.time.cols[i]] * times[,i] / b.vals[,bird.sp.of.time.cols[i]])
}
colnames(dists) <- gsub("log.time.pred[", "dist.pred[", colnames(dists), fixed = T)


head(dists)

###
### Figures
###

tree.sp.names <- c('Aidia','Pipturus','Premna','Ficus t.',
                   'Ficus p.','Triphasia','Melanolepis',
                   'Psychotria','Eugenia p.',
                   'Planchonella','Carica','Morinda',
                   'Aglaia','Eugenia r.','Meiogyne')
tree.sp.names <- c('Aidia\ncochinchinensis','Pipturus\nargenteus','Premna\nserratifolia','Ficus\ntinctoria',
                   'Ficus\nprolixa','Triphasia\ntrifolia','Melanolepis\nmultiglandulosa',
                   'Psychotria\nmariana','Eugenia\npalumbis',
                   'Planchonella\nobovata','Carica\npapaya','Morinda\ncitrifolia',
                   'Aglaia\nmarianensis','Eugenia\nreinwardtiana','Meiogyne\ncylindrocarpa')
bird.sp.names <- c("starling", "fruit dove", "g. white-eye",
                   "b. white-eye", "ground dove")

xmax <- 500 # max meters to put on plot (hard to see any differences when actual max is used)
ymax <- 0.03 # max y (probability density) to plot. Making this consistent helps see
# differences across bird species







op <- par()

setwd(figures.wd)
pdf("dispersal probability density.pdf", width = 5.25*(1+1/6), height = 14)

par(#mfrow=c(19,6),
  oma = c(0,0,0,0),
  mar = c(0,0,0,0),
  pin = c(0.7,0.7))

bird.alpha <- 20
bird.colors <- c(rgb(27,158,119, bird.alpha, maxColorValue = 255),
                 rgb(217,95,2, bird.alpha, maxColorValue = 255),
                 rgb(117,112,179, bird.alpha, maxColorValue = 255),
                 rgb(231,41,138, bird.alpha, maxColorValue = 255),
                 rgb(102,166,30, bird.alpha, maxColorValue = 255))

bird.colors.dark <- c(rgb(27,158,119, maxColorValue = 255),
                      rgb(217,95,2, maxColorValue = 255),
                      rgb(117,112,179, maxColorValue = 255),
                      rgb(231,41,138, maxColorValue = 255),
                      rgb(102,166,30, maxColorValue = 255))


bird.quadrants <- list()

bird.quadrants[[5]] <- matrix(c(0.8,1.0,
                                0.6,0.8,
                                0.4,0.6,
                                0.2,0.4,
                                0.0,0.2), byrow = T, ncol = 2)

bird.quadrants[[4]] <- matrix(c(0.75,1.0,
                                0.5,0.75,
                                0.25,0.5,
                                0,0.25), byrow = T, ncol = 2)

bird.quadrants[[3]] <- matrix(c(0.666,1.0,
                                0.333,0.666,
                                0,0.333), byrow = T, ncol = 2)

bird.quadrants[[2]] <- matrix(c(0.5,1.0,
                                0,0.5), byrow = T, ncol = 2)

bird.quadrants[[1]] <- matrix(c(0,1.0), byrow = T, ncol = 2)


my.add.pie <- function (z, x = 0, y = 0, labels = names(z), radius = 1, edges = 200,
                        clockwise = TRUE, init.angle = 90, density = NULL, angle = 45,
                        col = NULL, border = NULL, lty = NULL, label.dist = 1.1, lwd = 1,
                        ...)
{
  if (!is.numeric(z) || any(is.na(z) | z < 0))
    stop("'z' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(z))
  else labels <- as.graphicsAnnot(labels)
  z <- c(0, cumsum(z)/sum(z))
  dz <- diff(z)
  nz <- length(dz)
  asp <- get.asp()
  if (is.null(col))
    col <- if (is.null(density))
      c("#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B",
        "#9E67AB", "#CE7058", "#D77FB4")
  else par("fg")
  if (!is.null(col))
    col <- rep_len(col, nz)
  if (!is.null(border))
    border <- rep_len(border, nz)
  if (!is.null(lty))
    lty <- rep_len(lty, nz)
  angle <- rep(angle, nz)
  if (!is.null(density))
    density <- rep_len(density, nz)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = asp * radius * cos(t2p) + x, y = radius * sin(t2p) +
           y)
  }
  for (i in 1L:nz) {
    n <- max(2, floor(edges * dz[i]))
    P <- t2xy(seq.int(z[i], z[i + 1], length.out = n))
    polygon(c(P$x, 0 + x), c(P$y, 0 + y), density = density[i],
            angle = angle[i], border = border[i], col = col[i],
            lty = lty[i], lwd = lwd)
    P <- t2xy(mean(z[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      text(label.dist * (P$x - x) + x, label.dist * (P$y -
                                                       y) + y, labels[i], xpd = TRUE, adj = ifelse(P$x -
                                                                                                     x < 0, 1, 0), ...)
    }
  }
}

laymat <- matrix(rep(1:(7*19), times = rep(c(2,1,1,1,1,1,1), times = 19)),
                 nrow = 19, byrow = T)
layout(mat = laymat)

plot.new()
for(i in 1:length(gut_bird_spp)){
  plot.new()
  text(0,0,
       c("Micronesian\nStarling", "Mariana\nFruit Dove", "Golden\nWhite-eye",
         "Bridled\nWhite-eye", "White-throated\nGround Dove")[i],
       pos = 4,
       srt = 45, cex = 1.1, xpd = T)
}
plot.new()
# text(0,0, "all birds", pos = 4,
#     srt = 45, cex = 1.4, xpd = T)
text(0.5,0.5, "dispersal\nwithin\n250 m\nradius",# pos = 1,
     cex = 1.2, xpd = T)
# text(c(0,1), c(0,0),
#      c("-250","250"), xpd = T)

# Order alphabetically without groups that share common dispersers
new.tree.species.order <- c("ficus t.",
                            "aidia", "pipturus",
                            "papaya", "ficus p.", "melano", "pouteria", "premna", "psychotria",
                            "eugenia p.", "guamia", "triphasia",
                            "aglaia", "eugenia r.", "morinda")


new.tree.species.order.num <- match(new.tree.species.order, gut_tree_spp)

mist.500 <- c()
mafd.500 <- c()

for(j in new.tree.species.order.num){ #1:length(gut_tree_spp)){

  for(i in 1:length(gut_bird_spp)){
    if(i == 1){
      plot.new()
      text(1.75, 0.5, tree.sp.names[j], font = 3, pos = 2, xpd = T, cex = 1.3)
    }

    bt.col <- paste("dist.pred[",i,",",j,"]", sep = "")
    if(bt.col %in% colnames(dists)){

      plot(density(dists[,bt.col], from = 0, to = max(dists[,bt.col])),
           frame = T,
           col = bird.colors.dark[i],
           main = "",
           xaxt = "n",
           yaxt = "n",
           xlab = "",
           ylab = "",
           xlim = c(0, xmax),
           ylim = c(0, ymax),
           xpd = F)
      axis(1, at = c(0, 100, 200, 300, 400, 500), labels = NA,
           tck = 0.025,
           #lwd = 0.5,
           tick = T)


      text(xmax*1.1, ymax * 0.9, paste(round(median(dists[,bt.col])), "m", sep = " "), pos = 2)
      text(xmax*1.1, ymax * 0.75, paste(round(quantile(dists[,bt.col], 0.99)), "m", sep = " "), pos = 2)

      tail.val <- length(which(dists[,bt.col] > xmax)) / length(which(dists[,bt.col] > 0))
      if(i == 1){
        mist.500 <- c(mist.500, tail.val)
      }
      if(i == 2){
        mafd.500 <- c(mafd.500, tail.val)
      }

      if(tail.val > 0){
        points(x = xmax,
               y = tail.val,
               col = bird.colors.dark[i],
               pch = 16, cex = 1, xpd = F)
      } else{
        points(x = xmax,
               y = tail.val,
               col = bird.colors.dark[i],
               pch = 1, cex = 1, xpd = F)
      }
    } else{
      plot.new()
    }
  }

  t.cols <- which(colnames(dists) %in% paste("dist.pred[",1:5,",",j, "]", sep = ""))
  bird.sps <- as.numeric(unlist(lapply(strsplit(gsub("dist.pred[", "", colnames(dists)[t.cols], fixed = T),","), function(x) x[1])))


  n.samples <- dim(dists)[1] * 5/length(bird.sps)

  all.dists <- as.vector(dists[sample(1:dim(dists)[1], size = n.samples, replace = T), t.cols])


  angles <- c()
  for(i in 1:length(bird.sps)){
    angles <- c(angles, 2*pi*runif(n.samples,
                                   min = bird.quadrants[[length(bird.sps)]][i,1],
                                   max = bird.quadrants[[length(bird.sps)]][i,2]))
  }
  angles <- angles-pi

  plot(x = all.dists * cos(angles),
       y = all.dists * sin(angles),
       cex = 0.2,
       pch = 16,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       frame.plot = F,
       ylab = "",
       xlim = c(-250,250),
       ylim = c(-250,250),
       col = ifelse(all.dists > 250, NA, bird.colors[rep(bird.sps, each = n.samples)]))

  if(length(bird.sps) > 1){
    my.add.pie(rep(1,length(bird.sps)), col = NA, radius = 250, labels = NA,
               init.angle = 180,
               lwd = 0.02, 0)
  } else{
    draw.circle(0,0,250, lwd = 0.02)
  }


}

par(op)


par(new = TRUE)
par(fig = c(0.30, 0.95, 0, 0.45),#c(0.5, 1, 0, 0.25),
    mar = c(10,10,10,10),
    pin=c(1,1),
    xpd = T)

plot(NA, #density(dists[,"dist.pred[1,1]"], from = 0, to = max(dists[,"dist.pred[1,1]"])),
     frame = T,
     main = "",
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     xlim = c(0, xmax),
     ylim = c(0,0.03))
text(10, ymax / 2.25,
     pos = 2,
     "Dispersal
     probability", xpd = T)
text(xmax / 2, -0.0075,
     pos = 1,
     "Distance (m)", xpd = T)

text(0, 0, "0", pos = 2)
text(0, ymax, ymax, pos = 2)
text(0, 0, "0", pos = 1)
text(xmax, 0, xmax, pos = 1)

segments(x0 = -18,
         x1 = 518,
         y0 = 0,
         col = "grey", lwd = 0.1)

lines(density(dists[,"dist.pred[1,1]"],
              from = 0,
              to = max(dists[,"dist.pred[1,1]"])),
      xpd = F)



axis(1, at = c(0, 100, 200, 300, 400, 500), labels = NA,
     tck = 0.025,
     #lwd = 0.5,
     tick = T)

tail.val <- length(which(dists[,"dist.pred[1,1]"] > xmax)) / length(which(dists[,"dist.pred[1,1]"] > 0))

text(xmax*1.1, tail.val*0.475,
     cex = 0.75,
     pos = 4,
     paste("Probability of\ndispersal\n>", xmax, " m", sep = ""), xpd = T)

text(xmax*1.1, ymax * 0.9, paste(round(median(dists[,"dist.pred[1,1]"])), "m", sep = " "), pos = 2, cex = 0.75)
text(xmax*1.1, ymax * 0.75, paste(round(quantile(dists[,"dist.pred[1,1]"], 0.99)), "m", sep = " "), pos = 2, cex = 0.75)


text(xmax*1.1, ymax * 0.9, "median", pos = 4, cex = 0.75, xpd = T)
text(xmax*1.1, ymax * 0.75, "99th percentile", pos = 4, cex = 0.75, xpd = T)

text(xmax*1.125, c(ymax * 0.9, ymax * 0.75, tail.val),  "-")

points(x = xmax,
       y = tail.val,
       pch = 16, cex = 1, xpd = F)
dev.off()
par(op)

setwd(top.wd)


















times.med <- apply(times, 2, median) / 60
times.025 <- apply(times, 2, function(x) quantile(x, 0.025)) / 60
times.975 <- apply(times, 2, function(x) quantile(x, 0.975)) / 60

op <- par()

setwd(figures.wd)
pdf("gut passage time.pdf", width = 6, height = 4.5)



par(mfrow = c(3,7),
    mar = c(1.1,1.1,1.1,0.1))

med.matrix <- matrix(NA, nrow = 15, ncol = 5)

for(j in new.tree.species.order.num){ #1:length(gut_tree_spp)){

  if(j %in% new.tree.species.order.num[c(1,13)]){
    plot.new()
  }
  if(j == 10){
    plot(NA,
         xlim = c(0.5,6.5),
         ylim = c(2,500),
         frame = F,
         log = "y",
         xaxt = "n",
         yaxt = "n",
         las = 1,
         xlab = "",
         ylab = "")
    mtext(side = 2, line = -3, text = "Gut passage time\n(minutes)", cex = 0.8)
  }

  plot(NA,
       xlim = c(0.5,6.5),
       ylim = c(2,500),
       frame = F,
       log = "y",
       xaxt = "n",
       yaxt = "n",
       las = 1,
       xlab = "",
       ylab = "")
  text(0, 700, tree.sp.names[j], font = 3, xpd = T, pos = 4, cex = 0.82)
  if(j %in% new.tree.species.order.num[c(1,7,13)]){
    axis(2, at = 2^seq(from = 1, to = 9, by = 2), las = 1)
  } else{
    axis(2, at = 2^seq(from = 1, to = 9, by = 2), las = 1,
         labels = rep("",5))
  }



  for(i in 1:length(gut_bird_spp)){

    bt.col <- paste("log.time.pred[",i,",",j,"]", sep = "")
    if(bt.col %in% colnames(times)){

      points(i, times.med[bt.col], pch = 20+i)
      arrows(x0 = i,
             y0 = times.025[bt.col],
             y1 = times.975[bt.col],
             code = 3,
             angle = 90,
             length = 0.02)
      med.matrix[j,i] <- times.med[bt.col]
    }
  }
}




par(new = TRUE)
par(fig = c(0.5, 0.9, 0, 0.35),#c(0.5, 1, 0, 0.25),
    mar = c(10,10,10,10),
    pin=c(1,1),
    xpd = T)

if(j == tail(new.tree.species.order.num,1)){
  plot(NA,
       xlim = c(0.5,6.5),
       ylim = c(2,500),
       frame = F,
       log = "y",
       xaxt = "n",
       yaxt = "n",
       las = 1,
       xlab = "",
       ylab = "")
  points(x=rep(0, 5), y = 2^seq(from = 8, to = 2, length.out = 5),
         pch=c(21:25), xpd = T)
  text(x=rep(0, 5), y = 2^seq(from = 7.9, to = 1.9, length.out = 5),
       labels = c("Micronesian Starling", "Mariana Fruit Dove", "Golden White-eye",
                  "Bridled White-eye", "White-throated Ground Dove"),
       cex = 0.95,
       xpd = T,
       pos = 4)
}





dev.off()

par(op)

setwd(top.wd)