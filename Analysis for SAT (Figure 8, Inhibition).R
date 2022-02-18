
dane1 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_1.csv", row.names=1)
dane2 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_2.csv", row.names=1)
dane3 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_3.csv", row.names=1)
dane4 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_4.csv", row.names=1)
dane5 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_5.csv", row.names=1)
dane6 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_6.csv", row.names=1)
dane7 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_7.csv", row.names=1)
dane8 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/SAT_with_inhibition_8.csv", row.names=1)

dane <- merge(dane1, dane2, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane3, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane4, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane5, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane6, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane7, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane8, all.x = TRUE, all.y = TRUE)

data_size <- 40000

Select_contrast <- function(level, low_bound, high_bound) {

  selection <- numeric(data_size)
  if(level == 1) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 8.247)}
  if(level == 2) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 8.4585)}
  if(level == 3) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 8.7315)}
  if(level == 4) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 9.081)}
  if(level == 5) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 9.5295)}

return(selection)  
}

#selection <- Select_contrast(5, 0,3500)
#sum(selection)

Calculate_mean_inh <- function(selection, num_of_bins, width) {
  
  inh_lev <- numeric(num_of_bins)

  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    inh_lev[i] <- mean(dane$inhibition_level[selection & interval])
  }
  return(inh_lev)  
}

Calculate_d_prime <- function(selection, num_of_bins, errs, width) {
  
  correct <- numeric(num_of_bins)
  error <- numeric(num_of_bins)
  num_of_cors <- numeric(num_of_bins)
  num_of_errs <- numeric(num_of_bins)
  proc_correct <- numeric(num_of_bins)
  d_prime <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    correct[i] <- sum(selection & interval & dane$dec_50 == 1)
    error[i] <- sum(selection & interval & dane$dec_50 == -1)
    proc_correct[i] <- correct[i]/(error[i] + correct[i])
    num_of_cors[i] <- sum(selection & interval & (dane$dec_50 == 1))
    num_of_errs[i] <- sum(selection & interval & (dane$dec_50 == -1))
    if(num_of_errs[i] < errs) {proc_correct[i] <- NaN} 
    d_prime[i] <- 2*qnorm(proc_correct[i],mean = 0)
  }
  return(c(d_prime))  
  
}

Calculate_mean_RT <- function(selection, num_of_bins, errs, width) {
  
  mean_RT <- numeric(num_of_bins)
  num_of_cors <- numeric(num_of_bins)
  num_of_errs <- numeric(num_of_bins)

  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    mean_RT[i] <- mean(dane$time_50[selection & interval]) - 1000
    num_of_cors[i] <- sum(selection & interval & (dane$dec_50 == 1))
    num_of_errs[i] <- sum(selection & interval & (dane$dec_50 == -1))
    if(num_of_errs[i] < errs) {mean_RT[i] <- NaN} 
  }
  return(c(mean_RT))  
  
}

Calculate_stdev_RT <- function(selection, num_of_bins, errs, width) {
  
  stdev_RT <- numeric(num_of_bins)
  num_of_cors <- numeric(num_of_bins)
  num_of_errs <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    stdev_RT[i] <- sd(dane$time_50[selection & interval])
    num_of_cors[i] <- sum(selection & interval & (dane$dec_50 == 1))
    num_of_errs[i] <- sum(selection & interval & (dane$dec_50 == -1))
    if(num_of_errs[i] < errs) {stdev_RT[i] <- NaN} 
  }
  return(c(stdev_RT))  
  
}

Calculate_skew <- function(selection, num_of_bins, errs, width) {
  
  skew <- numeric(num_of_bins)
  num_of_cors <- numeric(num_of_bins)
  num_of_errs <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    skew[i] <- skewness(dane$time_50[selection & interval])
    num_of_cors[i] <- sum(selection & interval & (dane$dec_50 == 1))
    num_of_errs[i] <- sum(selection & interval & (dane$dec_50 == -1))
    if(num_of_errs[i] < errs) {skew[i] <- NaN} 
  }
  return(c(skew))  
  
}

Calculate_RT_diff <- function(selection, num_of_bins, errs, width) {
  
  correct_RT <- numeric(num_of_bins)
  error_RT <- numeric(num_of_bins)
  num_of_cors <- numeric(num_of_bins)
  num_of_errs <- numeric(num_of_bins)
  diff_RT <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    correct_RT[i] <- mean(dane$time_50[selection & interval & dane$dec_50 == 1])
    error_RT[i] <- mean(dane$time_50[selection & interval & dane$dec_50 == -1])
    diff_RT[i] <- error_RT[i] - correct_RT[i]
    num_of_cors[i] <- sum(selection & interval & (dane$dec_50 == 1))
    num_of_errs[i] <- sum(selection & interval & (dane$dec_50 == -1))
    if(num_of_errs[i] < errs) {diff_RT[i] <- NaN} 
  }
  return(c(diff_RT))  
}

#--------------------------
# Mean RT - not used
#--------------------------

opis <- 1.2
osie <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

low_bound <- 1000
high_bound <- 3000

width <- 0.02
errs <- 0
num_of_bins <- 7
l <- 1

selection <- Select_contrast(l, low_bound, high_bound)
removed <- 8000 - sum(selection)

inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
sterr <- (Calculate_stdev_RT(selection, num_of_bins, errs, width))^(1/2)
mean_RT <- Calculate_mean_RT(selection, num_of_bins, errs, width)

plot(inh_levels, mean_RT, type = "l", xlim = c(0.2, 0.1), ylim = c(400, 800), 
     xlab = "Inhibition", ylab = "Mean RT (ms)", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)
for(j in 1:7) {
  lines(c(inh_levels[j], inh_levels[j]), c(mean_RT[j] - sterr[j], mean_RT[j] + sterr[j]), lwd = 1, col = kolor[l])
}
points(inh_levels, mean_RT - sterr, cex = 1.5, pch = 45, col = kolor[l])
points(inh_levels, mean_RT + sterr, cex = 1.5, pch = 45, col = kolor[l])

for(l in 2:5) {
  selection <- Select_contrast(l, low_bound, high_bound)
  removed <- removed + 8000 - sum(selection)
  
  inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
  sterr <- (Calculate_stdev_RT(selection, num_of_bins, errs, width))^(1/2)
  mean_RT <- Calculate_mean_RT(selection, num_of_bins, errs, width)
  
  lines(inh_levels, mean_RT, type = "l", col = kolor[l], lwd = l)
  for(j in 1:7) {
    lines(c(inh_levels[j], inh_levels[j]), c(mean_RT[j] - sterr[j], mean_RT[j] + sterr[j]), lwd = 1, col = kolor[l])
  }
  points(inh_levels, mean_RT - sterr, cex = 1.5, pch = 45, col = kolor[l])
  points(inh_levels, mean_RT + sterr, cex = 1.5, pch = 45, col = kolor[l])
  
}

legend(0.2, 800, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
#abline(v = 0, lty = 2)

#-------------
# in the text

removed

mean_RT
sterr/(sum(selection)^(0.5))

selection <- Select_contrast(3, low_bound, high_bound)
sterr <- (Calculate_stdev_RT(selection, num_of_bins, errs, width))^(1/2)
mean_RT <- Calculate_mean_RT(selection, num_of_bins, errs, width)

mean_RT
sterr/(sum(selection)^(0.5))


#--------------------------
# Figure 8a: D-prime
#--------------------------

opis <- 1.2
osie <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

low_bound <- 1000
high_bound <- 3000

width <- 0.02
errs <- 0
num_of_bins <- 7
l <- 1

selection <- Select_contrast(l, low_bound, high_bound)
srednia <- Calculate_mean_RT(selection, num_of_bins, errs, width)
d_prime <- Calculate_d_prime(selection, num_of_bins, errs, width)

plot(srednia, d_prime, type = "l", xlim = c(200, 700), ylim = c(0,5),
     xlab = "Mean RT (ms)", ylab = "d'", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_contrast(l, low_bound, high_bound)
  srednia <- Calculate_mean_RT(selection, num_of_bins, errs, width)
  d_prime <- Calculate_d_prime(selection, num_of_bins, errs, width)
  
  lines(srednia, d_prime, type = "l", col = kolor[l], lwd = l)
}

legend(200, 5, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")

#---------------
# in the text

d_prime
srednia

selection <- Select_contrast(3, low_bound, high_bound)
srednia <- Calculate_mean_RT(selection, num_of_bins, errs, width)
d_prime <- Calculate_d_prime(selection, num_of_bins, errs, width)

d_prime
srednia



#------------------------------
# Figure 8b: RT difference

opis <- 1.2
osie <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

low_bound <- 1000
high_bound <- 3000

width <- 0.02
errs <- 50
num_of_bins <- 7
l <- 1

selection <- Select_contrast(l, low_bound, high_bound)
inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
RT_diff <- Calculate_RT_diff(selection, num_of_bins, errs, width)

plot(inh_levels, RT_diff, type = "l", ylim = c(-100, 200), xlim = c(0.2, 0.1),
     xlab = "Inhibition", ylab = "RT(E) - RT(C)", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_contrast(l, low_bound, high_bound)
  inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
  RT_diff <- Calculate_RT_diff(selection, num_of_bins, errs, width)
  
  lines(inh_levels, RT_diff, type = "l", col = kolor[l], lwd = l)
}

legend(0.2, 200, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
abline(h = 0, lty = 2)
#abline(v = 0, lty = 2)


#-----------------
# T-test

l <- 4
step <- 0.1/num_of_bins
i <- 5
interval <- (dane$inhibition_level > (0.1 + (i-1)*step)) & (dane$inhibition_level < (0.1 + i*step))
selection <- Select_contrast(l, low_bound, high_bound)
t.test(dane$time_50[dane$dec_50 != 0 & selection & interval] ~ dane$dec_50[dane$dec_50 != 0 & selection & interval])

inh_levels

#------------------------------
# Figure 8c: SD to mean ratio
#------------------------------

opis <- 1.2
osie <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

low_bound <- 1000
high_bound <- 3500

width <- 0.02
errs <- 0
num_of_bins <- 7
l <- 1

selection <- Select_contrast(l, low_bound, high_bound)
inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
stdev <- Calculate_stdev_RT(selection, num_of_bins, errs, width)
mean_RT <- Calculate_mean_RT(selection, num_of_bins, errs, width)
ratio <- stdev/mean_RT

plot(inh_levels, ratio, type = "l", xlim = c(0.2, 0.1), ylim = c(0.3,0.6),
     xlab = "Inhibition", ylab = "SD / Mean RT", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_contrast(l, low_bound, high_bound)
  inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
  stdev <- Calculate_stdev_RT(selection, num_of_bins, errs, width)
  mean_RT <- Calculate_mean_RT(selection, num_of_bins, errs, width)
  ratio <- stdev/mean_RT
  
  lines(inh_levels, ratio, type = "l", col = kolor[l], lwd = l)
}

legend(0.2, 0.45, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
#abline(v = 0, lty = 2)



#--------------------
# Skewness

library(e1071)


opis <- 1.2
osie <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

low_bound <- 1000
high_bound <- 3500

width <- 0.02
errs <- 0
num_of_bins <- 7
l <- 1

selection <- Select_contrast(l, low_bound, high_bound)
inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
skew <- Calculate_skew(selection, num_of_bins, errs, width)

plot(inh_levels, skew, type = "l", xlim = c(0.2, 0.1), ylim = c(1,3), 
     xlab = "Inhibition", ylab = "RT Skewness", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_contrast(l, low_bound, high_bound)
  inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
  skew <- Calculate_skew(selection, num_of_bins, errs, width)
  
  lines(inh_levels, skew, type = "l", col = kolor[l], lwd = l)
}

legend(0.2, 3, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
#abline(v = 0, lty = 2)

#-----------------------------------
# Inhibition levels vs. baseline firing rate

# levels_for_baseFR/5
# from data: 0.1137982 0.1215495 0.1358243 0.1501126 0.1640973 0.1785127 0.1864704
# simulated: 0.1138 0.1215 0.1358 0.1501 0.1641 0.1785 0.1865

#-----------------------------
# Figure 7a: BASELINE FIRING RATE vs. INHIBITION
#-----------------------------

opis <- 1.2
osie <- 1.2

dane <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/Dist_FR_Spon_vs_Inhibition.csv", row.names=1)

data_size <- length(dane$inhibition_level)
data_size

plot(dane$inhibition_level + rnorm(data_size, mean = 0, sd = 0.002), 
     dane$firing_rate_A, type = "p", xlim = c(0.2, 0.1), ylim = c(0,40), pch = 16, cex = 0.2, 
     xlab = "Inhibition", ylab = "Baseline Firing Rate (Hz)", 
     col = "light grey", cex.lab = opis, cex.axis = osie)
points(dane$inhibition_level + rnorm(data_size, mean = 0, sd = 0.002), dane$firing_rate_B,
       pch = 16, cex = 0.2, col = "light grey")

# simulated: 0.1138 0.1215 0.1358 0.1501 0.1641 0.1785 0.1865

inh_levels <- numeric(7)
inh_levels[1] <- 0.1138
inh_levels[2] <- 0.1215
inh_levels[3] <- 0.1358
inh_levels[4] <- 0.1501
inh_levels[5] <- 0.1641
inh_levels[6] <- 0.1785
inh_levels[7] <- 0.1865

mean_baseline <- numeric(7)
for(l in 1:7) {
  mean_baseline[l] <- 1/2*(mean(dane$firing_rate_A[dane$inhibition_level == inh_levels[l]]) + mean(dane$firing_rate_B[dane$inhibition_level == inh_levels[l]]))
}

sd_baseline <- numeric(7)
for(l in 1:7) {
  sd_baseline[l] <- 1/2*(sd(dane$firing_rate_A[dane$inhibition_level == inh_levels[l]]) + sd(dane$firing_rate_B[dane$inhibition_level == inh_levels[l]]))
}


lines(inh_levels, mean_baseline, type = "l", col = "black")
points(inh_levels, mean_baseline, pch = 19, cex = 1.5, col = "black")
for(l in 1:7) {
  lines(c(inh_levels[l], inh_levels[l]), c(mean_baseline[l] - sd_baseline[l], mean_baseline[l] + sd_baseline[l]), lwd = 1)
}
points(inh_levels, mean_baseline - sd_baseline, cex = 1.5, pch = 45, col = "black")
points(inh_levels, mean_baseline + sd_baseline, cex = 1.5, pch = 45, col = "black")

legend(0.14, 40, pch = c(19, NaN), lwd = c(1, 1), c("Mean", "SD"), cex = opis, bty = "n")

#-----------------------------
# Statistics

mean_baseline
sd_baseline

baseline <- c(dane$firing_rate_A, dane$firing_rate_B)
inhibition_v <- c(dane$inhibition_level, dane$inhibition_level)

equation <- baseline ~ inhibition_v
fit <- lm(equation)
summary(fit)
anova(fit)

#----------------------------
# DISPLACEMENT OF THE ATTRACTOR SETS
# not used
#----------------------------

osie <- 1.2
opis <- 1.2

dane1 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/Method/Dist_FR_A.csv", row.names=1)
dane2 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/Dist_FR_A_high_contrast.csv", row.names=1)

kolor1 <- "red"
kolor2 <- "orange"

plot(dane1[1:500, ]$firing_rate_A, dane1[1:500, ]$firing_rate_B, 
     ylim = c(0,70), xlim = c(0,70),
     type = "p",  
     xlab = "Firing Rate in A (Hz)", ylab = "Firing Rate in B (Hz)", 
     #xlab = "", ylab = "", 
     pch = 16, cex = 0.3, cex.axis = osie, cex.lab = opis, col = kolor1)

points(dane2$firing_rate_A, dane2$firing_rate_B, 
       pch = 16, cex = 0.3, col = kolor2)

legend(0, 70, legend = c("No inputs (inh: 0.13)", "Contrast: 13.53% (inh: 0.18)"), cex = opis, pch = 16, col=c(kolor1, kolor2), 
       bty = "n")



#----------------------------
# STABILITY vs. INHIBITION
# not used
#----------------------------

dane <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inhibition/Stability_vs_inhibition.csv", row.names=1)


Calculate_stable_runs <- function(selection, num_of_bins, width) {
  
  stable_runs <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$inhibition_level > (0.1 + (i-1/2)*step - width)) & (dane$inhibition_level < (0.1 + (i-1/2)*step + width))
    stable_runs[i] <- sum(dane$dec_50[selection & interval] == 0)/sum(selection & interval)
  }
  return(stable_runs)  
}


opis <- 1.2
osie <- 1.2

kolor <- "dark grey"

low_bound <- 0
high_bound <- 5000

width <- 0.01
num_of_bins <- 10

selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound)
inh_levels <- Calculate_mean_inh(selection, num_of_bins, width)
stable_runs <- Calculate_stable_runs(selection, num_of_bins, width)

plot(0.13 - inh_levels, stable_runs, type = "l", xlim = c(0.03, -0.07), ylim = c(0,1),
     xlab = "Difference vs. Baseline Inhibition Level", ylab = "% Stable Runs (over 5s)", 
     col = kolor, lwd = 2, cex.lab = opis, cex.axis = osie)


