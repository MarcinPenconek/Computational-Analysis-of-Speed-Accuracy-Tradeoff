
dane1 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Alpha/SAT_with_Random_Alpha_1.csv", row.names=1)
dane2 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Alpha/SAT_with_Random_Alpha_2.csv", row.names=1)
dane3 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Alpha/SAT_with_Random_Alpha_3.csv", row.names=1)

dane <- merge(dane1, dane2, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane3, all.x = TRUE, all.y = TRUE)

data_size <- 15000

Select_contrast_50 <- function(level, low_bound, high_bound) {
  
  selection <- numeric(data_size)
  if(level == 1) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 8.247)}
  if(level == 2) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 8.4585)}
  if(level == 3) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 8.7315)}
  if(level == 4) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 9.081)}
  if(level == 5) {selection <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == 9.5295)}
  
  return(selection)  
}

Select_contrast_45 <- function(level, low_bound, high_bound) {
  
  selection <- numeric(data_size)
  if(level == 1) {selection <- (dane$time_45 <= high_bound & dane$time_45 >= low_bound & dane$n_A == 8.247)}
  if(level == 2) {selection <- (dane$time_45 <= high_bound & dane$time_45 >= low_bound & dane$n_A == 8.4585)}
  if(level == 3) {selection <- (dane$time_45 <= high_bound & dane$time_45 >= low_bound & dane$n_A == 8.7315)}
  if(level == 4) {selection <- (dane$time_45 <= high_bound & dane$time_45 >= low_bound & dane$n_A == 9.081)}
  if(level == 5) {selection <- (dane$time_45 <= high_bound & dane$time_45 >= low_bound & dane$n_A == 9.5295)}
  
  return(selection)  
}

Select_contrast_40 <- function(level, low_bound, high_bound) {
  
  selection <- numeric(data_size)
  if(level == 1) {selection <- (dane$time_40 <= high_bound & dane$time_40 >= low_bound & dane$n_A == 8.247)}
  if(level == 2) {selection <- (dane$time_40 <= high_bound & dane$time_40 >= low_bound & dane$n_A == 8.4585)}
  if(level == 3) {selection <- (dane$time_40 <= high_bound & dane$time_40 >= low_bound & dane$n_A == 8.7315)}
  if(level == 4) {selection <- (dane$time_40 <= high_bound & dane$time_40 >= low_bound & dane$n_A == 9.081)}
  if(level == 5) {selection <- (dane$time_40 <= high_bound & dane$time_40 >= low_bound & dane$n_A == 9.5295)}
  
  return(selection)  
}

Select_contrast_35 <- function(level, low_bound, high_bound) {
  
  selection <- numeric(data_size)
  if(level == 1) {selection <- (dane$time_35 <= high_bound & dane$time_35 >= low_bound & dane$n_A == 8.247)}
  if(level == 2) {selection <- (dane$time_35 <= high_bound & dane$time_35 >= low_bound & dane$n_A == 8.4585)}
  if(level == 3) {selection <- (dane$time_35 <= high_bound & dane$time_35 >= low_bound & dane$n_A == 8.7315)}
  if(level == 4) {selection <- (dane$time_35 <= high_bound & dane$time_35 >= low_bound & dane$n_A == 9.081)}
  if(level == 5) {selection <- (dane$time_35 <= high_bound & dane$time_35 >= low_bound & dane$n_A == 9.5295)}
  
  return(selection)  
}


Calculate_mean_amplitude <- function(selection, num_of_bins, width) {
  
  amp_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    amp_lev[i] <- mean(dane$amplitude[selection & interval])
  }
  return(amp_lev)  
}

Calculate_mean_RT_50Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- mean(dane$time_50[selection & interval] - 1000)
  }
  return(RT_lev)  
}

Calculate_se_RT_50Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- sd(dane$time_50[selection & interval] - 1000)/(sum(selection & interval)^0.5)
  }
  return(RT_lev)  
}


Calculate_mean_RT_45Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- mean(dane$time_45[selection & interval] - 1000)
  }
  return(RT_lev)  
}

Calculate_se_RT_45Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- sd(dane$time_45[selection & interval] - 1000)/(sum(selection & interval)^0.5)
  }
  return(RT_lev)  
}

Calculate_mean_RT_40Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- mean(dane$time_40[selection & interval] - 1000)
  }
  return(RT_lev)  
}

Calculate_se_RT_40Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- sd(dane$time_40[selection & interval] - 1000)/(sum(selection & interval)^0.5)
  }
  return(RT_lev)  
}

Calculate_mean_RT_35Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- mean(dane$time_35[selection & interval] - 1000)
  }
  return(RT_lev)  
}

Calculate_se_RT_35Hz <- function(selection, num_of_bins, width) {
  
  RT_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    RT_lev[i] <- sd(dane$time_35[selection & interval] - 1000)/(sum(selection & interval)^0.5)
  }
  return(RT_lev)  
}

#-------------------------------------
# FIGURE 9f: Mean RT vs. Amplitude 
# frequency: 35, 40, 45 and 50Hz
# contrast level l = 1...5

opis <- 1.2
osie <- 1.2

low_bound <- 1000
high_bound <- 3500
non_decision <- 0

width <- 0.02
num_of_bins <- 12

data_size <- 15000

lev <- numeric(5)
lev[1] <- 8.247
lev[2] <- 8.4585
lev[3] <- 8.7315
lev[4] <- 9.081
lev[5] <- 9.5295

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"


l <- 3

selection <- Select_contrast_50(l, low_bound, high_bound)
amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
sterr <- Calculate_se_RT_50Hz(selection, num_of_bins, width)
mean_RT <- Calculate_mean_RT_50Hz(selection, num_of_bins, width)

sterr
mean_RT
amplitude_levels


plot(amplitude_levels[1:10], mean_RT[1:10], type = "l", xlim = c(0, 0.1), ylim = c(200, 1100), 
     xlab = "Amplitude (alpha = 8-12Hz)", ylab = "Mean RT (ms)", 
     col = kolor[l], lwd = l, cex.lab = opis, cex.axis = osie)
#for(j in 1:num_of_bins) {
#  lines(c(amplitude_levels[j], amplitude_levels[j]), c(mean_RT[j] - sterr[j], mean_RT[j] + sterr[j]), lwd = 1, col = kolor[l])
#}
#points(amplitude_levels, mean_RT - sterr, cex = 1.5, pch = 45, col = kolor[l])
#points(amplitude_levels, mean_RT + sterr, cex = 1.5, pch = 45, col = kolor[l])

selection <- Select_contrast_45(l, low_bound, high_bound)
amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
sterr <- Calculate_se_RT_45Hz(selection, num_of_bins, width)
mean_RT <- Calculate_mean_RT_45Hz(selection, num_of_bins, width)

sterr
mean_RT

lines(amplitude_levels, mean_RT, type = "l", col = kolor[l], lwd = l, lty = 3)

selection <- Select_contrast_40(l, low_bound, high_bound)
amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
sterr <- Calculate_se_RT_40Hz(selection, num_of_bins, width)
mean_RT <- Calculate_mean_RT_40Hz(selection, num_of_bins, width)

sterr
mean_RT
amplitude_levels

lines(amplitude_levels, mean_RT, type = "l", col = kolor[l], lwd = l, lty = 2)

selection <- Select_contrast_35(l, low_bound, high_bound)
amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
sterr <- Calculate_se_RT_35Hz(selection, num_of_bins, width)
mean_RT <- Calculate_mean_RT_35Hz(selection, num_of_bins, width)

sterr
mean_RT

lines(amplitude_levels, mean_RT, type = "l", col = kolor[l], lwd = l, lty = 5)
legend(0, 1100, col = kolor[l], lwd = l, lty = c(1, 3, 20, 5), c("Threshold 50Hz", "45Hz", "40Hz", "35Hz"), cex = opis, bty = "n")

#---------------------
# statistics


Select_contrast_correct <- function(level) {
  
  selection <- numeric(data_size)
  if(level == 1) {selection <- dane$n_A == 8.247}
  if(level == 2) {selection <- dane$n_A == 8.4585}
  if(level == 3) {selection <- dane$n_A == 8.7315}
  if(level == 4) {selection <- dane$n_A == 9.081}
  if(level == 5) {selection <- dane$n_A == 9.5295}
  
  return(selection)  
}

l <- 3
selection <- Select_contrast_correct(l)
of_interest <- dane$dec_50 != 0 & dane$time_50 < 3500 & dane$amplitude < 0.06

equation <- dane$time_50[selection & of_interest] ~ dane$amplitude[selection & of_interest]
model <- lm(equation)
summary(model)
anova(model)

l <- 3
selection <- Select_contrast_correct(l)
of_interest <- dane$dec_40 != 0 & dane$time_40 < 3500 & dane$amplitude < 0.06

equation <- dane$time_40[selection & of_interest] ~ dane$amplitude[selection & of_interest]
model <- lm(equation)
summary(model)
anova(model)

#---------------------
# FIGURE: Accuracy vs. Amplitude

Select_contrast_correct <- function(level) {
  
  selection <- numeric(data_size)
  if(level == 1) {selection <- dane$n_A == 8.247}
  if(level == 2) {selection <- dane$n_A == 8.4585}
  if(level == 3) {selection <- dane$n_A == 8.7315}
  if(level == 4) {selection <- dane$n_A == 9.081}
  if(level == 5) {selection <- dane$n_A == 9.5295}
  
  return(selection)  
}

Calculate_correct_50Hz <- function(selection, num_of_bins, width, low_bound, high_bound) {
  
  correct_lev <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > ((i-1/2)*step - width)) & (dane$amplitude < ((i-1/2)*step + width))
    correct_decisions <- dane$dec_50 == 1 & dane$time_50 >= low_bound & dane$time_50 <= high_bound
    correct_lev[i] <- mean(correct_decisions[selection & interval])
  }
  return(correct_lev)  
}


#-------------------
# Figure 9e

low_bound <- 1000
high_bound <- 3500

width <- 0.02
num_of_bins <- 10

l <- 1

selection <- Select_contrast_correct(l)
amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
correct <- Calculate_correct_50Hz(selection, num_of_bins, width, low_bound, high_bound)

correct

plot(amplitude_levels, correct, type = "l", xlim = c(0, 0.1), ylim = c(0.5, 1), 
     xlab = "Amplitude (alpha = 8-12Hz)", ylab = "Correct Decisions (%)", 
     col = kolor[l], lwd = l, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_contrast_correct(l)
  amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
  correct <- Calculate_correct_50Hz(selection, num_of_bins, width, low_bound, high_bound)
  
  correct
  
  lines(amplitude_levels, correct, type = "l", col = kolor[l], lwd = l)
}

legend(0, 0.75, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")

#-------------------------------------------
# Statistics

l <- 5
selection <- Select_contrast_correct(l)

range <- dane$amplitude < 0.02
sum(selection & range)

correct_vs_all <- (dane$dec_50 == 1 & dane$time_50 <= 3500)
mean(correct_vs_all[selection & range])
sd(correct_vs_all[selection & range])/(sum(selection & range)^0.5)

#-------------
# peak values

l <- 5
selection <- Select_contrast_correct(l)

range <- dane$amplitude > 0.05 & dane$amplitude < 0.06
sum(selection & range)

correct_vs_all <- (dane$dec_50 == 1 & dane$time_50 <= 3500)
mean(correct_vs_all[selection & range])
sd(correct_vs_all[selection & range])
sd(correct_vs_all[selection & range])/(sum(selection & range)^0.5)

#----------------

#install.packages("oglmx")
library(oglmx)

l <- 1
selection <- Select_contrast_correct(l)

range <- dane$amplitude < 0.06
sum(range)
correct_vs_all <- (dane$dec_50 == 1 & dane$time_50 >= 1000 & dane$time_50 <= 3500)

equation <- correct_vs_all[range & selection] ~ dane$amplitude[range & selection]
model <- probit.reg(equation)
summary(model)

#----------------
l <- 3
selection <- Select_contrast_correct(l)

range <- dane$amplitude > 0.08
sum(selection & range)

correct_vs_all <- (dane$dec_50 == 1 & dane$time_50 <= 3500)
mean(correct_vs_all[selection & range])
sd(correct_vs_all[selection & range])/(sum(selection & range)^0.5)

error_vs_all <- (dane$dec_50 == -1 & dane$time_50 <= 3500)
mean(error_vs_all[selection & range])
sum(error_vs_all[selection & range])
no_dec_vs_all <- (dane$dec_50 == 0 & dane$time_50 <= 3500)
mean(no_dec_vs_all[selection & range])
sd(no_dec_vs_all[selection & range])/(sum(selection & range)^0.5)

#--------------------------------------------
# not used


Select_observations <- function(lev, low_bound, high_bound) {
  selection <- matrix(1, nrow = 12, ncol = data_size)
  
  selection[1, ] <- (dane$time_10 <= high_bound & dane$time_10 >= low_bound & dane$n_A == lev)
  selection[2, ] <- (dane$time_15 <= high_bound & dane$time_15 >= low_bound & dane$n_A == lev)
  selection[3, ] <- (dane$time_20 <= high_bound & dane$time_20 >= low_bound & dane$n_A == lev)
  selection[4, ] <- (dane$time_25 <= high_bound & dane$time_25 >= low_bound & dane$n_A == lev)
  selection[5, ] <- (dane$time_30 <= high_bound & dane$time_30 >= low_bound & dane$n_A == lev)
  selection[6, ] <- (dane$time_35 <= high_bound & dane$time_35 >= low_bound & dane$n_A == lev)
  selection[7, ] <- (dane$time_40 <= high_bound & dane$time_40 >= low_bound & dane$n_A == lev)
  selection[8, ] <- (dane$time_45 <= high_bound & dane$time_45 >= low_bound & dane$n_A == lev)
  selection[9, ] <- (dane$time_50 <= high_bound & dane$time_50 >= low_bound & dane$n_A == lev)
  selection[10, ] <- (dane$time_55 <= high_bound & dane$time_55 >= low_bound & dane$n_A == lev)
  selection[11, ] <- (dane$time_60 <= high_bound & dane$time_60 >= low_bound & dane$n_A == lev)
  selection[12, ] <- (dane$time_65 <= high_bound & dane$time_65 >= low_bound & dane$n_A == lev)
  
  for(i in 1:12) {
    print(mean(selection[i, ]))
  }
  
  return(selection)
}

Select_levels <- function() {
  cases <- numeric(12)
  cases <- c(0,0,1,1,1,1,1,1,1,0,0,0)
  
  return(cases)
}

Calculate_correct <- function(selection, cases) {
  A <- numeric(12)
  
  A[1] <- mean(dane$dec_10 == 1 & selection[1, ] == 1)
  A[2] <- mean(dane$dec_15 == 1 & selection[2, ] == 1)
  A[3] <- mean(dane$dec_20 == 1 & selection[3, ] == 1)
  A[4] <- mean(dane$dec_25 == 1 & selection[4, ] == 1)
  A[5] <- mean(dane$dec_30 == 1 & selection[5, ] == 1)
  A[6] <- mean(dane$dec_35 == 1 & selection[6, ] == 1)
  A[7] <- mean(dane$dec_40 == 1 & selection[7, ] == 1)
  A[8] <- mean(dane$dec_45 == 1 & selection[8, ] == 1)
  A[9] <- mean(dane$dec_50 == 1 & selection[9, ] == 1)
  A[10] <- mean(dane$dec_55 == 1 & selection[10, ] == 1)
  A[11] <- mean(dane$dec_60 == 1 & selection[11, ] == 1)
  A[12] <- mean(dane$dec_65 == 1 & selection[12, ] == 1)
  
  B <- numeric(12)
  
  B[1] <- mean(dane$dec_10 == - 1 & selection[1, ] == 1)
  B[2] <- mean(dane$dec_15 == - 1 & selection[2, ] == 1)
  B[3] <- mean(dane$dec_20 == - 1 & selection[3, ] == 1)
  B[4] <- mean(dane$dec_25 == - 1 & selection[4, ] == 1)
  B[5] <- mean(dane$dec_30 == - 1 & selection[5, ] == 1)
  B[6] <- mean(dane$dec_35 == - 1 & selection[6, ] == 1)
  B[7] <- mean(dane$dec_40 == - 1 & selection[7, ] == 1)
  B[8] <- mean(dane$dec_45 == - 1 & selection[8, ] == 1)
  B[9] <- mean(dane$dec_50 == - 1 & selection[9, ] == 1)
  B[10] <- mean(dane$dec_55 == - 1 & selection[10, ] == 1)
  B[11] <- mean(dane$dec_60 == - 1 & selection[11, ] == 1)
  B[12] <- mean(dane$dec_65 == - 1 & selection[12, ] == 1)
  
  correct <- A/(A+B)
  
  return(correct[cases == 1])
}

Calculate_RT_diff <- function(selection, cases) {
  
  C20 <- mean(dane$time_10[dane$dec_10 == 1 & selection[1, ]])
  C25 <- mean(dane$time_15[dane$dec_15 == 1 & selection[2, ]])
  C30 <- mean(dane$time_20[dane$dec_20 == 1 & selection[3, ]])
  C35 <- mean(dane$time_25[dane$dec_25 == 1 & selection[4, ]])
  C40 <- mean(dane$time_30[dane$dec_30 == 1 & selection[5, ]])
  C45 <- mean(dane$time_35[dane$dec_35 == 1 & selection[6, ]])
  C50 <- mean(dane$time_40[dane$dec_40 == 1 & selection[7, ]])
  C55 <- mean(dane$time_45[dane$dec_45 == 1 & selection[8, ]])
  C60 <- mean(dane$time_50[dane$dec_50 == 1 & selection[9, ]])
  C65 <- mean(dane$time_55[dane$dec_55 == 1 & selection[10, ]])
  C70 <- mean(dane$time_60[dane$dec_60 == 1 & selection[11, ]])
  C75 <- mean(dane$time_65[dane$dec_65 == 1 & selection[12, ]])
  
  E20 <- mean(dane$time_10[dane$dec_10 == -1 & selection[1, ]])
  E25 <- mean(dane$time_15[dane$dec_15 == -1 & selection[2, ]])
  E30 <- mean(dane$time_20[dane$dec_20 == -1 & selection[3, ]])
  E35 <- mean(dane$time_25[dane$dec_25 == -1 & selection[4, ]])
  E40 <- mean(dane$time_30[dane$dec_30 == -1 & selection[5, ]])
  E45 <- mean(dane$time_35[dane$dec_35 == -1 & selection[6, ]])
  E50 <- mean(dane$time_40[dane$dec_40 == -1 & selection[7, ]])
  E55 <- mean(dane$time_45[dane$dec_45 == -1 & selection[8, ]])
  E60 <- mean(dane$time_50[dane$dec_50 == -1 & selection[9, ]])
  E65 <- mean(dane$time_55[dane$dec_55 == -1 & selection[10, ]])
  E70 <- mean(dane$time_60[dane$dec_60 == -1 & selection[11, ]])
  E75 <- mean(dane$time_65[dane$dec_65 == -1 & selection[12, ]])
  
  RT_error_minus_correct <- numeric(12)
  
  RT_error_minus_correct[1] <- E20 - C20
  RT_error_minus_correct[2] <- E25 - C25
  RT_error_minus_correct[3] <- E30 - C30
  RT_error_minus_correct[4] <- E35 - C35
  RT_error_minus_correct[5] <- E40 - C40
  RT_error_minus_correct[6] <- E45 - C45
  RT_error_minus_correct[7] <- E50 - C50
  RT_error_minus_correct[8] <- E55 - C55
  RT_error_minus_correct[9] <- E60 - C60
  RT_error_minus_correct[10] <- E65 - C65
  RT_error_minus_correct[11] <- E70 - C70
  RT_error_minus_correct[12] <- E75 - C75
  
  return(RT_error_minus_correct[cases == 1])  
}

Calculate_RT_mean <- function(selection, cases, non_decision) {
  
  srednia <- numeric(12)
  
  srednia[1] <- mean(dane$time_10[dane$dec_10 != 0 & selection[1, ]] - 1000 + non_decision)
  srednia[2] <- mean(dane$time_15[dane$dec_15 != 0 & selection[2, ]] - 1000 + non_decision)
  srednia[3] <- mean(dane$time_20[dane$dec_20 != 0 & selection[3, ]] - 1000 + non_decision)
  srednia[4] <- mean(dane$time_25[dane$dec_25 != 0 & selection[4, ]] - 1000 + non_decision)
  srednia[5] <- mean(dane$time_30[dane$dec_30 != 0 & selection[5, ]] - 1000 + non_decision)
  srednia[6] <- mean(dane$time_35[dane$dec_35 != 0 & selection[6, ]] - 1000 + non_decision)
  srednia[7] <- mean(dane$time_40[dane$dec_40 != 0 & selection[7, ]] - 1000 + non_decision)
  srednia[8] <- mean(dane$time_45[dane$dec_45 != 0 & selection[8, ]] - 1000 + non_decision)
  srednia[9] <- mean(dane$time_50[dane$dec_50 != 0 & selection[9, ]] - 1000 + non_decision)
  srednia[10] <- mean(dane$time_55[dane$dec_55 != 0 & selection[10, ]] - 1000 + non_decision)
  srednia[11] <- mean(dane$time_60[dane$dec_60 != 0 & selection[11, ]] - 1000 + non_decision)
  srednia[12] <- mean(dane$time_65[dane$dec_65 != 0 & selection[12, ]] - 1000 + non_decision)
  
  return(srednia[cases == 1])
}


Calculate_RT_sd <- function(selection, cases, non_decision) {
  
  stdev <- numeric(12)
  
  stdev[1] <- sd(dane$time_10[dane$dec_10 != 0 & selection[1, ]] - 1000 + non_decision)
  stdev[2] <- sd(dane$time_15[dane$dec_15 != 0 & selection[2, ]] - 1000 + non_decision)
  stdev[3] <- sd(dane$time_20[dane$dec_20 != 0 & selection[3, ]] - 1000 + non_decision)
  stdev[4] <- sd(dane$time_25[dane$dec_25 != 0 & selection[4, ]] - 1000 + non_decision)
  stdev[5] <- sd(dane$time_30[dane$dec_30 != 0 & selection[5, ]] - 1000 + non_decision)
  stdev[6] <- sd(dane$time_35[dane$dec_35 != 0 & selection[6, ]] - 1000 + non_decision)
  stdev[7] <- sd(dane$time_40[dane$dec_40 != 0 & selection[7, ]] - 1000 + non_decision)
  stdev[8] <- sd(dane$time_45[dane$dec_45 != 0 & selection[8, ]] - 1000 + non_decision)
  stdev[9] <- sd(dane$time_50[dane$dec_50 != 0 & selection[9, ]] - 1000 + non_decision)
  stdev[10] <- sd(dane$time_55[dane$dec_55 != 0 & selection[10, ]] - 1000 + non_decision)
  stdev[11] <- sd(dane$time_60[dane$dec_60 != 0 & selection[11, ]] - 1000 + non_decision)
  stdev[12] <- sd(dane$time_65[dane$dec_65 != 0 & selection[12, ]] - 1000 + non_decision)
  
  return(stdev[cases == 1])
}


#install.packages("e1071")
library(e1071)

Calculate_RT_skew <- function(selection, cases, non_decision) {
  
  skew <- numeric(12)
  
  skew[1] <- skewness(dane$time_10[dane$dec_10 != 0 & selection[1, ]] - 1000 + non_decision)
  skew[2] <- skewness(dane$time_15[dane$dec_15 != 0 & selection[2, ]] - 1000 + non_decision)
  skew[3] <- skewness(dane$time_20[dane$dec_20 != 0 & selection[3, ]] - 1000 + non_decision)
  skew[4] <- skewness(dane$time_25[dane$dec_25 != 0 & selection[4, ]] - 1000 + non_decision)
  skew[5] <- skewness(dane$time_30[dane$dec_30 != 0 & selection[5, ]] - 1000 + non_decision)
  skew[6] <- skewness(dane$time_35[dane$dec_35 != 0 & selection[6, ]] - 1000 + non_decision)
  skew[7] <- skewness(dane$time_40[dane$dec_40 != 0 & selection[7, ]] - 1000 + non_decision)
  skew[8] <- skewness(dane$time_45[dane$dec_45 != 0 & selection[8, ]] - 1000 + non_decision)
  skew[9] <- skewness(dane$time_50[dane$dec_50 != 0 & selection[9, ]] - 1000 + non_decision)
  skew[10] <- skewness(dane$time_55[dane$dec_55 != 0 & selection[10, ]] - 1000 + non_decision)
  skew[11] <- skewness(dane$time_60[dane$dec_60 != 0 & selection[11, ]] - 1000 + non_decision)
  skew[12] <- skewness(dane$time_65[dane$dec_65 != 0 & selection[12, ]] - 1000 + non_decision)
  
  return(skew[cases == 1])
}


Calculate_FR_threshold <- function(selection, cases) {
  
  FR <- numeric(12)
  
  FR[1] <- mean(dane$bound_10[dane$dec_10 != 0 & selection[1, ]])
  FR[2] <- mean(dane$bound_15[dane$dec_15 != 0 & selection[2, ]])
  FR[3] <- mean(dane$bound_20[dane$dec_20 != 0 & selection[3, ]])
  FR[4] <- mean(dane$bound_25[dane$dec_25 != 0 & selection[4, ]])
  FR[5] <- mean(dane$bound_30[dane$dec_30 != 0 & selection[5, ]])
  FR[6] <- mean(dane$bound_35[dane$dec_35 != 0 & selection[6, ]])
  FR[7] <- mean(dane$bound_40[dane$dec_40 != 0 & selection[7, ]])
  FR[8] <- mean(dane$bound_45[dane$dec_45 != 0 & selection[8, ]])
  FR[9] <- mean(dane$bound_50[dane$dec_50 != 0 & selection[9, ]])
  FR[10] <- mean(dane$bound_55[dane$dec_55 != 0 & selection[10, ]])
  FR[11] <- mean(dane$bound_60[dane$dec_60 != 0 & selection[11, ]])
  FR[12] <- mean(dane$bound_65[dane$dec_65 != 0 & selection[12, ]])
  
  return(FR[cases == 1])
}

#----------------------------------------------
# FIGURE: RT vs. d-prime

opis <- 1.2
osie <- 1.2

low_bound <- 1000
high_bound <- 3000
non_decision <- 0

data_size <- 15000

lev <- numeric(5)
lev[1] <- 8.247
lev[2] <- 8.4585
lev[3] <- 8.7315
lev[4] <- 9.081
lev[5] <- 9.5295

#selection <- Select_observations(lev[5], low_bound, high_bound)
#sum(selection)

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

selection <- Select_observations(lev[1], low_bound, high_bound)
cases <- Select_levels()
correct <- Calculate_correct(selection, cases)
d_prime <- 2*qnorm(correct,mean = 0)
srednia <- Calculate_RT_mean(selection, cases, non_decision)

plot(srednia, d_prime, type = "l", xlim = c(0, 650), ylim = c(0, 5), 
     xlab = "Mean RT (ms)", ylab = "d'", col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)
for(l in 2:5) {
  selection <- Select_observations(lev[l], low_bound, high_bound)
  cases <- Select_levels()
  correct <- Calculate_correct(selection, cases)
  d_prime <- 2*qnorm(correct,mean = 0)
  srednia <- Calculate_RT_mean(selection, cases, non_decision)

  lines(srednia, d_prime, type = "l", col = kolor[l], lwd = l)
}


legend(0, 5, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")


#-------------------------------------
# FIGURE: RT error - RT correct

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

selection <- Select_observations(lev[1], low_bound, high_bound)
cases <- Select_levels()
RT_error_minus_correct <- Calculate_RT_diff(selection, cases)
RT_error_minus_correct
#all_freq_levels <- c(10,15,20,25,30,35,40,45,50,55,60,65)
#freq_levels <- all_freq_levels[cases == 1] - 1.5
freq_levels <- Calculate_FR_threshold(selection, cases) - 2.5

plot(freq_levels, RT_error_minus_correct, type = "l", ylim = c(-30, 100), xlim = c(20, 50),
     xlab = "Firing Rate Threshold (Hz)", ylab = "RT(E) - RT(C)", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_observations(lev[l], low_bound, high_bound)
  cases <- Select_levels()
  RT_error_minus_correct <- Calculate_RT_diff(selection, cases)
  RT_error_minus_correct
  freq_levels <- Calculate_FR_threshold(selection, cases) - 3 + 0.5*l
  lines(freq_levels, RT_error_minus_correct, type = "l", col = kolor[l], lwd = l)
}


legend(20, 100, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
abline(h = 0, lty = 2)


#---------------------
# t-test

high_bound <- 3000
l <- 5
selection <- Select_observations(lev[l], low_bound, high_bound)

t.test(dane$time_20[dane$dec_20 != 0 & selection[3, ]] ~ dane$dec_20[dane$dec_20 != 0 & selection[3, ]])
t.test(dane$time_25[dane$dec_25 != 0 & selection[4, ]] ~ dane$dec_25[dane$dec_25 != 0 & selection[4, ]])
t.test(dane$time_30[dane$dec_30 != 0 & selection[5, ]] ~ dane$dec_30[dane$dec_30 != 0 & selection[5, ]])
t.test(dane$time_35[dane$dec_35 != 0 & selection[6, ]] ~ dane$dec_35[dane$dec_35 != 0 & selection[6, ]])
t.test(dane$time_40[dane$dec_40 != 0 & selection[7, ]] ~ dane$dec_40[dane$dec_40 != 0 & selection[7, ]])
t.test(dane$time_45[dane$dec_45 != 0 & selection[8, ]] ~ dane$dec_45[dane$dec_45 != 0 & selection[8, ]])
t.test(dane$time_50[dane$dec_50 != 0 & selection[9, ]] ~ dane$dec_50[dane$dec_50 != 0 & selection[9, ]])

# Results of t-test at 0.01 for freq 20-50Hz
# Contrast 1: N, N, N, N, N, Y, Y
# Contrast 2: Y, Y, Y, N, N, N, Y
# Contrast 3: Y, Y, Y, Y, N, N, Y
# Contrast 4: Y, Y, Y, Y, N, N, Y
# Contrast 5: Y, Y, Y, Y, N, N, N



#----------------------------------------------
# FIGURE: sd/mean ratio


opis <- 1.2
osie <- 1.2

low_bound <- 1000
high_bound <- 3500
non_decision <- 0

data_size <- 15000

lev <- numeric(5)
lev[1] <- 8.247
lev[2] <- 8.4585
lev[3] <- 8.7315
lev[4] <- 9.081
lev[5] <- 9.5295

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

selection <- Select_observations(lev[1], low_bound, high_bound)
cases <- Select_levels()
srednia <- Calculate_RT_mean(selection, cases, non_decision)
stdev <- Calculate_RT_sd(selection, cases, non_decision)
stdev_mean_ratio <- stdev/srednia

freq_levels <- Calculate_FR_threshold(selection, cases) - 2.5

plot(freq_levels, stdev_mean_ratio, type = "l", xlim = c(20, 50), ylim = c(0.37, 0.57), 
     xlab = "Firing Rate Threshold (Hz)", ylab = "SD / Mean RT", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_observations(lev[l], low_bound, high_bound)
  cases <- Select_levels()
  srednia <- Calculate_RT_mean(selection, cases, non_decision)
  stdev <- Calculate_RT_sd(selection, cases, non_decision)
  stdev_mean_ratio <- stdev/srednia
  
  freq_levels <- Calculate_FR_threshold(selection, cases) - 3 + 0.5*l
  lines(freq_levels, stdev_mean_ratio, type = "l", col = kolor[l], lwd = l)
}


legend(25, 0.57, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")



#----------------------------------------------
# FIGURE: skewness


opis <- 1.2
osie <- 1.2

low_bound <- 1000
high_bound <- 5000
non_decision <- 0

data_size <- 15000

lev <- numeric(5)
lev[1] <- 8.247
lev[2] <- 8.4585
lev[3] <- 8.7315
lev[4] <- 9.081
lev[5] <- 9.5295

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

selection <- Select_observations(lev[1], low_bound, high_bound)
cases <- Select_levels()
skew <- Calculate_RT_skew(selection, cases, non_decision)

freq_levels <- Calculate_FR_threshold(selection, cases) - 2.5

plot(freq_levels, skew, type = "l", xlim = c(20, 50), ylim = c(0, 3),
     xlab = "Firing Rate Threshold (Hz)", ylab = "RT Skewness", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)

for(l in 2:5) {
  selection <- Select_observations(lev[l], low_bound, high_bound)
  cases <- Select_levels()
  skew <- Calculate_RT_skew(selection, cases, non_decision)

  freq_levels <- Calculate_FR_threshold(selection, cases) - 3 + 0.5*l
  lines(freq_levels, skew, type = "l", col = kolor[l], lwd = l)
}


legend(20, 3, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")

#---------------------------------------
# FIGURE: Mean RT vs Amplitude of Alpha Waves
# Threshold 50Hz
#---------------------------------------


Calculate_sd_RT_vs_amplitude <- function(selection, num_of_bins, errs, width) {
  
  sd_RT <- numeric(num_of_bins)
  num_of_cors <- numeric(num_of_bins)
  num_of_errs <- numeric(num_of_bins)
  
  step <- 0.1/num_of_bins
  for(i in 1:num_of_bins) { 
    interval <- (dane$amplitude > (i-1/2)*step - width) & (dane$amplitude < (i-1/2)*step + width)
    sd_RT[i] <- sd(dane$time_50[selection & interval])
    num_of_cors[i] <- sum(selection & interval & (dane$dec_50 == 1))
    num_of_errs[i] <- sum(selection & interval & (dane$dec_50 == -1))
    if(num_of_errs[i] < errs) {mean_RT[i] <- NaN} 
  }
  return(c(sd_RT))  
}



#--------------------------
# Mean RT
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
high_bound <- 3500

width <- 0.01
errs <- 0
num_of_bins <- 15
l <- 1

selection <- Select_contrast(l, low_bound, high_bound)

amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
sterr <- (Calculate_sd_RT_vs_amplitude(selection, num_of_bins, errs, width))/(sum(selection)^0.5)
mean_RT <- Calculate_mean_RT_vs_amplitude(selection, num_of_bins, errs, width)

sterr
mean_RT

plot(amplitude_levels, mean_RT, type = "l", xlim = c(0, 0.1), ylim = c(400, 1200), 
     xlab = "Amplitude (alpha = 10Hz)", ylab = "Mean RT (ms)", 
     col = kolor[1], lwd = 1, cex.lab = opis, cex.axis = osie)
for(j in 1:num_of_bins) {
  lines(c(amplitude_levels[j], amplitude_levels[j]), c(mean_RT[j] - sterr[j], mean_RT[j] + sterr[j]), lwd = 1, col = kolor[l])
}
points(amplitude_levels, mean_RT - sterr, cex = 1.5, pch = 45, col = kolor[l])
points(amplitude_levels, mean_RT + sterr, cex = 1.5, pch = 45, col = kolor[l])

for(l in 2:5) {
  selection <- Select_contrast(l, low_bound, high_bound)
  
  amplitude_levels <- Calculate_mean_amplitude(selection, num_of_bins, width)
  sterr <- (Calculate_sd_RT_vs_amplitude(selection, num_of_bins, errs, width))/(sum(selection)^0.5)
  mean_RT <- Calculate_mean_RT_vs_amplitude(selection, num_of_bins, errs, width)
  
  lines(amplitude_levels, mean_RT, type = "l", col = kolor[l], lwd = l)
  for(j in 1:num_of_bins) {
    lines(c(amplitude_levels[j], amplitude_levels[j]), c(mean_RT[j] - sterr[j], mean_RT[j] + sterr[j]), lwd = 1, col = kolor[l])
  }
  points(amplitude_levels, mean_RT - sterr, cex = 1.5, pch = 45, col = kolor[l])
  points(amplitude_levels, mean_RT + sterr, cex = 1.5, pch = 45, col = kolor[l])
  
}

legend(0.2, 800, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
#abline(v = 0, lty = 2)

l <- 1
selection <- Select_contrast(l, low_bound, high_bound)
equation <- dane$time_50[selection] ~ dane$amplitude[selection]
fit <- lm(equation)
summary(fit)

