#-------------------------------------------------------------#
#                  Neural Network Decison Model               #
#                        Marcin Penconek                      #
#-------------------------------------------------------------#

# Dual Model reporting population-average firing rate

# model with Poisson inputs 
# speed-accuracy trade-off

# density: 0.55 i 0.36

# fixed firing 
# 0.07
# 0.006



#----------------------------------------------
# SECTION 1: Functions
#----------------------------------------------

Grouping_variable = function() {
  
  size <- n_groups*n_size + N_size
  ggroups <- c(rep(1:n_groups, each = n_size, len = n_groups*n_size), rep(n_groups + 1, len = N_size))
  groups <- matrix(0, nrow = size, ncol = n_groups + 2)
  groups[ ,1] <- ggroups
  for(i in 1:size) {
    groups[i, ggroups[i] + 1] <- 1
  }
  return(groups)
}


#----------------------------------------------
# Defining network connections


Network_definition = function() {
  
  network_connections <- matrix(FALSE, nrow = size, ncol = size)
  u <- runif(size*size, min = 0, max = 1)
  
  for(i in 1:size) {
    for(j in 1:size) {
      if((groups[i, 1] == n_groups + 1) & (groups[j, 1] == n_groups + 1)) {
        if(u[size*(i-1) + j] < N_density) {
          network_connections[i, j] <- TRUE
        }
      } else if(groups[i, 1] == groups[j, 1]) {
        if(u[size*(i-1) + j] < n_density) {
          network_connections[i, j] <- TRUE
        }
      } else {
        if(u[size*(i-1) + j] < n_to_N_density) {
          network_connections[i, j] <- TRUE
        }
      }
    }
  }
  return(network_connections)
}


#----------------------------------------------
# Initiating the network with random inputs

Network_init = function() {
  
  prob_s <- runif(size, min = 0, max = 1)
  t <- rexp(size, rate = r_default)
  s <- logical(size)
  
  for(i in 1:size) {
    if(prob_s[i] < inhibition_level) {
      s[i] <- TRUE
      t[i] <- rexp(1, rate = r_active)
    }
    else {
      s[i] <- FALSE
    }
  }
  state <- matrix(0, nrow = size, ncol = 2)
  state[ ,1] <- s
  state[ ,2] <- t
  
  return(state)
}


#-----------------------------------------------
# Defining stimulus (non-random) - n_A, n_B can be non-integer

Stimulus_definition_non_random = function() {
  
  stimulus <- c(rep(n_A, len = to_A), numeric(n_size - to_A), 
                rep(n_B, len = to_B), numeric(n_size - to_B), numeric(size - 2*n_size))
  
  return(stimulus)
}


#-----------------------------------------------
# Defining stimulus (time-random) - VERSION #4 (n_A, n_B can be non-integer)

Stimulus_definition_t_random = function(duration, input_A, input_B) {
  
  stimulus <- matrix(0, nrow = duration, ncol = size)
  
  connect_A <- rep(1, len = to_A)
  connect_B <- rep(1, len = to_B)
  
  r_A <- rpois(duration, lambda = input_A)
  r_B <- rpois(duration, lambda = input_B)
  
  for(k in 1:duration) {
    stimulus[k, ] <- c(connect_A*r_A[k], numeric(n_size - to_A), 
                       connect_B*r_B[k], numeric(n_size - to_B), numeric(size - 2*n_size))
  }
  return(stimulus)
}

#--------------------------
# Firing rate analysis based on time bins (roughly 30 ms)

Firing_Rate_Analysis <- function(sim, beg_index, end_index) {
  
  A_count <- sum(sim[beg_index:end_index, 6])
  B_count <- sum(sim[beg_index:end_index, 7])
  C_count <- sum(sim[beg_index:end_index, 8])
  time_bin <- (sim[end_index, 1] - sim[beg_index, 1])/1000
  
  return(c((A_count/n_size)/time_bin, (B_count/n_size)/time_bin, (C_count/N_size)/time_bin, time_bin))
}



#----------------------------------------
# Network evolution function

Firing_Rate_Evolution = function() {
  
  current_time <- 0
  current_stimulus <- numeric(size)
  stimulus_finished <- 0
  active_groups <- numeric(n_groups + 1)

  s <- state[ ,1]
  t <- state[ ,2]
  k <- 1

  simulation <- matrix(0, nrow = number_of_iterations, ncol = (1 + n_groups + 1 + 4))
  firing_rate <- matrix(0, nrow = number_of_bins, ncol = (1 + n_groups + 2))

  for(j in 1:number_of_iterations) {
    
    if(current_time > event_starts & current_time < event_ends) {
      current_stimulus <- stimulus[k,]
      k <- round((current_time - event_starts)/fixed_for) + 1
    }
    
    if(stimulus_finished == 0 & current_time > event_ends) {
      current_stimulus <- numeric(size)
      stimulus_finished <- 1
    }

    # Defining the neuron with the lowest t
    # Changing its state based on the activity of other neurons it is connected to
    
    current_time <- min(t)
    current_inhibition <- inhibition_level + amplitude*sin(current_time/period + phase)
    
    for(i in 1:size) if(t[i] == current_time) { Min_i <- i }
    
    for(i in 1:(n_groups + 1)) {
      active_groups[i] <- sum(s*groups[ ,i + 1])
    }
    
    simulation[j, ] <- c(current_time, active_groups, Min_i, 0, 0, 0)

    active <- sum(s)    
    cutoff <- (active/size)^2/current_inhibition
    
    active_connections <- sum(s*network_connections[Min_i, ])
    all_connections <- sum(network_connections[Min_i, ])
    
    change_factor <- (current_stimulus[Min_i] + active_connections)/(current_stimulus[Min_i] + all_connections)
    
    if(change_factor > cutoff) {
      s[Min_i] <- 1
      t[Min_i] <- current_time + 1/r_active
    } else {
      s[Min_i] <- 0
      t[Min_i] <- current_time + rexp(1, rate = r_default)
    }
    
    simulation[j, 6] <- s[Min_i]*(Min_i <= n_size)
    simulation[j, 7] <- s[Min_i]*(Min_i > n_size & Min_i <= 2*n_size)
    simulation[j, 8] <- s[Min_i]*(Min_i > 2*n_size)

    m <- round(j/time_bin)
    if(j/time_bin == m) {
      firing_rate[m, ] <- c(current_time, Firing_Rate_Analysis(simulation, j - time_bin + 1, j))
    }
      
  }
  return(firing_rate)
}




#----------------------------------------------
# SECTION 2.PARAMETERS
#----------------------------------------------

# Network Parameters

N_size <- 800
n_groups <- 2
n_size <- 100
size <- n_groups*n_size + N_size

n_density <- 0.55
n_to_N_density <- 0.36
N_density <- 0.36


#----------------------------------------------
# Global parameters

inhibition_level <- 0.13
r_default <- 0.006
r_active <- 0.07

#----------------------------------------------
# Alpha waves

phase <- 0
amplitude <- 0.05
alpha_freq <- 10
period <- (1000/alpha_freq)/(2*3.14)


#----------------------------------------------
# Simulation parameters

number_of_iterations <- 50000
# time_bin = #iterations in 30ms
time_bin <- 440
number_of_bins <- floor(number_of_iterations/time_bin)

threshold_level <- 75
threshold_time <- 500

to_A <- 50
to_B <- 50
fixed_for <- 30

#contrast levels: 4.98%, 6.39%, 8.21%, 10.54%, 13.53%
#n_A levels: 7.8735, 7.97925, 8.11575, 8.2905, 8.51475
#modified (doubled contrasts) n_A levels: 8.247, 8.4585, 8.7315, 9.081, 9.5295

#dot motion coherence levels: 3.2%, 6.4%, 12.8%, 25.6%, 51.2%
#n_A levels: 7.74, 7.98, 8.46, 9.42, 11.34

n_A <- 8.46
n_B <- 15 - n_A

event_starts <- 1000
event_ends <- 2000

duration <- round((event_ends - event_starts)/fixed_for + 1)

#-----------------------------------
# Parameters for Frequency Analysis

#updates_per_50ms <- 1730
#freq_step <- 420

updates_per_30ms <- 440
freq_step <- 145


#----------------------------------------------
# GRAPHIC PARAMETERS
#----------------------------------------------

osie <- 1.2
opis <- 1.2
grubosc <- 3

#-----------------------------------------------
# FIGURE 9a: Decision formation
#-----------------------------------------------

Firing_rate_plot <- function(x, y1, y2, y3, dec_1, dec_2) {
  
  #  if(decision[1] == 1) {
  #    decision_text <- c("A won: ", decision[2])
  #  } else if(decision[1] == -1) {
  #    decision_text <- c("B won: ", decision[2])
  #  } else {
  #    decision_text <- c("", "")
  #  }
  
  #  if(dec_1 == 0) {
  #    dectext <- "" 
  #  } else { 
  #    dectext <- round(dec_2, 0)
  #  }
  
  
  plot(x, y1, ylim = c(0,70), xlim = c(0,2000), type = "l",  
       xlab = "Time (ms)", ylab = "Firing Rate (Hz)", 
       #xlab = "", ylab = "", 
       cex.main = 1.5, cex.axis = osie, cex.lab = opis, col = "pink", lwd = grubosc)
  lines(x, y2, type = "l", col = "light blue", lwd=grubosc)
  
  #  legend(dec_2-200, 70, legend = dectext, cex = opis, bty="n")
  #  if(dec_1 != 0) {abline(v = dec_2)}
  #  if(n_A + n_B >0) {
  #    abline(v = event_starts - 1000, lwd=1, lty=2)
  #    abline(v = event_ends - 1000, lwd=1, lty=2)
  #  } else {
  #    legend(-1000, 70, legend = c("A", "B"), cex = opis, lwd=c(3,3), col=c("pink","light blue"), bty = "n")
  #  }
}

#dot motion coherence levels: 3.2%, 6.4%, 12.8%, 25.6%, 51.2%
#n_A levels: 7.74, 7.98, 8.46, 9.42, 11.34

n_A <- 11.34
n_B <- 15 - n_A
dec_thres <- 50

kolor1 <- "pink"
kolor2 <- "light blue"
grubosc <- 2

set.seed(1234)

for(i in 1:10) {
 
  if(i == 10) { 
    kolor1 = "red"
    kolor2 = "blue"
    grubosc = 3}

  groups <- Grouping_variable()
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  firing_rate <- Firing_Rate_Evolution()

#  threshold_reached <- 0
#  for(j in 1:number_of_bins) { 
#    if(threshold_reached > 0) { firing_rate[j,1] <- NA }
#    if(firing_rate[j,2] >= dec_thres & threshold_reached == 0) { 
#      threshold_reached <- firing_rate[j, 1]
#      firing_rate[j,2] <- dec_thres}
#  }
  if(i == 1) {
    plot(firing_rate[ , 1] - 1000, firing_rate[ , 2], ylim = c(0,70), xlim = c(0,1000), 
         type = "l", lwd=grubosc,   
         xlab = "Time from stimulus onset (ms)", ylab = "Firing Rate (Hz)", 
         #xlab = "", ylab = "", 
         cex.axis = osie, cex.lab = opis, col = kolor1)
    lines(firing_rate[ , 1] - 1000, firing_rate[ , 3], type = "l", lwd=grubosc, col = kolor2)
    abline(h = dec_thres, lwd=1, lty=2)
  } else { 
    lines(firing_rate[ , 1] - 1000, firing_rate[ , 2], type ="l", lwd=grubosc, col = kolor1) 
    lines(firing_rate[ , 1] - 1000, firing_rate[ , 3], type = "l", lwd=grubosc, col = kolor2)
    }

}

legend(0, 70, legend = c("A", "B"), cex = opis, lwd=c(3,3), col=c("red","blue"), bty = "n")
legend(550, 40, legend = "Alpha = 10Hz", cex = opis, bty = "n")

#---------------------------------
# not used


#-----------------------------------------------
# FIGURE 2: Phase space
#-----------------------------------------------


n_A <- 11.34
n_B <- 15 - n_A

kolor1 = "grey"
kolor2 = "grey"
kolor3 = "grey"
grubosc = 2

set.seed(1234)

for(i in 1:20) {
  
  if(i == 20) { 
    kolor1 = "black"
    kolor2 = "black"
    kolor3 = "black"
    grubosc = 3}
  
  dec_time <- 0
  
  groups <- Grouping_variable()
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  firing_rate <- Firing_Rate_Evolution()
  
  if(i == 1) {
    plot(firing_rate[ , 2], firing_rate[ , 3], 
         ylim = c(0,70), xlim = c(0,70),
         type = "l", lwd= grubosc,   
         xlab = "Firing Rate in A (Hz)", ylab = "Firing Rate in B (Hz)", 
         #xlab = "", ylab = "", 
         cex.axis = osie, cex.lab = opis, col = kolor1)
    lines(firing_rate[ , 2][firing_rate[,1] >= 1000], firing_rate[ , 3][firing_rate[,1] >= 1000],
          type = "l", lwd=grubosc, col = kolor2)
    lines(firing_rate[ , 2][firing_rate[,1] >= dec_time], firing_rate[ , 3][firing_rate[,1] >= dec_time],
          type = "l", lwd=grubosc, col = kolor3)
  } else { 
    lines(firing_rate[ , 2], firing_rate[ , 3], type ="l", lwd=grubosc, col = kolor1) 
    lines(firing_rate[ , 2][firing_rate[,1] >= 1000], firing_rate[ , 3][firing_rate[,1] >= 1000],
          type = "l", lwd=grubosc, col = kolor2)
    lines(firing_rate[ , 2][firing_rate[,1] >= dec_time], firing_rate[ , 3][firing_rate[,1] >= dec_time],
          type = "l", lwd=grubosc, col = kolor3)
  }
  legend(0, 70, legend = c("Decision Trajectory", "20 Trajectories"), cex = opis, lwd=c(3,3), col=c("black","grey"), bty = "n")
  
}

#------------------------------------------------
# FIGURE 3: Parameters controlling RT
#------------------------------------------------


Figure_3_plot <- function(x, y1, y2, y3, decision) {
  
  #  if(decision[1] == 1) {
  #    decision_text <- c("A won: ", decision[2])
  #  } else if(decision[1] == -1) {
  #    decision_text <- c("B won: ", decision[2])
  #  } else {
  #    decision_text <- c("", "")
  #  }
  
  if(decision[1] == 0) {
    dectext <- "" 
  } else { 
    dectext <- round(decision[2], 0)
  }
  
  
  plot(x, y3, ylim = c(0,70), xlim = c(-1000,2000), type = "l",  
       xlab = "Time (ms)", ylab = "", 
       #xlab = "", ylab = "", 
       cex.main = 1.5, cex.axis = osie, cex.lab = opis, col = "grey", lwd = 2)
  lines(x, y1, type = "l", col = "red", lwd=grubosc)
  lines(x, y2, type = "l", col = "blue", lwd=grubosc)
  #legend(decision[2]-200, 70, legend = dectext, cex = 1.2, bty="n")
  if(decision[1] != 0) {abline(v = decision[2])}
  #if(n_A + n_B >0) {
  #  abline(v = event_starts, lwd=1, lty=2)
  #  abline(v = event_ends, lwd=1, lty=2)
  #} else {
  #legend(-1000, 70, legend = c("A", "B"), cex = opis, lwd=c(3,3), col=c("red","blue"), bty = "n")
  #}
}

#-------------------------------------
# Model

set.seed(1234)

for(i in 1:1) {
  
  groups <- Grouping_variable()
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  firing_rate <- Firing_Rate_Evolution()

  Figure_3_plot(firing_rate[ ,1] - 1000, firing_rate[ ,2], firing_rate[ ,3], firing_rate[ ,4], c(0,0))
  
}


#----------------------------------------
# Analysis of the Temporal Distribution of Firing Rates
# Firing Rates when in the decision state

number_of_iterations <- 32000
number_of_bins <- floor(number_of_iterations/time_bin)

n_A <- 0
n_B <- 15
event_starts <- 0
event_ends <- 1000

time <- 0
firing_rate_A <- 0
firing_rate_B <- 0
firing_rate_C <- 0

dane <- data.frame(time, firing_rate_A, firing_rate_B, firing_rate_C)

groups <- Grouping_variable()
next_i <- 0

for(i in 1:500) {
  
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  firing_rate <- Firing_Rate_Evolution()
    
  if(firing_rate[35,3] >= 50 & firing_rate[number_of_bins,3] >= 50) {
    print("TAK")
    for(k in 1:(number_of_bins - 35)) {
      dane[next_i + k, ] <- firing_rate[34 + k,1:4]
    }
    next_i <- next_i + number_of_bins - 35
    }
  
  write.csv2(dane,'Dist_FR_B.csv')
  
}



#---------------------------------------------
# Firing Rates when in the spontaneous state

number_of_iterations <- 32000
number_of_bins <- floor(number_of_iterations/time_bin)

n_A <- 0
n_B <- 0
event_starts <- 0
event_ends <- 1000

time <- 0
firing_rate_A <- 0
firing_rate_B <- 0
firing_rate_C <- 0

dane <- data.frame(time, firing_rate_A, firing_rate_B, firing_rate_C)

groups <- Grouping_variable()
next_i <- 0

for(i in 1:500) {
  
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  firing_rate <- Firing_Rate_Evolution()
  
  if(firing_rate[number_of_bins,2] < 30 & firing_rate[number_of_bins,3] < 30) {
    
    for(k in 1:(number_of_bins - 35)) {
      dane[next_i + k, ] <- firing_rate[34 + k,1:4]
    }
    next_i <- next_i + number_of_bins - 35
  }
  
  write.csv2(dane,'Dist_FR_Spontaneous.csv')
  
}

#----------------------------
# ATTRACTOR SETS
#----------------------------

dane1 <- read.csv2("~/Studia/NeuroNetwork Model II/SAT paper/Dist_FR_A.csv", row.names=1)
dane2 <- read.csv2("~/Studia/NeuroNetwork Model II/SAT paper/Dist_FR_B.csv", row.names=1)
dane3 <- read.csv2("~/Studia/NeuroNetwork Model II/SAT paper/Dist_FR_Spontaneous.csv", row.names=1)

kolor1 <- "pink"
kolor2 <- "light blue"
kolor3 <- "grey"

plot(dane3[1:10000, ]$firing_rate_A, dane3[1:10000, ]$firing_rate_B, 
     ylim = c(0,70), xlim = c(0,70),
     type = "p",  
     xlab = "Firing Rate in A (Hz)", ylab = "Firing Rate in B (Hz)", 
     #xlab = "", ylab = "", 
     pch = 16, cex = 0.3, cex.axis = osie, cex.lab = opis, col = kolor3)

points(dane1[1:10000, ]$firing_rate_A, dane1[1:10000, ]$firing_rate_B, 
      pch = 16, cex = 0.3, col = kolor1)

points(dane2[1:10000, ]$firing_rate_A, dane2[1:10000, ]$firing_rate_B, 
      pch = 16, cex = 0.3, col = kolor2)

kolor1 <- "red"
kolor2 <- "blue"
kolor3 <- "black"

points(dane3[1:500, ]$firing_rate_A, dane3[1:500, ]$firing_rate_B, 
       pch = 16, cex = 0.3, col = kolor3)

points(dane1[1:500, ]$firing_rate_A, dane1[1:500, ]$firing_rate_B, 
       pch = 16, cex = 0.3, col = kolor1)

points(dane2[1:500, ]$firing_rate_A, dane2[1:500, ]$firing_rate_B, 
       pch = 16, cex = 0.3, col = kolor2)

legend(25, 70, legend = c("Decision State A", "Decision State B", "Spontaneous State"), cex = opis, pch = 16, col=c("red", "blue", "black"), 
       bty = "n")


# N = 43697 data points

summary(dane3$firing_rate_A)
quantile(dane3$firing_rate_A, probs = seq(0, 1, 0.05))
summary(dane3$firing_rate_B)
quantile(dane3$firing_rate_B, probs = seq(0, 1, 0.05))
(sd(dane3$firing_rate_A)+sd(dane3$firing_rate_B))/2

# Firing Rate in the spontaneous state: 4-20Hz (mean: 10.3, sd = 5.2)

summary(dane1$firing_rate_A)
quantile(dane1$firing_rate_A, probs = seq(0, 1, 0.05))
summary(dane2$firing_rate_B)
quantile(dane2$firing_rate_B, probs = seq(0, 1, 0.05))
sd(dane1$firing_rate_A)
sd(dane2$firing_rate_B)
(sd(dane1$firing_rate_A)+sd(dane2$firing_rate_B))/2

# Firing Rate in the winning attractor: 50-65Hz (mean: 59.7, sd = 4.7)

summary(dane2$firing_rate_A)
quantile(dane2$firing_rate_A, probs = seq(0, 1, 0.05))
summary(dane1$firing_rate_B)
quantile(dane1$firing_rate_B, probs = seq(0, 1, 0.05))
sd(dane2$firing_rate_A)
sd(dane1$firing_rate_B)
(sd(dane2$firing_rate_A)+sd(dane1$firing_rate_B))/2

# Firing Rate in the dominated pool: 0-5Hz (mean: 2, sd = 1.4)
