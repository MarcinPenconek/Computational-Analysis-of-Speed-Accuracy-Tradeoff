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



Freq_1000 <- function() {
  
  i <- 1
  while((i < number_of_bins - 1) & (firing_rate[i, 1] < event_starts)) {
    i <- i + 1
  }
  return(c(firing_rate[i, 2], firing_rate[i, 3]))
}

Which_won <- function(bound) {

  dec <- 0
  time_stop <- 0
  FR_stop <- 0
  
  i <- 1
  while((i < number_of_bins - 1) & firing_rate[i, 2] < bound) {
    i <- i + 1
  }
  if(i < number_of_bins - 1) {
    dec <- 1
    time_stop <- firing_rate[i, 1]
    FR_stop <- firing_rate[i, 2]
  }
  return(c(dec, time_stop, FR_stop))
  
}

#----------------------------------------------
# GRAPHIC PARAMETERS
#----------------------------------------------

osie <- 1.2
opis <- 1.2
grubosc <- 3


#----------------------------------------
# Ramping with alpha waves
# and random contrast/n_A levels: from 8.247 to 9.5295
# GENERATING DATABASE

#n_A_lev <- numeric(5)
#n_A_lev[1] <- 8.247
#n_A_lev[2] <- 8.4585
#n_A_lev[3] <- 8.7315
#n_A_lev[4] <- 9.081
#n_A_lev[5] <- 9.5295

number_of_iterations <- 65000
time_bin <- 440
number_of_bins <- floor(number_of_iterations/time_bin)

event_starts <- 1000
event_ends <- 3500
duration <- round((event_ends - event_starts)/fixed_for + 1)

FRA <- 0
FRB <- 0

dec_50 <- 0
time_50 <- 0
WFRA <- 0

dane <- data.frame(n_A, FRA, FRB, dec_50, time_50, WFRA, amplitude, phase, alpha_freq)

groups <- Grouping_variable()

runs_per_amplitude <- 500

for(j in 1:runs_per_amplitude) {
  n_A <- runif(1, min = 8.247, max = 9.5295)
  n_B <- 15 - n_A
  
  set_random <- round(runif(1, min = 1, max = 500))
  phase <- runif(1, min = 0, max = 2*3.14)
  alpha_freq <- runif(1, min = 8, max = 12)
  
  for(k in 1:6) {
    
    amplitude <- 0.02*(k-1)
    set.seed(set_random)
    
    network_connections <- Network_definition()
    state <- Network_init()
    s <- state[ ,1]
    t <- state[, 2]
    
    stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
    firing_rate <- Firing_Rate_Evolution()

    FFR <- Freq_1000()
    FRA <- FFR[1]
    FRB <- FFR[2]
    won <- Which_won(50)
    dec_50 <- won[1]
    time_50 <- won[2]
    WFRA <- won[3]
    
    dane[6*(j-1) + k, ] <- c(n_A, FRA, FRB, dec_50, time_50, WFRA, amplitude, phase, alpha_freq)
  }
  
}


#write.csv2(dane,'Ramping_with_Random_Alpha.csv')

#----------------------------------
# FIGURE 9d: Ramping
#----------------------------------

opis <- 1.2
osie <- 1.2

dane <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Alpha/Ramping_with_Random_Alpha.csv")

#dane <- merge(dane1, dane2, all.x = TRUE, all.y = TRUE)

#input <- (dane$n_A > 0)

selection <- dane$dec_50 == 1
sum(selection)
slope <- 100*(dane$WFRA - dane$FRA)/(dane$time_50 - 1000)

data_size <- sum(selection)
data_size

plot(dane$amplitude[selection] + rnorm(data_size, 0, 0.002), slope[selection], 
     type = "p", xlim = c(-0.005, 0.105), ylim = c(0,15), pch = 16, cex = 0.5, 
     xlab = "Amplitude (alpha = 8-12Hz)", ylab = "Integration Speed (Hz/100ms)", 
     col = "light grey", cex.lab = opis, cex.axis = osie)

lines(c(0, 0.02, 0.04, 0.06, 0.08), 
      c(mean(slope[selection & dane$amplitude == 0]), 
        mean(slope[selection & dane$amplitude == 0.02]),
        mean(slope[selection & dane$amplitude == 0.04]),
        mean(slope[selection & dane$amplitude == 0.06]),
        mean(slope[selection & dane$amplitude == 0.08])),
      type = "l", col = "black")

points(c(0, 0.02, 0.04, 0.06, 0.08), 
      c(mean(slope[selection & dane$amplitude == 0]), 
        mean(slope[selection & dane$amplitude == 0.02]),
        mean(slope[selection & dane$amplitude == 0.04]),
        mean(slope[selection & dane$amplitude == 0.06]),
        mean(slope[selection & dane$amplitude == 0.08])),
       pch = 19, cex = 1.5, col = "black")

for(k in 1:5) {
  lines(c(0.02*(k-1), 0.02*(k-1)), c(mean(slope[selection & dane$amplitude == 0.02*(k-1)]) - sd(slope[selection & dane$amplitude == 0.02*(k-1)]), 
                   mean(slope[selection & dane$amplitude == 0.02*(k-1)]) + sd(slope[selection & dane$amplitude == 0.02*(k-1)])), lwd = 1)
                   
  points(0.02*(k-1), mean(slope[selection & dane$amplitude == 0.02*(k-1)]) - sd(slope[selection & dane$amplitude == 0.02*(k-1)]),
         cex = 1.5, pch = 45, col = "black")  
  points(0.02*(k-1), mean(slope[selection & dane$amplitude == 0.02*(k-1)]) + sd(slope[selection & dane$amplitude == 0.02*(k-1)]),
         cex = 1.5, pch = 45, col = "black")  
}

legend(0.06, 15, pch = c(19, NaN), lwd = c(1, 1), c("Mean", "SD"), cex = opis, bty = "n")


#-----------------------
# statistics

sum(selection & dane$amplitude == 0)
sum(selection & dane$amplitude == 0.02)
sum(selection & dane$amplitude == 0.04)
sum(selection & dane$amplitude == 0.06)
sum(selection & dane$amplitude == 0.08)
sum(selection & dane$amplitude == 0.1)

k <- 5
mean(slope[selection & dane$amplitude == 0.02*(k-1)])
sd(slope[selection & dane$amplitude == 0.02*(k-1)])


equation <- slope[selection] ~ dane$amplitude[selection]
fit <- lm(equation)
summary(fit)
anova(fit)


plot(fitted.values(fit), residuals(fit))
abline(h = 0)

