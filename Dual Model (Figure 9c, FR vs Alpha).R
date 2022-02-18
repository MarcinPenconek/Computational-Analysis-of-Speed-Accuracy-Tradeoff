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

updates_per_30ms <- 440
freq_step <- 145


#----------------------------------------------
# GRAPHIC PARAMETERS
#----------------------------------------------

osie <- 1.2
opis <- 1.2
grubosc <- 3


#---------------------------------------------
# Baseline Firing Rates vs. Alpha Waves

number_of_iterations <- 32000
time_bin <- 1430 #updates per 100ms
number_of_bins <- floor(number_of_iterations/time_bin)

n_A <- 0
n_B <- 0
event_starts <- 0
event_ends <- 1000

phase <- 0
amplitude <- 0.05
alpha_freq <- 10
period <- (1000/alpha_freq)/(2*3.14)

time <- 0
firing_rate_A <- 0
firing_rate_B <- 0
firing_rate_C <- 0


dane <- data.frame(time, firing_rate_A, firing_rate_B, firing_rate_C, amplitude, phase, alpha_freq)

groups <- Grouping_variable()
next_i <- 0

for(l in 1:6) {
  
  if(l == 1) {amplitude <- 0}
  if(l == 2) {amplitude <- 0.02}
  if(l == 3) {amplitude <- 0.04}
  if(l == 4) {amplitude <- 0.06}
  if(l == 5) {amplitude <- 0.08}
  if(l == 6) {amplitude <- 0.1}

  for(i in 1:200) {
    
    phase <- runif(1, min = 0, max = 2*3.14)
    alpha_freq <- runif(1, min = 8, max = 12)

    network_connections <- Network_definition()
    state <- Network_init()
    s <- state[ ,1]
    t <- state[, 2]
    
    stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
    firing_rate <- Firing_Rate_Evolution()
    
    if(firing_rate[number_of_bins,2] < 30 & firing_rate[number_of_bins,3] < 30) {
      
      for(k in 1:(number_of_bins - 10)) {
        dane[next_i + k, ] <- c(firing_rate[9 + k,1:4], amplitude, phase, alpha_freq)
      }
      next_i <- next_i + number_of_bins - 10
    }
  }
  #write.csv2(dane,'Dist_FR_Spon_vs_Random_Alpha_1.csv')
  
}


#-----------------------------
# Figure 9c: BASELINE FIRING RATE vs. ALPHA
#-----------------------------

opis <- 1.2
osie <- 1.2

dane1 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Alpha/Dist_FR_Spon_vs_Random_Alpha_1.csv")
dane2 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Alpha/Dist_FR_Spon_vs_Random_Alpha_2.csv")
dane <- merge(dane1, dane2, all.x = TRUE, all.y = TRUE)

data_size <- length(dane$amplitude)
data_size

plot(dane$amplitude + rnorm(data_size, mean = 0, sd = 0.002), 
     dane$firing_rate_A, type = "p", xlim = c(-0.005, 0.105), ylim = c(0, 20), pch = 16, cex = 0.2, 
     xlab = "Amplitude (alpha: 8-12Hz)", ylab = "Baseline Firing Rate (Hz)", 
     col = "light grey", cex.lab = opis, cex.axis = osie)
points(dane$amplitude + rnorm(data_size, mean = 0, sd = 0.002), dane$firing_rate_B,
       pch = 16, cex = 0.2, col = "light grey")


mean_AB <- numeric(6)
for(l in 1:6) {
  mean_AB[l] <- 1/2*(mean(dane$firing_rate_A[dane$amplitude == 0.02*(l - 1)]) + mean(dane$firing_rate_B[dane$amplitude == 0.02*(l - 1)]))
}

sd_AB <- numeric(6)
for(l in 1:6) {
  sd_AB[l] <- 1/2*(sd(dane$firing_rate_A[dane$amplitude == 0.02*(l - 1)]) + sd(dane$firing_rate_B[dane$amplitude == 0.02*(l - 1)]))
}

mean_AB
sd_AB

lines(c(0, 0.02, 0.04, 0.06, 0.08, 0.1), mean_AB, type = "l", col = "black")

#plot(c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), mean_AB, type = "l", 
#     xlim = c(0, 0.07), ylim = c(9, 11), pch = 16, cex = 0.2, 
#     xlab = "Amplitude (alpha: 10Hz)", ylab = "Baseline Firing Rate (Hz)", 
#     cex.lab = opis, cex.axis = osie)

points(c(0, 0.02, 0.04, 0.06, 0.08, 0.1), mean_AB, pch = 19, cex = 1.5, col = "black")
#abline(h = 10, lty = 2)

#mean_C <- numeric(6)
#for(l in 1:6) {
#  mean_C[l] <- mean(dane$firing_rate_C[dane$amplitude == 0.02*(l - 1)])
#}

#sd_C <- numeric(8)
#for(l in 1:8) {
#  sd_C[l] <- sd(dane$firing_rate_C[dane$amplitude == 0.01*(l - 1)])
#}

#mean_C
#lines(c(0, 0.02, 0.04, 0.06, 0.08, 0.1), mean_C, type = "l", col = "black")
#points(c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07), mean_C, pch = 18, cex = 1.5, col = "black")

for(l in 1:6) {
  lines(c(0.02*(l - 1), 0.02*(l - 1)), c(mean_AB[l] - sd_AB[l], mean_AB[l] + sd_AB[l]), lwd = 1)
}
points(c(0, 0.02, 0.04, 0.06, 0.08, 0.1), mean_AB - sd_AB, cex = 1.5, pch = 45, col = "black")
points(c(0, 0.02, 0.04, 0.06, 0.08, 0.1), mean_AB + sd_AB, cex = 1.5, pch = 45, col = "black")

legend(0.06, 20, pch = c(19, NaN), lwd = c(1, 1), c("Mean", "SD"), cex = opis, bty = "n")
#abline(h = 10, lty = 2)

#-----------------------------
# Statistics

mean_AB
sd_AB

baseline <- c(dane$firing_rate_A, dane$firing_rate_B)
ampli_v <- c(dane$amplitude, dane$amplitude)

equation <- baseline ~ ampli_v
fit <- lm(equation)
summary(fit)
anova(fit)


