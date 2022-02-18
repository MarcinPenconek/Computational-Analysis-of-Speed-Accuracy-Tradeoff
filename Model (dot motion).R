#-------------------------------------------------------------#
#                  Neural Network Decison Model               #
#                        Marcin Penconek                      #
#-------------------------------------------------------------#


# Different decision function

# model with the Poisson inputs for dot motion experiment
# density: 0.55 i 0.36

# fixed firing 
# 0.07
# 0.006

# decision function based on frequency threshold



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
# Defining stimulus (time-random)

Stimulus_definition_t_random = function(duration) {
  
  stimulus <- matrix(0, nrow = duration, ncol = size)
  
  connect_A <- rep(1, len = to_A)
  connect_B <- rep(1, len = to_B)
  
  r_A <- rpois(duration, lambda = n_A)
  r_B <- rpois(duration, lambda = n_B)
    
  for(k in 1:duration) {
    stimulus[k, ] <- c(connect_A*r_A[k], numeric(n_size - to_A), 
                       connect_B*r_B[k], numeric(n_size - to_B), numeric(size - 2*n_size))
  }
  return(stimulus)
}


#----------------------------------------
# Network evolution function

Network_evolution = function() {
  
  current_time <- 0
  current_stimulus <- numeric(size)
  stimulus_finished <- 0
  
  s <- state[ ,1]
  t <- state[ ,2]
  k <- 1

  simulation <- matrix(0, nrow = number_of_iterations, ncol = (1 + n_groups + 1 + 4))

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
    for(i in 1:size) if(t[i] == current_time) { Min_i <- i }
    
    active_groups <- numeric(n_groups + 1)
    for(i in 1:(n_groups + 1)) {
      active_groups[i] <- sum(s*groups[ ,i + 1])
    }
    
    simulation[j, ] <- c(current_time, active_groups, Min_i, 0, 0, 0)

    active <- sum(s)    
    cutoff <- (active/size)^2/inhibition_level
    
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
    simulation[j, 6] <- s[Min_i]*(Min_i <= 100)
    simulation[j, 7] <- s[Min_i]*((Min_i > 100) & (Min_i <= 200))
    simulation[j, 8] <- s[Min_i]*(Min_i > 200)
  }
  return(simulation)
}

#----------------------------------------------
# Visualization
  
Population_plot <- function(x, y1, y2, y3, dec, dec_time) {
  

  if(dec == 0) {
    dectext <- "" 
  } else { 
    dectext <- round(dec_time, 0)
  }
  
  
  plot(x, y3, ylim = c(0,100), xlim = c(0,6000), type = "l",  
       #xlab = "Time (ms)", ylab = "% active", 
       xlab = "", ylab = "", 
       main = "a                                                                                                                                     ", 
       cex.main = 1.5, cex.axis = 1.2)
  lines(x, y1, type = "l", col = "red", lwd=1.5)
  lines(x, y2, type = "l", col = "blue", lwd=1.5)
  legend(dec_time-200, 105, legend = dectext, cex = 1.2, bty="n")
  if(dec != 0) {abline(v = dec_time)}
  if(n_A + n_B >0) {
    abline(v = event_starts, lwd=1, lty=2)
    abline(v = event_ends, lwd=1, lty=2)
  } else {
      legend(200, 100, legend = c("A", "B"), cex = 1.5, lwd=c(3,3), col=c("red","blue"), bty = "n")
    }
}


#----------------------------------------------
# Analysis of a decision

Freq_Analysis <- function(beg_index, end_index) {
  
  A_count <- sum(simulation[beg_index:end_index, 6])
  B_count <- sum(simulation[beg_index:end_index, 7])
  time_bin <- (simulation[end_index, 1] - simulation[beg_index, 1])/1000
  
  return(c((A_count/100)/time_bin, (B_count/100)/time_bin, time_bin))
}

#------------------------------------
# Alternative decision function for reaction time paradigm

Decision_made_attractor <- function() {
  decision <- 0
  decision_time <- 0
  current_time <- 0
  current_max <- 0
  ongoing_max <- 0
  
  start_i <- 1
  while(simulation[start_i, 1] < event_starts) {
    start_i <- start_i + 1
  }
  
  num_of_bins <- round((number_of_iterations - start_i - updates_per_30ms)/freq_step)
  
  for(k in 1:num_of_bins) {
    
    freq <- Freq_Analysis(start_i + k*freq_step - updates_per_30ms, start_i + k*freq_step)
    current_max <- max(freq[1], freq[2])
    current_time <- simulation[start_i + k*freq_step, 1]
    
    if(decision != 0 & ongoing_max < threshold_level*threshold_time) {
      ongoing_max <- (ongoing_max + current_max*(current_time - previous_time))
      if(ongoing_max < threshold_level*(current_time - decision_time)) {
        decision <- 0
      }
      previous_time <- current_time
    }
    
    if(decision == 0) {
      if(current_max >= threshold_level) {
        decision_time <- current_time
        previous_time <- current_time
        if(current_max == freq[1]) { 
          decision <- 1
        } else if(current_max == freq[2]) {
          decision <- -1
        }
      }
    }
  }
  
  if(ongoing_max >= threshold_level*threshold_time) {
    return(c(decision, decision_time))
  } else {
    return(c(0, 0))
  }  
  
}

#----------------------------------
# Decision function for the reaction-time paradigm
# bound = 52.5Hz

Which_wins_when <- function(bound) {
  
  bound_level <- 0
  bound_time <- 0
  pool_A <- 0
  pool_B <- 0
  which_wins <- 0
  
  b <- start_b
  freq <- Freq_Analysis(b - updates_per_30ms, b)
  
  while(max(freq[1], freq[2]) < bound & b + freq_step < number_of_iterations) {
    b <- b + freq_step
    freq <- Freq_Analysis(b - updates_per_30ms, b)
  }
  
  if(b + freq_step + 1 < number_of_iterations) { 
    bound_level <- max(freq[1], freq[2])
    bound_time <- simulation[b, 1]
    pool_A <- freq[1]
    pool_B <- freq[2]
    which_wins <- sign(freq[1] - freq[2])
  }
  
  return(c(bound_level, bound_time, pool_A, pool_B, which_wins))
}

#------------------------------------
# Decision function for randomly delayed fixed duration paradigm (delay 500 to 1500ms after end of stimulus)

Which_pool_wins <- function() {

  dec_delay <- runif(1, min = 500, max = 1500)
  
  delayed_A <- 0
  delayed_B <- 0
  which_wins <- 0
  
  if(event_ends + dec_delay < simulation[number_of_iterations,1]) {
    
    b <- 14*(event_ends + dec_delay)
    
    if(i < number_of_iterations) {
      while(simulation[b, 1] < (event_ends + dec_delay)) {
        b <- b + 1
      }
    }

    if(i < number_of_iterations) {
      freq <- Freq_Analysis(b - updates_per_30ms, b)
      delayed_A <- freq[1]
      delayed_B <- freq[2]
      which_wins <- sign(freq[1] - freq[2])
    }
  }

  
  return(c((event_ends + dec_delay), delayed_A, delayed_B, which_wins))
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
# Simulation parameters

number_of_iterations <- 45000

threshold_level <- 52.5
threshold_time <- 500
time_bins <- 30

to_A <- 50
to_B <- 50

n_A <- 15
n_B <- 0

#n_A <- 7.74
#n_B <- 15 - n_A

event_starts <- 1000
event_ends <- 3000

fixed_for <- time_bins
duration <- round((event_ends - event_starts)/fixed_for + 1)

start_b <- round(event_starts*14)

#-----------------------------------
# Parameters for Frequency Analysis

updates_per_30ms <- 440
freq_step <- 145


#------------------------------------------------
# SECTION 3: Running the model
#------------------------------------------------


for(i in 1:1) {
  
  groups <- Grouping_variable()
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration)
  simulation <- Network_evolution()
  decision <- Decision_made_attractor()

  Population_plot(simulation[ ,1], simulation[ ,2], simulation[ ,3], simulation[ ,4]*(100/N_size), decision[1], decision[2])

}



#------------------------------------------------
# SECTION 4: Generating a sample of simulations
#------------------------------------------------

# in the reaction-time paradigm stimulus duration is 2000ms (starting at 1000ms); number of iterations = 45000
# in the fixed-duration (working memory) paradigm stimulus duration is 1000ms; change the decision function
# number of runs per coherence level: 500; number of iterations = 55000

delayed_time <- 0
delayed_A <- 0
delayed_B <- 0
delayed_dec <- 0

bound_th <- 0
time_th <- 0
poolA_th <- 0
poolB_th <- 0
dec_th <- 0

over_time_A <- 0
over_time_B <- 0

dane <- data.frame(N_size, n_groups, n_size, n_density, n_to_N_density, N_density, 
                   inhibition_level, r_default, r_active, number_of_iterations, threshold_level, time_bins, 
                   to_A, to_B, n_A, n_B, over_time_A, over_time_B, event_starts, event_ends, 
                   bound_th, time_th, poolA_th, poolB_th, dec_th)

groups <- Grouping_variable()

runs_per_coherence <- 5

for(l in 1:6) {
  
  #dot motion coherence levels: 3.2%, 6.4%, 12.8%, 25.6%, 51.2%
  #n_A levels: 7.5, 7.74, 7.98, 8.46, 9.42, 11.34
  
  if(l == 1 ) {n_A <- 7.5}
  if(l == 2 ) {n_A <- 7.74}
  if(l == 3 ) {n_A <- 7.98}
  if(l == 4 ) {n_A <- 8.46}
  if(l == 5 ) {n_A <- 9.42}
  if(l == 6 ) {n_A <- 11.34}
  
  n_B <- 15 - n_A
  print("Next level started")
  
  for(i in 1:runs_per_coherence) {

    network_connections <- Network_definition()
    state <- Network_init()
    s <- state[ ,1]
    t <- state[, 2]
    
    stimulus <- Stimulus_definition_t_random(duration)
    
    over_time_A <- mean(stimulus[ ,1:100])
    over_time_B <- mean(stimulus[ ,101:200])
    
    simulation <- Network_evolution()
    win_th <- Which_wins_when(threshold_level)
    #decision_att <- Decision_made_attractor()
    #delayed <- Which_pool_wins()
    
    dane[i + (l-1)*runs_per_coherence, ] <- c(N_size, n_groups, n_size, n_density, n_to_N_density, N_density, 
                   inhibition_level, r_default, r_active, number_of_iterations, threshold_level, time_bins,
                   to_A, to_B, n_A, n_B, over_time_A, over_time_B, event_starts, event_ends, 
                   win_th)
    
  }
}

#write.csv2(dane,'RT_dot_motion_4.csv')

View(dane)
getwd() 

#-----------------------------------------------



