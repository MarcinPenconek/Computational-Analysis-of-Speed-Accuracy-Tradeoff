#-------------------------------------------------------------#
#                  Neural Network Decison Model               #
#                        Marcin Penconek                      #
#-------------------------------------------------------------#


# Model for SAT Analysis with Firing Rate Threshold

# model with random inputs speed-accuracy trade-off
# density: 0.55 i 0.36
# iterations: 100000

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




#----------------------------------------
# Network evolution function

Network_evolution = function() {
  
  current_time <- 0
  current_stimulus <- numeric(size)
  stimulus_finished <- 0
  active_groups <- numeric(n_groups + 1)

  s <- state[ ,1]
  t <- state[ ,2]
  k <- 1

  simulation <- matrix(0, nrow = number_of_iterations, ncol = (1 + n_groups + 1 + 3))

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
    
    for(i in 1:(n_groups + 1)) {
      active_groups[i] <- sum(s*groups[ ,i + 1])
    }
    
    simulation[j, ] <- c(current_time, active_groups, Min_i, 0, 0)

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
    simulation[j, 7] <- s[Min_i]*(Min_i > 100 & Min_i <= 200)
    
  }
  return(simulation)
}

#----------------------------------------------
# Visualization
  
Activation_plot <- function(x, y1, y2, y3, decision) {
  
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
  
  
  plot(x, y3, ylim = c(0,100), xlim = c(0,6000), type = "l",  
       #xlab = "Time (ms)", ylab = "% active", 
       xlab = "", ylab = "", 
       main = "a                                                                                                                                     ", 
       cex.main = 1.5, cex.axis = 1.2)
  lines(x, y1, type = "l", col = "red", lwd=1.5)
  lines(x, y2, type = "l", col = "blue", lwd=1.5)
  legend(decision[2]-200, 105, legend = dectext, cex = 1.2, bty="n")
  if(decision[1] != 0) {abline(v = decision[2])}
  if(n_A + n_B >0) {
    abline(v = event_starts, lwd=1, lty=2)
    abline(v = event_ends, lwd=1, lty=2)
  } else {
      legend(200, 100, legend = c("A", "B"), cex = 1.5, lwd=c(3,3), col=c("red","blue"), bty = "n")
    }
}



Freq_Analysis <- function(beg_index, end_index) {
  
  A_count <- sum(simulation[beg_index:end_index, 6])
  B_count <- sum(simulation[beg_index:end_index, 7])
  time_bin <- (simulation[end_index, 1] - simulation[beg_index, 1])/1000
  
  return(c((A_count/100)/time_bin, (B_count/100)/time_bin, time_bin))
}


Which_wins_when <- function(bound) {
  
  bound_level <- 0
  bound_time <- 0
  pool_A <- 0
  pool_B <- 0
  which_wins <- 0
  
  i <- 1
  while(simulation[i, 1] < event_starts + 30) {
    i <- i + 1
  }
  freq <- Freq_Analysis(i - updates_per_30ms, i)
  
  while(max(freq[1], freq[2]) < bound & i + freq_step < number_of_iterations) {
    i <- i + freq_step
    freq <- Freq_Analysis(i - updates_per_30ms, i)
  }
  
  if(i + freq_step + 1 < number_of_iterations) { 
    bound_level <- max(freq[1], freq[2])
    bound_time <- simulation[i, 1]
    pool_A <- freq[1]
    pool_B <- freq[2]
    which_wins <- sign(freq[1] - freq[2])
  }
  
  return(c(bound_level, bound_time, pool_A, pool_B, which_wins))
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

number_of_iterations <- 55000

to_A <- 50
to_B <- 50
fixed_for <- 30

#contrast levels: 4.98%, 6.39%, 8.21%, 10.54%, 13.53%
#n_A levels: 7.8735, 7.97925, 8.11575, 8.2905, 8.51475
#modified (doubled contrasts) n_A levels: 8.247, 8.4585, 8.7315, 9.081, 9.5295


n_A <- 7
n_B <- 12 - n_A

event_starts <- 1000
event_ends <- 3500

duration <- round((event_ends - event_starts)/fixed_for + 1)

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
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  simulation <- Network_evolution()

  Activation_plot(simulation[ ,1], simulation[ ,2], simulation[ ,3], simulation[ ,4]*(100/N_size), c(0,0))

}



#------------------------------------------------
# SECTION 4: Generating a sample of simulations
#------------------------------------------------


over_time_A <- 0
over_time_B <- 0

bound_10 <- 0
bound_15 <- 0
bound_20 <- 0
bound_25 <- 0
bound_30 <- 0
bound_35 <- 0
bound_40 <- 0
bound_45 <- 0
bound_50 <- 0
bound_55 <- 0
bound_60 <- 0
bound_65 <- 0

time_10 <- 0
time_15 <- 0
time_20 <- 0
time_25 <- 0
time_30 <- 0
time_35 <- 0
time_40 <- 0
time_45 <- 0
time_50 <- 0
time_55 <- 0
time_60 <- 0
time_65 <- 0

poolA_10 <- 0
poolA_15 <- 0
poolA_20 <- 0
poolA_25 <- 0
poolA_30 <- 0
poolA_35 <- 0
poolA_40 <- 0
poolA_45 <- 0
poolA_50 <- 0
poolA_55 <- 0
poolA_60 <- 0
poolA_65 <- 0

poolB_10 <- 0
poolB_15 <- 0
poolB_20 <- 0
poolB_25 <- 0
poolB_30 <- 0
poolB_35 <- 0
poolB_40 <- 0
poolB_45 <- 0
poolB_50 <- 0
poolB_55 <- 0
poolB_60 <- 0
poolB_65 <- 0

dec_10 <- 0
dec_15 <- 0
dec_20 <- 0
dec_25 <- 0
dec_30 <- 0
dec_35 <- 0
dec_40 <- 0
dec_45 <- 0
dec_50 <- 0
dec_55 <- 0
dec_60 <- 0
dec_65 <- 0

dane <- data.frame(N_size, n_groups, n_size, n_density, n_to_N_density, N_density, 
          inhibition_level, r_default, r_active, number_of_iterations, 
          to_A, to_B, fixed_for, n_A, n_B, over_time_A, over_time_B, event_starts, event_ends, 
          bound_10, time_10, poolA_10, poolB_10, dec_10, bound_15, time_15, poolA_15, poolB_15, dec_15,
          bound_20, time_20, poolA_20, poolB_20, dec_20, bound_25, time_25, poolA_25, poolB_25, dec_25, 
          bound_30, time_30, poolA_30, poolB_30, dec_30, bound_35, time_35, poolA_35, poolB_35, dec_35, 
          bound_40, time_40, poolA_40, poolB_40, dec_40, bound_45, time_45, poolA_45, poolB_45, dec_45, 
          bound_50, time_50, poolA_50, poolB_50, dec_50, bound_55, time_55, poolA_55, poolB_55, dec_55,
          bound_60, time_60, poolA_60, poolB_60, dec_60, bound_65, time_65, poolA_65, poolB_65, dec_65)

groups <- Grouping_variable()

runs_per_contrast <- 500
lambda0 <- 15

for(l in 1:5) {
  
  #contrast levels: 4.98%, 6.39%, 8.21%, 10.54%, 13.53%
  #n_A levels: 7.8735, 7.97925, 8.11575, 8.2905, 8.51475
  #modified (doubled contrasts) n_A levels: 8.247, 8.4585, 8.7315, 9.081, 9.5295
  #n_A levels: 9.8964, 10.1502, 10.4778, 10.8972, 11.4354
  
  if(l == 1 ) {n_A <- lambda0*(1/2+4.98/100)}
  if(l == 2 ) {n_A <- lambda0*(1/2+6.39/100)}
  if(l == 3 ) {n_A <- lambda0*(1/2+8.21/100)}
  if(l == 4 ) {n_A <- lambda0*(1/2+10.54/100)}
  if(l == 5 ) {n_A <- lambda0*(1/2+13.53/100)}
  
  n_B <- lambda0 - n_A
  print("Next level started")

  for(i in 1:runs_per_contrast) {
    
    network_connections <- Network_definition()
    state <- Network_init()
    s <- state[ ,1]
    t <- state[, 2]
    
    stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
    
    over_time_A <- mean(stimulus[ ,1:100])
    over_time_B <- mean(stimulus[ ,101:200])
    
    simulation <- Network_evolution()

    win_10 <- Which_wins_when(10)
    win_15 <- Which_wins_when(15)
    win_20 <- Which_wins_when(20)
    win_25 <- Which_wins_when(25)
    win_30 <- Which_wins_when(30)
    win_35 <- Which_wins_when(35)
    win_40 <- Which_wins_when(40)
    win_45 <- Which_wins_when(45)
    win_50 <- Which_wins_when(50)
    win_55 <- Which_wins_when(55)
    win_60 <- Which_wins_when(60)
    win_65 <- Which_wins_when(65)
    
    dane[i + (l-1)*runs_per_contrast, ] <- c(N_size, n_groups, n_size, n_density, n_to_N_density, N_density, 
                   inhibition_level, r_default, r_active, number_of_iterations,
                   to_A, to_B, fixed_for, n_A, n_B, over_time_A, over_time_B, event_starts, event_ends, 
                   win_10, win_15, win_20, win_25, win_30, win_35, win_40, win_45, win_50, win_55, 
                   win_60, win_65)
    
  }
}


write.csv2(dane,'SAT_with_FR_Threshold_1.csv')




View(dane)
getwd() 

#-----------------------------------------------

