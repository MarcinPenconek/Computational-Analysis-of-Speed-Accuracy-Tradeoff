#-------------------------------------------------------------#
#                  Neural Network Decison Model               #
#                        Marcin Penconek                      #
#-------------------------------------------------------------#


# Different decision function

# model for initial conditions analysis
# density: 0.55 i 0.36
# iterations: 100000

# fixed firing 
# 0.07
# 0.006

# different decision functions



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
# Random initial conditions


Network_init_random = function() {
  
  inh_1 <- runif(1, min = 0, max = 1)
  inh_2 <- runif(1, min = 0, max = 1)
  
  prob_s <- runif(size, min = 0, max = 1)
  t <- rexp(size, rate = r_default) 
  s <- logical(size)
  
  for(i in (2*n_size+1):size) {
    if(prob_s[i] < inhibition_level) {
      s[i] <- TRUE
      t[i] <- rexp(1, rate = r_active)
    }
    else {
      s[i] <- FALSE
    }
  }
  
  for(i in 1:n_size) {
    if(prob_s[i] < inh_1) {
      s[i] <- TRUE
      t[i] <- rexp(1, rate = r_active)
    }
    else {
      s[i] <- FALSE
    }
  }
  
  for(i in (n_size+1):(2*n_size)) {
    if(prob_s[i] < inh_2) {
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
    simulation[j, 7] <- s[Min_i]*(Min_i > 100 & Min_i <= 200)
    simulation[j, 8] <- s[Min_i]*(Min_i > 200)
    
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



Firing_Rate_Analysis <- function(sim, beg_index, end_index) {
  
  A_count <- sum(sim[beg_index:end_index, 6])
  B_count <- sum(sim[beg_index:end_index, 7])
  C_count <- sum(sim[beg_index:end_index, 8])
  time_bin <- (sim[end_index, 1] - sim[beg_index, 1])/1000
  
  return(c((A_count/n_size)/time_bin, (B_count/n_size)/time_bin, (C_count/N_size)/time_bin, time_bin))
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

number_of_iterations <- 15000
threshold_level <- 50

to_A <- 50
to_B <- 50
fixed_for <- 30

#contrast levels: 4.98%, 6.39%, 8.21%, 10.54%, 13.53%
#n_A levels: 7.8735, 7.97925, 8.11575, 8.2905, 8.51475
#modified (doubled contrasts) n_A levels: 8.247, 8.4585, 8.7315, 9.081, 9.5295


n_A <- 0
n_B <- 0

event_starts <- 0
event_ends <- 0

duration <- round((event_ends - event_starts)/fixed_for + 1)

#-----------------------------------
# Parameters for Frequency Analysis

updates_per_30ms <- 440
freq_step <- 145




#------------------------------------------------
# SECTION 4: Generating a sample of simulations
#------------------------------------------------
# GENERATING DATABASE

number_of_simulations <- 500

over_time_A <- 0
over_time_B <- 0

end_state <- NaN
initial_time <- 50
end_time <- 1000

initial_FR_A <- 0
initial_FR_B <- 0
initial_FR_C <- 0
initial_bin <- 0
end_FR_A <- 0
end_FR_B <- 0
end_FR_C <- 0
end_bin <- 0

dane <- data.frame(N_size, n_groups, n_size, n_density, n_to_N_density, N_density, 
          inhibition_level, r_default, r_active, number_of_iterations, threshold_level, 
          to_A, to_B, fixed_for, n_A, n_B, over_time_A, over_time_B, event_starts, event_ends, 
          initial_time, initial_FR_A, initial_FR_B, initial_FR_C, initial_bin, 
          end_time, end_FR_A, end_FR_B, end_FR_C, end_bin, end_state)

groups <- Grouping_variable()

for(i in 1:number_of_simulations) {
  
  network_connections <- Network_definition()
  state <- Network_init_random()
  s <- state[ ,1]
  t <- state[, 2]
  
  stimulus <- Stimulus_definition_t_random(duration, n_A, n_B)
  simulation <- Network_evolution()

  initial_not_done <- 1
  end_not_done <- 1
  
  for(j in 1:number_of_iterations) {
    if(simulation[j,1] > initial_time & initial_not_done) {
      initial_FR <- Firing_Rate_Analysis(simulation, j - updates_per_30ms, j)
      initial_not_done <- 0
    }
  }  
  
  for(j in 1:number_of_iterations) {
    if(simulation[j,1] > end_time & end_not_done) {
      end_FR <- Firing_Rate_Analysis(simulation, j - updates_per_30ms, j)
      end_not_done <- 0
    }
  }
  
  end_state <- NaN

  if(end_FR[1] > 40 & end_FR[2] < 40) { end_state <- 1 }
  if(end_FR[1] < 40 & end_FR[2] > 40) { end_state <- -1 }
  if(end_FR[1] < 40 & end_FR[2] < 40) { end_state <- 0 }
  
  dane[i, ] <- c(N_size, n_groups, n_size, n_density, n_to_N_density, N_density, 
                 inhibition_level, r_default, r_active, number_of_iterations, threshold_level,
                 to_A, to_B, fixed_for, n_A, n_B, over_time_A, over_time_B, event_starts, event_ends, 
                 initial_time, initial_FR, end_time, end_FR, end_state)

}


#write.csv2(dane,'Basin_of_attraction.csv')

#----------------------------------------------
# GRAPHIC PARAMETERS
#----------------------------------------------

osie <- 1.2
opis <- 1.2
grubosc <- 3


#---------------------------
# Plot: Basin of Attraction

dane <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/Method/Basin_of_attraction.csv", row.names=1)

kolor <- "black"

for(k in 1:500) {
  if(dane$end_state[k] == 1) { kolor[k] <- "red" }
  if(dane$end_state[k] == -1) { kolor[k] <- "blue" }
  if(dane$end_state[k] == 0) { kolor[k] <- "light grey" }
}


plot(dane$initial_FR_A, dane$initial_FR_B, 
     ylim = c(0,70), xlim = c(0,70),
     type = "p",  
     xlab = "Initial Firing Rate in A (Hz)", ylab = "Initial Firing Rate in B (Hz)", 
     #xlab = "", ylab = "", 
     pch = 16, cex = 1, cex.axis = osie, cex.lab = opis, col = kolor)

legend(25, 70, legend = c("Convergence to A", "Convergance to B", "Conv. to Spont. State"), cex = opis, pch = 16, col=c("red","blue","light grey"), bty = "n")




