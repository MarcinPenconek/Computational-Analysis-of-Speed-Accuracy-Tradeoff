#-------------------------------------------------------------#
#                  Neural Network Decison Model               #
#                        Marcin Penconek                      #
#-------------------------------------------------------------#


# Fixed firing model
# density: 0.55 i 0.35
# rates: 0.07 i 0.006



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
# Defining stimulus

Stimulus_definition_non_random = function() {
  
  stimulus <- c(rep(n_A, len = to_A), numeric(n_size - to_A), 
                rep(n_B, len = to_B), numeric(n_size - to_B), numeric(size - 2*n_size))
  
  return(stimulus)
}




#----------------------------------------
# Network evolution function

Network_evolution = function() {
  
  current_time <- 0
  current_stimulus <- numeric(size)
  stimulus_started <- 0
  stimulus_finished <- 0

  s <- state[ ,1]
  t <- state[ ,2]

  simulation <- matrix(0, nrow = number_of_iterations, ncol = (5 + n_groups + 1))
  reported_everactive <- NA
  
  for(j in 1:number_of_iterations) {
    
    if(current_time < when_start) {
      everactive <- s
    } else {
      reported_everactive <- sum(everactive)
    }
    
    if(stimulus_started == 0 & current_time > event_starts & current_time < event_ends) {
      current_stimulus <- stimulus
      stimulus_started <- 1
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
    
    active_groups <- numeric(n_groups + 1)
    for(i in 1:(n_groups + 1)) {
      active_groups[i] <- sum(s*groups[ ,i + 1])
    }
    
    simulation[j, ] <- c(current_time, Min_i, s[Min_i], 0, reported_everactive, active_groups)

    active <- sum(s)    
    cutoff <- (active/size)^2/current_inhibition
    
    active_connections <- sum(s*network_connections[Min_i, ])
    all_connections <- sum(network_connections[Min_i, ])
    
    change_factor <- (current_stimulus[Min_i] + active_connections)/(current_stimulus[Min_i] + all_connections)
    
    if(change_factor > cutoff) {
      s[Min_i] <- 1
      everactive[Min_i] <- 1
      simulation[j, 4] <- 1
      t[Min_i] <- current_time + 1/r_active
    } else {
      s[Min_i] <- 0
      t[Min_i] <- current_time + rexp(1, rate = r_default)
    }
    
  }
  return(simulation)
}

#----------------------------------------------
# Visualization
  
Population_plot <- function(x, y1, y2, y3, decision, ...) {
  
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
       xlab = "", ylab = "% Active",
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  lines(x, y1, type = "l", col = "red", lwd=1.5)
  lines(x, y2, type = "l", col = "blue", lwd=1.5)
  legend(decision[2]-200, 105, legend = dectext, cex = 1.2, bty="n")
  if(decision[1] != 0) {abline(v = decision[2])}
  if(n_A + n_B >0) {
    abline(v = event_starts, lwd=1, lty=2)
    abline(v = event_ends, lwd=1, lty=2)
  } else {
    legend(200, 100, legend = c("A", "B"), cex = 2, lwd=c(3,3), col=c("red","blue"), bty = "n")
  }
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
amplitude <- 0.1
alpha_freq <- 10
period <- (1000/alpha_freq)/(2*3.14)

#----------------------------------------------
# Simulation parameters

number_of_iterations <- 40000
threshold_level <- 75

to_A <- 50
to_B <- 50

n_A <- 10
n_B <- 0

event_starts <- 1000
event_ends <- 2000


#------------------------------------------------
# SECTION 3: Running the model
#------------------------------------------------


for(i in 1:1) {
  
  groups <- Grouping_variable()
  network_connections <- Network_definition()
  state <- Network_init()
  s <- state[ ,1]
  t <- state[, 2]
  
  when_start <- runif(1, min = 0, max = 5000)
  stimulus <- Stimulus_definition_non_random()
  simulation <- Network_evolution()
  decision <- c(0,0)

  
  Population_plot(simulation[ ,1], simulation[ ,6], simulation[ ,7], simulation[ ,8]*(100/N_size), decision)
  
#  if(i == 1) {
#    plot(simulation[, 1], simulation[, 5], ylim = c(0,1000), xlim = c(0, 6000), type = "l")
#  } else {
#        if(decision[1] != 0) 
#          {
#          lines(simulation[, 1], simulation[, 5], type = "l")
#        }
#  }
  
}

#plot(log(simulation[,6]), log(simulation[,7]), type="l")


#----------------------------------------------
# GRAPHIC PARAMETERS
#----------------------------------------------

osie <- 1.2
opis <- 1.2

#-----------------------------------------------
# FIGURE: Scatter plot
#-----------------------------------------------

set.seed(1234)

n_A <- 0
n_B <- 0

groups <- Grouping_variable()
network_connections <- Network_definition()
state <- Network_init()
s <- state[ ,1]
t <- state[, 2]

when_start <- runif(1, min = 0, max = 5000)
stimulus <- Stimulus_definition_non_random()
simulation <- Network_evolution()

plot(simulation[,1][simulation[,1] > 1000 & simulation[,1] < 2000 & simulation[,4] == 1 & simulation[,2] >= 0  & simulation[,2] <= 1000] - 1000, 
     simulation[,2][simulation[,1] > 1000 & simulation[,1] < 2000 & simulation[,4] == 1 & simulation[,2] >= 0 & simulation[,2] <= 1000], 
     pch = 124, col = "black", cex = 0.3, xlab = "Time (ms)", ylab = "Neurons", cex.axis = osie, cex.lab = opis)

