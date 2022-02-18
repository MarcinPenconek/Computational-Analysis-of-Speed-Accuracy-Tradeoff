
#---------------------------------------------------
# ANALYSIS OF DOT MOTION DISCRIMINATION EXPERIMENT
#---------------------------------------------------

# Reaction Time Paradigm
# Data: RT_results_dot_motion.csv

dane1 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/Method/RT_dot_motion_1.csv", row.names=1)
dane2 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/Method/RT_dot_motion_2.csv", row.names=1)
dane3 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/Method/RT_dot_motion_3.csv", row.names=1)
dane4 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/Method/RT_dot_motion_4.csv", row.names=1)

dane <- merge(dane1, dane2, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane3, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane4, all.x = TRUE, all.y = TRUE)

selection <- (dane$dec_th != 0 & dane$time_th >= 1100 & dane$time_th <= 3000)
correct <- (dane$dec_th == 1)
inputs <- dane$n_A
labels <- c(1, 3.2, 6.4, 12.8, 25.6, 51.2)

#----------------------------------------------
# GRAPHIC PARAMETERS
#----------------------------------------------

osie <- 1.2
opis <- 1.2

#----------------------------------------------
# Accuracy

#dot motion coherence levels: 3.2%, 6.4%, 12.8%, 25.6%, 51.2%
#n_A levels: 7.5, 7.74, 7.98, 8.46, 9.42, 11.34

coherence <- (inputs == 7.5) + 3.2*(inputs == 7.74) + 6.4*(inputs == 7.98) + 12.8*(inputs == 8.46) + 25.6*(inputs == 9.42) + 51.2*(inputs == 11.34)


proc_correct <- c(mean(100*correct[selection & coherence == 1]), mean(100*correct[selection & coherence == 3.2]), 
                  mean(100*correct[selection & coherence == 6.4]), mean(100*correct[selection & coherence == 12.8]), 
                  mean(100*correct[selection & coherence == 25.6]), mean(100*correct[selection & coherence == 51.2]))

proc_correct

plot(labels, proc_correct, log = "x", xlim = c(1, 100), ylim = c(50,100), xlab = "Coherence %", ylab = "% Correct", 
       pch = 19, cex.axis = osie, cex.lab = opis, cex = 1.5)

t <- 1:100
lines(t, 100*(1 - 1/2*exp(-((t/12.5)^1.18))))
lines(t, 100*(1 - 1/2*exp(-((t/6)^1.7))), lty = 2)
lines(t, 100*(1 - 1/2*exp(-((t/15)^1.1))), lty = 2)

legend(9, 65, legend = c("Model", "Weibull Function", "Experimental Data"), pch = c(19, NaN, NaN), lty = c(NaN,1,2), cex = 1, bty = "n")

sd_pr_correct <- c(sd(100*correct[selection & coherence == 1]), sd(100*correct[selection & coherence == 3.2]), 
                   sd(100*correct[selection & coherence == 6.4]), sd(100*correct[selection & coherence == 12.8]), 
                   sd(100*correct[selection & coherence == 25.6]), sd(100*correct[selection & coherence == 51.2]))
n_pr_correct <- 12000*c(mean(selection & coherence == 1), mean(selection & coherence == 3.2), 
                  mean(selection & coherence == 6.4), mean(selection & coherence == 12.8), 
                  mean(selection & coherence == 25.6), mean(selection & coherence == 51.2))
(12000-sum(n_pr_correct))/12000

sd_pr_correct
se <- sd_pr_correct/n_pr_correct^0.5
se

#-----------------------

f <- function(x, a, b) {
  1 - 1/2*exp(-(x/b)^a)
}

correct_selected <- correct[selection]
coherence_selected <- coherence[selection]

fit <- nls(correct_selected ~ f(coherence_selected,a,b), start=list(a=1.2, b=10))
summary(fit)

#--------------------------------
# Reaction Time


RT <- dane$time_th - 1000

RT_correct <- c(mean(RT[selection & coherence == 1]), mean(RT[correct & selection & coherence == 3.2]), 
                mean(RT[correct & selection & coherence == 6.4]), mean(RT[correct & selection & coherence == 12.8]), 
                mean(RT[correct & selection & coherence == 25.6]), mean(RT[correct & selection & coherence == 51.2]))

RT_error <- c(NaN, mean(RT[!correct & selection & coherence == 3.2]), 
              mean(RT[!correct & selection & coherence == 6.4]), mean(RT[!correct & selection & coherence == 12.8]), 
              mean(RT[!correct & selection & coherence == 25.6]), mean(RT[!correct & selection & coherence == 51.2]))

sd_correct <- c(sd(RT[selection & coherence == 1]), sd(RT[correct & selection & coherence == 3.2]), 
                sd(RT[correct & selection & coherence == 6.4]), sd(RT[correct & selection & coherence == 12.8]), 
                sd(RT[correct & selection & coherence == 25.6]), sd(RT[correct & selection & coherence == 51.2]))

sd_error <- c(NaN, sd(RT[!correct & selection & coherence == 3.2]), 
              sd(RT[!correct & selection & coherence == 6.4]), sd(RT[!correct & selection & coherence == 12.8]), 
              sd(RT[!correct & selection & coherence == 25.6]), sd(RT[!correct & selection & coherence == 51.2]))

sd_correct
sd_error

n_sd_correct <- 12000*c(mean(selection & coherence == 1), mean(correct & selection & coherence == 3.2), 
                  mean(correct & selection & coherence == 6.4), mean(correct & selection & coherence == 12.8), 
                  mean(correct & selection & coherence == 25.6), mean(correct & selection & coherence == 51.2))
n_sd_correct

n_sd_error <- 12000*c(NaN, mean(!correct & selection & coherence == 3.2), 
                        mean(!correct & selection & coherence == 6.4), mean(!correct & selection & coherence == 12.8), 
                        mean(!correct & selection & coherence == 25.6), mean(!correct & selection & coherence == 51.2))
n_sd_error

se_correct <- sd_correct/n_sd_correct^0.5
se_error <- sd_error/n_sd_error^0.5
se_correct
se_error


RT_exp_correct <- c(825.8, 806.4, 758.4, 674.9, 541.7, 432.2)
RT_exp_error <- c(NaN, 844.5, 831.3, 829.9, 736, NaN)

RT_error[6] <- NaN

plot(labels, RT_error, log = "x", type = "l", xlim = c(1, 100), ylim = c(0, 1000), xlab = "Coherence %", ylab = "Reaction Time (ms)", 
     pch = 19, cex.axis = osie, cex.lab = opis, cex = 1.5)
points(labels, RT_correct, type = "l", cex = 1.5)
points(labels, RT_exp_correct, type = "l", lty = 2)
points(labels, RT_exp_error, type = "l", lty = 2)

#points(labels, RT_exp_correct, pch = 15, col = "grey", cex = 2)
#points(labels, RT_exp_error, pch = 0, col = "grey", cex = 2)
points(labels, RT_correct, pch = 19, cex = 1.5)
points(labels, RT_error, pch = 1, cex = 1.5)
legend(4, 300, legend = c("Model: Correct Decisions", "Model: Error Decisions", "Experimental Data"), pch = c(19, 1, NaN), lty = c(1,1,2), cex = 1, bty = "n")


RT_correct
RT_error
RT_all <- c(RT_correct, RT_error)

RT_exp <- c(RT_exp_correct, RT_exp_error)

RT_diff <- abs(RT_exp - RT_all)/RT_exp
RT_diff

mean(RT_diff, na.rm = TRUE)
max(RT_diff, na.rm = TRUE)


