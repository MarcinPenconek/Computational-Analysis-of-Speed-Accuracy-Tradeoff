
dane1 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inputs/SAT_with_inputs_1.csv", row.names=1)
dane2 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inputs/SAT_with_inputs_2.csv", row.names=1)
dane3 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inputs/SAT_with_inputs_3.csv", row.names=1)
dane4 <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inputs/SAT_with_inputs_4.csv", row.names=1)

dane_18inh <- read.csv2("~/Marcin/Studia/NeuroNetwork Model II/SAT paper/SAT Analysis/Inputs/SAT_with_inputs_inh_18.csv", row.names=1)

dane <- merge(dane1, dane2, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane3, all.x = TRUE, all.y = TRUE)
dane <- merge(dane, dane4, all.x = TRUE, all.y = TRUE)

data_size <- 20000

Search_space <- function(a, b, diam, ers) {
  d <- diam
  numb_of_errors <- 0
  while(numb_of_errors < ers) {
    d <- d + diam
    selection <- (2*dane$over_time_A > (a - d) & 2*dane$over_time_A < (a + d) & 2*dane$over_time_B > (b - d) & 2*dane$over_time_B < (b + d))
    correct <- (dane$dec_50 * (dane$over_time_A - dane$over_time_B) > 0)
    error <- (dane$dec_50 * (dane$over_time_A - dane$over_time_B) < 0)
    numb_of_errors <- sum(error & selection)
  }
  RT_correct <- mean(dane$time_50[correct & selection])
  RT_error <- mean(dane$time_50[error & selection])
  numb_of_correct <- sum(correct & selection)
  numb_of_errors <- sum(error & selection)
  return(c(RT_correct, RT_error, numb_of_correct, numb_of_errors, d))
}


wynik <- Search_space(7,8,0.01,10)
wynik

#---------------------------------------
# FIGURE: Errors
# not used


osie <- 1.2
opis <- 1.2

x <- 0:20
y <- 0:20
plot(x, y, xlim = c(0,20), ylim = c(0,20), type ="l", lty = 1, xlab = "Input to A", ylab = "Input to B", 
     cex.lab = opis, cex.axis = osie)

kolor <- numeric(data_size)
wielkosc <- numeric(data_size) + 0.2
for(i in 1:data_size) {
  if(dane$dec_50[i] == sign(dane$n_B[i] - dane$n_A[i])) { 
    kolor[i] <- "orange"
    wielkosc[i] <- 0.5
    }
  if(dane$dec_50[i] == sign(dane$n_A[i] - dane$n_B[i])) { kolor[i] <- "grey"}
  if(dane$dec_50[i] == 0) { kolor[i] <- "white"}
}

points(dane$n_A, dane$n_B, pch = 16, col = kolor, cex = wielkosc)


#---------------------------------------
# FIGURE 6b: Cases with RT on errors shorter than RT on correct

osie <- 1.2
opis <- 1.2

set.seed(1234)

diam <- 0.2
x <- 0:20
y <- 0:20
plot(x, y, xlim = c(0,20), ylim = c(0,20), type ="l", lty = 1, xlab = "Input to A", ylab = "Input to B", 
     cex.lab = opis, cex.axis = osie)

for(i in 1:50000) {
  kolor <- "orange"
  a <- runif(1, min = 0, max = 20)
  b <- runif(1, min = 0, max = 20)
  if((20-a)^2 + (20-b)^2 < 50) { points(a, b, pch = 16, cex = 1, col = kolor) }
}

ile_kropek <- 0
ile_bledow <- 0

diams <- numeric(1000)
for(i in 1:1000) {
  kolor <- "light grey"
  a <- runif(1, min = 0, max = 20)
  b <- runif(1, min = 0, max = 20)
  wynik <- Search_space(a, b, diam, 50)
  diams[i] <- wynik[5]
  if(wynik[5] <= 2.5) {
    roznica <- wynik[2] - wynik[1]
    if(roznica < 0) { 
      kolor <- "red"
      ile_bledow <- ile_bledow + 1
      }
    if(roznica > 0) { kolor <- "light grey"}
    ile_kropek <- ile_kropek + 1
    points(a, b, pch = 16, col = kolor)
  }
} 
legend(10, 5, col = c("red", "orange"), pch = 16, legend = c("RT(E) < RT(C)", "Overstimulation"), cex = 1.2, bty = "n")

x <- 0:20
y <- 0:20
lines(x, y, lty = 1)




#-------------------------------
# Figure 5a: Saturation


Calculate_correct <- function(slope, tolerance, point, interval, numb) {
  
  correct <- (dane$dec_50 * (dane$over_time_A - dane$over_time_B) > 0)
  error <- (dane$dec_50 * (dane$over_time_A - dane$over_time_B) < 0)
  bin_selection <- (((2*dane$over_time_A + 2*dane$over_time_B) > (point - interval)) & ((2*dane$over_time_A + 2*dane$over_time_B) < (point + interval)))
  slope_selection_A <- (dane$over_time_A/dane$over_time_B > (slope - tolerance)) & (dane$over_time_A/dane$over_time_B < (slope + tolerance))
  slope_selection_B <- (dane$over_time_B/dane$over_time_A > (slope - tolerance)) & (dane$over_time_B/dane$over_time_A < (slope + tolerance))
  slope_selection <- slope_selection_A | slope_selection_B
  all_relevant <- sum(correct[bin_selection & slope_selection]) + sum(error[bin_selection & slope_selection])
  proc_correct <- sum(correct[bin_selection & slope_selection])/all_relevant
  
  if(all_relevant > numb) {
    return(c(proc_correct, all_relevant))
  } else {
    return(c(NaN, all_relevant))
  }
  
}


wynik_correct <- Calculate_correct(1.5, 0.1, 7, 1, 30)
wynik_correct

#---------------------------------------
# FIGURE: Overstimulation/Saturation - version with contrasts (not used)

tolerance <- 0.25
interval <- 2
lambda0 <- 15

osie <- 1.2
opis <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

l <- 1

if(l == 1 ) {n_A <- lambda0*(1/2+4.98/100)}
if(l == 2 ) {n_A <- lambda0*(1/2+6.39/100)}
if(l == 3 ) {n_A <- lambda0*(1/2+8.21/100)}
if(l == 4 ) {n_A <- lambda0*(1/2+10.54/100)}
if(l == 5 ) {n_A <- lambda0*(1/2+13.53/100)}

n_B <- lambda0 - n_A
slope <- n_A/n_B
slope

lamb <- 1:40
correct_lamb <- numeric(40)
for(i in 1:40) {
  correct_lamb[i] <- Calculate_correct(slope, tolerance, i, interval, 100)[1]
}

plot(lamb, correct_lamb, type = "l", xlim = c(0,40), ylim = c(0.45,1), lwd = l, col = kolor[l], 
     xlab = "Sum of Inputs to A and B", ylab = "% Correct", cex.lab = opis, cex.axis = osie)

for(j in 2:5) {

  l <- j
  if(l == 1 ) {n_A <- lambda0*(1/2+4.98/100)}
  if(l == 2 ) {n_A <- lambda0*(1/2+6.39/100)}
  if(l == 3 ) {n_A <- lambda0*(1/2+8.21/100)}
  if(l == 4 ) {n_A <- lambda0*(1/2+10.54/100)}
  if(l == 5 ) {n_A <- lambda0*(1/2+13.53/100)}
  n_B <- lambda0 - n_A
  slope <- n_A/n_B
  
  lamb <- 1:40
  correct_lamb <- numeric(40)
  for(i in 1:40) {
    correct_lamb[i] <- Calculate_correct(slope, tolerance, i, interval, 100)[1]
  }
  lines(lamb, correct_lamb, lwd = l, col = kolor[l])
  
}

legend(0, 0.7, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
abline(v = 30, lty = 2)



#---------------------------------------
# FIGURE 5a: Overstimulation/Saturation - broad-band version


interval <- 2
lambda0 <- 15

osie <- 1.2
opis <- 1.2

kolor <- numeric(2)
kolor[1] <- "grey"
kolor[2] <- "dark grey"

slope <- 1.25
tolerance <- 0.5

lamb <- 1:40
correct_lamb <- numeric(40)
for(i in 1:40) {
  correct_lamb[i] <- Calculate_correct(slope, tolerance, i, interval, 100)[1]
}

plot(lamb, correct_lamb, type = "l", xlim = c(0,40), ylim = c(0.5,1), lwd = 2, col = kolor[1], 
     xlab = "Sum of Inputs to A and B", ylab = "% Correct", cex.lab = opis, cex.axis = osie)

slope <- 1.75
tolerance <- 0.5

lamb <- 1:40
correct_lamb <- numeric(40)
for(i in 1:40) {
  correct_lamb[i] <- Calculate_correct(slope, tolerance, i, interval, 100)[1]
}
lines(lamb, correct_lamb, lwd = 3, col = kolor[2])


legend(0, 0.63, col = kolor, lwd = c(2, 4), c("Low Contrast", "High Contrast"), cex = opis, bty = "n")
abline(v = 30, lty = 2)


#---------------------------------------
# FIGURE: Comparison of Errors
# not used

osie <- 1.2
opis <- 1.2

data_size <- 20000

x <- 0:20
y <- 0:20
plot(x, y, xlim = c(0,20), ylim = c(0,20), type ="l", lty = 1, xlab = "Input to A", ylab = "Input to B", 
     cex.lab = opis, cex.axis = osie)

kolor <- numeric(data_size)
wielkosc <- numeric(data_size) + 0.2
for(i in 1:data_size) {
  if(dane$dec_50[i] == sign(dane$n_B[i] - dane$n_A[i])) { 
    kolor[i] <- "orange"
    wielkosc[i] <- 0.5
  }
  if(dane$dec_50[i] == sign(dane$n_A[i] - dane$n_B[i])) { kolor[i] <- "light grey"}
  if(dane$dec_50[i] == 0) { kolor[i] <- "white"}
}

points(dane$n_A[dane$n_A > dane$n_B], dane$n_B[dane$n_A > dane$n_B], pch = 16, col = kolor[dane$n_A > dane$n_B], cex = wielkosc[dane$n_A > dane$n_B])

data_size <- 10000

kolor <- numeric(data_size)
wielkosc <- numeric(data_size) + 0.2
for(i in 1:data_size) {
  if(dane_18inh$dec_50[i] == sign(dane_18inh$n_B[i] - dane_18inh$n_A[i])) { 
    kolor[i] <- "orange"
    wielkosc[i] <- 0.5
  }
  if(dane_18inh$dec_50[i] == sign(dane_18inh$n_A[i] - dane_18inh$n_B[i])) { kolor[i] <- "light grey"}
  if(dane_18inh$dec_50[i] == 0) { kolor[i] <- "white"}
}

points(dane_18inh$n_A[dane_18inh$n_A < dane_18inh$n_B], dane_18inh$n_B[dane_18inh$n_A < dane_18inh$n_B], pch = 16, col = kolor[dane_18inh$n_A < dane_18inh$n_B], cex = wielkosc[dane_18inh$n_A < dane_18inh$n_B])

legend(0, 20, col = "orange", pch = 16, "Errors", cex = opis, bty = "n")


#---------------------------------------
# FIGURE: Overstimulation/Saturation - version with contrasts

data_size <- 10000
dane <- dane_18inh

tolerance <- 0.25
interval <- 2
lambda0 <- 15

osie <- 1.2
opis <- 1.2

kolor <- numeric(5)
kolor[1] <- 5
kolor[2] <- "light blue"
kolor[3] <- 4
kolor[4] <- "blue"
kolor[5] <- "dark blue"

l <- 1

if(l == 1 ) {n_A <- lambda0*(1/2+4.98/100)}
if(l == 2 ) {n_A <- lambda0*(1/2+6.39/100)}
if(l == 3 ) {n_A <- lambda0*(1/2+8.21/100)}
if(l == 4 ) {n_A <- lambda0*(1/2+10.54/100)}
if(l == 5 ) {n_A <- lambda0*(1/2+13.53/100)}

n_B <- lambda0 - n_A
slope <- n_A/n_B
slope

lamb <- 1:40
correct_lamb <- numeric(40)
for(i in 1:40) {
  correct_lamb[i] <- Calculate_correct(slope, tolerance, i, interval, 100)[1]
}

plot(lamb, correct_lamb, type = "l", xlim = c(0,40), ylim = c(0.45,1), lwd = l, col = kolor[l], 
     xlab = "Sum of Inputs to A and B", ylab = "% Correct", cex.lab = opis, cex.axis = osie)

for(j in 2:5) {
  
  l <- j
  if(l == 1 ) {n_A <- lambda0*(1/2+4.98/100)}
  if(l == 2 ) {n_A <- lambda0*(1/2+6.39/100)}
  if(l == 3 ) {n_A <- lambda0*(1/2+8.21/100)}
  if(l == 4 ) {n_A <- lambda0*(1/2+10.54/100)}
  if(l == 5 ) {n_A <- lambda0*(1/2+13.53/100)}
  n_B <- lambda0 - n_A
  slope <- n_A/n_B
  
  lamb <- 1:40
  correct_lamb <- numeric(40)
  for(i in 1:40) {
    correct_lamb[i] <- Calculate_correct(slope, tolerance, i, interval, 100)[1]
  }
  lines(lamb, correct_lamb, lwd = l, col = kolor[l])
  
}

legend(0, 0.7, col = kolor, lwd = c(1, 2, 3, 4, 5), c("Contrast 4.98%", "6.39%", "8.21%", "10.54%", "13.53%"), cex = opis, bty = "n")
abline(v = 30, lty = 2)

