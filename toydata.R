require(ggplot2)

generate_data <- function(p0 = 0.2, p1 = 0.5, mu0 = 0, mu1 = 10, length = 100, 
                          X0 = 0, Y0 = mu0) {
  X = rep(0, length + 1)
  Y = rep(0, length + 1)
  
  X[1] = X0
  Y[1] = Y0
  
  for(i in 2:length) {
    # process model
    if (X[i-1] == 0) {
      X[i] = rbinom(1, 1, p0)
    } else {
      X[i] = rbinom(1, 1, p1)
    }
    
    # measure model
    if(X[i] == 0) {
      Y[i] = rnorm(1, mu0, 1)
    } else {
      Y[i] = rnorm(1, mu1, 1)
    }
  }
  
  list(X = X, Y = Y)
}