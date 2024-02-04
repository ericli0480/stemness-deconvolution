################################################################################
##  Purpose: Implement Method 4 in real data analysis using Shiny
##    -- Weighted least squares with constraints and quadratic programming 
##    -- Cross-validation to select h 
##  Functions to be used for Method 4.
##  ---------------------------------------------------------------------------
##  Platform, Bugs and Side-Effects:
##    -- Version 18
##    -- R x64 4.3.1 under Mac OS X
##    -- May get the error that matrix D is not positive definite. This is an 
##       issue with the numerical values in the dataset. In these cases our
##       model will not work without adjusting these values slightly.
##  ---------------------------------------------------------------------------
##  Author: Eric Li
##          Elkins High School
##          Missouri City, TX 77459
##          Phone: 281-302-9126
##          Email: ericjli0480@gmail.com
################################################################################

################################################################################
##
##                  Functions for Method 4 with a pre-specified h
##
################################################################################

# Shiny app libraries
library(shiny)
library(shinyFiles)
library(bslib)
library(thematic)

# Method 4 libraries
library(ggplot2)
library(quadprog)
library(ggstar)

# Kernel function: Gaussian kernel
kernel.function <- function(x) {
  1 / sqrt(2 * pi) * exp(-x * x / 2)
}

# Estimate stemness of all cell types at a single spot
method4.single.location <- function(x0, y0, h, T.data, K) {
  # x0 and y0 are the spatial coordinates of the spot of interest.
  # h is a positive scalar, the bandwidth we are using.
  # T.data is the training data in a data.frame format (e.g., dat1)
  # K is the total number of cell types
  
  # Extract cell weights and stemness
  W.mat <- as.matrix(T.data[, paste("W", 1:K, sep = "")])
  Y <- as.matrix(T.data[, "Y.orig"])
  
  # Find distance between the spot of interest and all other spots.
  # NOTE: distance as calculated here is a vector, whose length is the 
  # number of spots in the training data
  distance <- sqrt((T.data$x.coord - x0) ** 2 + (T.data$y.coord - y0) ** 2)
  
  # Calculate the kernel weights of the spots in the training data.
  kw <- kernel.function(distance / h) / h
  
  # Constrained optimization to estimate the stemness of all cell types 
  # at this spot
  Omega.W.mat <- kw * W.mat
  d.vec <- as.numeric(t(Y) %*% Omega.W.mat)
  D.mat <- t(W.mat) %*% Omega.W.mat
  A.mat <- cbind(diag(K), -diag(K))
  b0.vec <- c(rep(0, K), rep(-1, K))
  out <- solve.QP(
    Dmat = D.mat,
    dvec = d.vec,
    Amat = A.mat,
    bvec = b0.vec,
    meq = 0,
    factorized = FALSE
  )
  
  this.alpha.hat <- out$solution
  this.alpha.hat.unconstrained <- out$unconstrained.solution
  return(list(A = this.alpha.hat,
              B = this.alpha.hat.unconstrained))
}

# Run Method 4 at all the spots with a pre-specified h
method4.run <- function(h, dat1, K, M4progress) {
  # h is a positive scalar containing the bandwidth
  # dat1 contains weights, stemnesses, and spatial x, y coordinates
  # K is the number of cell types
  # M4progress is a progress indicator for the Shiny app
  
  N <- nrow(dat1) # Total number of cells
  
  alpha <- matrix(, nrow = N, ncol = K) 
  alpha.unconstrained <- matrix(, nrow = N, ncol = K)
  
  # Loop through all the spots, estimate stemnesses for each
  for (i in 1:N) {
    # Set detail for progress bar without incrementing it
    M4progress$inc(0, detail = paste("Running at spot ", i, sep = ""))
    
    # Coordinates of our target spot
    x0 <- dat1$x.coord[i]
    y0 <- dat1$y.coord[i]
    
    # Run method 4
    out <-
      method4.single.location(
        x0 = x0,
        y0 = y0,
        h = h,
        T.data = dat1,
        K = K
      )
    
    # Save constrained and unconstrained solutions
    alpha[i,] <- out$A
    alpha.unconstrained[i,] <- out$B
    
    # Increment progress indicator
    M4progress$inc(1/N)
  }
  
  # Return both constrained and unconstrained solutions
  return(list(alpha = alpha, alpha.unconstrained = alpha.unconstrained))
}


################################################################################
##
##                       Leave-one-out Cross-Validation to select h
##
################################################################################

# Calculate the cross-validated loss function for
# a given positive scalar bandwidth h
CV <- function(h, dat1 = dat1, K) {
  # The arguments are the same as in previous functions for method 4
  
  N <- nrow(dat1)
  CV.loss <- 0 # Return value
  
  for (i in 1:N) {
    # Find the training and validation datasets
    V.data <- dat1[i, , drop = F]
    T.data <- dat1[-i,]
    
    # Coordinates of the validation dataset, which is a single spot
    x0 <- V.data$x.coord[1]
    y0 <- V.data$y.coord[1]
    
    # Run Method 4 on the single spot in the validation dataset
    method4.return <- method4.single.location(
      x0 = x0,
      y0 = y0,
      h = h,
      T.data = T.data,
      K = K
    )
    
    # Save the constrained and unconstrained solutions
    this.alpha.hat <- method4.return$A
    this.alpha.hat.unconstrained <- method4.return$B # Not used here
    
    # Reconstruct the stemness at this spot using the values of alpha
    # we just found. Compare it with the true stemness.
    this.W <- as.numeric(V.data[, paste("W", 1:K, sep = "")])
    this.Y.hat <- sum(this.W * this.alpha.hat)
    CV.loss <- CV.loss + (V.data$Y.orig[1] - this.Y.hat) ** 2
  }
  
  return(CV.loss / N) # Average cross-validation loss over all spots
}

# Plot the cross-validation loss (U-shaped curves)
CV.plot <- function(lo, hi, num.pts, dat1, K, CVprogress) {
  # lo and hi are the ranges of h for which we will plot the CV loss.
  # num.pts is the number of points to be plotted.
  # dat1: the input data.
  # K: the number of cell types
  # CVprogress: the progress indicator object for the Shiny app.
  
  CV.loss <- NULL # The values of the cross-validation losses
  h.series <-
    seq(lo, hi, length = num.pts) # h values we will calculate CV loss for
  
  # Calculate cross-validation loss for each value of h
  for (h in h.series) {
    # Set detail for progress bar without incrementing it
    CVprogress$inc(0, detail = paste("Calculating CV.loss(", 
                                     round(h, digits = 4), 
                                     ")", 
                                     sep = ""))
    
    this.CV.loss <- CV(h = h, dat1 = dat1, K = K)
    CV.loss <- c(CV.loss, this.CV.loss)
    
    # Increment progress bar
    CVprogress$inc(1/num.pts)
  }
  
  # The value of h with lowest CV loss
  # The match function finds its index
  h.best <- h.series[match(min(CV.loss), CV.loss)]
  
  # Generate the cross-validation plot
  CV.df <- data.frame(h.series = h.series, CV.loss = CV.loss)
  CVplotOutput <- ggplot(CV.df, aes(x = h.series, y = CV.loss)) +
    geom_point() +
    geom_vline(xintercept = h.best, color = "red") + # Indicate h.best
    xlab("Bandwidth h") + ylab("Cross-validation Loss") +
    labs(title = paste("Cross-Validation (Best: h = ", h.best, ")", sep = ""))
  
  CVplotOutput
}