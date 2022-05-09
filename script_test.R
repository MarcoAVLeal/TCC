library(manipulate)


# A naive implementation of the Nadaraya-Watson estimator
nw <- function(x, X, Y, h, K = dnorm) {
  
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  
  # Matrix of size n x length(x) (rbind() is called for ensuring a matrix
  # output if x is a scalar)
  Kx <- rbind(sapply(X, function(Xi) K((x - Xi) / h) / h))
  
  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  
  # Means at x ("drop" to drop the matrix attributes)
  drop(W %*% Y)
  
}

# Generate some data to test the implementation
set.seed(12345)
n <- 100
eps <- rnorm(n, sd = 2)
m <- function(x) x^2 * cos(x)
# m <- function(x) x - x^2 # Works equally well for other regression function
X <- rnorm(n, sd = 2)
Y <- m(X) + eps
x_grid <- seq(-10, 10, l = 500)



manipulate::manipulate({
  
  # Plot data
  plot(X, Y)
  rug(X, side = 1); rug(Y, side = 2)
  lines(x_grid, m(x_grid), col = 1)
  lines(x_grid, nw(x = x_grid, X = X, Y = Y, h = h), col = 2)
  legend("topright", legend = c("True regression", "Nadaraya-Watson"),
         lwd = 2, col = 1:2)
  
}, h = manipulate::slider(min = 0.01, max = 10, initial = 0.5, step = 0.01))


