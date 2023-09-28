
n0 <- 20
n1 <- 40
rho_n <- n1 / n0
x0 <- rnorm(n=n0)
x1 <- rnorm(1, n=n1)
x <- c(x0, x1)
y0 <- rep(0, n0)
y1 <- rep(1, n1)
y <- c(y0, y1)

mod <- glm(y ~ x, family=binomial())
mod_coef <- mod$coefficients

dr_loglike <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sum1 <- sum(alpha + x1 * beta)
  sum2 <- sum(log(1 + rho_n * exp(alpha + x*beta)))
  -(sum1 - sum2)
}

optim_result <- optim(c(0, 0), dr_loglike, method="BFGS")
dr_coef <- optim_result$par
dr_coef
dr_coef[1] + log(rho_n)
mod_coef

# Run glm and then add log(n1 / n0) to intercept

x_min <- -4
x_max <- 5
x_values <- seq(x_min, x_max, 0.01)
dr_ratios <- exp(dr_coef[1] + x_values * dr_coef[2])
plot(NULL, type="n", xlim=c(x_min, x_max), ylim=c(0, 2),
     xlab="x", ylab="Density")
lines(density(x0), col=1)
lines(density(x1), col=2)
lines(x_values, dr_ratios, col=3)
legend("topleft", legend=c("Y=0", "Y=1", "Ratio"), col=1:3, lty=1)


##### Now do this with sim data
data <- generate_data(t_max=20, dof=10, n_internal=100, n_external=100, ar_param=0.5)
n0 <- sum(data$is_external)
n1 <- sum(data$is_internal)
rho_n <- n1 / n0

mod <- glm(is_internal ~ x2,
           family=binomial(), data=data)
X_all <- model.matrix(mod)
X <- as.matrix(X_all[,-1])
X0 <- as.matrix(X[data$is_external,])
X1 <- as.matrix(X[data$is_internal,])
mod_coef <- mod$coefficients
mod_coef[1] <- mod_coef[1] - log(rho_n)

dr_loglike <- function(theta) {
  alpha <- theta[1]
  beta <- theta[-1]
  sum1 <- sum(alpha + c(X1 %*% beta))
  sum2 <- sum(log(1 + rho_n * exp(X_all %*% theta)))
  -(sum1 - sum2)
}

p <- ncol(X_all)
optim_result <- optim(rep(0, p), dr_loglike, method="BFGS")
dr_coef <- optim_result$par
dr_coef - mod_coef
data$raw_tilt_ratios <- c(exp(X_all %*% mod_coef))

# This is wrong now
for (j in seq(p-1)) {
  x_j <- X[,j]
  x_min <- min(x_j)
  x_max <- max(x_j)
  x_values <- seq(x_min, x_max, 0.1)
  dr_ratios <- exp(mod_coef[1] + x_values * mod_coef[1+j])
  plot(NULL, type="n", xlim=c(x_min, x_max), ylim=c(0, 3),
       xlab="x", ylab="Density")
  density0 <- density(x0, from=x_min, to=x_max)
  density1 <- density(x1, from=x_min, to=x_max)
  lines(density0, col=1)
  lines(density1, col=2)
  points(x_j, data$raw_tilt_ratios, col=3)
  lines(density0$x, density1$y / density0$y)
  legend("topleft", legend=c("Y=0", "Y=1", "Ratio"), col=1:3, lty=1)
}

plot(data$x1[data$is_internal], data$x2[data$is_internal])
points(data$x1[data$is_external], data$x2[data$is_external], col=2)
boxplot(data$raw_tilt_ratios[data$is_internal])
boxplot(data$raw_tilt_ratios[data$is_external])
hist(data$raw_tilt_ratios)
mean(data$raw_tilt_ratios[data$is_internal])
mean(data$raw_tilt_ratios[data$is_external])


