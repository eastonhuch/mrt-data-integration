dat <- data.frame(y=rnorm(100))
mean_est_fun <- function(data) {
  print(data)
  y <- data$y
  function(theta) {
    y - theta
  }
}
fitted_model <- m_estimate(
  estFUN = mean_est_fun, 
  data   = dat,
  root_control = setup_root_control(start = c(1)))
coef(fitted_model)
mean(dat$y)
vcov(fitted_model)
var(dat$y) / nrow(dat)
dat