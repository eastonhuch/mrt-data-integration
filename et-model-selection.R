set.seed(1)
data <- generate_data(n_internal=400, n_external=400)

# Bad
tilt_mod1 <- glm(
  is_internal ~ 1,
  family=binomial(), data=data)
print(AIC(tilt_mod1))

# Worse
tilt_mod2 <- glm(
  is_internal ~ x1,
  family=binomial(), data=data)
print(AIC(tilt_mod2))

# Better
tilt_mod3 <- glm(
  is_internal ~ x2,
  family=binomial(), data=data)
print(AIC(tilt_mod3))

# Worse
tilt_mod4 <- glm(
  is_internal ~ x1 + x2,
  family=binomial(), data=data)
print(AIC(tilt_mod4))

# Much better
tilt_mod5 <- glm(
  is_internal ~ x1*x2,
  family=binomial(), data=data)
print(AIC(tilt_mod5))

# Much worse
tilt_mod6 <- glm(
  is_internal ~ bs(x1, df=2, degree=2) + bs(x2, df=2, degree=2),
  family=binomial(), data=data)
print(AIC(tilt_mod6))

# Still bad
tilt_mod7 <- glm(
  is_internal ~ bs(x1, df=3, degree=2) + bs(x2, df=3, degree=2),
  family=binomial(), data=data)
print(AIC(tilt_mod7))

# Really good
tilt_mod8 <- glm(
  is_internal ~ bs(x1, df=3, degree=2)*I(bs(x2, df=3, degree=2)),
  family=binomial(), data=data)
print(AIC(tilt_mod8))

# Only slightly better
tilt_mod9 <- glm(
  is_internal ~ bs(x1, df=3, degree=2)*I(bs(x2, df=3, degree=2)),
  family=binomial(), data=data)
print(AIC(tilt_mod9))

# Starting to get worse
tilt_mod10 <- glm(
  is_internal ~ bs(x1, df=4, degree=2)*I(bs(x2, df=4, degree=2)),
  family=binomial(), data=data)
print(AIC(tilt_mod10))
