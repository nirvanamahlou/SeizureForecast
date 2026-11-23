pkgs <- c("Kendall", "moments", "tseries","urca","forecast","aTSA","seastests")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
values <- c(
  112,118,132,129,121,135,148,148,136,119,104,118,   # 1389
  115,126,141,135,125,149,170,170,158,133,114,140,   # 1390
  145,150,178,163,172,178,199,199,184,162,146,166,   # 1391
  171,180,193,181,183,218,230,242,209,191,172,194,   # 1392
  196,196,236,235,229,243,264,272,237,211,180,201,   # 1393
  204,188,235,227,234,264,302,293,259,229,203,229,  # 1394
  242,233,267,269,270,315,364,347,312,274,237,278,  # 1395
  284,277,317,313,318,374,413,405,355,306,271,306,   # 1396
  315,301,356,348,355,422,465,467,404,347,305,336,   # 1397
  340,318,362,348,363,435,491,505,404,359,310,337,   # 1398
  360,342,406,396,420,472,548,559,463,407,362,405,   # 1399
  417,391,419,461,472,535,622,606,508,461,390,432    # 1400
)


ts <- ts(values, start = c(1389, 1), frequency = 12)

plot(ts_data,
     xlab = "year",
     ylab = "number of patients",
     main = "سری زمانی تعداد بیماران تشنج حاد (ماهانه)",
     col  = "steelblue",
     lwd  = 2)
grid()

months   <- cycle(ts)
kruskal  <- kruskal.test(ts ~ as.factor(months))
print(kruskal)

fit_stl <- stl(ts, s.window="periodic")
seasonal  <- fit_stl$time.series[,"seasonal"]
kruskal.test(seasonal  ~ as.factor(cycle(ts)))

mk_res <- MannKendall(ts)
cat("Mann–Kendall:\n"); print(mk_res)

stationary.test(ts)

t      <- 1:length(ts)
df     <- data.frame(y = as.numeric(ts), t = t)
fit_lin  <- lm(y ~ t,          data = df)
fit_quad <- lm(y ~ t + I(t^2), data = df)

cat("\n--- Linear model summary ---\n")
print(summary(fit_lin))

cat("\n--- Quadratic model summary ---\n")
print(summary(fit_quad))

stl  <- stl(ts, s.window = "periodic")
remainder <- stl_res$time.series[ , "remainder"]

cat("\nLjung–Box on STL remainder:\n")
print(Box.test(remainder, lag = 20, type = "Ljung-Box", fitdf = 0))

y      <- as.numeric(ts)
mu     <- mean(y)
sigma2 <- var(y)
acfs   <- acf(y,  plot = FALSE)$acf   
pacfs  <- pacf(y, plot = FALSE)$acf   
sk     <- skewness(y)
ku     <- kurtosis(y)

cat("\nDescriptive Statistics:\n")
cat(sprintf("  Mean:     %.2f\n", mu))
cat(sprintf("  Variance: %.2f\n", sigma2))
cat(sprintf("  ACF1:     %.3f\n", acfs[2]))
cat(sprintf("  PACF1:    %.3f\n", pacfs[1]))
cat(sprintf("  Skewness: %.3f\n", sk))
cat(sprintf("  Kurtosis: %.3f\n", ku))

acf(ts,
    main = "ACF ",
    xlab = "Lag",
    ylab = "ACF")


pacf(ts,
     main = "PACF ",
     xlab = "Lag",
     ylab = "PACF")


adf_drift <- ur.df(ts, type = "drift", selectlags = "AIC")
summary(adf_drift)

adf_trend <- ur.df(ts, type = "trend", selectlags = "AIC")
summary(adf_trend)
 
fit<- auto.arima(ts) 
summary(fit)
checkresiduals(fit)     

cat("Mean of residuals =", mean(residuals(fit)), "\n")
shapiro.test(residuals(fit))

lag_lb <- 24
lb <- Box.test(residuals(fit), lag = lag_lb,
               type = "Ljung-Box",
               fitdf = length(fit$coef))  
print(lb)


coefs <- coef(fit)
ar_ns <- coefs[grep("^ar", names(coefs))]     # φ₁, φ₂, …
ma_ns <- coefs[grep("^ma", names(coefs))]     # θ₁, θ₂, …


roots_ar_ns <- polyroot(c(1, -ar_ns))
roots_ma_ns <- polyroot(c(1,  ma_ns))

cat("Non‐seasonal AR roots (Modulus):\n")
print(Mod(roots_ar_ns))
cat("Non‐seasonal MA roots (Modulus):\n")
print(Mod(roots_ma_ns))


sar <- coefs[grep("^sar", names(coefs))]    
sma <- coefs[grep("^sma", names(coefs))]     

roots_sar <- polyroot(c(1, -sar))
roots_sma <- polyroot(c(1,  sma))

cat("Seasonal AR roots (Modulus):\n")
print(Mod(roots_sar))
cat("Seasonal MA roots (Modulus):\n")
print(Mod(roots_sma))





n.ahead <- 12
pred <- predict(fit, n.ahead = n.ahead)

ts_pred <- ts(pred$pred,
              start = c(end(ts)[1], end(ts)[2] + 1),
              frequency = frequency(ts))

ts.plot(ts, ts_pred,
        col = c("black", "red"),
        lty = c(1, 2),
        main = " predict",
        xlab = "year",
        ylab = "number of patients")
legend("topleft",
       legend = c("original","forecast"),
       col = c("black","red"),
       lty = c(1,2))
grid()


final <- Arima(ts,
                   order            = c(2,1,1),
                   seasonal         = c(0,1,0),
                   include.constant = FALSE)


fc <- forecast::forecast(final, h = 12, level = c(80, 95))

print(fc)
plot(fc, main="predict", xlab="year", ylab="number of patients",shadecols = c("lightblue","lightgray"))
grid()
