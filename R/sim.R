# simulate data
ts = arima.sim(n = 63, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488), d = 1),
               sd = sqrt(0.1796))

ts = arima.sim(n = 63, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
               sd = sqrt(0.1796))

library(TSA)
ts.sim <- arima.sim(list(order = c(1, 1, 2), ar = 0.7, ma = c(0.5, -0.2)), n = 1000)
eacf(diff(ts.sim))
ts.plot(ts.sim)

arima(ts.sim, order = c(1, 1, 2))
