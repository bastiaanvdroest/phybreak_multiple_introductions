y <- c()
p <- list(
  gen.init = 1e-4,
  gen.sample = 1,
  gen.growth = 18.87829,
  gen.culling = 100
)

a <- (1-p$gen.init)/p$gen.init
r <- p$gen.growth
S <- p$gen.sample
C <- p$gen.culling

for (t in seq(0,50,by=0.1)){
  if (t < 10){
    y <- c(y, 1/(1+a*exp(-r*t)))
  } else if (t >=10 & t < 35){
    y <- c(y, S*1/(1+a*exp(-r*t)))
  } else if (t >= 35){
    y <- c(y, S*1/(1+a*exp(-r*35)) * exp(-C*(t-35)))
  }
}
plot(seq(0,50,by=0.1),y,xlab="Days",ylab="Transmission probability",main="Infectiousness Function")
plot(y[350:450])
