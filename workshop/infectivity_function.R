### Define extra data needed
#culling.times = 

### Define parameters
parameters <- list("trans.init" = 1e-4,
                   "trans.culling" = 5,
                   "trans.growth" = 1,
                   "trans.sample" = 1)

### Define which parameters should be estimated
helpers <- list("est.tG" = TRUE,
                "est.tS" = F)

### Make room in sampleslot for estimated parameters
samplers <- list("tG" = c(),
                 "tS" = c())

### Define the update steps
update_tG <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  p <- pbe1$p
  
  ### change to proposal state
  p$trans.growth <- exp(log(p$trans.growth) + rnorm(1, 0, 0.1))
  ### update proposal environment
  copy2pbe1("p", le)
  ### calculate likelihood
  propose_pbe("mG")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikgen - pbe0$logLikgen
  
  ### accept
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("mG")
  }
}

update_tS <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  p <- pbe1$p
  v <- pbe1$v
  
  p$trans.sample <- exp(log(p$trans.sample) + rnorm(1, 0, 0.1))
  # 
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate likelihood
  propose_pbe("mG")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikgen - pbe0$logLikgen
  
  ### accept
  if (runif(1) < exp(logaccprob) & p$trans.sample <= 1) {
    accept_pbe("mG")
  }
}


updaters <- list(update_tG, update_tS)

### Generation time interval
infect_function <- function(time, inftimes, le, nodetimes, 
                            host, log = FALSE,
                            test.arguments = FALSE){
  
  p <- le$p
  v <- le$v
  
  if (is.null(culling.times)) {
    stop("culling times of hosts must be provided")
  } else if(class(culling.times[1]) == "Date"){
    culling.times <- as.numeric(difftime(culling.times, reference_date))
  }
  
  if(test.arguments) return()
  
  if(is.null(p$trans.init))
    stop("initial fraction infected is missing")
  if(is.null(p$trans.growth))
    stop("growth factor of infectiousness is missing")
  if(is.null(p$trans.sample))
    stop("reduction factor after first positive sample is missing")
  if(is.null(p$trans.culling))
    stop("decay factor after culling is missing")
  
  a <- (1-p$trans.init)/p$trans.init
  r <- p$trans.growth
  S <- p$trans.sample
  C <- p$trans.culling
  
  # Calculate normalization factor by calculating mean AUC of infectiousness function
  AUCs <- unlist(lapply(1:length(v$inftimes), function(i){
    samtime = v$nodetimes[i] - v$inftimes[i]
    cultime = culling.times[i] - v$inftimes[i]
    if (r*samtime < 100){
      probs = sum((log(a+exp(r*samtime)) - log(a+1)) / r,
                  S * ( log(a+exp(r*cultime)) - log(a+exp(r*samtime)) ) / r,
                  (S / (1 + a*exp(-r*cultime))) / C)
    } else {
      probs = sum((r*samtime - log(a+1)) / r,
                  S * ( r*(cultime - samtime) ) / r,
                  S / C)
    }
    return(probs)
  }))
  norm_factor <- 1/mean(AUCs)
  
  # Use culling times of infectors is rest of calculations
  cultimes <- culling.times[match(inftimes, v$inftimes)]
  
  samtimes <- as.numeric(nodetimes - inftimes)
  cultimes <- as.numeric(cultimes - inftimes)
  hosttimes <- as.numeric(time - inftimes)
  
  if (length(hosttimes) == 0){
    probs = 1
  } else {
  
    if(is.null(host)){
      if(length(hosttimes) != length(nodetimes)){
        probs <- 0.1
        j <- 1
      } else {
        probs <- c()
        j <- 0
      }
      for (i in 1:length(samtimes)){
        if(hosttimes[i+j] < 0)
          probs <- c(probs, 0)
        else if(hosttimes[i+j] < samtimes[i])
          probs <- c(probs, 1/(1+a*exp(-r*hosttimes[i+j])))
        else if(hosttimes[i+j] >= samtimes[i] & hosttimes[i+j] < cultimes[i])
          probs <- c(probs, S/(1+a*exp(-r*hosttimes[i+j])))
        else if(hosttimes[i+j] >= cultimes[i] & hosttimes[i+j] < cultimes[i] + 5)
          probs <- c(probs, S/(1+a*exp(-r*cultimes[i])) * exp(-C*(hosttimes[i+j]-cultimes[i])))
        else 
          probs <- c(probs, 0)
      }
    }
  }
  
  if(log)
    return(log(probs*norm_factor))
  else
    return(probs*norm_factor)
}
