time <- Sys.time()

library(phybreak)
library(parallel)
library(Rdsm)

# Get arguments from the command line
args = commandArgs(trailingOnly=TRUE)
id = args[1]
dir = args[2]
intro = as.numeric(args[3])
n = as.numeric(args[4])
wh = as.numeric(args[5])

# Simulate outbreak
sim <- sim_phybreak(obsize = n, introductions = intro, wh.history = wh)
names(sim$sequences) <- 1:n

# Create phybreak object of simulated sequences and sampling times
phy <- phybreak(sim, prior.mu.sd = 5e-5, prior.sample.mean.sd = 0.5,
                prior.wh.history.mean = 50, prior.wh.history.shape = 3)

s <- sample_phybreak(phy, nsample= 2.5e4, nchains = 3)
save(sim, s, file = sprintf("path/to/save/results.RData", dir, intro, id))

print(sprintf("Duration of run: %s", Sys.time() - time))