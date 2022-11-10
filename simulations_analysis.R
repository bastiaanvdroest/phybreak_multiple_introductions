library(ggplot2)
library(dplyr)
library(scales)
library(phybreak)
library(reshape2)

source("~/Documents/phybreak/R/thin_phybreak.R")

#### Functions ####

errorClassification <- function(sim.inf, s.inf){
  sim.inf.numeric <- sapply(sim.inf, function(x){
    if (x == "index") return(0)
    else return(as.numeric(strsplit(x, "\\.")[[1]][2]))
  })
  sim.trees <- sapply(1:length(sim.inf.numeric), function(i){
    return(tail(phybreak:::.ptr(sim.inf.numeric, i),1))
  })
  names(sim.trees) <- names(sim.inf)

  error.matrix <- matrix(data = 0, nrow = 1, ncol = 6)
  colnames(error.matrix) <- c(paste0(rep(LETTERS[1:2], each = 3), 1:3))#, "merge")
  false.inf <- which(sim.inf != s.inf)
  for (i in false.inf){
    ## Both infectors are not an index: Type 1 error
    if (sim.inf[i] != "index" & s.inf[i] != "index") {
      ## In same cluster: Type A1 error, else Type B1 error
      if (sim.trees[sim.inf[i]] == sim.trees[s.inf[i]]) {
        error.matrix[1,1] <- error.matrix[1,1] + 1
      } else {
        error.matrix[1,4] <- error.matrix[1,4] + 1
      }
      next()
    }
    
    ## If simulated infector is index: Type 2 error
    if (sim.inf[i] == "index") {
      ## In same cluster: Type A2 error, else Type B2 error (merge)
      if (sim.trees[names(sim.inf)[i]] == sim.trees[s.inf[i]]) {
        error.matrix[1,2] <- error.matrix[1,2] + 1
      } else {
        error.matrix[1,5] <- error.matrix[1,5] + 1
      }
      
       next()
    }
    
    ## If estimated infector is index: Type 3 error
    if (s.inf[i] == "index") {
      ## In same cluster: Type A3 error, else Type B3 error (split)
      if (sim.trees[sim.inf[i]] == sim.trees[names(sim.inf)[i]]) {
        error.matrix[1,3] <- error.matrix[1,3] + 1
      } else {
        error.matrix[1,6] <- error.matrix[1,6] + 1
      }
      next()
    }
  }
  return(error.matrix)
}

colorGenerator <- function(n){
  return(hue_pal()(n))
}

#### Data loading ####
# Path to folder with RData files of simulations
datadir = "path/to/folder"

files = list.files(datadir)

################################################################################

d <- do.call(rbind,lapply(files, function(f){
  
  #obsize = strsplit(strsplit(f, "intro")[[1]][1], "obsize")[[1]][2]
  obsize = strsplit(strsplit(f, "i")[[1]][1], "s")[[1]][2]
  numb = strsplit(strsplit(f, "_")[[1]][3], ".RData")[[1]][1]
  #histtime = strsplit(strsplit(f, "t")[[1]][2], "_")[[1]][1]
  #histtime = strsplit(strsplit(f, "t")[[1]][2], "mu")[[1]][1]
  #histtime = strsplit(strsplit(f, "t")[[1]][2], "ob")[[1]][1]
  #histtime = strsplit(strsplit(f, "t")[[1]][2], "p")[[1]][1]
  #sintro = strsplit(strsplit(f, "intro")[[1]][2], "hist")[[1]][1]
  sintro = strsplit(strsplit(f, "i")[[1]][2], "_")[[1]][1]
  #obtime = strsplit(strsplit(f, "obtime")[[1]][2], "\\.")[[1]][1]
  prior = strsplit(f, "_")[[1]][2]
  coalrate = strsplit(strsplit(f, "_")[[1]][3], "wh")[[1]][2]
  #nkeep = 2.5e4
  
  #if (grepl("mupr", sintro)) sintro <- strsplit(sintro, "mupr")[[1]][1]

  intro <- tryCatch({
    load(paste0(datadir,f))
    #s <- thin.phybreak(s, nkeep = nkeep)
    #s <- thin.phybreak(s, nkeep = nkeep)
    #s <- thin.phybreak(s2, nkeep = nkeep)
    #score <- scoreClustering(s, 1)
    
    intro = ifelse( (quantile(s$s$introductions, 0.025) <= as.numeric(sintro) & 
                       quantile(s$s$introductions, 0.975) >= as.numeric(sintro)),
                    1, 0)
    #simLUCA = transphylo2phybreak(sim)$v$nodetimes[as.numeric(obsize)+1]
    #LUCA = sapply(1:nkeep, function(i) min(s$s$nodetimes[,i])) 
    #seccoal = transphylo2phybreak(sim)$v$nodetimes[order(transphylo2phybreak(sim)$v$nodetimes)][2]
    mu = ifelse( (confidenceInterval(s$s$mu, percentile = T)[1] <= 1e-4 & 
                    confidenceInterval(s$s$mu, percentile = T)[2] >= 1e-4),
                 1, 0)
    mG = ifelse( (confidenceInterval(s$s$mG, percentile = T)[1] <= 1 & 
                    confidenceInterval(s$s$mG, percentile = T)[2] >= 1),
                 1, 0)
    mS = ifelse( (confidenceInterval(s$s$mS, percentile = T)[1] <= 1 & 
                    confidenceInterval(s$s$mS, percentile = T)[2] >= 1),
                 1, 0)
    wh.hist = ifelse( (confidenceInterval(s$s$wh.h, percentile = T)[1] <= 50 & 
                         confidenceInterval(s$s$wh.h, percentile = T)[2] >= 50),
                      1, 0)
    wh.slope = ifelse( (confidenceInterval(s$s$wh.s, percentile = T)[1] <= 1 & 
                          confidenceInterval(s$s$wh.s, percentile = T)[2] >= 1),
                       1, 0)
    
    infectors = transtree(s)$infector
    trueinf = sum(sim$sim.infectors == infectors)
    
    index = which(infectors == "index")
    trueindex = which(sim$sim.infectors == "index")
    splits = sum(is.na(match(index,  trueindex)))
    merges = sum(is.na(match(trueindex, index)))
    
    data.frame(intro = intro,
               wh.h = wh.hist,
               mu = mu,
               mG = mG,
               mS = mS,
               wh.s = wh.slope,
               infector = trueinf,
               splits = splits,
               merges = merges)},
    error = function(err){
      print(err)
      data.frame(extseq = NULL, intro = NULL, score = NULL)}
  ) 
  
  
  
  if (nrow(intro)>0)
    return(data.frame(obsize = obsize,
                    number = numb,
                    sim.intro = as.numeric(sintro),
                    prior = prior,
                    coalrate = coalrate,
                    intro))
                    # mu = median(s$s$mu),
                    # mG = median(s$s$mG),
                    # mS = median(s$s$mS),
  else
    return(NULL)

}))

d %>% group_by(prior) %>% summarise(intro = sum(intro), wh.h = sum(wh.h), mu = sum(wh.h), mG = sum(mG), mS = sum(mS),
                                        wh.s = sum(wh.s), infectors = mean(infector), splits = mean(splits), merges = mean(merges))

d_sim.intro <- d %>% group_by(sim.intro, coalrate) %>% summarise(intro = sum(intro), wh.h = sum(wh.h), mu = sum(wh.h), mG = sum(mG), mS = sum(mS),
                                    wh.s = sum(wh.s), infectors = mean(infector), splits = mean(splits), merges = mean(merges))

d_sim.intro$coalrate <- factor(d_sim.intro$coalrate, levels = c(10, 50, 250))

## visualize number of true infectors, splits or merges in groups of coalrate
p1 <- ggplot(d_sim.intro, aes(x = coalrate, y = infectors)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  facet_grid(~ sim.intro) +
  theme_bw()

p2 <- ggplot(d_sim.intro, aes(x = coalrate, y = splits)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  facet_grid(~ sim.intro) +
  theme_bw()
  
p3 <- ggplot(d_sim.intro, aes(x = coalrate, y = merges)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  facet_grid(~ sim.intro) +
  theme_bw()

## plot infectors, splits and merges in one figure
cowplot::plot_grid(p1, cowplot::plot_grid(p2, p3, nrow = 2))

ess <- do.call(rbind,lapply(files, function(f){
  
  obsize = strsplit(strsplit(f, "i")[[1]][1], "s")[[1]][2]
  numb = strsplit(strsplit(f, "_")[[1]][4], ".RData")[[1]][1]
  sintro = strsplit(strsplit(f, "i")[[1]][2], "_")[[1]][1]
  coalrate = strsplit(strsplit(f, "_")[[1]][3], "wh")[[1]][2]
  prior = strsplit(f, "_")[[1]][2]
  nkeep = 2.5e4
  
  intro <- tryCatch({
    load(paste0(datadir,f))
    #s <- thin.phybreak(s, nkeep = nkeep)
    m <- get_mcmc(s)
    #s <- thin.phybreak(s, nkeep = nkeep)
    #s <- thin.phybreak(s2, nkeep = nkeep)
    #score <- scoreClustering(s, 1)
    intro <- coda::effectiveSize(m[,"introductions"])
    mg <- coda::effectiveSize(m[,"mG"])
    ms <- coda::effectiveSize(m[,"mS"])
    r <- coda::effectiveSize(m[,"wh.slope"])
    rh <- coda::effectiveSize(m[,"wh.history"])
    mu <- coda::effectiveSize(m[,"mu"])
    
    data.frame(intro = intro,
               r = r,
               rh = rh,
               mu = mu,
               mg = mg,
               ms = ms)},
    error = function(err){
      print(err)
      data.frame(extseq = NULL, intro = NULL, score = NULL)}
  ) 
  
  
  
  if (nrow(intro)>0)
    return(data.frame(obsize = obsize,
                      number = numb,
                      sim.intro = as.numeric(sintro),
                      intro))
  # mu = median(s$s$mu),
  # mG = median(s$s$mG),
  # mS = median(s$s$mS),
  else
    return(NULL)
  
}))

ess %>% group_by(sim.intro) %>% summarise(intro = mean(intro), mu = mean(mu), r = mean(r), mG = mean(mg), mS = mean(ms),
                                          rhist = mean(rh))
# 
# d <- d %>% filter(history != "Error")

################################################################################
d <- do.call(rbind,lapply(files, function(f){
  print(f)
  obsize = strsplit(strsplit(f, "i")[[1]][1], "s")[[1]][2]
  numb = strsplit(strsplit(f, "_")[[1]][4], ".RData")[[1]][1]
  sintro = strsplit(strsplit(f, "i")[[1]][2], "_")[[1]][1]
  coalrate = strsplit(strsplit(f, "_")[[1]][3], "wh")[[1]][2]
  prior = strsplit(f, "_")[[1]][2]
  
  # Use for analsyis of mink simulations
  # obsize = 63
  # numb = strsplit(strsplit(f, "_")[[1]][5], ".RData")[[1]][1]
  # sintro = strsplit(strsplit(f, "_")[[1]][4], "_")[[1]][1]
  # coalrate = 20
  # prior = NA
  
  intro <- tryCatch({
    load(paste0(datadir,f))
    logLik.tree <- logLik.phybreak(phybreak(sim, use.tree = T, wh.history = as.numeric(coalrate)))
    
    intro = -(as.numeric(sintro) - sum(transtree(s)$infector == "index"))
    mu = mean(s$s$mu) / 1e-4
    mG = mean(s$s$mG)
    mS = mean(s$s$mS)
    wh.hist =  median(s$s$wh.h) /as.numeric(coalrate)
    wh.slope = median(s$s$wh.s) / 1
    print(intro)
    ir = mean(s$s$ir) / (as.numeric(sintro) / max(sim$sample.times))
    logLik =  mean(s$s$logLik) - logLik.tree
    
    infectors = transtree(s)$infector
    trueinf = sum(sim$sim.infectors == infectors)
    
    falseinf = sim$sample.hosts[which(sim$sim.infectors != infectors)]
    if(length(falseinf) > 0){
      inf.sets = infectorsets(s, falseinf, minsupport = 1/21)
      perc.inf = sum(sapply(falseinf, function(x) {
        if(sim$sim.infectors[x] %in% inf.sets[[x]]$infector) return(1)
        else return(0)
      }))
    } else {
      perc.inf = 0
    }
    
    index = which(infectors == "index")
    trueindex = which(sim$sim.infectors == "index")
    splits = ifelse(length(index)==length(trueindex), 0, sum(is.na(match(index,  trueindex))))
    merges = sum(is.na(match(trueindex, index)))
    
    data.frame(intro = intro,
               wh.h = wh.hist,
               mu = mu,
               mG = mG,
               mS = mS,
               ir = ir,
               wh.s = wh.slope,
               infector = trueinf,
               perc.infector = trueinf + perc.inf,
               splits = splits,
               merges = merges,
               logLik = logLik)},
    error = function(err){
      print(err)
      print(f)
      data.frame(extseq = NULL, intro = NULL, score = NULL)}
  ) 
  
  
  
  if (nrow(intro)>0) 
    return(data.frame(obsize = obsize,
                      number = numb,
                      sim.intro = as.numeric(sintro),
                      prior = prior,
                      coalrate = coalrate,
                      intro))

  else
    return(NULL)
  
}))

d_sim.intro <- d %>% group_by(sim.intro, coalrate) %>% 
  summarise(intro = mean(sim.intro) + mean(intro), wh.h = mean(wh.h), mu = mean(mu), mG = mean(mG), mS = mean(mS), ir = mean(ir), wh.s = mean(wh.s), 
            infectors = mean(infector), perc.infector = mean(perc.infector), splits = mean(splits), merges = mean(merges), logLik = mean(logLik))
d_sim.intro$coalrate <- factor(1/as.numeric(d_sim.intro$coalrate), levels = c(1/10,1/50,1/250))
# For mink simulations:
# d_sim.intro$coalrate <- factor(1/as.numeric(d_sim.intro$coalrate), levels = 1/20)

## Bar plot of number of introductions
di <- ggplot(d_sim.intro, aes(x = coalrate)) +
  geom_bar(aes(y=intro), stat = "identity") +
  geom_hline(aes(yintercept = sim.intro), linetype = "dashed") +
  facet_grid(~ sim.intro, scales = "free_y") +
  theme_bw() +
  labs(x = "Coalescent rate of history host", y = "Estimated median\nnumber of introductions") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Bar plot of number of true infectors and number of infectors in 95% CI
ds <- ggplot(d_sim.intro, aes(x = coalrate)) +
  geom_bar(aes(y=infectors/20 * 100), stat = "identity") +
  geom_bar(aes(y=perc.infector / 20 * 100), stat = "identity", color = "grey50", alpha = 0.3) +
  geom_hline(yintercept = 95, linetype = "dashed") +
  facet_grid(~ sim.intro) + 
  theme_bw() +
  labs(x="Coalescent rate of history host", y = "Correctly identified\ninfectors (%)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


################################################################################

d.error <- do.call(rbind,lapply(files, function(f){
  
  obsize = strsplit(strsplit(f, "i")[[1]][1], "s")[[1]][2]
  numb = strsplit(strsplit(f, "_")[[1]][4], ".RData")[[1]][1]
  sintro = strsplit(strsplit(f, "i")[[1]][2], "_")[[1]][1]
  coalrate = strsplit(strsplit(f, "_")[[1]][3], "wh")[[1]][2]
  prior = strsplit(f, "_")[[1]][2]
  
  # For mink simulations:
  # obsize = 63
  # numb = strsplit(strsplit(f, "_")[[1]][5], ".RData")[[1]][1]
  # sintro = strsplit(strsplit(f, "_")[[1]][4], "_")[[1]][1]
  # coalrate = 20
  # prior = NA
  
  intro <- tryCatch({
    load(paste0(datadir,f))
    infectors = transtree(s)$infector
    trueinf = sum(sim$sim.infectors == infectors)
    
    df <- data.frame(errorClassification(sim$sim.infectors, infectors),
    trueinf = trueinf)
    
    data.frame(df, total = rowSums(df))},
    error = function(err){
      print(err)
      data.frame(extseq = NULL, intro = NULL, score = NULL)}
  ) 
  
  if (nrow(intro)>0)
    return(data.frame(obsize = obsize,
                      number = numb,
                      sim.intro = as.numeric(sintro),
                      prior = prior,
                      coalrate = coalrate,
                      intro))
  else
    return(NULL)
  
}))

d.error.summary <- d.error %>% group_by(sim.intro, coalrate) %>% 
  summarise(SCC = mean(A1), SHC = mean(A2), SCH = mean(A3),
            MCC = mean(B1), MHC = mean(B2), MCH = mean(B3),
            merge = mean(merge), trueinf = mean(trueinf))
d.error.summary <- data.frame(d.error.summary)

d.error.plot <- melt(d.error.summary, id.vars = c("sim.intro", "coalrate"))
d.error.plot$coalrate <- factor(1/as.numeric(d.error.plot$coalrate), levels = c(1/10,1/50,1/250))
d.error.plot$value <- d.error.plot$value/20*100
de <- ggplot(d.error.plot, aes(x = coalrate, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  facet_grid(~ sim.intro) +
  scale_fill_manual(labels = c("S:C2C", "S:H2C", "S:C2H",
                               "M:C2C", "M:H2C", "M:C2H", 
                               "No error"),
                    values = c(SCC = "dark red", SHC = "red", SCH = "salmon",
                               MCC = "dark blue", MHC = "royalblue", MCH = "sky blue",
                               trueinf = "dark grey")) +
  labs(fill = "Error type", y = "Identified infectors (%)", x = "Coalescent rate of history host")+
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

cpg <- cowplot::plot_grid(di, ds, nrow = 2, labels = c("A", "B"), label_size = 30, align = "v")  
cowplot::plot_grid(cpg, de, ncol = 2, labels = c("", "C"),  label_size = 30, rel_widths = c(3/7,4/7))

ggsave("~/Documents/phd_files/nertsen/artikel/figures/error_typing_coalrate_sm.eps", device=cairo_ps, width = 17.8, height = 9)