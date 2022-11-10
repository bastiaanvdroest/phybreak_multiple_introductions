library(readxl)
library(phybreak)

options(stringsAsFactors = FALSE)

# phydata: Data frame with sample name, host name, sample date and sequence
# phydata_culling: Data frame including culling dates

#### Make phydata object ####
phydata <- read.csv(phydata, header = T)
phydata_culling <- as.data.frame(read_excel(phydata_culling))
coord <- as.data.frame(read_excel(coord))

phydata$Time <- as.Date(phydata$Time,"%d-%m-%Y")

### Match locations with farms
excel <- as.data.frame(read_excel(overview, sheet = "Overzicht_Bedrijven"))
clusters <- excel
excel <- do.call(rbind,apply(excel, 1, function(x){
  x <- as.data.frame(t(x))
  lat <- coord$latitude[which(coord$postcode == x$`Postcode:`)]
  lat <- lat/10^(nchar(lat)-2)
  lon <- coord$longitude[which(coord$postcode == x$`Postcode:`)]
  lon <- lon/10^(nchar(lon)-1)
  
  farm <- unlist(strsplit(as.character(x$`Farm:`), "NB"))
  if(farm[2]=="1A") farm[2] <- "1"
  farm <- paste("NB-EMC-",farm[2],sep="")
  return(data.frame(farm=farm,
                    latitude=lat,
                    longitude=lon))
}))
phydata$Latitude <- excel$latitude[match(phydata$Host, excel$farm)]
phydata$Longitude <- excel$longitude[match(phydata$Host, excel$farm)]

### Fill in missing diagnostic date and add culling date

phydata <- phydata[!is.na(phydata$Host),]
phydata <- do.call(rbind,apply(phydata, 1, function(x){
  x <- data.frame(t(x))
  if (is.na(x$Time)){
    host <- unlist(strsplit(as.character(x$Host), "-"))
    host <- paste(host[c(1,3)],collapse="")
    if (host == "NBun") {
      x$Culling <- NA
      return(data.frame(x))
    } else {
      t <- phydata_culling$Diagnose[which(phydata_culling$NB == host)]
      x$Time <- unlist(strsplit(as.character(t), " "))[1]
      x$Culling <- phydata_culling$Ruimingsdatum[which(phydata_culling$NB == host)]
    }
  } else {
    host <- unlist(strsplit(x$Host, "-"))
    host <- paste(host[c(1,3)],collapse="")
    if (host == "NB1") host <- "NB1a"
    else if (host == "NBLosloop") host <- "NB69"
    t <- phydata_culling$Diagnose[which(phydata_culling$NB == host)]
    x$Time <- unlist(strsplit(as.character(t), " "))[1]
    x$Culling <- phydata_culling$Ruimingsdatum[which(phydata_culling$NB == host)]
  }
  return(x)
}))

phydata <- phydata[!is.na(phydata$Time),]
phydata$Time <- as.Date(phydata$Time)

phydata <- phydata[!is.na(phydata$Culling),]

seqs <- do.call(rbind, strsplit(phydata$Sequence, ""))
rownames(seqs) <- phydata$Sample

### Remove N's and dels from sequences
Ns <- unlist(lapply(phydata$Sequence, function(x){
  ns <- as.numeric(gregexpr("N", x)[[1]])
  #dels <- as.numeric(gregexpr("-", x)[[1]])
  return(sort(c(ns[ns>0])))#, dels[dels>0])))
}))

Ns <- sort(unique(Ns))
seqs2 <- seqs[,-Ns]

culling.times <- phydata[,c(2,7)][!duplicated(phydata$Host),]
phydata <- phybreakdata(sequences = seqs2,
                        sample.times = phydata$Time, 
                        sample.names = phydata$Sample, 
                        culling.times = phydata$culling,
                        host.names = phydata$Host)

culling.times <- as.Date(culling.times[match(phydata$sample.hosts[1:61], culling.times$Host),][,2])
assign("culling.times", get("culling.times", environment()), phybreak:::userenv)

# Remove everything from global environment to save space
rm(list=setdiff(ls(), c("phydata", "clusters","culling.times")))

#### Sample phybreak ####

phyb <- phybreak(phydata, wh.history = 20,
                 prior.sample.mean.mean = 10,
                 prior.wh.history.mean = 20, prior.wh.history.shape = 3,
                 prior.intro.rate = 6/190, prior.intro.rate.shape = 3,
                 trans.model = "user", infectivity_file = "~/Documents/phd_files/nertsen/infectivity_function.R") 

s <- sample_phybreak(phyb, nsample = 2e4, nchains = 1, parallel=F)

## Color farms according to clusters from Lu et al. (2021)
library(scales)
colorGenerator <- function(n){
  return(hue_pal()(n))
}

colors <- colorGenerator(5)
excel <- read_excel(clusterinfofile)
clusters <- do.call(rbind,lapply(unique(s$d$hostnames), function(x){
  if(x == "External") {
    cl <- "unknown"
  } else {
    x <- unlist(strsplit(x, "-"))[3]
    if(x == "NB1") x = "NB1A"
    cl <- excel$`Mink cluster`[which(excel$`FarmID`==x)]
    if (is.na(cl)) cl <- ""
  }
  print(cl)
  if (cl == "A/D") color <- colors[1]

  if (cl == "A") color <- colors[1]
  if (cl == "B") color <- colors[2]
  if (cl == "C") color <- colors[3]
  if (cl == "D") color <- colors[4]
  if (cl == "E") color <- colors[5]
  if (cl == "F") color <- colors[6]
  if (cl == "G") color <- colors[7]
  if (cl == "unknown") color <- "black"

  return(data.frame("cluster"=cl, "color"=color))
}))

# Plot the transmission tree of the maximum parent credibility (MPC) tree
plotTrans(s, plot.which = "mpc", arrow.lwd = 3, arrow.col = brewer.pal(5,"Reds"), 
          label.col = clusters$color, label.cex = 1.1, arrow.length = 0.15, 
          axis.cex = 1.2, title.cex = 1.2)

# Plot the phylogenetic tree of the MPC tree
plotPhylo(s, plot.which = "mpc")

# Plot the phylogenetic and transmission tree of the MPC tree
plotPhyloTrans(s, plot.which = "mpc", tree.col = clusters$color)