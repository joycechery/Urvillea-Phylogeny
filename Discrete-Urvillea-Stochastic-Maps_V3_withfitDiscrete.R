#function to computer posterior probabilities at the nodes (L49-55) by liam revell from http://blog.phytools.org/2012/12/plotting-node-piecharts-on-top-of.html
#function to plot simmaps and node posterior probabilities (L63-76--"Plot_simmap") by Dr. Michael May (UC Berkeley)

#this code plots stochastic character maps along the branches of the phylogeny and piecharts at the nodes
library(phytools)
library(grDevices)
library(geiger)

setwd("~/Desktop/")
tree <-read.tree("Tree.tre")

#tree$tip.label
tree <- reroot(tree, 1, interactive = TRUE)

# upload data 
data  <- read.csv("Data.csv", row.names = 1) 
datum <- as.data.frame(cbind(rownames(data), data[,1]), row.names = FALSE)
datum <- datum[complete.cases(datum), ]

# species
species <- datum$V1

# reformat data
datas <- datum$V2
names(datas) <- species

# drop tips
pruned_tree <- keep.tip(tree, species)
pruned_tree=ladderize(pruned_tree)


####character transition model
fitER<-fitDiscrete(pruned_tree,datas,model= "ER")
fitSYM<-fitDiscrete(pruned_tree,datas,model= "SYM")
fitARD<-fitDiscrete(pruned_tree,datas,model= "ARD")

#Compare the AIC scores from the model outputs below, picking the model with the lowest AIC
fitER$opt$aic
fitSYM$opt$aic
fitARD$opt$aic

# stochastic mapping, use the model with the lowest AIC score
simmap.trees <- make.simmap(pruned_tree, datas, model = "ER", pi = "estimated", nsim=1000)

# set colors & character states
col_vec <- c("dodgerblue4", "orangered3","gold","maroon3", "orangered2", "steelblue3")
states<- c("Ontogeny_1", "Ontogeny_2", "Ontogeny_3",  "Ontogeny_4", "Ontogeny_5", "Ontogeny_6") 

# function to compute the node states
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

XX<-sapply(simmap.trees,foo)
pies<-t(apply(XX,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=1000))

#generate summary of stochastic maps with pies of posterior at nodes..
#see below (line 72) for plot_simmap function

plot_simmap(time_tree = simmap.trees[[1]], 
            tree = simmap.trees[[1]], 
            simmaps = simmap.trees, 
            states = states,
            show.tip.label = T,
            lwd = 8,
            label.cex = .8,
            label.offset=.015,
            colors = col_vec, edge.width=0.1, nt=10001)

add.simmap.legend(colors=col_vec, prompt=FALSE,x=105,y=180, fsize =.9)
nodelabels(pie=pies,cex=0.8,piecol=col_vec, lwd=1)
legend(x=0,y=1, legend = states, col = col_vec, pch = 20, yjust = 0, bty = 'n', cex =.7, pt.cex = 4)


