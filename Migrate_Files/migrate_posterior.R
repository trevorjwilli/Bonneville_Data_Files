#!/usr/bin/env Rscript --vanilla

####################
#-Posterior Graphs-#
####################

# This code uses the bayes file from Migrate to recreate the posterior graphs
# using ggplot2

# To run this code, copy to the directory with the bayesfile and run 
# 'Rscript migrate_posterior.R bayesfile'

args = commandArgs(trailingOnly = TRUE)

library(dplyr)
library(ggplot2)
library(gridExtra)

data <- read.delim(args[1] ,sep = " ",header = FALSE, comment.char = "#")
colnames(data) <- c("Locus", "Parameter", "HPC50", "HPC95", "Value", "Count", "Frequency", "Cum_Freq", "PriorFreq")

params <- scan(args[1], what = "list", sep = "\n")
params <- params[grep("#@", params)]
params <- unlist(strsplit(params, split= "[[:blank:]]+"))
params <- matrix(params, ncol = 3, byrow = T)

plotlist <- list()
for(i in 1:nrow(params)) {
  param <- data %>% filter(Parameter == as.numeric(params[i,2]), Locus == length(unique(data$Locus)))
  param$HPC <- factor(ifelse(param$HPC50 == 1 & param$HPC95 == 1, "50", ifelse(param$HPC95 == 1 & param$HPC50 == 0, "95", "100")), levels = c("100", "95", "50"))

  graph <- ggplot(param, aes(x = Value, y = Frequency)) +
    geom_bar(stat = "identity", aes(fill = HPC, col = HPC)) +
    scale_fill_manual(values = c("gray", "gray55", "black")) +
    scale_color_manual(values = c("gray", "gray55", "black")) +
    geom_hline(yintercept = param[1,"PriorFreq"], col = "red") +
    xlab(params[i,3]) +
    theme_classic() +
    theme(legend.position = "none")
  plotlist[[i]] <- graph
}

grid.arrange(grobs = plotlist, ncol = 2)

