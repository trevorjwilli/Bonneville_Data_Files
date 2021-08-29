###########################
#-Migrate_probabilities.R-#
###########################

prob_table <- function(dir) {
  
  getval <- function(x) {
    vals <- x[grep('All\\s+-\\d+', x, ignore.case = FALSE)]
    vals <- unlist(strsplit(vals, '\\s+'))
    val <- as.numeric(vals[4])
    val
  }
  
  p_mod <- function(vec) {
    pval <- exp(vec - max(vec))/sum(exp(vec - max(vec)))
    bf <- vec - max(vec)
    out <- data.frame(logmL = vec, LBF = bf, prob = pval)
    out
  }
  
  infiles <- list.files(dir, full.names = TRUE)
  print(infiles)
  models <- list.files(dir)
  infiles2 <- lapply(infiles, readLines)
  modvals <- lapply(infiles2, getval)
  
  out <- p_mod(unlist(modvals))
  out <- data.frame(models, out)
  out
}

aliciae <- prob_table("outfiles_between/outfiles_aliciae")
aliciae <- data.frame(species = rep('Lepidomeda aliciae', 4), aliciae)

ardens <- prob_table("outfiles_between/outfiles_ardens")
ardens <- data.frame(species = rep('Catostomus ardens', 4), ardens)

atraria <- prob_table("outfiles_between/outfiles_atraria")
atraria <- data.frame(species = rep('Gila atraria', 4), atraria)

bairdii <- prob_table("outfiles_between/outfiles_bairdii")
bairdii <- data.frame(species = rep('Cottus bairdii', 4), bairdii)

balteatus <- prob_table("outfiles_between/outfiles_balteatus")
balteatus <- data.frame(species = rep('Richardsonius balteatus', 4), balteatus)

osculus <- prob_table("outfiles_between/outfiles_osculus")
osculus <- data.frame(species = rep('Rhinichthys osculus', 4), osculus)

phlegethontis <- prob_table("outfiles_between/outfiles_phlegethontis")
phlegethontis <- data.frame(species = rep('Iotichthys phlegethontis', 4), phlegethontis)

platyrhynchus <- prob_table("outfiles_between/outfiles_platyrhynchus")
platyrhynchus <- data.frame(species = rep('Catostomus platyrhynchus', 4), platyrhynchus)

alldata <- rbind(ardens, platyrhynchus, bairdii, atraria, phlegethontis, aliciae, osculus, balteatus)

write.csv(alldata, "modsel_between.csv")

williamsoni <- prob_table("williamsoni/outfiles_williamsoni")
williamsoni
