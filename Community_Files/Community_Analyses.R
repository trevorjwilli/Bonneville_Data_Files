#devtools::install_github("trevorjwilli/CommUtilsR")
library(vegan)
library(sf)
library(dplyr)
library(CommUtilsR)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(cowplot)
library(bipartite)
library(adespatial)

cscore_scale <- function(x) {
  x <- ifelse(x > 0, 1, 0)
  score <- mean(designdist(t(x), method = "((A-J)*(B-J))/(A*B)", terms = "binary", name = "C-score"))
  score
}

data <- read.csv("Occurrence.csv")
loc <- read.csv("Localities.csv", stringsAsFactors = F)

occ <- sxsmat(data, c("Species", "site"))

occ <- occ[,-which(colSums(occ) == 1)]

loc <- loc %>% filter(site %in% rownames(occ))

occ <- occ[order(rownames(occ)),]
rownames(loc) <- loc$site
loc <- loc[order(rownames(loc)),]

##### Null Model #####

### Full Matrix ###

set.seed(352)
null <- oecosimu(comm = occ, nestfun = cscore_scale, method = "curveball", nsimul = 5000, burnin = 500, alternative = "two.sided")
null

# GSL #

GSL_sites <- loc$site[which(loc$basin == "GSL")]

occ_gsl <- occ[which(rownames(occ) %in% GSL_sites),]
occ_gsl <- occ_gsl[,-which(colSums(occ_gsl) <= 1)]

set.seed(82)
null_gsl <- oecosimu(comm = occ_gsl, nestfun = cscore_scale, method = "curveball", nsimul = 5000, burnin = 500, alternative = "two.sided")
null_gsl

# GSLD #

GSLD_sites <- loc$site[which(loc$basin == "GSLD")]

occ_gsld <- occ[which(rownames(occ) %in% GSLD_sites),]
occ_gsld <- occ_gsld[,-which(colSums(occ_gsld) <= 1)]

set.seed(64)
null_gsld <- oecosimu(comm = occ_gsld, nestfun = cscore_scale, method = "curveball", nsimul = 5000, burnin = 500, alternative = "two.sided")
null_gsld

# SEV #

SEV_sites <- loc$site[which(loc$basin == "SEV")]

occ_sev <- occ[which(rownames(occ) %in% SEV_sites),]
occ_sev <- occ_sev[,-which(colSums(occ_sev) <= 1)]

set.seed(9872)
null_sev <- oecosimu(comm = occ_sev, nestfun = cscore_scale, method = "curveball", nsimul = 5000, burnin = 500, alternative = "two.sided")
null_sev

##### Probabilistic #####

library(cooccur)

cooc <- cooccur(t(occ), spp_names = TRUE)
cooc
plot(cooc)

test <- plot(cooc)
test <- test$data

levs <- c("G. atraria", "R. cataractae", "R. osculus", "C. ardens", "C. bairdii", "L. aliciae", "C. platyrhynchus", "R. balteatus")
test$sp1 <- gsub('(\\w)\\w+ (\\w+)', '\\1. \\2', test$X1)
test$sp2 <- gsub('(\\w)\\w+ (\\w+)', '\\1. \\2', test$X2)

test$sp1 <- factor(test$sp1, levels = levs)
test$sp2 <- factor(test$sp2, levels = levs)

ggplot(test, aes(x = sp1, y = sp2, fill = as.factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(labels = c("Negative", "Random", "Positive"), values = c("gold2", "grey", "lightskyblue2")) +
  scale_x_discrete(guide = guide_axis(angle = 45), 
                   labels = c(expression(italic("R. cataractae")),
                              expression(italic("R. osculus")),
                              expression(italic("C. ardens")),
                              expression(italic("C. bairdii")),
                              expression(italic("L. aliciae")),
                              expression(italic("C. platyrhynchus")),
                              expression(italic("R. balteatus")))) +
  scale_y_discrete(labels = c(expression(italic("G. atraria")),
                              expression(italic("R. cataractae")),
                              expression(italic("R. osculus")),
                              expression(italic("C. ardens")),
                              expression(italic("C. bairdii")),
                              expression(italic("L. aliciae")),
                              expression(italic("C. platyrhynchus")),
                              expression(italic("R. balteatus")))) +
  theme_minimal() +
  theme(axis.line.x = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10))

##### Nestedness #####

set.seed(232)
nested <- oecosimu(comm = occ, nestfun = nestednodf, method = "curveball", nsimul = 5000, burnin = 500, alternative = "two.sided")
nested

arra <- nestedrank(occ, return.matrix = T)
arra
nestmat <- arra$nested.matrix
locmatch <- loc[order(match(loc$site, rownames(nestmat))),]
locmatch
colnames(nestmat) <- gsub('(\\w)\\w+ (\\w+)', '\\1. \\2', colnames(nestmat))
heatmap(nestmat, Rowv = NA, Colv = NA, col = c("white", "grey"),
        distfun = NA, hclustfun = NA, reorderfun = NA, revC = T, 
        scale = "none", labRow = locmatch$basin, labCol = as.expression(lapply(colnames(nestmat), function(x) bquote(italic(.(x))))), margins = c(8,5))

##### PERMANOVA #####

bray <- vegdist(occ, binary = T)

set.seed(6)
pnova <- adonis2(bray~basin, data = loc)
pnova

##### NMDS #####

set.seed(87078)
nmds <- metaMDS(occ)
nmds
nmds$points

scrs <- as.data.frame(scores(nmds, display = "sites"))
scrs <- cbind(scrs,loc)

spscrs <- as.data.frame(scores(nmds, display = "species"))
spscrs$sp <- gsub('(\\w)\\w+ (\\w+)', '\\1. \\2', rownames(spscrs))

hull <- scrs %>%
  group_by(basin) %>%
  slice(chull(NMDS1, NMDS2))
scrs

ggplot(data = scrs,mapping = aes(x = NMDS1, y = NMDS2)) +
  geom_point(shape = 21, size = 2, aes(fill = basin)) +
  geom_polygon(data = hull, alpha = 0.5, aes(fill = basin)) +
  scale_fill_jco(name = "Basin") +
  geom_segment(data = spscrs, aes(x=0, xend = NMDS1, y=0, yend=NMDS2),
  arrow = arrow(length = unit(0.25, "cm")), alpha = 0.5) +
  geom_text_repel(data = spscrs, aes(x = NMDS1, y = NMDS2, label = sp, fontface = "italic"), size = 3) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12))

### Spatial ###

set.seed(8282)
cand <- listw.candidates(loc[,4:5], nb = c('del', 'gab', 'mst'), weights = c('binary', 'flin'))
sel_w <- listw.select(decostand(occ, "hellinger"), cand, MEM.autocor = "positive")
spat <- sel_w$best$MEM.select
sel_w
mod1 <- rda(decostand(occ, "hellinger") ~ ., data = spat)
anova(mod1)
RsquareAdj(mod1)

loc_gsl <- loc %>% filter(basin == "GSL")
loc_gsld <- loc %>% filter(basin == "GSLD")
loc_sev <- loc %>% filter(basin == "SEV")

set.seed(2098)
cand <- listw.candidates(loc_gsl[,4:5], nb = c('del', 'gab', 'mst'), weights = c('binary', 'flin'))
sel_w <- listw.select(decostand(occ_gsl, "hellinger"), cand, MEM.autocor = "positive")
spat <- sel_w$best$MEM.select
sel_w
modgsl <- rda(decostand(occ_gsl, "hellinger") ~ ., data = spat)
anova(modgsl)
RsquareAdj(modgsl)

set.seed(138)
cand <- listw.candidates(loc_gsld[,4:5], nb = c('del', 'gab', 'mst'), weights = c('binary', 'flin'))
sel_w <- listw.select(decostand(occ_gsld, "hellinger"), cand, MEM.autocor = "positive")
spat <- sel_w$best$MEM.select
sel_w
modgsld <- rda(decostand(occ_gsld, "hellinger") ~ ., data = spat) # Will have ERROR since spatial predictors not significant

set.seed(33)
cand <- listw.candidates(loc_sev[,4:5], nb = c('del', 'gab', 'mst'), weights = c('binary', 'flin'))
sel_w <- listw.select(decostand(occ_sev, "hellinger"), cand, MEM.autocor = "positive")
spat <- sel_w$best$MEM.select
sel_w
modsev <- rda(decostand(occ_sev, "hellinger") ~ ., data = spat)
