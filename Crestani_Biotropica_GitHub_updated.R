##R scripts used for data processing, statistical analyses, and figure generation 
##for the study "Morphological traits and the degree of frugivory in avian fruit consumer"


#####################Bird species morphologies and their functional roles######


morpho.rel <- read.table("morpho_rel_dados_25.csv",sep= ";", header = T, dec=",")

summary(morpho.rel)

############# Distribution of the data#############

shapiro.test(morpho.rel$Frug)  

# Histogram
hist(morpho.rel$Frug, breaks = 20, main = "Histogram of the response variable")

# Q-Q plot
qqnorm(morpho.rel$Frug)
qqline(morpho.rel$Frug, col = "red")

######################Correlation Spearman for non-parametrics data
###################### using raw values#############################

library(Hmisc)  # Required for rcorr()
dados_sel <- morpho.rel[, 2:8]
resultado <- rcorr(as.matrix(dados_sel), type = "spearman")  

pvals <- resultado$P[lower.tri(resultado$P)]

# Benjamini-Hochberg (FDR) - less conservative than Bonferroni
pval_fdr <- p.adjust(pvals, method = "BH")

print(pval_fdr)

library(corrplot)
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# p-values

rc <- rcorr(as.matrix(dados_sel), type = "pearson")
r_mat <- rc$r
p_mat <- rc$P

# FDR-adjusted p-values
p_vals <- p_mat[lower.tri(p_mat)]
p_fdr <- p.adjust(p_vals, method = "BH")


p_fdr_mat <- matrix(NA, nrow = ncol(dados_sel), ncol = ncol(dados_sel))
rownames(p_fdr_mat) <- colnames(dados_sel)
colnames(p_fdr_mat) <- colnames(dados_sel)
p_fdr_mat[lower.tri(p_fdr_mat)] <- p_fdr
p_fdr_mat[upper.tri(p_fdr_mat)] <- t(p_fdr_mat)[upper.tri(p_fdr_mat)]


r_long <- melt(r_mat, varnames = c("Var1", "Var2"), value.name = "cor")
p_long <- melt(p_fdr_mat, varnames = c("Var1", "Var2"), value.name = "p_fdr")


df_corr <- left_join(r_long, p_long, by = c("Var1", "Var2")) %>%
  mutate(
    sig = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01 ~ "**",
      p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(round(cor, 2), sig)
  )


df_corr <- df_corr %>% filter(as.character(Var1) != as.character(Var2))


df_corr$Var1 <- factor(df_corr$Var1, levels = colnames(dados_sel))
df_corr$Var2 <- factor(df_corr$Var2, levels = colnames(dados_sel))


ggplot(df_corr, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 4) +
  scale_fill_viridis(option = "viridis", limits = c(-1, 1),begin = 0.2,end = 0.8, name = "Correlation (r)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  coord_fixed() +
  labs(title = "Trait correlations using FDR", x = "", y = "")


# ============================================================================
# ALLOMETRIC RELATIONSHIPS (Tables S3 and S4)
# ============================================================================

morpho_all <- read.table(
  "morpho_dados_25.csv", sep = ";", header = TRUE,
  dec = ",", stringsAsFactors = FALSE
)

required_traits <- c("Bm", "Bl", "Bew", "Gizz", "Int", "Liv")
stopifnot(all(required_traits %in% names(morpho_all)))

# Fit the five published allometric regressions and inspect model diagnostics.
# A separate model is fitted for each response variable.
inspect_allometry <- function(model, title) {
  cat("\n--- ", title, " ---\n", sep = "")
  print(summary(model))
  plot(fitted(model), resid(model), pch = 16,
       xlab = "Fitted values", ylab = "Residuals", main = title)
  abline(h = 0, lty = 2)
  qqnorm(resid(model), main = paste("Normal Q-Q:", title))
  qqline(resid(model))
  print(shapiro.test(resid(model)))
}

# ---------------------------------------------------------------------------
# 1. Allometric relationships: all bird species (Table S3)
# Log-transformed responses: beak width, intestine length, and liver mass.
# ---------------------------------------------------------------------------

stopifnot(all(morpho_all$Bew > 0, na.rm = TRUE),
          all(morpho_all$Int > 0, na.rm = TRUE),
          all(morpho_all$Liv > 0, na.rm = TRUE))

allometry_all <- list(
  body_length = lm(Bl ~ Bm, data = morpho_all),
  beak_width = lm(log(Bew) ~ Bm, data = morpho_all),
  gizzard_area = lm(Gizz ~ Bm, data = morpho_all),
  intestine_length = lm(log(Int) ~ Bm, data = morpho_all),
  liver_mass = lm(log(Liv) ~ Bm, data = morpho_all)
)

for (trait in names(allometry_all)) {
  inspect_allometry(allometry_all[[trait]],
                    paste("All species:", trait))
}

# ---------------------------------------------------------------------------
# 2. Allometric relationships: excluding Euphonia species (Table S4)
# -----------------------------------------------------------------

stopifnot("species" %in% names(morpho_all))
morpho_no_euphonia <- morpho_all[
  !grepl("^Euphonia([_ ]|$)", as.character(morpho_all$species),
         ignore.case = TRUE),
  , drop = FALSE
]

stopifnot(nrow(morpho_no_euphonia) < nrow(morpho_all))
stopifnot(all(morpho_no_euphonia$Bew > 0, na.rm = TRUE),
          all(morpho_no_euphonia$Liv > 0, na.rm = TRUE))

allometry_no_euphonia <- list(
  body_length = lm(Bl ~ Bm, data = morpho_no_euphonia),
  beak_width = lm(log(Bew) ~ Bm, data = morpho_no_euphonia),
  gizzard_area = lm(Gizz ~ Bm, data = morpho_no_euphonia),
  intestine_length = lm(Int ~ Bm, data = morpho_no_euphonia),
  liver_mass = lm(log(Liv) ~ Bm, data = morpho_no_euphonia)
)

for (trait in names(allometry_no_euphonia)) {
  inspect_allometry(allometry_no_euphonia[[trait]],
                    paste("Excluding Euphonia:", trait))
}



# ============================================================
# 3. Beta regression: morphological predictors of frugivory
# Reference: Table S5
# ============================================================
morpho_glm <- read.table("morpho_rel_dados_25.csv", sep = ";",
                         header = TRUE, dec = ",")

stopifnot("species" %in% names(morpho_glm))
stopifnot(any(grepl("^Euphonia([_ ]|$)",
                   as.character(morpho_glm$species),
                   ignore.case = TRUE)))

# Convert percentage frugivory to proportions and adjust boundary
# values (0 or 1) to the open interval required by beta regression.
morpho_glm$Frug_prop <- morpho_glm$Frug / 100
n_glm <- nrow(morpho_glm)
morpho_glm$Frug_prop <- (morpho_glm$Frug_prop * (n_glm - 1) + 0.5) / n_glm

# Log-transform predictors, retaining biologically meaningful zeros.
log_constant <- 0.001
morpho_glm$log_Gizz <- log(morpho_glm$Gizz + log_constant)
morpho_glm$log_Liv  <- log(morpho_glm$Liv  + log_constant)
morpho_glm$log_Bl   <- log(morpho_glm$Bl   + log_constant)
morpho_glm$log_Bew  <- log(morpho_glm$Bew  + log_constant)
morpho_glm$log_Bm   <- log(morpho_glm$Bm   + log_constant)
morpho_glm$log_Int  <- log(morpho_glm$Int  + log_constant)

# Check collinearity using the same predictor set as the beta model.
glm_formula <- Frug_prop ~ log_Gizz + log_Liv + log_Bl +
  log_Bew + log_Bm + log_Int
vif_model <- lm(glm_formula, data = morpho_glm)
car::vif(vif_model)

library(glmmTMB)
mod2 <- glmmTMB(glm_formula, family = beta_family(), data = morpho_glm)

# Residual diagnostics.
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod2, plot = TRUE)

# Coefficient estimates and Wald tests; compare against Table S5.
summary(mod2)
car::Anova(mod2)

# Optional effect plots (on the transformed predictor scale).
library(ggeffects)
pred_gizzard <- ggpredict(mod2, terms = "log_Gizz")
pred_liver   <- ggpredict(mod2, terms = "log_Liv")
plot(pred_gizzard)
plot(pred_liver)



# ============================================================
# 4. Bird morphology and fruit consumption: PCA and Procrustes
# References: Figures 2-3 and Tables S6-S10
# ============================================================
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)

# 4.1 Interaction space: PCA (Tables S6-S7; Figure 2)
inter <- read.table("interactions_25.csv", sep = ";", header = TRUE)
stopifnot("species" %in% names(inter))
# Standardize bird identifiers to match the phylogenetic tree.
inter$species <- gsub(" +", "_", trimws(inter$species))
stopifnot(!anyDuplicated(inter$species))
rownames(inter) <- as.character(inter$species)
stopifnot(ncol(inter) == 10L, all(vapply(inter[-1], is.numeric, logical(1))))
inter.std <- decostand(inter[, 2:10], method = "standardize")
inter.pca <- rda(inter.std, scale = FALSE)
summary(inter.pca)

inter.eig <- inter.pca$CA$eig
inter.explained <- 100 * inter.eig / sum(inter.eig)
inter.cumulative <- cumsum(inter.explained)
print(data.frame(axis = seq_along(inter.eig),
                 explained_percent = inter.explained,
                 cumulative_percent = inter.cumulative))

inter.sites <- as.data.frame(scores(inter.pca, display = "sites", scaling = 2))
inter.sites$species <- inter$species
names(inter.sites)[1:2] <- c("PC1", "PC2")
inter.variables <- as.data.frame(scores(inter.pca, display = "species", scaling = 2))
inter.variables$variable <- rownames(inter.variables)
names(inter.variables)[1:2] <- c("PC1", "PC2")

plot_interactions <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_segment(data = inter.variables,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_text_repel(data = inter.variables,
                  aes(x = PC1, y = PC2, label = variable),
                  max.overlaps = Inf) +
  geom_point(data = inter.sites, aes(x = PC1, y = PC2),
             size = 3) +
  geom_text_repel(data = inter.sites,
                  aes(x = PC1, y = PC2, label = species),
                  max.overlaps = Inf) +
  labs(x = sprintf("PC1 (%.1f%%)", inter.explained[1]),
       y = sprintf("PC2 (%.1f%%)", inter.explained[2])) +
  theme_minimal()
print(plot_interactions)

# 4.2 Morphological trait space: PCA (Tables S8-S9; Figure 3)
morpho.pca.data <- read.table("morpho_dados_25.csv",
                              sep = ";", header = TRUE, dec = ",")
# Use the same bird identifiers as the interaction matrix and bird tree.
morpho.pca.data$species <- gsub(" +", "_", trimws(morpho.pca.data$species))
morpho_sp <- morpho.pca.data %>%
  group_by(species) %>%
  summarise(across(c(Bm, Bl, Bew, Gizz, Int, Liv),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")
morpho.vars <- morpho_sp[, c("Bm", "Bl", "Bew", "Gizz", "Int", "Liv")]
morpho.std <- decostand(morpho.vars, method = "standardize")
rownames(morpho.std) <- as.character(morpho_sp$species)
morpho.pca <- rda(morpho.std, scale = FALSE)
summary(morpho.pca)

morpho.eig <- morpho.pca$CA$eig
morpho.explained <- 100 * morpho.eig / sum(morpho.eig)
morpho.cumulative <- cumsum(morpho.explained)
print(data.frame(axis = seq_along(morpho.eig),
                 explained_percent = morpho.explained,
                 cumulative_percent = morpho.cumulative))

morpho.sites <- as.data.frame(scores(morpho.pca, display = "sites", scaling = 2))
morpho.sites$species <- morpho_sp$species
names(morpho.sites)[1:2] <- c("PC1", "PC2")
morpho.variables <- as.data.frame(scores(morpho.pca, display = "species", scaling = 2))
morpho.variables$variable <- rownames(morpho.variables)
names(morpho.variables)[1:2] <- c("PC1", "PC2")

plot_morphology <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_segment(data = morpho.variables,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_text_repel(data = morpho.variables,
                  aes(x = PC1, y = PC2, label = variable),
                  max.overlaps = Inf) +
  geom_point(data = morpho.sites, aes(x = PC1, y = PC2), size = 3) +
  geom_text_repel(data = morpho.sites,
                  aes(x = PC1, y = PC2, label = species),
                  max.overlaps = Inf) +
  labs(x = sprintf("PC1 (%.1f%%)", morpho.explained[1]),
       y = sprintf("PC2 (%.1f%%)", morpho.explained[2])) +
  theme_minimal()
print(plot_morphology)

# 4.3 Procrustes / PROTEST (Table S10)
stopifnot(setequal(as.character(inter$species),
                   as.character(morpho_sp$species)))
morpho_order <- match(as.character(inter$species),
                      as.character(morpho_sp$species))
morpho.std.aligned <- morpho.std[morpho_order, , drop = FALSE]
stopifnot(identical(as.character(inter$species),
                    rownames(morpho.std.aligned)))
morpho.pca.aligned <- rda(morpho.std.aligned, scale = FALSE)

set.seed(1001)
concord <- procrustes(inter.pca, morpho.pca.aligned)
procrustes_test <- protest(inter.pca, morpho.pca.aligned)
print(concord)
print(procrustes_test)

# ============================================================
# Phylogenetic analyses 
# ============================================================

#crestani_v01
library(geiger)
library(picante)
library(phytools)
library(ape)
library(caper)
library(plyr)
library(reshape)
library(reshape2)
library(data.table)
library(dplyr)
library(Matrix)
library(lme4)#
library(tibble)
library(phangorn)
library(plotly)
library(phylobase)
library(ape)
library(phyloextra)#
library(phylotools)
library(adephylo)#
library(geiger)
library(Hmisc)
library(rphast)#
library("readxl")
library(phylosignal)
require(gdata)
library(here)
library(ggpubr)#
library(ggrepel)#points names into ggplot
library(ggtree) #to plot tree in ggplot #
library('TreeTools')#
library(ggplot2)
library(lattice)
library(phylolm)

setwd("/media/usuario/Elements/Erison/extra/Ana")
#consensus tree of the 1k bird tree taken by birdtree.org
#phylobird_all<-read.nexus("crestani_08_2025.nex")#hacket all:9993 OUTs x1000
#phylobird<-averageTree(phylobird,method="quadratic.path.difference") #criando a arvore de consenso
#phylobird
#vcv<-vcv(phylobird)
#writeNexus(phylobird,"bird_phy_v05.nex")
phylobird<-read.nexus("bird_phy_v05.nex")
#write.csv(vcv,"vcv_bird_plant.csv")
plotTree(phylobird)
is.rooted(phylobird)
order(phylobird$tip.label)
order(trait$species)
#signal

#trait<-read_excel("crestani_morpho_08_25.xlsx")
trait <- read.csv("morpho_dados_25.csv", 
                  stringsAsFactors = FALSE, 
                  sep = "\t", 
                  dec = ",")
trait$species<-stringr::str_replace(trait$species, ' ', '_') #trocando espaço por "_"
traitm <- trait %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


trait$species
rownames(traitm)<-traitm$species
name.check(phylobird,traitm$species)

#trait$ord<-order(phylobird$tip.label) #Organizando espécies na mesma ordem da arvore filogenetica
#trait
#trait<-arrange(trait,ord)
#trait


#separando traits e ligando aos nomes das especies
#gizzard
giz<-c(traitm$Gizz)
giz<-setNames(giz,traitm$species)
name.check(phylobird,giz)
# intestine length
int<-c(traitm$Int)
int<-setNames(int,traitm$species)
#skin
#skin<-c(trait2$skin)
#skin<-setNames(skin,trait$species)
#liver mass
liv<-c(traitm$Liv)
liv<-setNames(liv,traitm$species)
#surf
#surf<-c(trait$Surf)
#surf<-setNames(surf,trait$species)
#frug
frug<-c(traitm$Frug)
frug<-setNames(frug,traitm$species)
#body length
body_l<-c(traitm$Bl)
body_l<-setNames(body_l,traitm$species)
#body mass
bew<-c(traitm$Bew)
bew<-setNames(bew,traitm$species)

#body mass
body_m<-c(traitm$Bm)
body_m<-setNames(body_m,traitm$species)


################# Traits signal###################
trait<-data.frame(traitm)
rownames(trait)<-traitm$species
trait<-trait[,-1]
#phy4_t<-phylo4d(phylobird,tip.data=trait) #by gizzard area
#dotplot.phylo4d(phy4_t)

#gizzard
phy4_giz<-phylo4d(phylobird,tip.data=giz) #by gizzard area
phyloSignal(phy4_giz)#Not significant
dotplot.phylo4d(phy4_giz,trait.labels="Gizzard area mm²",tree.ladderize=T)

#intestine
phy4_int<-phylo4d(phylobird,tip.data=int) #by intestine length
phyloSignal(phy4_int)#Not significant
dotplot.phylo4d(phy4_int)

#liver
phy4_liv<-phylo4d(phylobird,tip.data=liv) #by liver mass
phyloSignal(phy4_liv)# NOT significant
dotplot.phylo4d(phy4_liv)

#phy4_skin<-phylo4d(phylobird,tip.data=skin) #by skin
#phyloSignal(phy4_skin)#NOT significant
#dotplot.phylo4d(phy4_skin)

#frug
phy4_frug<-phylo4d(phylobird,tip.data=frug) #by frugivory degree
phyloSignal(phy4_frug)#significant
dotplot.phylo4d(phy4_frug,trait = names(tdata(phy4_frug)))

#body_l
phy4_body_l<-phylo4d(phylobird,tip.data=body_l) #by body length
phyloSignal(phy4_body_l)# significant
dotplot.phylo4d(phy4_body_l,trait = names(tdata(phy4_body_l)))

#body_m
phy4_body_m<-phylo4d(phylobird,tip.data=body_m) #by body mass
phyloSignal(phy4_body_m)# significant
dotplot.phylo4d(phy4_body_m,trait = names(tdata(phy4_body_m)))

#surf
#phy4_surf<-phylo4d(phylobird,tip.data=surf) #by surf
#phyloSignal(phy4_surf)#NOT significant
#dotplot.phylo4d(phy4_surf,trait = names(tdata(phy4_surf)))

#bew
phy4_bew<-phylo4d(phylobird,tip.data=bew) #by beak width
phyloSignal(phy4_bew)#significant
dotplot.phylo4d(phy4_bew,trait = names(tdata(phy4_bew)))


#significant are
#gizzard diameter, intestine length, frugivory degree and body mass/length and beak width
#traits variance are correlated with species phylogeny
#######Dotplot treee
rownames(trait)<-trait$species
trait_d<-cbind(trait$Bm,trait$Bl,trait$Frug,trait$Bew)
colnames(trait_d)<-c("body mass","body length","frugivory degree","beak width")
rownames(trait_d)<-trait$species
phy4_t<-phylo4d(phylobird,tip.data=trait_d)
dotplot.phylo4d(phy4_t, dot.cex=4, tip.cex=2,trait.cex=2)



png("dotplot_signal_v06.png",2100,1200)
dotplot.phylo4d(phy4_t, dot.cex=4, tip.cex=2.8,trait.cex=2)
dev.off()




##PGLS
#How is the relationship between frugivory and morphological traits varies with
#phylogeny as co-variable?

#trait2<-read.csv("morpho_rel_media_v02.csv",sep=",", head=T)#
bird_specie<-trait$species
trait<-data.frame(trait)
trait$frug_log<-log10((trait$Frug+1)/100)

trait$species <- rownames(trait)
 gizf <- try(comparative.data(phy = phylobird, 
               data = trait, 
               names.col = "species",
               vcv = TRUE))

summary.pgls(pgls(frug_log~Gizz,data=gizf,lambda = "ML"))# NOT significant
summary.pgls(pgls(frug_log~Int,data=gizf,lambda = "ML")) #NOT significant
summary.pgls(pgls(frug_log~Liv,data=gizf,lambda = "ML"))#NOT significant
#summary.pgls(pgls(frug~skin,data=gizf)) #significant
#summary.pgls(pgls(frug~surf,data=gizf)) #NOT significant
summary.pgls(pgls(frug~Bl,data=gizf,lambda = "ML")) #NOT significant
summary.pgls(pgls(frug_log~Bm,data=gizf,lambda = "ML")) #NOT significant

summary.pgls(pgls(frug_log~Bew,data=gizf,lambda = "ML")) # significant
#the relationship between frugivory degree gizzard area has correlation with phylogenetic
#distance between the species.
#The same to frugivory degree and liver mass#### NO more
plot(gizf$data)

#plotting significant relationships
giz_g <- ggplot(data = trait, aes(x = Bew, y = frug_log)) +
  geom_smooth(method = "glm", 
              method.args = list(family = gaussian),  # Assumindo dados contínuos
              color = "black", 
              se = FALSE) +
  geom_point(size = 5) +
  labs(y = "Frugivory degree (log)", x = "Beak width (mm)") +
  geom_label_repel(aes(label = bird_specie),
                   box.padding = 1.6,
                   point.padding = 0,
                   segment.color = 'grey50') +
  theme_classic() +  # Usando tema padrão primeiro
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 30))
giz_g

ggsave("phylo_giz_frug_v5.jpeg",height=8,width=14,dpi=300)


#liv_g<-ggplot(data=trait, aes(x=liv,y=frug))+
#  geom_smooth(method=lm, col="black")+
#  geom_point(size=5)+
#  scale_color_brewer(palette = "RdYlBu")+
#  labs(y="Frugivory degree (%)",x=" Liver mass (mg)")+
#  geom_label_repel(aes(label = bird_specie),
#                   box.padding   = 1.6,
#                   point.padding = 0,
#                   segment.color = 'grey50') +
#  theme_classic2()+annotate("text",x=0.04,y=0.7,size=7, label="b=1320.3 p=0.04 r²=0.41")+
#  theme(axis.title = element_text(size =20),
#        axis.text = element_text(size = 20),legend.text=element_text(size=30))
#liv_g
#
#ggsave("phylo_liv_frug.jpeg",height=8,width=14,dpi=300)

max_depth <- max(node.depth.edgelength(phylobird))

 print(max_depth)

 tree_g <- ggtree(phylobird) +
   geom_tiplab(
     size = 6,
     offset = 0.1,
     hjust = 0,
     align = FALSE
   ) +
   xlim(0, max_depth * 1.3) +  # Adiciona 30% de margem
   theme_tree()

ggsave("birdtree_v5.jpeg",height=8,width=14,dpi=300)


ar1<-ggarrange(tree_g,giz_g,nrow=2,ncol=1,labels="AUTO",hjust = -7,
               legend="none")
png("ana_v6.png",1600,1000)
ar1
dev.off()



ps_k<-phylosig(phylobird,frug,method="lambda",test=T)
print(ps_k)
plot(ps_k)

##Install packages and libraries
#library(devtools)
#install_github("hmorlon/PANDA",ref="Benoit_phylosignal", dependencies = TRUE)
#devtools::install_github("jinyizju/V.PhyloMaker")
#suppressPackageStartupMessages(library("V.PhyloMaker")) 
suppressPackageStartupMessages(library(RPANDA))


# Load packages
library(ape)
library(tidyverse)

# Read the phylogenetic trees already supplied with the repository.
# plant_tree.nex contains Newick-formatted text despite its .nex extension.
phylo.plant <- ape::read.tree("plant_tree.nex")
phylo.bird <- ape::read.nexus("bird_phy_v05.nex")

# Read the bird-by-plant consumption-rate matrix.
net <- read.csv("interactions_25.csv", h = TRUE, row.names = 1,
                sep = ";", check.names = FALSE)
rownames(net) <- gsub(" +", "_", trimws(rownames(net)))

# Check that the new plant genus names and bird identifiers match the trees.
if (!setequal(colnames(net), phylo.plant$tip.label)) {
  stop("Plant names in interactions_25.csv do not match plant_tree.nex.")
}
if (!setequal(rownames(net), phylo.bird$tip.label)) {
  stop("Bird identifiers differ between the network rows and bird tree tips.")
}
net <- net[phylo.bird$tip.label, phylo.plant$tip.label, drop = FALSE]

plot(phylo.plant, cex = 0.9, main = "Phylogenetic Tree of Plant Species")
plot(phylo.bird, cex = 0.9, main = "Phylogenetic Tree of Bird Species")



#Check data
network <-  as.matrix(net) # int. matrix with plants in columns and birds in rows
network

# compute phylogenetic signals in species interactions
output.nbpart1 <- phylosignal_network(network, tree_A = phylo.plant, tree_B = phylo.bird, method = "Jaccard_weighted", correlation = "Pearson", nperm=10000, permutation = "nbpartners") 
print(output.nbpart1)
# nb_A   nb_B  mantel_cor_A   pvalue_A    mantel_cor_B    pvalue_B     
# 9      16      -0.0112673     0.46         0.1264960       0.13 
       


#method = "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities 
#permutations	indicates which permutations must be performed to evaluate the significance of the Mantel correlation: "nbpartners" keeps constant the number of partners per species and shuffling at random their identity

#**interpretation: closely related plant species tend not to interact with similar birds. On the other hand, closely related birds tended to interact with similar plants (p-value B was marginally significant = 0.05)**#

#Number of species in guild A = "nb_A"
#Number of species in guild B = "nb_B"
#Correlation for guild A (in this case plants) = "mantel_cor_A"
# its associated upper p-value = "pvalue_upper_A"    
#***if "pvalue_upper_A"<0.05, there is a significant phylogenetic signal in guild A.
#Correlation for guild B (in this case birds) = "mantel_cor_B"
# its associated upper p-value = "pvalue_upper_B"    
#if "pvalue_upper_B"<0.05, 
#***there is a significant phylogenetic signal in guild A.***