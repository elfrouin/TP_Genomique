# à faire avant la séance de mercredi
source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")

###############################################################################
####               1. Libraries et répertoire de travail                   ####
###############################################################################

# libraries
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)

# set working directory 
# @@ à compléter/modifier @@
setwd("~/.../QIIME_analysis/Analysis/")


###############################################################################
####                   2. Création de l'objet phyloseq                     ####
###############################################################################
# chargement des fichhiers dans R studio
otu="otu_table_filtered.json"
map="Scripts_R/mapfile_phyloseq.csv"
tree="rep_phylo.tre"

treefile  <- read_tree(tree)
a         <- import_biom(otu, treefile, parseFunction = parse_taxonomy_greengenes)

# Creation de l'objet phyloseq :
mapfile   <- import_qiime_sample_data(map)
finalobject <- merge_phyloseq(a,mapfile)

# Visualiser l'objet phyloseq
finalobject
# vous pouvez accéder aux informations sur otu_table, sample_data, tax_table, phy_tree
otu_table(finalobject)
sample_data(finalobject)
# ~~ Question 1 ~~
# Faire de même pour tax_table et phy_tree



###############################################################################
####          4. Construction d'un arbre phylogénétique                    ####
###############################################################################
# construction d'un arbre (level =Phylum)
data_tree <- tax_glom(finalobject,"Phylum")
plot_tree(data_tree,color="Step",size="abundance",label.tips="Phylum",ladderize = "TRUE", nodelabf = nodeplotblank,base.spacing = 0.025,plot.margin = 0.4, text.size = 3)

# ~~ Question 2 ~~
# tester la construction d'un arbre au rang de la Class 
### liste des rangs disponibles
rank_names(finalobject)


# construction d'un arbre sur un sous ensemble  (ici niveau phylum = Actinobacteria)
subact <- subset_taxa(finalobject, Phylum=="Actinobacteria")
keeptaxa <- names(sort(taxa_sums(subact), TRUE))
final <-  prune_taxa(keeptaxa, subact)
# ~~ Question 3 ~~
# Combien il y a-t-il d'OTUs dans l'objet phyloseq final ?
plot_tree(final,color="Step",size="abundance",label.tips="Family",ladderize = "TRUE", nodelabf = nodeplotblank,base.spacing = 0.025,plot.margin = 0.1, text.size = 2)

# pour limiter la representation de l'arbre aux OTUs essentiels (par ex les 20 premiers)
keeptaxa <- names(sort(taxa_sums(subact), TRUE)[1:20])
final <-  prune_taxa(keeptaxa, subact)
plot_tree(final,color="Step",size="abundance",label.tips="Family",ladderize = "TRUE", nodelabf = nodeplotblank,base.spacing = 0.025,plot.margin = 0.1, text.size = 2)

###############################################################################
####    5. Heatmap de l'abondance des Firmicutes par échantillons          ####
###############################################################################
table_firmicutes <- subset_taxa(finalobject, Phylum=="Firmicutes")
table_firmicutes <- prune_taxa(names(sort(taxa_sums(table_firmicutes),TRUE)[1:200]), table_firmicutes)
plot_heatmap(table_firmicutes,taxa.label = "Family")
