###############################################################################
####               1. Libraries et répertoire de travail                   ####
###############################################################################

# libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)


# set working directory 
# @@ à compléter/modifier @@
setwd("~/.../QIIME_analysis/Analysis/")



###############################################################################
####                 2. Import et formatage des données                    ####
###############################################################################

# import data 
sample_phylum= read.table('taxa_summary/otu_table_filtered_L2.txt', skip = 1, comment.char= '', header = T,sep ='\t', dec='.', row.names = 1)

# visualize data
head(sample_phylum)
# ~~ Question 1 ~~
# Que pensez vous des titres des colonnes ?

# remane columns
names(sample_phylum) = gsub("\\..*", "", names(sample_phylum))
names(sample_phylum) = gsub("X", "D", names(sample_phylum))
#order columns
sample_phylum =sample_phylum[ , order(as.numeric(gsub("D", "", names(sample_phylum))))]

# visualize data
head(sample_phylum)


###############################################################################
####                       3. Barplot avec ggplot2                         ####
###############################################################################

# add a column ID o the DATAframe
sample_phylum$id = rownames(sample_phylum)

# format data 
data = melt(sample_phylum)
# ~~ Question 2 ~~
# à quoi ser la commande melt ?

# barplot
ggplot(data, aes(fill = data[,1], y = value,x=variable)) +
       geom_bar(stat='identity') +
       scale_fill_brewer(palette="Spectral")+
       xlab("Days") +
       ylab("Relative abundances of the bacterial phyla") +
       theme(axis.text.x = element_text(angle = 90, hjust = 1))+
       ggtitle("OTU-based community composition in the gut microbiota")


###############################################################################
####                             4. Exercices                              ####
###############################################################################
# ~~ Question 3 ~~
# Réaliser un barplot présentant tous les classes (et non plus les phyla ...)

# import data 
sample_class= read.table('taxa_summary/otu_table_filtered_L3.txt', skip = 1, comment.char= '', header = T,sep ='\t', dec='.', row.names = 1)


###############################################################################
####                 5. Heatmap des OTUs par échantillons                  ####
###############################################################################

sample_otus= read.table('otu_table_filtered.tsv', skip = 1, comment.char  = '', header = T,sep ='\t', dec='.', row.names = 1)
names(sample_otus) = gsub("\\..*", "", names(sample_otus))
names(sample_otus) = gsub("X", "D", names(sample_otus))
#order columns
sample_otus =sample_otus[ , order(as.numeric(gsub("D", "", names(sample_otus))))]

# transform the dataframe to a matrix
mat = as.matrix(sample_otus)

mat_prop= mat/rowSums(mat)
# ~~ Question 4 ~~
# à quoi correspond mat_prop ?

# ~~ Question 5 ~~
# Faire un heatmap à l'aide de la matrice mat_prop et de la commande pheatmap 
## pour plus de clareté
# 1. ne pas clusteriser les colonnes (cluster_cols = F)
# 2. ne pas afficher les labels des lignes (show_rownames = F) 
# 3. palette de couleur grise : c('#FFFFFF','#E3E3E3','#C6C6C6','#AAAAAA','#8E8E8E','#717171','#555555','#393939','#1C1C1C','#000000')
## aller voir l'aide ! 

