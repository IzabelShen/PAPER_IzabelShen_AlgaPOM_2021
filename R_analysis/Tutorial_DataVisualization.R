###########################
### R scripts for plotting and data visualization 
### Data from 
### "Quality of phytoplankton deposition structures bacterial communities at the water-sediment interface"
### By D Izabel-Shen, S Albert, M Winder, H Farnelid and F Nascimento
### 
###Prepared 9 February 2021
### Author: Dandan Izabel-Shen, Stockholm University; dand.shen <at> gmail.com
###########################

# Before you start
# Make sure you are using the latest version of R (and Rstudio)
#The following packages (and their dependencies) are needed to run the whole analysis 

#vegan
#ggplot2
#plyr
#gplots
#phyloseq
#gridExtra
#dichromat
#ggpubr
#VennDiagram
#RColorBrewer

##########################################################
#####Plotting within-sample (Alpha) diveristy water sample
#########################################################
# load R libraries for this section 
library(ggplot2)
library(gridExtra)

#read the file 
## Reading dataset 
Comm_Diversity <-read.csv("/InputFiles/Alpha.csv", header=T, row.names = 1)

#select the data for species richness and evenness separately
Richness <-Comm_Diversity[1:29,]
Evenness <-Comm_Diversity[30:58,]

#set different levels of microcosms as factors in each dataset
Richness$Microcosms <-factor(Richness$Microcosms, levels =(unique(Richness$Microcosms)), ordered=TRUE)
Evenness$Microcosms <-factor(Evenness$Microcosms, levels =(unique(Evenness$Microcosms)), ordered=TRUE)

#set color code as factors in each dataset
Richness$Color <-as.character(Richness$Color)
Evenness$Color <-as.character(Evenness$Color)


# plotting using ggplot2

#select the data for x and y axes
p1 <-ggplot(Richness, aes(x=Microcosms,y=Values, fill=Microcosms)) + 

#set boxplot function 
geom_boxplot(fatten=0.70, alpha=0.80)+ 

#set scattering points to show the deviation of biological replicates in this case
geom_point(pch=21, position=position_jitterdodge(),alpha=0.45, size=2.5) + 

#manual set the color scales to fill the each box
scale_fill_manual(values=c("Control"="deepskyblue3","100D"="tan4", "80D_20C"="orange1", "50D_50C"="gold1", "20D_80C"="olivedrab2", "100C"="darkgreen")) + 

#use 'theme' function to set specific requirement for plotting 
#set a theme with white background and black gridlines
theme_bw()+ 

#remove the minor gridlines for x and y axes
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 

#set title for the figure 
labs(title="A                              Richness") + 

#set color, size, and fonts for the title of x and y axes
theme(axis.title.y=element_text(colour = 'black',size = 14, family="serif")) + 
theme(axis.title.x=element_text(colour = 'black', size = 14,family="serif" )) + 

#set color, size and fonts for the text of x and y axes
theme(axis.text.y=element_text(colour = 'black', size = 14, family="serif")) + 
theme(axis.text.x=element_text(colour = 'black', size = 14, family="serif")) +  

#set the text parameters for legends
theme(legend.text=element_text(size=14, family="serif")) + 
theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
theme(legend.background=element_rect(fill=NA)) + 
theme(legend.title=element_text(color="black", size=13, face='bold', family="serif")) + 

#set the shape and size of the symbols for legend
guides(color=guide_legend(override.aes=list(shape=15, size=3))) + 
guides(shape=guide_legend(override.aes=list(size=2.5))) + 
theme(legend.key.size=unit("0.5","cm")) 

#the same function and computing notes used to plot evenness
p2<-ggplot(Evenness, aes(x=Microcosms,y=Values, fill=Microcosms)) + 
geom_boxplot(fatten=0.70, alpha=0.80)+ 
geom_point(pch=21, position=position_jitterdodge(),alpha=0.45, size=2.5) + 
scale_fill_manual(values=c("Control"="deepskyblue3","100D"="tan4", "80D_20C"="orange1", "50D_50C"="gold1", "20D_80C"="olivedrab2", "100C"="darkgreen")) + 
theme_bw()+ 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
labs(title="B                              Evenness") + 
theme(axis.title.y=element_text(colour = 'black',size = 14, family="serif")) + 
theme(axis.title.x=element_text(colour = 'black', size = 14, family="serif")) + 
theme(axis.text.y=element_text(colour = 'black', size = 14, family="serif")) + 
theme(plot.title=element_text(colour= 'black', size=11, face='bold', family="serif")) + 
theme(axis.text.x=element_text(colour = 'black', size = 14, family="serif")) +  
theme(legend.text=element_text(size=14, family="serif")) + 
theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
theme(legend.background=element_rect(fill=NA)) + 
theme(legend.title=element_text(color="black", size=13, face='bold', family="serif")) + 
guides(color=guide_legend(override.aes=list(shape=15, size=3))) + 
guides(shape=guide_legend(override.aes=list(size=2.5))) + 
theme(legend.key.size=unit("0.5","cm")) 

#visualize the plots
p1  # richness
p2  #evenness

#### put p1 and p2 together in a figure and save as PDF, use the command below
#for example
#set width and height for the output 
pdf("/Your directory/Figure.pdf", width =14, height=10)  
#use 'nrow' if you want to display two panels in a row, or use 'ncol' if you want to display two panel in column.
grid.arrange(p1,p2, nrow=2, widths = c(1.6,1.6), heights=c(2,2))  
dev.off() # to clear off figure visualization in the R window 


##########################################################
#####Plotting between-sample (Beta) diveristy water sample
##########################################################
# load R libraries for this section 
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(plyr)  ### to plot the hulls

#read a phyloseq object containing ASV abundance, taxonomic affiliation and sample info
waterFinal <-readRDS("/InputFiles/waterFinal.rds")

#Transform the counttable to relative abundance
Relat <-transform_sample_counts(waterFinal, function(OTU)  OTU/sum(OTU)*100)
#ceheck if the relative abundances are displayed in the dataset
otu_table(Relat)

#####################
#######Bray-Curtis dissimilarity matrix
#####################
nmds1<-ordinate(waterFinal, "NMDS", "bray")

#check the stress level for the NMDS ordination
nmds1 $ stress 

#create a dataframe including sample info and the matrix info extracted from Bray Curtis
Mymeta = data.frame (Microcosms=c("Control","Control","Control", "Control",	"Control","100D",	"100D", "100D","100D", "80D_20C", "80D_20C", "80D_20C","80D_20C", "80D_20C", "50D_50C","50D_50C",	"50D_50C"	,"50D_50C",	"50D_50C",	"20D_80C", "20D_80C","20D_80C", "20D_80C" ,"20D_80C" ,"100C", "100C",	"100C" , "100C" ,"100C"))
NMDS = data.frame (NMDS1=nmds1$points[,1], NMDS2=nmds1$points[,2], Microcosms=Mymeta$Microcosms)
NMDS$Microcosms <-factor(NMDS$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C", "20D_80C", "100C"))

###Cretae a dataset for polygon plotting
find_hull <-function(NMDS) NMDS[chull(NMDS$NMDS1,NMDS$NMDS2),]
hulls <-ddply(NMDS, "Microcosms",find_hull) 

#plotting Bray-Curtis based NMDS ordination (computing notes for each line please see Alpha diversity)
p3 <-ggplot(data=NMDS, aes(NMDS1, NMDS2, fill=Microcosms, color=Microcosms)) + 
geom_point(size=2.5) + 
geom_polygon(data = hulls, aes(x=NMDS1, y=NMDS2), alpha=0.3) + 
scale_fill_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen")) + 
scale_color_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen"))  + 
theme_bw()+ 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+ 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
labs(title="C                      Bray-Curtis NMDS") + 
theme(axis.title.y=element_text(colour = 'black',size = 14,family='serif')) + 
theme(axis.title.x=element_text(colour = 'black',size = 14, family='serif')) + 
theme(axis.text.y=element_text(colour = 'black', size = 14,family='serif' )) + 
theme(axis.text.x=element_text(colour = 'black', size = 14, family='serif')) + 
theme(legend.text=element_text(size=13,family="serif")) + 
theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
theme(legend.background=element_rect(fill=NA)) + 
theme(legend.title=element_text(color="black", size=14, face='bold', family='serif')) + 
guides(color=guide_legend(override.aes=list(shape=15, size=3))) + 
guides(shape=guide_legend(override.aes=list(size=2.5))) + 
theme(legend.key.size=unit("0.5","cm")) 

#visualize the plot
p3

#if you want to save the figure as PDF 
#### put p1 and p2 together in a figure for manuscript and save as PDF
#for example
#set width and height for the output 
#pdf("/Your directory/BrayNMDS.pdf", width = 7.5, height=5.8)
#dev.off()# to clear off figure visualization in the R window 


#####################
#######UniFrac dissimilarity matrix
#####################

#For UniFrac, a dissimilarity matrix taking phylogenetic breadth of samples into account
#to get calculate phylogenetic breadth, sequences of all ASVs are required
#calculation of phylogenetic breadth was carried out in QIIM2

# Install miniconda3 in your Linxus or Unix terminal
#activate QIIM2 environment in the miniconda3
source activate qiime2-2019.10

# import fasta files into qiime2 format 
qiime tools import \
  --input-path /InputFiles/WaterSeqs.fna  \
  --output-path /Your directory/WaterSeqs.qza \  # change to your directory
  --type 'FeatureData[Sequence]'

# 1_The first step to build phylogenetic tree is to alignment the sequences 
qiime alignment mafft\
 --i-sequences /Your directory/WaterSeqs.qza \
 --o-alignment /Your directory/aligned-rep-seps.qza 
 
##Time to run : 1 minute

# 2_Masking sites
qiime alignment mask \
  --i-alignment /Your directory/aligned-rep-seps.qza \
  --o-masked-alignment /Your directory/masked-aligned-rep-seps.qza 
  
# 3_Creating a tree
qiime phylogeny fasttree \
  --i-alignment /Your directory/masked-aligned-rep-seps.qza \
  --o-tree /Your directory/unrooted-tree.qza
  
# 4_Midpoiting rooting
#Fastree creates an unrooted tree.

qiime phylogeny midpoint-root \
  --i-tree /Your directory/unrooted-tree.qza \
  --o-rooted-tree /Your directory/rooted-tree.qza

qiime tools export \
 --input-path /Your directory/rooted-tree.qza \
 --output-path /Your directory/epxorted-tree

########Until here, a phylogenetic tree file is created using QIIM2

#back to R analysis
#add phylogenetic tree data compoent to the phyloseq object 
tree <-read_tree("/Input/tree.nwk")
newPhyloObject =merge_phyloseq(waterFinal,  tree)
newPhyloObject ## check if the phyloseq object is right 
#a phyloseq object reads as follows
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3692 taxa and 29 samples ]
#sample_data() Sample Data:       [ 29 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 3692 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 3692 tips and 3002 internal nodes ]
#refseq()      DNAStringSet:      [ 3692 reference sequences ]

#save a phyloseq object
saveRDS (newPhyloObject, "/InputFiles/newPhyloObject.rds") # change to your directory

#set random number generator (RNG) seed, for perfect reporudcibility
set.seed(711)

#calculate UniFrac using the new phyloseq object 
#weighted UniFrac, matrix weights on abundant taxa
uniFrac_weighted <-UniFrac(newPhyloObject, weighted = TRUE, normalized = TRUE, parallel = FALSE, fast = TRUE)

#unweighted UniFrac, matrix considers both presence-absence of taxa
uniFrac_unweighted <-UniFrac(newPhyloObject, weighted=FALSE, normalized = TRUE, parallel = FALSE, fast = TRUE)  #default weighted=FALSE

#a quick glimpse of the two matrices
uniFrac_weighted
uniFrac_unweighted

#construct UniFrac dissimilarity based NMDS ordinations
nmds_Weight <-metaMDS(uniFrac_weighted)  # weighted distance first 
nmds_unweight <-metaMDS(uniFrac_unweighted)
  
#check the stress levels of the ordinations 
nmds_Weight$stress 
nmds_unweight$stress 

 
########Weighted UniFrac-distance NMDS######## 

#create a dataframe including sample info and distance matrix data point for weighted UniFrac 
Mymeta = data.frame (Microcosms=c("Control","Control","Control", "Control",	"Control","100D",	"100D", "100D","100D", "80D_20C", "80D_20C", "80D_20C","80D_20C", "80D_20C", "50D_50C","50D_50C",	"50D_50C"	,"50D_50C",	"50D_50C",	"20D_80C", "20D_80C","20D_80C", "20D_80C" ,"20D_80C" ,"100C", "100C",	"100C" , "100C" ,"100C"), Replicate= c("Control_1", "Control_2", "Control_3",	"Control_4",	"Cotrol_5",	"100D_2",	"100D_3", "100D_4", "100D_5", "80D/20C_1", "80D/20C_2",	"80D/20C_3",	"80D/20C_4",	"80D/20C_5"	,"50D/50C_1",	"50D/50C_2",	"50D/50C_3",	"50D/50C_4",	"50D/50C_5",	"20D/80C_1", "20D/80C_2", "20D/80C_3"	,"20D/80C_4",	"20D/80C_5", "100C_1",	"100C_2",	"100C_3",	"100C_4",	"100C_5"))
NMDS_UniFrac_weight = data.frame (NMDS1=nmds_Weight$points[,1], NMDS2=nmds_Weight$points[,2], Microcosms=Mymeta$Microcosms)
NMDS_UniFrac_weight$Microcosms <-factor(NMDS_UniFrac_weight$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C", "20D_80C", "100C"))

#Create a dataset for polygon of weight UniFrac
find_hull_1 <-function(NMDS_UniFrac_weight) NMDS_UniFrac_weight[chull(NMDS_UniFrac_weight$NMDS1,NMDS_UniFrac_weight$NMDS2),]
hulls_1 <-ddply(NMDS_UniFrac_weight, "Microcosms",find_hull_1) 

#plot UniFrac weighted (notes for each function see Alpha diversity plotting)
p4 <-ggplot(data=NMDS_UniFrac_weight, aes(NMDS1, NMDS2, fill=Microcosms, color=Microcosms)) + 
geom_point(size=3.3) + geom_polygon(data = hulls_1, aes(x=NMDS1, y=NMDS2), alpha=0.3) + 
scale_fill_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen")) + 
scale_color_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen"))  + 
theme_bw()+ 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+ 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
labs(title="Weighted-UniFrac NMDS") + theme(axis.title.y=element_text(colour = 'black',size = 11)) + 
theme(axis.title.x=element_text(colour = 'black',size = 11)) + 
theme(axis.text.y=element_text(colour = 'black', size = 11)) + 
theme(plot.title=element_text(colour= 'black', size=11, face='bold')) + 
theme(axis.text.x=element_text(colour = 'black', size = 11)) + 
theme(legend.title=element_text(color="black", size=12, face='bold', family="serif")) + 
theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
theme(legend.background=element_rect(fill=NA))  + 
guides(color=guide_legend(override.aes=list(shape=15, size=3))) + 
guides(shape=guide_legend(override.aes=list(size=2.5))) + 
theme(legend.key.size=unit("0.5","cm"))


######## UnWeighted UniFrac-distance NMDS########
#similar notes as for Weighted UniFrac 
Mymeta = data.frame (Microcosms=c("Control","Control","Control", "Control",	"Control","100D",	"100D", "100D","100D", "80D_20C", "80D_20C", "80D_20C","80D_20C", "80D_20C", "50D_50C","50D_50C",	"50D_50C"	,"50D_50C",	"50D_50C",	"20D_80C", "20D_80C","20D_80C", "20D_80C" ,"20D_80C" ,"100C", "100C",	"100C" , "100C" ,"100C"), Replicate= c("Control_1", "Control_2", "Control_3",	"Control_4",	"Cotrol_5",	"100D_2",	"100D_3", "100D_4", "100D_5", "80D/20C_1", "80D/20C_2",	"80D/20C_3",	"80D/20C_4",	"80D/20C_5"	,"50D/50C_1",	"50D/50C_2",	"50D/50C_3",	"50D/50C_4",	"50D/50C_5",	"20D/80C_1", "20D/80C_2", "20D/80C_3"	,"20D/80C_4",	"20D/80C_5", "100C_1",	"100C_2",	"100C_3",	"100C_4",	"100C_5"))
NMDS_UniFrac_Unweight = data.frame (NMDS1=nmds_unweight$points[,1], NMDS2=nmds_unweight$points[,2], Microcosms=Mymeta$Microcosms)
NMDS_UniFrac_Unweight$Microcosms <-factor(NMDS_UniFrac_weight$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C", "20D_80C", "100C"))
find_hull_2 <-function(NMDS_UniFrac_Unweight) NMDS_UniFrac_Unweight[chull(NMDS_UniFrac_Unweight$NMDS1,NMDS_UniFrac_Unweight$NMDS2),]
hulls_2 <-ddply(NMDS_UniFrac_Unweight, "Microcosms",find_hull_2) 

#plot UniFrac weighted
p5 <-ggplot(data=NMDS_UniFrac_Unweight, aes(NMDS1, NMDS2, fill=Microcosms, color=Microcosms)) + 
geom_point(size=2.5) + 
geom_polygon(data = hulls_2, aes(x=NMDS1, y=NMDS2), alpha=0.3) + 
scale_fill_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen")) + 
scale_color_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen"))  + 
theme_bw()+ 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+ 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
labs(title="D                       UniFrac-unweighted NMDS") + 
theme(axis.title.y=element_text(colour = 'black',size = 14, family="serif")) + 
theme(axis.title.x=element_text(colour = 'black',size = 14, family="serif")) + 
theme(axis.text.y=element_text(colour = 'black', size = 14, family="serif")) + 
theme(plot.title=element_text(colour= 'black', size=14, face='bold', family="serif")) + 
theme(axis.text.x=element_text(colour = 'black', size = 14,family="serif" )) + 
theme(legend.text=element_text(size=14, family="serif")) + 
theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
theme(legend.background=element_rect(fill=NA)) + 
theme(legend.title=element_text(color="black", size=13, face='bold', family="serif")) + 
guides(color=guide_legend(override.aes=list(shape=15, size=3))) + 
guides(shape=guide_legend(override.aes=list(size=2.5))) + 
theme(legend.key.size=unit("0.5","cm")) 

#put p1,p2, p3, p5 panels together to create a figure, like Figure 1 in the manuscript
pdf("/InputFiles/Figure1.pdf", width =14, height=10) #change to your directory
grid.arrange(p1, p3, p2, p5, ncol=2, nrow=2, widths = c(1.6,1.6), heights=c(2,2))
dev.off()



#########################################################################
#####Create a phylogenetic tree containing the 100 ASVs
#####that a significant change in relative abundances between treatments 
#########################################################################
###load libraries for this section
library(phyloseq)
library(ape)  # for write DNA sequences in a file

#load phyloseq object
newPhyloObject <-readRDS ("/InputFiles/newPhyloObject.rds")

#load dataframe containing 100 ASVs 
sigtab1 <-read.csv("/InputFiles/sigtab1.csv", header=T, row.names=1)

#select the row names of the dataframe
keepASVs <-rownames(sigtab1)

#use 'prune' to subset the phyloseq only containing the 100 ASVs
signObject <-prune_taxa(keepASVs, newPhyloObject)

#get the fasta sequences of these 100 ASVs
Seqs100 <-refseq(signObject)

#save the fasta sequences in a file
write.dna(Seqs100,"/InputFiles/Seqs100.fasta",format="fasta") #CHANGE to your own folder path 


#use QIIME2 to construct a phylogenetic tree 
#NOTE: Remove the gaps of the fasta file using BBEdit and change the files extension to 'fna' used as import file in QIIME2
#install miniconda3 in Unix or Linxus terminal and install QIIME2 in your conda environment  

#execute terminal#
#activate conda environment 
source activate qiime2-2019.10

#### import fasta files into qiime2 format 
qiime tools import \
  --input-path /InputFiles/Seqs100.fna  \
  --output-path /Your directory/Seqs100.qza \
  --type 'FeatureData[Sequence]'

#### 1_The first step to build phylogenetic tree is to alignment the sequences 
qiime alignment mafft\
 --i-sequences /Your directory/Seqs100.qza \
 --o-alignment /Your directory/aligned-rep-seps.qza 
 
##Time to run : 1 minute

##### 2_Masking sites
qiime alignment mask \
  --i-alignment /Your directory/aligned-rep-seps.qza \
  --o-masked-alignment /Your directory/masked-aligned-rep-seps.qza 
  
#### 3_Creating a tree
qiime phylogeny fasttree \
  --i-alignment /Your directory/masked-aligned-rep-seps.qza \
  --o-tree /Your directory/unrooted-tree.qza
  
#### 4_Midpoiting rooting
#Fastree creates an unrooted tree.

qiime phylogeny midpoint-root \
  --i-tree /Your directory/unrooted-tree.qza \
  --o-rooted-tree /Your directory/rooted-tree.qza

qiime tools export \
 --input-path /Your directory/rooted-tree.qza \
 --output-path /Your directory/epxorted-100tree
 
 ##a tree file is in the folder 'exported-100tree'in your local directory
 #this tree file is called 'Smalltree.nwk'

#To visualize the phylogenetic tree
#import the tree files 'Smalltree.nwk' into iTOL tool (https://itol.embl.de/)



###########################################################################
#####Plotting heatmaps with the relative abundance of 100 ASVs across samples
###########################################################################

###load libraries for this section
library("gplots"); packageVersion("gplots")  # for plotting heatmap
library("dichromat"); packageVersion("dichromat") # for creating color palette 

#load dataframe
otu100 <-read.csv("/InputFiles/otu100_noCTR.csv", header=T, row.names=1)
#transpose the dataframe
totu100 <-t(otu100) ## ASVs in columns, sample in row


#pearson correlation dissimilarity matrix
row_dist_algae <-as.dist(1-cor(t(totu100), method="pearson"))
col_dist_algae <-as.dist(1-cor(totu100, method="pearson"))

#Clustering 
col_hc_algae <-hclust(col_dist_algae, method="average")
row_hc_algae <-hclust(row_dist_algae, method="average")

#display as dendrogram
col_d_algae <-as.dendrogram(col_hc_algae)
row_d_algae <-as.dendrogram(row_hc_algae)

#make the color for heatmap
my_palette<-colorRampPalette(c("royalblue4","blue","white","pink","tomato"))((n=800))

#plot heatmap 
p6<-heatmap.2 (t(totu100), Rowv=as.dendrogram(col_hc_algae), col=my_palette, dendrogram="row", Colv=NA, density.info="none",trace="none", margins=c(5,18), key=TRUE,scale="column",keysize=0.9, srtCol = 48, cexCol=0.65, cexRow=0.48, lhei=c(0.5,5))
#visualize the panel
p6


#############################################################
#Bubbleplot for significant ASVs at the family level ########
#############################################################
Taxosign <-read.csv("/InputFiles/FamilySig.csv", header=T)

#reorder the family level for plotting
Taxosign$Family <-factor(Taxosign$Family, levels=c("Nitrincolaceae","Pseudomonadaceae","Methylophagaceae","Hyphomonadaceae","Rhodobacteraceae","Burkholderiaceae","Rhodocyclaceae","Sporichthyaceae", "Corynebacteriales order", "Ilumatobacteraceae","Flavobacteriaceae", "Crocinitomicaceae","JGI 0000069-P22 order", "Candidatus Kaiserbacteria order", "Candidatus Pacebacteria order"))

#reorder the treatment ID for plotting
Taxosign$Microcosms <-factor(Taxosign$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C","20D_80C", "100C"))

#Bubble plot
p7<-ggplot(Taxosign, aes(x=Microcosms , y=Family, size=Abundance, fill=Family)) + 
geom_point(color="black", alpha=0.80, shape=21) + 

#manual set the size of the bubble (or big cicles)
scale_size(range=c(1.5,22), name="Mean relative abundance (% of total reads)", breaks=c(0, 0.5, 2, 5, 12)) +

#manual set the color code for each taxa displayed
scale_fill_manual(values=c("Corynebacteriales order"= "#FF9900","Sporichthyaceae" = "#6699FF", "Hyphomonadaceae" = "#FFFF00", "Ilumatobacteraceae" = "bisque", "Rhodobacteraceae" = "#FF0000", "Burkholderiaceae" = "#660099", "Rhodocyclaceae" ="#666600", "Methylophagaceae"="#99FF99", "Pseudomonadaceae"="#663300", "Nitrincolaceae"="#FFCCFF", "Crocinitomicaceae"="#00CCFF","Flavobacteriaceae"="#006600","JGI 0000069-P22 order"="#FF33CC", "Candidatus Kaiserbacteria order"="#000000","Candidatus Pacebacteria order"="#999999"), guide=FALSE)+
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.background =element_rect(fill = 'white', colour = 'black')) + 
theme(axis.text.x =element_text(colour = 'black',size = 11, angle = 25,hjust = 1,family='serif' )) +
theme(axis.text.y=element_blank()) + 
theme(axis.title.x=element_text(colour = 'black',size = 12, family='serif')) + 
theme(axis.title.y=element_blank()) + 
labs(size="NA", y="Family level", x="Microcosms") + 
theme(legend.title=element_blank()) +  
theme(legend.text=element_text(size=10, family='serif')) +
theme(axis.ticks.y = element_blank())

#visualize the plot
p7


########################################################################
###########plot overall taxonomic composition the class and family level
########################################################################

######Class level  >1%  relative abundance collectively across samples
#remove the text angle=30, hjust=1
TaxoComp1 <-read.csv("/InputFiles/Class.csv", header=T)

#reorder taxa
TaxoComp1$Class <-factor(TaxoComp1$Class, levels=c("Acidimicrobiia", "Actinobacteria",	"Thermoleophilia",	"Alphaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria"	, "Campylobacteria"	,	"Bacteroidia",	"Verrucomicrobiae",	"Parcubacteria", "Gracilibacteria", "Other"))

#reorder treatment ID
TaxoComp1$Microcosms <-factor(TaxoComp1$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C","20D_80C", "100C"))

#set individual samples for plotting
TaxoComp1$Sample <-factor(TaxoComp1$Sample)

#plotting##
p8<-ggplot(TaxoComp1, aes(factor(Sample), Relative, fill=Class, order=as.numeric(Class))) + 
geom_bar(stat="identity", width=0.64) + 
scale_fill_manual(values=c("Acidimicrobiia"="midnightblue", "Actinobacteria"="steelblue3","Thermoleophilia"="lightskyblue2","Alphaproteobacteria"="darkred","Deltaproteobacteria"="khaki2","Gammaproteobacteria"="darkslategray", "Campylobacteria"= "forestgreen", "Bacteroidia"="darkseagreen2", "Verrucomicrobiae"="burlywood4", "Parcubacteria"="plum1","Gracilibacteria"= "darkmagenta", "Other"="lightgrey"), name="Class level")  + 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + 
theme(panel.background =element_rect(fil="white", color="black")) + 
theme(axis.text.x =element_blank()) + 
theme(axis.text.y=element_text(colour = 'black', size = 11, family='serif')) + 
theme(axis.title.x=element_text(colour = 'black',size = 11, family='serif')) + 
theme(axis.ticks.x = element_blank()) + 
theme(axis.title.y=element_text(colour = 'black',size = 11, family='serif')) + 
theme(plot.title=element_text(size=11, colour='black', face='bold', family='serif')) + 
theme(legend.title=element_text(size=10, family='serif')) +  
theme(legend.text=element_text(size=10, family='serif', face='italic')) + 
theme(legend.key.size=unit("0.38","cm")) + 

#use 'face_grid' function to separate treatments
facet_grid(. ~ Microcosms, scales="free", space="free_x") + 
labs(title="", y="Relative abundance (% of total sequence reads)", x="") + 
theme(strip.text=element_text(size=9,face='bold', family='serif')) + 
theme(strip.background=element_rect(size=1)) + 
theme(panel.margin=unit(0.18,"line")) + 
theme(strip.text.x=element_text(margin=margin(0.09,0,0.09,0,"cm"))) + 
theme(plot.margin=margin(0.03,0,0.05,0.05,"cm"))

#visualize plot
p8



######Family level  >3%  relative abundance collectively across samples
TaxoComp2 <-read.csv("/InputFiles/Family.csv", header=T)
TaxoComp2$Family <-factor(TaxoComp2$Family, levels=c("Mycobacteriaceae",	"Sporichthyaceae",	"Caulobacteraceae",	"Rhodobacteraceae",	"SAR11_Clade III",	"Bacteriovoracaceae",	"Moraxellaceae",	"Methylophagaceae",	"Pseudomonadaceae",	"Nitrincolaceae",	"Porticoccaceae",	"Saccharospirillaceae", "Colwelliaceae", "Burkholderiaceae",	"Arcobacteraceae", "Crocinitomicaceae", "NS11-12 marine group", "Flavobacteriaceae", "JGI 0000069-P22 order", "Other"))
TaxoComp2$Microcosms <-factor(TaxoComp2$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C","20D_80C", "100C"))
TaxoComp2$Sample <-factor(TaxoComp2$Sample)

#plotting
p9<-ggplot(TaxoComp2, aes(factor(Sample), Relative, fill=Family, order=as.numeric(Family))) + 
geom_bar(stat="identity", width=0.64) + 
scale_fill_manual(values=c("Mycobacteriaceae"="aquamarine4",	"Sporichthyaceae"="aquamarine1","Caulobacteraceae"="bisque4",	"Rhodobacteraceae"="bisque",	"SAR11_Clade III"="forestgreen",	"Bacteriovoracaceae"="tomato1",	"Moraxellaceae"="darkgoldenrod1",	"Methylophagaceae"="darkgoldenrod4",	"Pseudomonadaceae"="deepskyblue4",	"Nitrincolaceae"="deepskyblue1",	"Porticoccaceae"="darkred",	"Saccharospirillaceae"="darkorchid4", "Colwelliaceae"="darkorchid1", "Burkholderiaceae"="darksalmon",	"Arcobacteraceae"="darkslateblue", "Crocinitomicaceae"="deeppink1", "NS11-12 marine group"="gray80", "Flavobacteriaceae"="khaki1", "JGI 0000069-P22 order"="lightpink", "Other"="black"), name="Family level")  + 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + 
theme(panel.background =element_rect(fill="white", color="black")) + 
theme(axis.text.x =element_blank()) + 
theme(axis.text.y=element_text(colour = 'black', size = 11, family='serif')) + 
theme(axis.title.x=element_text(colour = 'black',size = 11, family='serif')) + 
theme(axis.ticks.x = element_blank()) + 
theme(axis.title.y=element_text(colour = 'black',size = 11, family='serif')) + 
theme(plot.title=element_text(size=11, colour='black', face='bold', family='serif')) + 
theme(legend.title=element_text(size=10, family='serif')) +  
theme(legend.text=element_text(size=10, family='serif', face='italic')) + 
theme(legend.key.size=unit("0.38","cm")) + 
facet_grid(. ~ Microcosms, scales="free", space="free_x") + 
labs(title="", y="Relative abundance (% of total sequence reads)", x="") + 
theme(strip.text=element_text(size=9,face='bold', family='serif')) + 
theme(strip.background=element_rect(size=1)) + 
theme(panel.margin=unit(0.18,"line")) + 
theme(strip.text.x=element_text(margin=margin(0.09,0,0.09,0,"cm"))) + 
theme(plot.margin=margin(0.03,0,0.05,0.05,"cm"))

#visualize 
p9


########################################################
###########Plotting Slurries taxonomic composition
########################################################

#load libraries for this section
library(phyloseq)
library(ggplot2)
library("gridExtra");packageVersion("gridExtra") #v 2.3

#load Slurries phyloseq
Slurries <-readRDS("InputFiles/Slurries.rds")
SlurriesTaxa <-read.csv("/InputFiles/SlurriesTaxa.csv", header=T, row.names=1)
Slurriestaxa <-rownames(SlurriesTaxa)
Slurries1 <-prune_taxa(Slurriestaxa, Slurries)

#plotting
p11<-plot_bar(Slurries1, fill="class", title="Class level") + theme(axis.title.y=element_text(colour = 'black',size = 9, family="serif")) + 
             theme(axis.title.x=element_text(colour = 'black', size = 12,family="serif" )) + 
             theme(axis.text.y=element_text(colour = 'black', size = 10, family="serif")) + 
             theme(plot.title=element_text(colour= 'black', size=12, family="serif")) + 
             theme(axis.text.x=element_text(colour = 'black', size = 10, family="serif")) +
             theme(legend.title=element_text(size=11, family='serif')) +  
            theme(legend.text=element_text(size=11, family='serif', face='italic')) +
            theme(legend.key.size=unit("0.38","cm"))


Alpha= subset_taxa(Slurries1, class=="Alphaproteobacteria")
p12 <-plot_bar(Alpha, "Sample", "Abundance", "family", title="Alphaproteobacteria") + theme(axis.title.y=element_text(colour = 'black',size = 9, family="serif")) + 
             theme(axis.title.x=element_text(colour = 'black', size = 12,family="serif" )) + 
             theme(axis.text.y=element_text(colour = 'black', size = 10, family="serif")) + 
             theme(plot.title=element_text(colour= 'black', size=12, family="serif")) + 
             theme(axis.text.x=element_text(colour = 'black', size = 10, family="serif")) +
             theme(legend.title=element_text(size=11, family='serif')) +  
            theme(legend.text=element_text(size=11, family='serif', face='italic')) +
            theme(legend.key.size=unit("0.38","cm"))



Gamma= subset_taxa(Slurries1, class=="Gammaproteobacteria")
p13 <-plot_bar(Gamma, "Sample", "Abundance", "family", title="Gammaproteobacteria") + theme(axis.title.y=element_text(colour = 'black',size = 9, family="serif")) + 
             theme(axis.title.x=element_text(colour = 'black', size = 12,family="serif" )) + 
             theme(axis.text.y=element_text(colour = 'black', size = 10, family="serif")) + 
             theme(plot.title=element_text(colour= 'black', size=12, family="serif")) + 
             theme(axis.text.x=element_text(colour = 'black', size = 10, family="serif")) +
             theme(legend.title=element_text(size=11, family='serif')) +  
            theme(legend.text=element_text(size=11, family='serif', face='italic')) +
            theme(legend.key.size=unit("0.38","cm"))


Camy= subset_taxa(Slurries1, class=="Campylobacteria")
p14 <-plot_bar(Camy, "Sample", "Abundance", "family", title="Campylobacteria") + theme(axis.title.y=element_text(colour = 'black',size = 9, family="serif")) + 
             theme(axis.title.x=element_text(colour = 'black', size = 12,family="serif" )) + 
             theme(axis.text.y=element_text(colour = 'black', size = 10, family="serif")) + 
             theme(plot.title=element_text(colour= 'black', size=12, family="serif")) + 
             theme(axis.text.x=element_text(colour = 'black', size = 10, family="serif")) +
             theme(legend.title=element_text(size=11, family='serif')) +  
            theme(legend.text=element_text(size=11, family='serif', face='italic')) +
            theme(legend.key.size=unit("0.38","cm"))


Graci= subset_taxa(Slurries1, class=="Gracilibacteria")
p15 <-plot_bar(Graci, "Sample", "Abundance", "order", title="Gracilibacteria") + theme(axis.title.y=element_text(colour = 'black',size = 9, family="serif")) + 
             theme(axis.title.x=element_text(colour = 'black', size = 12,family="serif" )) + 
             theme(axis.text.y=element_text(colour = 'black', size = 10, family="serif")) + 
             theme(plot.title=element_text(colour= 'black', size=12, family="serif")) + 
             theme(axis.text.x=element_text(colour = 'black', size = 10, family="serif")) +
             theme(legend.title=element_text(size=11, family='serif')) +  
            theme(legend.text=element_text(size=11, family='serif', face='italic')) +
            theme(legend.key.size=unit("0.38","cm"))



Bacteroidia= subset_taxa(Slurries1, class=="Bacteroidia")
p16 <-plot_bar(Bacteroidia, "Sample", "Abundance", "family", title="Bacteroidia") + theme(axis.title.y=element_text(colour = 'black',size = 9, family="serif")) + 
             theme(axis.title.x=element_text(colour = 'black', size = 12,family="serif" )) + 
             theme(axis.text.y=element_text(colour = 'black', size = 10, family="serif")) + 
             theme(plot.title=element_text(colour= 'black', size=12, family="serif")) + 
             theme(axis.text.x=element_text(colour = 'black', size = 10, family="serif")) +
             theme(legend.title=element_text(size=11, family='serif')) +  
            theme(legend.text=element_text(size=11, family='serif', face='italic')) +
            theme(legend.key.size=unit("0.38","cm"))

####put all panels together in a figure
pdf("/Your directory/FigSlurri.pdf", width =12, height=7.5) #Change to your directory 
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3, nrow=2)
dev.off()


###############################################################
#Plotting Pearson's correlation between NO3 and NH4 and PO43
###############################################################

#load libraries for this section
library("ggpubr")
library("gridExtra")

#load dataframe
NutrientNO_PO <-read.csv("/InputFiles/NOPOtestout.csv", header=T, row.names=1)
#one outline of NO4 near 40 umol /L was removed 

#plot
p17<-ggscatter(NutrientNO_PO, x = "NOX", y = "PO4",  
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson", size=3,
            ylab = (expression(paste("PO4 concentration", " ",  mu, "M"))), xlab = (expression(paste("NOx concentration"," ", mu, "M"))), cor.coef.size = 6, font.family="serif") + 
            labs(title="Pearson's correlation") + 
            theme(axis.title.y=element_text(colour = 'black',size =17.5, family="serif")) + 
            theme(axis.title.x=element_text(colour = 'black', size = 17.5,family="serif" )) + 
            theme(axis.text.y=element_text(colour = 'black', size = 17.5, family="serif")) + 
            theme(plot.title=element_text(colour= 'black', size=19, family="serif")) + 
            theme(axis.text.x=element_text(colour = 'black', size = 17.5, family="serif"))



##################################################
#Plotting bacterial abundance
#################################################
##load library for this section
library(ggplot2)

#load bacterial abundance data with mean and standard deviation
BAMean <-read.csv("/InputFiles/BAMean.csv", header=TRUE, row.names=1)

#reorder the treatment ID for plotting
BAMean$Treatment <-factor(BAMean$Treatment, levels=c("Control", "100D", "80D_20C", "50D_50C", "20D_80C", "100C"))

#reorder the time point for plotting
BAMean$Time <-factor(BAMean$Time, levels=c("Day7", "Day24"))

#plotting
P18 <-ggplot(BAMean, aes(x=Treatment, y=Mean, fill=Treatment)) + 
geom_bar(position="dodge", stat = "identity") +

#add the error bar for standard deviation 
geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(.9)) + 
#manual fill the color of the bars
scale_fill_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen"))+ 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + 
theme(panel.background =element_rect(fill = 'white', colour = 'black')) + 
theme(axis.text.x =element_text(colour = 'black',size = 8, angle = 25,hjust = 1)) +
theme(axis.text.y=element_text(colour = 'black', size = 8)) + 
theme(axis.title.x=element_text(colour = 'black',size = 9, face = 'bold')) + 
theme(axis.title.y=element_text(colour = 'black',size = 9, face = 'bold')) + 
labs(y="Bacterial abundance(cells/ml)", x="Control & treatments", color="Incubation") + 
theme(plot.title=element_text(size=10.5, face = "bold")) + 
theme(legend.title=element_text(size=8)) +  
theme(legend.text=element_text(size=8))  + 
facet_grid(. ~ Time, scales="free", space="free_x") + 
strip.text=element_text(size=9,face='bold')) + 
theme(strip.background=element_rect(size=1)) + 
theme(panel.margin=unit(0.18,"line")) + 
theme(strip.text.x=element_text(margin=margin(0.09,0,0.09,0,"cm")))



####################################################
#NMDS plot Bray-Curtis using Absolute abundance
####################################################
##load libraries for this section 
library(ggplot2)
library(plyr) # for hull

##load phyloseq object
AbsoPhyloObject <-readRDS("/InputFiles/AbsoPhyloObject.rds")

#run Bray-Curtis dissimilarity matrix in 
nmds1<-ordinate(AbsoPhyloObject, "NMDS", "bray")
nmds1 $ stress

#create dataframe for plot
Mymeta = data.frame (Microcosms=c("Control","Control","Control", "Control",	"Control","100D",	"100D", "100D","100D", "80D_20C", "80D_20C", "80D_20C","80D_20C", "80D_20C", "50D_50C","50D_50C",	"50D_50C"	,"50D_50C",	"50D_50C",	"20D_80C", "20D_80C","20D_80C", "20D_80C" ,"20D_80C" ,"100C", "100C",	"100C" , "100C" ,"100C"))
NMDS = data.frame (NMDS1=nmds1$points[,1], NMDS2=nmds1$points[,2], Microcosms=Mymeta$Microcosms)
NMDS$Microcosms <-factor(NMDS$Microcosms, levels=c("Control","100D", "80D_20C", "50D_50C", "20D_80C", "100C"))

##Cretae dataset for polygon
find_hull <-function(NMDS) NMDS[chull(NMDS$NMDS1,NMDS$NMDS2),]
hulls <-ddply(NMDS, "Microcosms",find_hull) 

#Bray-Curtis NMDS########################
p19 <-ggplot(data=NMDS, aes(NMDS1, NMDS2, fill=Microcosms, color=Microcosms)) + 
geom_point(size=2.5) + geom_polygon(data = hulls, aes(x=NMDS1, y=NMDS2), alpha=0.3) + 
scale_fill_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen")) + 
scale_color_manual(values=c("deepskyblue3","tan4","orange1","gold1", "olivedrab2", "darkgreen"))  + 
theme_bw()+ theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+ 
theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
labs(title="Bray-Curtis NMDS based on Absolute abundance") + 
theme(axis.title.y=element_text(colour = 'black',size = 14,family='serif')) + 
theme(axis.title.x=element_text(colour = 'black',size = 14, family='serif')) + 
theme(axis.text.y=element_text(colour = 'black', size = 14,family='serif' )) + 
theme(axis.text.x=element_text(colour = 'black', size = 14, family='serif')) + 
theme(legend.text=element_text(size=13,family="serif")) + 
theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
theme(legend.background=element_rect(fill=NA)) + 
theme(legend.title=element_text(color="black", size=14, face='bold', family='serif')) + 
guides(color=guide_legend(override.aes=list(shape=15, size=3))) + 
guides(shape=guide_legend(override.aes=list(size=2.5))) + 
theme(legend.key.size=unit("0.5","cm")) + 
theme(plot.title = element_text(colour='black', size = 18.5, face = 'bold', family="serif"))

#visualize
p19


#############################################
#Plotting overlap between Sediment and water
#############################################

#load libraries for this section
library('VennDiagram'); packageVersion('VennDiagram')
#'1.6.20'

#load water and sediment dataset
Water<-read.csv("/InputFiles/waterFinalcount.csv",header=T, row.names=1)
Sediment<-read.csv("/InputFiles/SedFinalcount.csv",header=T, row.names=1)

#check the number of ASVs
nrow(Water)
#3692
nrow(Sediment)
#6570

#analyze the overlap between the two
overlapp<-list(intersect(Water[,1], Sediment[,1]))

#check the number of ASVs overlapped
overlapp 
#547

#draw Venndiagram
p20 <-draw.pairwise.venn(area=3692,area2=6570,cross.area=547,
                          category = c("Water", "Sediment"),fill = c("deepskyblue3", "chocolate4"),lty = "blank",cex = 2,cat.cex = 2,cat.pos = c(140, 50),cat.dist = 0.09,cat.just = list(c(-1, -0.04), c(0.07, 1.8)),ext.pos = 10,ext.dist = -0.05,ext.length = 0.85,ext.line.lwd = 2,ext.line.lty = "dashed")



#############################################################
#Selection of Rare and abundant from the sediment-water overlap
#abundant (>0.1% relative abundance)
#rare (<0.1 % relative abundance) ###########
#############################################################

#load the table containing both sediment and water overlaps
WaterSed <-read.csv("/InputFiles/WaterSed.csv", row.names=1, header=T)
nrow(WaterSed)
#[1] 534

###Select abundant ASVs from sediment datasets
Sedabunt <-subset(WaterSed, Sed_T0>0.1 & Sed_CTR>0.1 & Sed_100D>0.1 & Sed_80D>0.1 & Sed_50D >0.1 & Sed_20D >0.1 & Sed_100C>0.1) 
nrow(Sedabunt)
#8

#save table
write.csv(Sedabunt, file="/Your directory/Sedabunt.csv")

#Select rare ASVs from sediment datasets
Sedrare <-subset(WaterSed, !rownames(WaterSed) %in% rownames(Sedabunt))
nrow(Sedrare)
#526
#save table
write.csv(Sedrare, file="/Your directory/Sedrare.csv")


###Select abundant ASVs from water datasets
Watabunt <-subset(WaterSed, Wat_CTR>0.1 & Wat_100D>0.1 & Wat_80D>0.1 & Wat_50D >0.1 & Wat_20D >0.1 & Wat_100C>0.1) 
nrow(Watabunt)
#19
write.csv(Watabunt, file="/Your directory/Watabunt.csv")

###Select rare ASVs from water datasets
Watrare <-subset(WaterSed, !rownames(WaterSed) %in% rownames(Watabunt))
nrow(Watrare)
#515
write.csv(Watabunt, file="/Your directory/Watrare.csv")



####################################################
#taxonomic composition of the overlap
####################################################
library("ggplots2")
library("RColorBrewer")

###load file 
#Class level  >1%  relative abundance collectively across samples
MergedClass <-read.csv("/InputFiles/MergedClass.csv", header=T)

#reorder the class level
MergedClass$Class <-factor(MergedClass$Class, levels=c("Acidimicrobiia", "Actinobacteria","Alphaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria","Bacteroidia","Campylobacteria", "Clostridia",  "Gracilibacteria", "Nitrospira", "Oxyphotobacteria", "Other"))

#reorder treatment ID
MergedClass$Microcosms <-factor(MergedClass$Microcosms, levels=c("Sed_T0", "Sed_Control", "Sed100D", "Sed80D_20C", "Sed50D_50C", "Sed20D_80C", "Sed100C", "Wat_Control","Wat100D", "Wat80D_20C", "Wat50D_50C","Wat20D_80C", "Wat100C"))

#reorder for sediment and water plot
MergedClass$Phase <-factor(MergedClass$Phase, levels=c("Sediment" ,"Water"))

#Set the colors
display.brewer.pal(n=12, name="Paired")

#plot barplot
p21<-ggplot(MergedClass, aes(factor(Microcosms), Abundance, fill=Class, order=as.numeric(Class))) + 
geom_bar(stat="identity", width=0.64) + 
scale_fill_manual(values=c(brewer.pal(n=12, name="Paired")), name="Class level")  + 
theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + 
theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + 
theme(panel.background =element_rect(fil="white", color="black")) + 
theme(axis.text.x =element_text(colour = 'black', size = 13, family='serif', angle=30, hjust=1)) + 
theme(axis.text.y=element_text(colour = 'black', size = 14, family='serif')) + 
theme(axis.title.x=element_text(colour = 'black',size = 14, family='serif')) + 
theme(axis.ticks.x = element_blank()) + 
theme(axis.title.y=element_text(colour = 'black',size = 14, family='serif')) + 
theme(plot.title=element_text(size=13, colour='black', face='bold', family='serif')) + 
theme(legend.title=element_text(size=13, family='serif')) +  
theme(legend.text=element_text(size=13, family='serif', face='italic')) + 
theme(legend.key.size=unit("0.30","cm")) + 
facet_grid(. ~ Phase, scales="free", space="free_x") + 
labs(title="", y="Relative abundance (% of total sequence reads)", x="") + 
theme(strip.text=element_text(size=14,face='bold', family='serif')) + 
theme(strip.background=element_rect(size=1)) + 
theme(panel.margin=unit(0.28,"line")) + 
theme(strip.text.x=element_text(margin=margin(0.09,0,0.09,0,"cm"))) + 
theme(plot.margin=margin(0.03,0,0.05,0.05,"cm"))


###Congratulations! You've scuccessfully recreate all of the figures in the paper. Hope you have fun!
###Feel free to contact me if you have any issues with scripts when performing data visualization!
######################################## END #######################################################
###For a taste of the full analysis pipeline and to recreate all of the figures in the paper, check out this public workflow:
#####
