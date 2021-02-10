###########################
### R scripts for plotting and data visualization 
### Data from 
### "Quality of phytoplankton deposition structures bacterial communities at the water-sediment interface"
### By D Izabel-Shen, S Albert, M Winder, H Farnelid and F Nascimento
### Prepared 9 February 2021
### Author: Dandan Izabel-Shen, Stockholm University; dand.shen <at> gmail.com

# Before you start
# Make sure you are using the latest version of R (and Rstudio)
#The following packages (and their dependencies) are needed to run the whole analysis 
#stats
#agricolae
#car
#picante
#vegan
#phyloseq
#DESeq2
#BiocManager



########################################################
# One way-ANOVA test 
#the effect of POM mixture on the bacterial abundance 
#on each time point (Day 7 and Day 24)
########################################################

#load libraries for this section 
library(stats)
library(agricolae)
library(car)

#load data
BA <-read.csv("/Users/Dandan/Desktop/WaterSediments/Water/BA/BA.csv", header=T, row.names=1)

#set time and treatment as factors for analysis
Time <-as.factor(BA$Time)
Treatment <-as.factor(BA$Treatment)

##Repeated measures ANOVA to test the time effects 
##Error(Sample) for random effect
BARepteat = aov(log ~ Treatment* Time + Error(Sample), data=BA)

#check the statistics output
summary(BARepteat)


#######subset Day 7 dataset#################
#subset Day7 data from the whole dataset
BA_Day7 <-BA[c(1:30), c(1:5)]

#set treatment as factor
Treatment <-as.factor(BA_Day7$Treatment)

#use linear model 
Fit1 <-lm(BA_Day7$log ~ Treatment, data=BA_Day7)

#run normality test and check if the residuls of linear model are normally distributed
shapiro.test(Fit1$residuals) # Analysis of Variance Table #
# Shapiro-Wilk normality test

##Since the normality assumption is met, proceed with one-way ANOVA test
Day7 <-anova(Fit1)

##check the ANOVA test results 
Day7

#check which level (treatment) the significance was detected, the post hoc 
Posthos1 <-HSD.test(lm(BA_Day7$log ~ Treatment, data=BA_Day7), 'Treatment', alpha=0.1)
Posthos1


#######subset Day 24 dataset#################
BA_Day24 <-BA[c(31:60), c(1:5)]
Treatment <-as.factor(BA_Day24$Treatment)
Fit2 <-lm(BA_Day24$log ~ Treatment, data=BA_Day24)

shapiro.test(Fit2$residuals) # Analysis of Variance Table #
###Shapiro-Wilk normality test

#test the homogeneity of samples
leveneTest(Fit2)

#Met ANOVA assumptions, proceed with analysis 
Day24 <-anova(Fit2)
Day24

#proceed with Post hoc test
Posthoc2 <-HSD.test(lm(BA_Day24$log ~ Treatment, data=BA_Day24), 'Treatment', alpha=0.1)
Posthoc2




########################################################
#One-way ANOVA with linear model on alpha diversity
#######################################################

#load library for this section 
library(stats)
library(car)

#load alpha diversity data
Alpha <-read.csv("/InputFiles/Alpha.csv",header=T, row.names = 1)

#select each diversity index from the dataset
Richness <-Alpha[c(1:29), (1:7)]
Evenness <-Alpha[c(30:58), (1:7)]

#set the factors for test
Microcosms_Richness <-as.factor(Richness$Microcosm)
Microcosms_Evenness <-as.factor(Evenness$Microcosm)

## test on Richness
R1<-lm(Richness$log ~Microcosms_Richness, data=Richness)

#normality test
shapiro.test(R1$residuals)
#homogeneity test 
leveneTest(R1)

# ANOVA assumption is met, proceed with anova test in linear model
anova(lm(Richness$log ~Microcosms_Richness, data=Richness))


## test on Evenness
E1<-lm(Evenness$log ~Microcosms_Evenness, data=Evenness)

#normality test
shapiro.test(p5$residuals)
#homogeneity test 
leveneTest(p5)

#ANOVA assumption is met,proceed with anova test in linear model
anova(lm(Evenness$log ~Microcosms_Evenness, data=Evenness))




#################################################
#One way ANOVA on nutrients analyses
#similar analysis as done on Bacterial abundance
#################################################

#load libraries for this section 
library(stats)
library(agricolae)
library(car)

#read dataset 
Nutrient <-read.csv("/InputFiles/Nutrient.csv", header=T, row.names=1)

#repeated measurement with time effects..
Time <-as.factor(Nutrient$Time)
Treatment <-as.factor(Nutrient$Treatment)

#Extract for NH4
NH <-Nutrient[1:60,]
NH$Time <-as.factor(NH$Time)
NH$Treatment <-as.factor(NH$Treatment)

#Extract for PO4
PO <-Nutrient[61:120,]
PO$Time <-as.factor(PO$Time)
PO$Treatment <-as.factor(PO$Treatment)

#Extract for NOx
NO <-Nutrient[121:180,]
NO$Time <-as.factor(NO$Time)
NO$Treatment <-as.factor(NO$Treatment)

##Run the repeat measurement ANOVA
#a mixed-effect model is NOT, as the design is balanced
#Error(Sample) for random effect_ NH
T1 = aov(Value ~ Treatment * Time + Error(Sample), data=NH)
#check results
summary(T1)

#Error(Sample) for random effect_ PO
T2 = aov(Value ~ Treatment * Time + Error(Sample), data=PO)
summary(T2)

#Error(Sample) for random effect_NO
T3 = aov(Value ~ Treatment * Time + Error(Sample), data=NO)
summary(T3)


#######
#####Test the differences on nutrient for Day 0 and day 24 separately
######
## subset Day T0dataset
#######################

########
##NO data
########
NO_T0 <-NO[c(1:30), c(1:6)]
#set treatment factor
Treatment <-as.factor(NO_T0$Treatment)
#linear model 
Fit6 <-lm(NO_T0$Value ~ Treatment, data=NO_T0)

#Shapiro-Wilk normality test
shapiro.test(Fit1$residuals) # Analysis of Variance Table #
#Homogeneity test
leveneTest(Fit1)

#ANOVA assumption is met, proceed with test
NO_T0test <-anova(Fit6)
HSD.test(lm(NO_T0$Value ~ Treatment, data=NO_T0), 'Treatment', alpha=0.1)

##########
###PO data
##########
PO_T0 <-PO[c(1:30), c(1:6)]
#set treatment factor
Treatment <-as.factor(PO_T0$Treatment)
#linear model 
Fit7 <-lm(PO_T0$Value ~ Treatment, data=PO_T0)

#Shapiro-Wilk normality test
shapiro.test(Fit7$residuals) # Analysis of Variance Table #
#Homogeneity test
leveneTest(Fit7)

#ANOVA assumption is met, proceed with test
PO_T0test <-anova(Fit7)
HSD.test(lm(PO_T0$Value ~ Treatment, data=PO_T0), 'Treatment', alpha=0.1)

#########
##NH data
########
NH_T0 <-NO[c(1:30), c(1:6)]
#set treatment factor
Treatment <-as.factor(NH_T0$Treatment)
#linear model 
Fit8 <-lm(NH_T0$Value ~ Treatment, data=NH_T0)

#Shapiro-Wilk normality test
shapiro.test(Fit8$residuals) # Analysis of Variance Table #
#Homogeneity test
leveneTest(Fit8)

#ANOVA assumption is met, proceed with test
NH_T0test <-anova(Fit8)
HSD.test(lm(NH_T0$Value ~ Treatment, data=NH_T0), 'Treatment', alpha=0.1)


#########################
## subset Tf (day 24) dataset
#######################

##########
##NO data
########
NO_Tf <-NO[c(31:60), c(1:6)]
Treatment <-as.factor(NO_Tf$Treatment)
Fit9 <-lm(NO_Tf$log ~ Treatment, data=NO_Tf)
shapiro.test(Fit2$residuals) # Analysis of Variance Table
leveneTest(Fit9)
NO_Tftest <-anova(Fit9)
HSD.test(lm(NO_Tf$log ~ Treatment, data=NO_Tf), 'Treatment', alpha=0.1)

##########
##PO data
########
PO_Tf <PO[c(31:60), c(1:6)]
Treatment <-as.factor(PO_Tf$Treatment)
Fit10 <-lm(PO_Tf$log ~ Treatment, data=PO_Tf)
shapiro.test(Fit10$residuals) # Analysis of Variance Table
leveneTest(Fit10)
PO_Tftest <-anova(Fit10)
HSD.test(lm(PO_Tf$log ~ Treatment, data=PO_Tf), 'Treatment', alpha=0.1)

##########
##NH data
########
NH_Tf <PO[c(31:60), c(1:6)]
Treatment <-as.factor(NH_Tf$Treatment)
Fit11 <-lm(NH_Tf$log ~ Treatment, data=NH_Tf)
shapiro.test(Fit11$residuals) # Analysis of Variance Table
leveneTest(Fit11)
NH_Tftest <-anova(Fit11)
HSD.test(lm(NH_Tf$log ~ Treatment, data=NH_Tf), 'Treatment', alpha=0.1)




#####################################################
# Beta community diversity analyses
#####################################################

#load libraries for this section 
library(vegan)
library(phyloseq)

#read phyloseq object 
newPhyloObject <-readRDS("/InputFiles/newPhyloObject.rds")

###Bray-Curtis distance matrix
newPhyloObjectRelat <-phyloseq::distance(transform_sample_counts(newPhyloObject, function(x) x/sum(x)*100), method="bray")
sampledf <- data.frame(sample_data(newPhyloObject))
#PERMANOVA
T1<-adonis(formula = newPhyloObjectRelat ~ Microcosm, data = sampledf) 

###UniFrac distance matrix
uniFrac_weighted <-UniFrac(newPhyloObject, weighted = TRUE, normalized = TRUE, parallel = FALSE, fast = TRUE)
uniFrac_unweighted <-UniFrac(newPhyloObject, weighted=FALSE, normalized = TRUE, parallel = FALSE, fast = TRUE)  #default weighted=FALSE

#PERMANOVA
Test2 <-adonis(uniFrac_weighted~ Microcosm, data = sampledf)
Test3 <-adonis(uniFrac_unweighted~ Microcosm, data = sampledf)

#Pairwise Mantel test for all three distance matrices
Mantel statistic based on Pearson's product-moment correlation 
##### Mantel test for all resemblances 
M1<- mantel(newPhyloObjectRelat, uniFrac_weighted, method="pearson", permutations=999, strata = NULL,
     na.rm = FALSE, parallel = getOption("mc.cores"))
M2<-mantel(newPhyloObjectRelat, uniFrac_unweighted, method="pearson", permutations=999, strata = NULL,
    na.rm = FALSE, parallel = getOption("mc.cores"))
M3<-mantel(uniFrac_weighted, uniFrac_unweighted, method="pearson", permutations=999, strata = NULL,
    na.rm = FALSE, parallel = getOption("mc.cores"))

#Return results 
M1
M2
M3



##########################################################
#####DESeq2 differential abundance analyses and plotting
##########################################################
#subset the phyloseq object 
#select treatments with sample ID diatom-enriched ASVs "100D", "80D_20C" vs. cyanobacteria-enriched ASVs "20D_80C", "100C"

#load libraries for this section 
library("DESeq2"); packageVersion("DESeq2") #v 1.26.0
library("phyloseq"); packageVersion("phyloseq") #v 1.30.0
library("ggplot2"); packageVersion("ggplot2") #v 3.3.2
library("BiocManager"); packageVersion("BiocManager") #v 1.30.10"

#Read phyloseq object
newPhyloObject <-readRDS("/InputFiles/newPhloObject.rds")

#select the data that for pairwise comparison
subset_DES <-subset_samples(newPhyloObject, Microcosm=="100D" | Microcosm=="80D_20C" | Microcosm=="20D_80C" | Microcosm=="100C")

#add the conditions DiatomPrefer CyanoPrefer to the levels to the metadata
conditions_meta <-read.csv("/InputFiles/conditions_meta.csv", header=T)
meta_DES <-sample_data(subset_DES)

#round the data to integer 
otu_table(subset_DES)<-round(otu_table(subset_DES))

#adding new variable to define formula of design function in DESeq2
meta_DES <-cbind(p1, conditions_meta) # 'Diatom_favored' as reference level, 'Cyano_favored'

#update the sample data for DESeq comparison
sample_data(subset_DES) <-met_DES
#check the updated sample data
sample_data(subset_DES)

#Convert to DESeq2 and call funciton##
alga_DES1 <-phyloseq_to_deseq2(subset_DES, ~conditions)

#Checking the levels of conditions and the order of the levels 
levels(alga_DES1$conditions)
## it should read "Cyano_favored"  "Diatom_favored"

#Setting Diatom_favored as reference group
alga_DES1$conditions <- relevel(alga_DES1$conditions, ref = "Diatom_favored")
levels(alga_DES1$conditions)   
#[1] "Diatom_favored" "Cyano_favored"

#Start differential abundance analysis
alga_DES1 = DESeq(alga_DES1, test="Wald", fitType="parametric")

#read the results of DESEq
res1 <-results(alga_DES1, cooksCutoff = FALSE)

#Benjamini-Hochberg prcedure with an adjusted alpha of 0.20
alpha=0.20

#extract adjusted P values smaller than alpha=0.2
sigtab1 = res1[which(res1$padj < alpha), ]

#calculate the number of ASVs met this requirement
sum(res1$padj < alpha, na.rm=TRUE )
#100 # 100 ASVs in this case

#rename taxonomic rank in the phyloseq object
colnames(tax_table(subset_DES)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")

#make a dataframe for DESeq output
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(subset_DES)[rownames(sigtab1), ], "matrix"))

#check the frist 50 rows of the dataframe
head(sigtab1)

#check the dimension of the dataframe
dim(sigtab1)
# 100 13  # 100ASVs and 13 columns for variables 

#save the output
write.csv(sigtab1, file="/Your directory/sigtab1.csv") # use your current directory

#save a table for DES results for all ASVs
res1 = cbind(as(res1, "data.frame"), as(tax_table(subset_DES)[rownames(res1), ], "matrix"))
write.csv(res1, file="/Your directory/res1.csv") #use your current directory


###check the clustering of the heatmaps 
##ANOSIM analysis test 
#whether there is a significant difference between two or more groups of sampling unites

#load library for this section
library(vegan)

#load dataframe
OTU100 <-read.csv("/InputFiles/otu100_noCTR.csv", header=T, row.names=1) 
### Rows are samples in this case (AVSs) and columns are response variables in this case are the relative abundance in each samples

##load grouping table 
ClusterTest <-read.csv("/InputFiles/ClustersTest.csv", header=T, row.names=1)

#select the grouping variable
Grouping <-as.factor(ClusterTest$Grouping)

#check the grouping levels
levels(Grouping)
#[1] "CyanoGr"  "DiatomGr" ### CyanoGr 72 ASVs, DiatomGr 28 ASVs

#Run ANOSIM 
anosim(OTU100, Grouping, permutations=999, distance="bray", strata=NULL)




###########################################
#Checking the phylogenetic relatedness
##########################################

#load library for this section
library(picante)

#########Using all taxa found in that sediments####
#read the tree files
Alltree <-read.tree("/InputFiles/tree.nwk")

Alltree.p <-as.phylo(Alltree)
All_dist <-cophenetic(Alltree.p)
Strategy<-read.table("/InputFiles/Strategy.txt", header=T, sep="\t", row.names=1)
destroyX = function(es) {
    f = es
    for (col in c(1:ncol(f))){ #for each column in dataframe
        if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
            colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
        }
    }
    assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

destroyX(Strategy)

#Net relatedness index
NRI <-ses.mpd (Strategy, All_dist, null.model="taxa.labels", abundance.weighted=FALSE, runs=999, iterations=1000)

#Nearest taxon index
NTI <-ses.mntd(Strategy, All_dist, null.model="taxa.labels",abundance.weighted=FALSE, runs=999,iterations=100)



###########################################
#Spearman's correlation
##########################################
#load libraries for this section
library(stats)

#load table containing relative abundance of bacterial families and nutrient data
FamilySig <-read.csv("/InputFiles/FamilySig.csv", header=T)

#Corynebacteriales order
FamilyCory <-FamilySig[1:6,]
FamilyCory1 <-cor.test(x=FamilyCory$NH4,y=FamilyCory$Abundance, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyCory2 <-cor.test(x=FamilyCory$PO4, y=FamilyCory$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyCory3 <-cor.test(x=FamilyCory$NOx, y=FamilyCory$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyCory1
FamilyCory2 
FamilyCory3

#Sporichthyaceae
FamilySpor<-FamilySig[7:12,]
FamilySpor1 <-cor.test(x=FamilySpor$NH4,y=FamilySpor$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilySpor2 <-cor.test(x=FamilySpor$PO4, y=FamilySpor$Abundance, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilySpor3 <-cor.test(x=FamilySpor$NOx, y=FamilySpor$Abundance, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilySpor1
FamilySpor2 
FamilySpor3

#Hyphomonadaceae
FamilyHyph<-FamilySig[13:18,]
FamilyHyph1 <-cor.test(x=FamilyHyph$NH4,y=FamilyHyph$Abundance, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyHyph2 <-cor.test(x=FamilyHyph$PO4, y=FamilyHyph$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyHyph3 <-cor.test(x=FamilyHyph$NOx, y=FamilyHyph$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyHyph1
FamilyHyph2 
FamilyHyph3

#Ilumatobacteraceae
FamilyIllum<-FamilySig[19:24,]
FamilyIllum1 <-cor.test(x=FamilyIllum$NH4,y=FamilyIllum$Abundance, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyIllum2 <-cor.test(x=FamilyIllum$PO4, y=FamilyIllum$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyIllum3 <-cor.test(x=FamilyIllum$NOx, y=FamilyIllum$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyIllum1
FamilyIllum2 
FamilyIllum3

#Rhodobacteraceae
FamilyRoh <-FamilySig[25:30,]
FamilyRoh1 <-cor.test(x=FamilyRoh$NH4,y=FamilyRoh$Abundance, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyRoh2 <-cor.test(x=FamilyRoh$PO4, y=FamilyRoh$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyRoh3 <-cor.test(x=FamilyRoh$NOx, y=FamilyRoh$Abundance, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyRoh1
FamilyRoh2 
FamilyRoh3

#Methylophagaceae
FamilyMethylo <-FamilySig[31:36,]
FamilyMet1 <-cor.test(x=FamilyMethylo$NH4, y=FamilyMethylo$Abundance,method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyMet2 <-cor.test(y=FamilyMethylo$Abundance, x=FamilyMethylo$PO4, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyMet3 <-cor.test(y=FamilyMethylo$Abundance, x=FamilyMethylo$NOx, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyMet1
FamilyMet2
FamilyMet3

#Nitrincolaceae
FamilyNitrin <-FamilySig[37:42,]
FamilyNi1 <-cor.test(x=FamilyNitrin$NH4, y=FamilyNitrin$Abundance,method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyNi2 <-cor.test(y=FamilyNitrin$Abundance, x=FamilyNitrin$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyNi3 <-cor.test(y=FamilyNitrin$Abundance, x=FamilyNitrin$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyNi1
FamilyNi2
FamilyNi3

#Pseudomonadaceae
FamilyPs <-FamilySig[43:48,]
FamilyPs1 <-cor.test(x=FamilyPs$NH4, y=FamilyPs$Abundance,method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyPs2 <-cor.test(y=FamilyPs$Abundance, x=FamilyPs$PO4, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyPs3 <-cor.test(y=FamilyPs$Abundance, x=FamilyPs$NOx, method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyPs1
FamilyPs2
FamilyPs3

#Burkholderiaceae
FamilyBur <-FamilySig[49:54,]
FamilyBur1 <-cor.test(x=FamilyBur$NH4, y=FamilyBur$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyBur2 <-cor.test(y=FamilyBur$Abundance, x=FamilyBur$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyBur3 <-cor.test(y=FamilyBur$Abundance, x=FamilyBur$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyBur1
FamilyBur2
FamilyBur3

#Rhodocyclaceae
FamilyRhodo <-FamilySig[55:60,]
FamilyRhodo1 <-cor.test(x=FamilyRhodo$NH4, y=FamilyRhodo$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyRhodo2 <-cor.test(y=FamilyRhodo$Abundance, x=FamilyRhodo$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyRhodo3 <-cor.test(y=FamilyRhodo$Abundance, x=FamilyRhodo$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyRhodo1
FamilyRhodo2
FamilyRhodo3

#Crocinitomicaceae
FamilyCro <-FamilySig[61:66,]
FamilyCro1 <-cor.test(x=FamilyCro$NH4, y=FamilyCro$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyCro2 <-cor.test(y=FamilyCro$Abundance, x=FamilyCro$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyCro3 <-cor.test(y=FamilyCro$Abundance, x=FamilyCro$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyCro1
FamilyCro2
FamilyCro3

#Flavobacteriaceae1
FamilyFlav <-FamilySig[67:72,]
FamilyFlav1 <-cor.test(x=FamilyFlav$NH4, y=FamilyFlav$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyFlav2 <-cor.test(y=FamilyFlav$Abundance, x=FamilyFlav$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyFlav3 <-cor.test(y=FamilyFlav$Abundance, x=FamilyFlav$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyFlav1
FamilyFlav2
FamilyFlav3

#JGI 0000069-P22 order
FamilyJP <-FamilySig[73:78,]
FamilyJP1 <-cor.test(x=FamilyJP$NH4, y=FamilyJP$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyJP2 <-cor.test(y=FamilyJP$Abundance, x=FamilyJP$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyJP3 <-cor.test(y=FamilyJP$Abundance, x=FamilyJP$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyJP1
FamilyJP2
FamilyJP3

#Candidatus Kaiserbacteria order
FamilyCa <-FamilySig[79:84,]
FamilyCa1 <-cor.test(x=FamilyCa$NH4, y=FamilyCa$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyCa2 <-cor.test(y=FamilyCa$Abundance, x=FamilyCa$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyCa3 <-cor.test(y=FamilyCa$Abundance, x=FamilyCa$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyCa1
FamilyCa2
FamilyCa3

#Candidatus Pacebacteria order
FamilyPa <-FamilySig[85:90,]
FamilyPa1 <-cor.test(x=FamilyPa$NH4, y=FamilyPa$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyPa2 <-cor.test(y=FamilyPa$Abundance, x=FamilyPa$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyPa3 <-cor.test(y=FamilyPa$Abundance, x=FamilyPa$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyPa1
FamilyPa2
FamilyPa3

#Flavobacteriaceae2
FamilyFlaAnother <-FamilySig[91:96,]
FamilyFlaAnother1 <-cor.test(x=FamilyPa$NH4, y=FamilyPa$Abundance,method="spearman", alternative="less", exact=FALSE, conf.level=0.95)
FamilyFlaAnother2 <-cor.test(y=FamilyPa$Abundance, x=FamilyPa$PO4, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyFlaAnother3 <-cor.test(y=FamilyPa$Abundance, x=FamilyPa$NOx, method="spearman", alternative="greater", exact=FALSE, conf.level=0.95)
FamilyFlaAnother1
FamilyFlaAnother2
FamilyFlaAnother3



####################################################
#Analysis using absolute abundance
###################################################

#load libraries for this section
library(stats)
library(vegan)
######linear model of regression analysis##
##test the correlation between relative vs. absolute abundances

#load relative abundance data and absolute abundance data
#load relative abundance table
RA <-read.csv("/Users/DandaN/Desktop/RARegnew.csv",row.names=1, header=T)
#transpose the table, so ASVs in column
tRA <-t(RA)

#load absolute abundance table (relative abundance x total cell counts in each sample for each ASV)
AA <-read.csv("/Users/DandaN/Desktop/AARegnew.csv",row.names=1, header=T)
tAA <-t(AA)

set.seed(125)
#create a dataframe to collect the output 
df<-data.frame(r1=rep(NA,100),r2=rep(NA,100),p=rep(NA,100))

#create a loop to run 100 linear regressions independently
for(i in 1:100){
    reg= summary(lm(tRA[,i] ~ tAA[,i])) 
    df[i,'r1']<-reg$r.squared
    df[i,'r2']<-reg$adj.r.squared
    df[i,'p']<-reg$coefficients[2,1]
    rownames(df)[i]<-colnames(tRA)[i]}
    }

#check the first 50 lines to verify the output
head(df)

#change the column name of the table
colnames(df) <-c("R2", "adj.R2", "P-value")
head(df)

########
##community analysis with absolute abundance
########
library(vegan)
library(ggplot2)
library(phyloseq)

#load phyloseq object for absolute abundance
AbsoPhyloObject <-readRDS ("/InputFiles/AbsoPhlyoObject.rds")

#PERMANOVA testing
#the effects of POM mixtures on community variation measured in terms of absolute abundance

#Bray-Curtis distance matrix using absolute abundance
AADistance <-phyloseq::distance (AbsoPhyloObject, method="bray")
sampledf <- data.frame(sample_data(AbsoPhyloObject))

##running adonis test
adonis(formula = AADistance ~ Microcosm, data = sampledf) 

########
##Mantel and PROTEST test
#########
#for correlation between beta-diversity 
#measured from relative abundance absolute abundance

#load phyloseq object with relative abundance
newPhyloObject <-readRDS("/InputFiles/newPhyloObject.rds")

#Bray-Curtis distance matrix using relative abundance
newPhyloObjectRelat <-phyloseq::distance(transform_sample_counts(newPhyloObject, function(x) x/sum(x)*100), method="bray")
sampledf <- data.frame(sample_data(newPhyloObject))

#Mantel test
mantel(newPhyloObjectRelat, AADistance, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))

#Procrustes rotation of two configurations and PROTEST
#Function procrustes rotates a configuration to maximum similarity with another configuration.
protest(newPhyloObjectRelat, AADistance, scale=FALSE, symmetric=T)

 
