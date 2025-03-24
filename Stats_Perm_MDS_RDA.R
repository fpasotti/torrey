# Species and environmental data dune ## for examples stackoveflow
dune2.spe <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/dune2.spe.txt', row.names = 1)
dune2.env <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/dune2.env.txt', row.names = 1)

# Species attributes
dune2.traits <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/dune2.traits.txt', row.names = 1)
dune2.ell <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/dune2.ell.txt', row.names = 1)

# Original dune dataset from vegan:
# If you need to install vegan package, use the following:
# install.packages ('vegan')
library (vegan)
require(dplyr)
data (dune) # matrix with species data (20 samples in rows and 30 species in columns)
data (dune.env)# matix of environmental variables (20 samples in rows and 5 environmental variables in columns)


############ LONG TERM COMPARISON MACROFAUNA ############
#samples from 2011 identified by me and Ann and Uli will be compared to samples Id by Nene --> For the Simper analysis I try to group together organisms that have not been Id to the 
#same level or with the same precision --> too high chance of error type I


##I will now be copying the process from Dune dataset ###It worked (look at how the dune dataset is organised!)

## D= Densities
## B= Biomass
TaxaB <- data.frame(TAXA_MACRO_BIOM)## load multivariate Biomass of Macrofauna from xcl imported dataset
TaxaD <- data.frame(MACRO_DENS_LONG2)##load multivariate Densities of Macrofauna from xcl imported dataset
factor <- data.frame(FACTORS_MACRO_BIOM) ## Load Factors for both (24 observations)
require(vegan)## we need vegan library
require(dplyr)## to reshape the df for SIMPER
factor$STATION <- factor(factor$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
#station_names <-mapvalues(factor$STATION, from = c("FARO", "ISLAD", "CREEK"), to = c("FARO", "ISLA D", "CREEK")) ## add the space in the factor ISLAD
levels(factor$STATION)[levels(factor$STATION)=="ISLAD"] <- "ISLA D"
factor$YEAR <-  factor(factor$YEAR, levels=c("2011","2016"))## define factor YEAR levels
factor$STATIONYEAR  <- with(factor, interaction(factor$STATION,  factor$YEAR), drop = TRUE, levels=c("CREEK/2011","ISLA D/2011","FARO/2011","CREEK/2016","ISLA D/2016","FARO/2016"))## DEfine interactions for PERMDISP later and fpr boxplot

factor$STATION##call the factors to see they have been defined properly
factor$YEAR
factor$STATIONYEAR
 str(factor)

## perform nMDS on the DENSITIES of Macrofauna LONG TERM
metaMDS_L_D_Ma <- metaMDS(TaxaD, distance="bray")## I prefer this type of distance matrix calculation for biological data
factor # dataframe with metadata for the plot with factors STATION YEAR and STATIONYEAR
###TRYING PLOTTING MDS below WITH ggplot2
#build a data frame with NMDS coordinates and metadata
MDS1_L_D_Ma = metaMDS_L_D_Ma$points[,1]
MDS2_L_D_Ma = metaMDS_L_D_Ma$points[,2]
NMDS_Long_D_Ma = data.frame(MDS1 = MDS1_L_D_Ma, MDS2 = MDS2_L_D_Ma, Station = FACTORS_MACRO_BIOM$STATION, YEAR = FACTORS_MACRO_BIOM$YEAR)

str(NMDS_Long_D_Ma)

NMDS_Long_D_Ma$StationYEAR <-with(NMDS_Long_D_Ma, interaction(NMDS_Long_D_Ma$Station,  NMDS_Long_D_Ma$YEAR), drop = TRUE, levels=c("CREEK/2011","ISLA D/2011","FARO/2011","CREEK/2016","ISLA D/2016","FARO/2016"))## DEfine interactions for PERMDISP later and fpr boxplot
NMDS_Long_D_Ma$YEAR <- factor(NMDS_Long_D_Ma$YEAR, levels = c('2011', '2016'))
NMDS_Long_D_Ma$Station <- factor(NMDS_Long_D_Ma$Station, levels = c('FARO', 'ISLAD', 'CREEK'))

ord<-ordiellipse(metaMDS_D, factor$STATIONYEAR, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red'),label = F)##draw="ploygon" allows to color fill the ellipses
ord


#### calculating the values to show the ellipses with a vegan package function
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame() # create the df with these values which have to be calculated for the interaction factor StationYEAR
for(g in levels(NMDS_Long_D_Ma$StationYEAR)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_Long_D_Ma[NMDS_Long_D_Ma$StationYEAR==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell <- data.frame() # create the df with these values which have to be calculated for the interaction factor StationYEAR
for(g in levels(NMDS_Long_D_Ma$StationYEAR)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_Long_D_Ma[NMDS_Long_D_Ma$StationYEAR==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

NMDS.mean=aggregate(NMDS_Long_D_Ma[,1:2],list(group=NMDS_Long_D_Ma$YEAR),mean)

### now ellipses can be added with function geom_path where group=g will determine the shape and hence the points represented by
#the ellipses and its coloration
#geom_path(data=df_ell, aes(x=MDS1, y=MDS2, color=group), size=1, linetype=2)+

### Try plotting with ggplot ###

str(NMDS_Long_D_Ma)
str(df_ell)
  
require(ggplot2)
plotG<- ggplot(NMDS_Long_D_Ma, aes(x=MDS1, y=MDS2), show.legend=NA) +
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2, color=group), size=1, linetype=2, show.legend = NA)+
  scale_color_manual(values = c('darkgreen', 'black', 'red','darkgreen', 'black', 'red')) +
  annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)+
  #scale_fill_manual(values=c('darkgreen', 'black', 'red'))+# to give the color values I wish for 
  geom_point (data= NMDS_Long_D_Ma, shape = shps[NMDS_Long_D_Ma$YEAR]) + # to ask for adding different shapes to the two levels of factor YEAR
  theme_bw() +
  labs(title = "Five years Analysis - Macrofauna densities")

plot(plotG)

str(df_ell)

### it works till here but yet I can only show the Interaction factor with the ellipse since they are calculated on its levels, but I could add
#the YEAR factor level at the center of the YEAR group ellipses (in this case it makes sense since the multivariate densities do show differences between years, otherwise it would become 
#confusing)


#### from Lara ########
#Make Metric MDS | k is the number of dim | du is the distance matrix
#fit<-cmdscale(du,eig=TRUE,k=2) 
#fit$points

#Save coordinates as dataframe
#mds<-as.data.frame(fit$points)
#mds

#Plot MDS 
#p1<-s.class(mds,fac=factor$Area,col=c("#CC66FF","#CC0000","#33FF00","#FFCC00","#0066CC","#33CCFF"),
            #grid=TRUE,cgrid=1)
p1<-s.class(NMDS_Long_D_Ma,fac=NMDS_Long_D_Ma$StationYEAR,col=c('red', 'darkgreen','black','red', 'darkgreen','black'),
            grid=TRUE,cgrid=0.5)
head(NMDS_Long_D_Ma)
str(NMDS_Long_D_Ma)


###### Using only Vegan #### perform nMDS on the DENSITIES of Macrofauna LONG TERM ####
metaMDS_L_D_Ma <- metaMDS(TaxaD, distance="bray")## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# empty plot
plot(metaMDS_D, type = 'n',main='Five years Analysis - 15 m sampling', subtitle='Macrofauna densities')##Non-metric fit, R^2=0.97 / Linear fit, R^2=0.83
# add points
points(metaMDS_D, col = cols[factor$STATION], pch = shps[factor$YEAR])
# add legend
legend('topleft', col=cols, legend=levels(factor$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor$YEAR), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_D, factor$STATIONYEAR, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

stressplot(metaMDS_D) ##YES!!!

TaxaD_2 <- TaxaD^0.25
##perform PERMANOVA with interaction factor on DENSITIES MACROFAUNA LONG TERM
adonis(TaxaD_2~factor$STATION*factor$YEAR, permutation = 999, method='bray')## the PERMANOVA results 
## give me a significant interaction factor STATION x YEAR p=***

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_TaxaD=vegdist(TaxaD_2, method="bray")## Bray curtis similarity matrix on fourth root transformed data
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_D=(betadisper(d_TaxaD,factor$STATIONYEAR, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_D)## for factor STATIONYEAR p=*** hence there are differences in the dispersions 
## and we must be carefull in considering the factor effect coz there is also a dispersion effect
## as different coz the differences may be due to the differences in dispersions

par(cex.axis=0.8)
boxplot(bd_D, cex=0.8,main="Macrofauna Multivar Dens, PermDisp STxYEAR p=***")
permutest(bd_D)






######### SIMPER Analysis ######### --> The SIMilarity PERcentages breakdown (SIMPER) procedure attempts to assess the average 
##percent contribution of individual variables (e.g. species) to the dissimilarity between stations
##in a Bray-Curtis dissimilarity matrix. This allows users to identify variables (e.g. species) that are 
##likely to be the major contributors to any difference between groups (factors) detected by methods such as 
##PERMANOVA. 
## I want to modify my dataset to group to higher taxon level to make the SIMPER more meaningful and less biased towards rare taxa but more looking at general taxon (functional) trends

require(dplyr)

TaxaDSim <- TaxaD %>% mutate(Amphipoda= Lyanissidae + Phoxocephalidae + Eusiridae + Oedicerotidae + Type_1 + Type_8 +Corophiodae)#I sum all the amphipods in one
TaxaDSim <-  select (TaxaDSim,-c(Lyanissidae, Phoxocephalidae, Eusiridae, Oedicerotidae, Type_1, Type_8, Corophiodae))#I get rid of the old columns
TaxaDSim <- TaxaDSim %>% mutate(Other_Bivalvia= Nucula_sp+ Thyasiridae)
TaxaDSim <-  select (TaxaDSim,-c(Nucula_sp, Thyasiridae))
TaxaDSim <- TaxaDSim %>% mutate(Other_Polychaeta= Polychaeta.sp1+ Spionidae + Ophelina + Travisia + Terbellidae+Orbinidae )
TaxaDSim <-  select (TaxaDSim,-c(Polychaeta.sp1, Spionidae , Ophelina , Travisia , Terbellidae,Orbinidae))
TaxaDSim <- TaxaDSim %>% mutate(Cumacea= Leuconidae+Bodotriidae+Nannastacidae)
TaxaDSim <-  select (TaxaDSim,-c(Leuconidae,Bodotriidae,Nannastacidae))
TaxaDSim <- TaxaDSim %>% mutate(Isopoda=Serolidae+Munnidae)
TaxaDSim <-  select (TaxaDSim,-c(Serolidae,Munnidae))
TaxaDSim_2 <- TaxaDSim^0.25

### SIMPER analysis Maxcrofauna Densities Long Term
simD= simper(TaxaDSim_2, group = factor$STATIONYEAR, permutations= 9999)## run the simper() on the transformed dataset and not on a previously made dist matrix cause simper() does use Bray curtis similarity matrix itself itself
simD
lapply(simD,FUN=function(TaxaDSim_2){TaxaDSim_2$overall})## to get the overall dissimilarity between groups

SUM <- summary(simD, digits=3, ordered=TRUE)
str(SUM)
for (i in seq_along(SUM)) {
  filename = paste(i, ".csv")
  write.csv(SUM[[i]], filename)
} ## to save into distinct csv files each of the SIMPER df produced by the Summary function of SIMPER


#### BIOMASS nMDS LONG TERM

TaxaB <- data.frame(TAXA_MACRO_BIOM)##load table from xls dataset ## factor was loaded earlier

factor$STATION
factor$YEAR
factor$STATIONYEAR


## nMDS
metaMDS_B <- metaMDS(TaxaB, distance="bray")## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# empty plot
plot(metaMDS_B, type = 'n',main='Macrofauna Biomass Five years Analysis')
# add points
points(metaMDS_B, col = cols[factor$STATION], pch = shps[factor$YEAR])
# add legend
legend('topleft', col=cols, legend=levels(factor$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor$YEAR), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_B, factor$STATIONYEAR, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

### alternative drawing
ordiellipse(metaMDS_B, factor$STATIONYEAR, display = "sites", kind = "se", conf=0.95,draw="line", col=c('green', 'black', 'red'),label = F)##draw="line" allows to only have the lines that delimit  the ellipses
stressplot(metaMDS_B) ##YES!!!

### Now some PERMANOVA MACROFAUNA BIOMASS LONG TERM
##adonis() partitions the dissimilarities between sites on the basis of
##certain factors - so it is like ANOVA but with multivariate responses
##using any dissimilarity and you can think of this analysis as looking at
##whether the distances between points within a group (factor, site etc)
##are greater than the between groups distances or distances to other
##groups.

##If adonis() shows a significant difference between sites, this
##difference may be due to differences in the spread or variance of the
##respective groups (measured by betadisper on the centroid of each group and then to calculate the squared deviations from these points ), just like in standard ANOVA. 
##to test for this.

##let's TRANSFORM the data fourth root
TaxaB_2 <- TaxaB^0.25
##perform PERMANOVA with interaction factor
adonis(TaxaB_2~factor$STATION*factor$YEAR, permutation = 9999, method='bray')## the PERMANOVA results give me a significant interaction factor STATION x YEAR


##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_TaxaB=vegdist(TaxaB_2, method="bray")
class(d_TaxaB)
as.matrix(d_TaxaB)
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_B=(betadisper(d_TaxaB,factor$STATION, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_B)## for factor STATION p=NS hence there are no differences in the dispersion of the variances hence interpretation of STATION effect
#is safe 

##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd2_B=(betadisper(d_TaxaB,factor$STATIONYEAR, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd2_B)## for factor STATIONYEAR p=*  differences at each STATION between YEARS can be due to differences in the variances dispersion
##as different coz the differences should not be due to differences in dispersions in the interaction factor
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.8)
boxplot(bd2_B, cex=0.8,main="MACRO Multivar Biom FIVE YEARS TREND, PermDisp p=*")

TukeyHSD.betadisper(TaxaB_2~factor$STATION*factor$YEAR, method='bonf') #### Can't find function!


##let's do a pairwise comparison to see what groups differ significantly (we have a=groups=2 and n=4)
##number of unique permutations=(2*4)!/(2!(4!)^2))=35 --> minimum significance level 0.028 --> it is ok
library(vegan)
library(lattice)
library(permute)
library(RVAideMemoire)
library(ecodist)
#d_TaxaB=distance(TaxaB_2,method="bray-curtis")

pairwise.perm.manova(d_TaxaB, factor$STATIONYEAR, test = "Wilks", nperm = 9999, 
                     progress = TRUE, p.method = "bonferroni") ##it cannot deal with low number of unique permutations!!!need a p(MC) which can't exist in R!


###about WHY NOT doing 2-way PERMDISP types:
##A while ago, I had asked a similar question to Marti Anderson. Here's is
#> > what she had replied:
 # > > 
  #> >         "I have actually moved away from doing 2-way PERMDISP types of
#> >         analyses. The reason for this is that there are inherent
#> >         difficulties in testing homogeneity of dispersions across the
#> >         main effects in cases where there is an interaction between the
#> >         two factors in their *locations*, when tested using PERMANOVA.
#> >         Think of it this way - suppose there is no differences in the
#> >         locations of samples in multivariate space among the levels of
#> >         factor A and in level 1 of factor B, but there are such
#> >         differences (once again, in *location*) in level 2 of factor B
#> >         (which would correspond to an interaction in PERMANOVA). Then,
#> >         clearly, this would be detected as a difference in dispersion,
#> >         even though it is entirely due to the locations of cells
#> >         (combinations of factors) in the two-way design. Thus, it seems
#> >         to me that for the two-way design, to test dispersions, one can:
#> >         
#> >         (i) test for differences in dispersion among the individual axb
#> >         cells (which is simply a one-way test where the individual cells
#> >         are treated as levels of a single factor) and
#> >         
#> >         (ii) test for an interaction in a two-way PERMANOVA and if this
#> >         is NS, then one can examine the test of dispersions (or
#> >         locations, for that matter) separately and independently for the
#> >         two factors."
#> 
 # > Thanks for posting this Etienne. I purposely wrote betadisper() to
#> accept a single grouping variable because I was struggling to see how
#> you would do the analysis with multiple factors and so handled the
#> simple, easy case first. It has been on my TODO list to revist this, but
#> Marti's comments are good enough for me to tick this one off the list,
#> and perhaps add a note to the help page for betadisper() summarising the
#> above points.
#> 
## Cheers,
##

## Simper Analysis MACRODFAUNA BIOMASS LONG TERM on the dissimilarity matrix we built before (Bray Curtis)
TaxaBSim <- TaxaB %>% mutate(Other_Bivalvia= Nucula_sp+ Thyasiridae)## I am making the dataset more meaningful by grouping at higher taxonomic level some rare taxa
TaxaBSim <-  select (TaxaBSim,-c(Nucula_sp, Thyasiridae))
TaxaBSim <- TaxaBSim %>% mutate(Other_Polychaeta= Polychaeta.sp_1+ Spionidae + Ophelina + Travisia + Terebellidae+Orbinidae )
TaxaBSim <-  select (TaxaBSim,-c(Polychaeta.sp_1, Spionidae , Ophelina , Travisia , Terebellidae,Orbinidae))
TaxaBSim <- TaxaBSim %>% mutate(Cumacea= Leuconidae+Bodotriidae+Nannastacidae)
TaxaBSim <-  select (TaxaBSim,-c(Leuconidae,Bodotriidae,Nannastacidae))

TaxaBSim_2 <- TaxaBSim^0.25


## simper
simB_LONG =simper(TaxaBSim_2, group = factor$STATIONYEAR, permutations= 9999)
SUM2 <- summary(simB_LONG)
lapply(simB_LONG,FUN=function(TaxaBSim_2){TaxaBSim_2$overall})## to get the overal dissimilarity between groups

str(SUM2)
for (i in seq_along(SUM2)) {
  filename = paste(i, ".csv")
  write.csv(SUM2[[i]], filename)
} ## to save into distinct csv files each of the SIMPER df produced by the Summary function of SIMPER


## nMDS on the DENSITIES of MACROFAUNA LONG TERM
TaxaD <- data.frame(MACRO_DENS_LONG)##load table from xls dataset
factor <- data.frame(FACTORS_MACRO_BIOM)
require(vegan)
factor$STATION <- factor(factor$STATION, levels=c("FARO","ISLAD","CREEK"))
levels(factor$STATION)[levels(factor$STATION)=="ISLAD"] <- "ISLA D"
factor$YEAR <-  factor(factor$YEAR, levels=c("2011","2016"))
factor$STATIONYEAR  <- with(factor, interaction(factor$STATION,  factor$YEAR), drop = TRUE)

TaxaD
factor
factor$STATION
factor$YEAR
factor$STATIONYEAR

##Diversity Indices
H_Macro_Dens <- diversity(TaxaD)
J_Macro_Dens <- H_Macro_Dens/log(specnumber(TaxaD))
H_Macro_Dens <- as.data.frame(H_Macro_Dens)
J <- as.data.frame(J)
k <- sample(nrow(TaxaD), 6)
K <- as.data.frame(k)
R <- renyi(TaxaD [k,])##We can really regard a site more diverse if all of its Renyi diversities are higher than 
##in another site.
plot(R)

Plot_R <- as.data.frame(c('FARO/2011 (23)', 'ISLAD/2011 (9)', 'CREEK/2011 (24)', 'FARO/2015 (22)', 'ISLAD/2015 (12)', 'CREEK/2015 (20)'))

 
library(vegan)

metaMDS_D <- metaMDS(TaxaD, distance="bray")## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# yoy first make an empty plot
plot(metaMDS_D, type = 'n',main='Macrofauna Densities 5 YEARS TREND')
# add points
points(metaMDS_D, col = cols[factor$STATION], pch = shps[factor$YEAR])
# add legend
legend('topleft', col=cols, legend=levels(factor$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor$YEAR), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_D, factor$STATIONYEAR, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses
##the stations are more homogeneous in 2016 compared to 2011 where they differ more in terms of community structure

stressplot(metaMDS_B) ##YES!!!

results4 <- envfit(metaMDS_D,ENV_Long,permutations=999, strata=NULL,choices=c(1,2),display="sites")
plot(results4)



### Now some PERMANOVA on the DENSITIES of macrofauna LONG TERM

##let's transform the data fourth root
TaxaD_2 <- TaxaD^0.25
##perform PERMANOVA with interaction factor
adonis(TaxaD_2~factor$STATION*factor$YEAR, permutation = 999, method='bray')


##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), to understand whether
##it is the DIESPERSION of your group from the CENTROID that determines the difference of the 
##CENTROID themselves 
## we make the distance matrix ourselves
d_TaxaD=vegdist(TaxaD, method="bray")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_D=(betadisper(d_TaxaD,factor$STATION, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_D)## for factor STATION p=NS hence there are no  differences in the dispersions of variances
par(cex.axis=0.8)
boxplot(bd_D, cex=0.8)

bd1_D=(betadisper(d_TaxaD,factor$YEAR, type = "centroid"))
permutest(bd_D)## for factor YEAR p=NS hence there are  differences in the dispersions of vairances
par(cex.axis=0.8)
boxplot(bd1_D, cex=0.8)
##Subsequently, we can conduct a PERMDISP on this distance matrix testing for differences in dispersion among the individual STATIONxYEAR
#cells (which is simply a one-way test where the individual
#>  samples are treated as levels of a single factor (STATIONYEAR)
bd2_D=(betadisper(d_TaxaD,factor$STATIONYEAR, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd2_D)## for factor STATIONYEAR p=* hence there are differences in the dispersion of variances  at each STATION between YEARS and we should be careful considering the interaction term

par(cex.axis=0.8)
boxplot(bd2_D, cex=0.8,main="MACRO Multivar Dens 5 YEARS TREND, PermDisp STxYEAR p=*")


##let's do a pairwise comparison to see what groups differ significantly (we have a=groups=2 and n=4)
##number of unique permutations=(2*4)!/(2!(4!)^2))=35 --> minimum significance level 0.028 --> it is ok
library(vegan)
library(RVAideMemoire)
pairwise.perm.manova(d_TaxaD, factor$STATIONYEAR, test = "Wilks", nperm = 999, 
                     progress = TRUE, p.method = "bonferroni")###only significant STATION effect Isla D differs from the others but 
##no effect of STATIONYEAR or YEAR --> very high variances in 2011 compared to 2016


####################### SEASONAL ANALYSIS MACROFAUNA ###############


Taxa_SB <- data.frame(SEASON_MACRO_BIOM_updated)## load multivariate Biomass of Macrofauna from xcl imported dataset
Taxa_SD <- data.frame(MACRO_Season_Dens_Yoldia1Col)##load multivariate Densities of Macrofauna from xcl imported dataset
factor_SD <- data.frame(SEASON_MACRO_DENS_factors)#load factors for Densities Seasonal
factor_SB <- data.frame(SEASON_MACRO_BIOM_factors)

### to make all columns as.numeric in case they erroneously are read as text
ix <- 1:36 # tot columns number is 36
MACRO_Season_Dens_Yoldia1Col[ix] <- lapply(MACRO_Season_Dens_Yoldia1Col[ix], as.numeric) 
str(MACRO_Season_Dens_Yoldia1Col)



factor_SB$STATION <- factor(factor_SB$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(factor_SB$STATION)[levels(factor_SB$STATION)=="ISLAD"] <- "ISLA D"
levels(factor_SB$DEPTH)[levels(factor_SB$DEPTH)=="9"] <- "9m"
levels(factor_SB$DEPTH)[levels(factor_SB$DEPTH)=="15"] <- "15m"
levels(factor_SD$DEPTH)[levels(factor_SD$DEPTH)=="9"] <- "9m"
levels(factor_SD$DEPTH)[levels(factor_SD$DEPTH)=="15"] <- "15m"

factor_SB$SEASON <-  factor(factor_SB$SEASON, levels=c("Summer15","Spring15", "Spring16"))## define factor SEASON levels
factor_SB$STATIONSEASON  <- with(factor_SB, interaction(factor_SB$STATION,  factor_SB$SEASON), drop = TRUE )## DEfine interactions for PERMDISP later and fpr boxplot

factor_SD$STATION <- factor(factor_SD$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(factor_SD$STATION)[levels(factor_SD$STATION)=="ISLAD"] <- "ISLA D"
factor_SD$SEASON <-  factor(factor_SD$SEASON, levels=c("Summer15","Spring15", "Spring16"))## define factor SEASON levels
factor_SD$STATIONSEASON  <- with(factor_SD, interaction(factor_SD$STATION,  factor_SD$SEASON), drop = TRUE)## DEfine interactions for PERMDISP later and fpr boxplot


##call the factors and check
factor_SD$STATION
factor_SB$STATION
factor_SB$SEASON
factor_SD$SEASON
factor_SD$STATIONSEASON
factor_SB$STATIONSEASON
factor_SD


## nMDS on the BIOMASS of Macrofauna SEASONAL
metaMDS_SB <- metaMDS(Taxa_SB, distance="bray")## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 19)
# empty plot
plot(metaMDS_SB, type = 'n',main='Macrofauna Biomass SEASONAL analysis')### Non-metric fit R^2=0.96, Linear fit, R^2=0.829
# add points
points(metaMDS_SB, col = cols[factor_SB$STATION], pch = shps[factor_SB$SEASON])
# add legend
legend('topleft', col=cols, legend=levels(factor_SB$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor_SB$SEASON), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_SB, factor_SB$STATIONSEASON, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses


stressplot(metaMDS_SB) ##YES!!!


##PERMANOVA##
#Now some PERMANOVA on the biomass of macrofauna##adonis relies on a long-understood phenomenon that allows one to partition 
## sums of squared deviations from a centroid in two different ways, he most widely recognized method, used, e.g., for ANOVA and MANOVA,
##is to first identify the relevant centroids and then to calculate the squared deviations from these points

##let's transform the data fourth root
Taxa_SB_2 <- Taxa_SB^0.25
##perform PERMANOVA with interaction factor
adonis(Taxa_SB_2~factor_SB$STATION*factor_SB$SEASON, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*



##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Taxa_SB_Yoldia=vegdist(Taxa_SB_Yoldia, method="bray")
bd__SB_Yoldia=(betadisper(d_Taxa_SB_Yoldia,factor_SB$STATION, type = "centroid"))
permutest(bd__SB_Yoldia)### Permdisp STATION p=*,  SEASON=NS, STATIONSEASON p=*
simSB_Yoldia =simper(Taxa_SB_Yoldia, group = factor_SB$STATIONSEASON, permutations= 9999)
simSB_Yoldia
summary(simSB_Yoldia)


##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Taxa_SB=vegdist(Taxa_SB_2, method="bray")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd__SB=(betadisper(d_Taxa_SB,factor_SB$SEASON, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd__SB)## for factor STATION p=*** hence there is a factor effect but ALSO a dispersion effect 


par(cex.axis=0.6)
boxplot(bd__SB, cex=0.6,main="MACRO Multivar Biom SEASONAL analysis, PermDisp p=***", names = c("FARO / SUM15","ISLA D / SUM15","CREEK / SUM15","FARO / SPR15","ISLA D / SPR15","CREEK / SPR15","FARO / SPR16","ISLA D / SPR16","CREEK / SPR16"),las=2)

###SIMPER analysis
simSB= simper(Taxa_SB_2, group = factor_SB$STATIONSEASON, permutations= 9999)## run the simper() on the transformed dataset and not on a previously made dist matrix cause simper() does use Bray curtis itself

SUM3 <- summary(simSB, digits=3, ordered=TRUE)
str(SUM3)
for (i in seq_along(SUM3)) {
  filename = paste(i, ".csv")
  write.csv(SUM3[[i]], filename)
} ## to save into distinct csv files each of the SIMPER df produced by the Summary function of SIMPER


## nMDS on the DENSITIES of Macrofauna Seasonal

metaMDS_SD <- metaMDS(Taxa_SD, distance="bray")## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 19)
# empty plot
plot(metaMDS_SD, type = 'n',main='Macrofauna Densities SEASONAL analysis')### Non-metric fit R^2=0.96, Linear fit, R^2=0.829
# add points
points(metaMDS_SD, col = cols[factor_SD$STATION], pch = shps[factor_SD$SEASON])
# add legend
legend('topleft', col=cols, legend=levels(factor_SD$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor_SD$SEASON), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_SD, factor_SD$STATIONSEASON, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses


stressplot(metaMDS_SD) ##YES!!!



### PERMANOVA on the biomass of macrofauna##adonis relies on a long-understood phenomenon that allows one to partition 
## sums of squared deviations from a centroid in two different ways, he most widely recognized method, used, e.g., for ANOVA and MANOVA,
##is to first identify the relevant centroids and then to calculate the squared deviations from these points

##let's transform the data fourth root
Taxa_SD_2 <- Taxa_SD^0.25
##perform PERMANOVA with interaction factor
adonis(Taxa_SD_2~factor_SD$STATION*factor_SD$SEASON, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=**

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Taxa_SD=vegdist(Taxa_SD_2, method="bray")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
# for factor STATION p=NS  hence  OK! there are no differences in the dispersions we can trust the PERMANOVA 

##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd2__SD=(betadisper(d_Taxa_SD,factor_SD$STATIONSEASON, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd2__SD)## for factor STATIONYEAR p=NS hence the centroid dispersion is not driving the significance we find in this factor
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd2__SD, cex=0.6,main="MACRO Multivar Dens SEASONAL analysis, PermDisp p=NS", names = c("FARO / SUM15","ISLA D / SUM15","CREEK / SUM15","FARO / SPR15","ISLA D / SPR15","CREEK / SPR15","FARO / SPR16","ISLA D / SPR16","CREEK / SPR16"),las=2)

###SIMPER analysis
Taxa_SD_2 <- Taxa_SD^0.25
simSD= simper(Taxa_SD_2, group = factor_SB$STATIONSEASON, permutations= 9999)## run the simper() on the transformed dataset and not on a previously made dist matrix cause simper() does use Bray curtis itself
#I used the factor of the Biomass df cause I did not use the grouped df I made by hand but the complete list
SUM4 <- summary(simSD, digits=3, ordered=TRUE)
str(SUM4)
for (i in seq_along(SUM3)) {
  filename = paste(i, ".csv")
  write.csv(SUM3[[i]], filename)
} ## to save into distinct csv files each of the SIMPER df produced by the Summary function of SIMPER








#########----------- YOLDIA case ######---------
### I want to make a subset of my dataset where I only look at Yoldia size classes biomass 
library(dplyr)
Taxa_SB_Yoldia <- subset(Taxa_SB_2[,1:4])
### I want to make an analysis of the Yoldia Size classes trends for biimass SEASONAL analysis
adonis(Taxa_SB_Yoldia~factor_SB$STATION*factor_SB$SEASON, permutation = 999, method='bray')##
##Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#factor_SB$STATION                   2   1.20385 0.60192 28.9267 0.49283  0.001 ***
#factor_SB$SEASON                    2   0.32608 0.16304  7.8353 0.13349  0.001 ***
#factor_SB$STATION:factor_SB$SEASON  4   0.22613 0.05653  2.7168 0.09257  0.021 *  
#Residuals                          33   0.68668 0.02081         0.28111           
#Total                              41   2.44274                 1.00000           
#---Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
pairwise.perm.manova(Taxa_SB_Yoldia, factor_SB$STATION, test = "Wilks", nperm = 999, 
                     progress = TRUE, p.method = "bonferroni")
### STATIONSEASON -> from the biomass of Yoldia 4 classes only Faro shows significant differences between Spring15 and Spring16. No other
##differences are found.
### SEASON -> As an overall there is no seasonal trend p> 0.05
## STATION -> all stations differ from one another p=**


str(Taxa_SB_Yoldia)
Taxa_SD# Densities now
Taxa_SD_Yoldia <- subset(Taxa_SD[,1:4])
Taxa_SD_Yoldia <- 
str(Taxa_SD_Yoldia)
str(Taxa_SD_Yoldia$Yoldia_big)
str(factor_SD)
adonis(Taxa_SD_Yoldia~factor_SB$STATION*factor_SB$SEASON, permutation = 999, method='bray')##
##Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#factor_SB$STATION                   2    2.2963 1.14816 12.4288 0.34094  0.001 ***
#factor_SB$SEASON                    2    0.8939 0.44697  4.8385 0.13273  0.001 ***
#factor_SB$STATION:factor_SB$SEASON  4    0.4965 0.12411  1.3435 0.07371  0.206    
#Residuals                          33    3.0485 0.09238         0.45262           
#Total                              41    6.7352                 1.00000           
#---Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.perm.manova(Taxa_SD_Yoldia, factor_SB$SEASON, test = "Wilks", nperm = 999, 
                     progress = TRUE, p.method = "bonferroni")
## factor STATION -> they all differ from eachother p=** / factor SEASON -> spring16 differs from both Spring15 p=** than Summer15p=*

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Taxa_SD_Yoldia=vegdist(Taxa_SD_Yoldia, method="bray")
bd__SD_Yoldia=(betadisper(d_Taxa_SD_Yoldia,factor_SB$STATIONSEASON, type = "centroid"))
permutest(bd__SD_Yoldia)### Permdisp STATION p=*,  SEASON=NS
simSB_Yoldia =simper(Taxa_SB_Yoldia, group = factor_SB$STATIONSEASON, permutations= 9999)
simSB_Yoldia
summary(simSB_Yoldia)

par(cex.axis=0.5)
boxplot(bd__SD_Yoldia, cex=0.4,main="Yoldia SC Densities SEASONAL analysis, PermDisp ST p=*/ SEA p=NS", las=2)


## Yoldia_big##

adonis(Taxa_SD_Yoldia$Yoldia_big~factor_SB$STATION*factor_SB$SEASON, permutation = 999, method='euclidean')##
##Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#factor_SB$STATION                   2   1428383  714192 11.4253 0.35036  0.001 ***
#factor_SB$SEASON                    2    430704  215352  3.4451 0.10565  0.048 *  
#factor_SB$STATION:factor_SB$SEASON  4    154946   38736  0.6197 0.03801  0.659    
#Residuals                          33   2062816   62510         0.50598           
#Total                              41   4076849                 1.00000           
#---Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Densities univariate ###
results= aov(Taxa_SD_Yoldia$Yoldia_big ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SD_Yoldia) ## STATIO p=***, SEASON p=**, STATIONSEASON p=*
Anova(results, type = "III")
TukeyHSD(results)## Isla D different from both Faro and Creek (p=***) - Spring16 different from summer 15 p=*

results= aov(Taxa_SD_Yoldia$Yoldia_med ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SD_Yoldia)## STATION p=*
Anova(results, type = "III")
TukeyHSD(results)## Faro differs from Creek (p=***) and Isla D (p=*)- Spring16 differs from Spring 15 p=**

results= aov(Taxa_SD_Yoldia$Yoldia_small ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SD_Yoldia)## STATION p=* but borderline
Anova(results, type = "III")
TukeyHSD(results)## Faro differs from Creek (p=*) and Isla D (p=*)/ Spring16 differs from Spring 15 p=*** / Creek spring 15 differs Creek spring16

results= aov(Taxa_SD_Yoldia$Yoldia_ss ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SD_Yoldia)## NS 
Anova(results, type = "III")

### Biomass ###
results= aov(Taxa_SB_Yoldia$Yoldia_big ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SB_Yoldia)## STATION p=***, SEASON p=**, STATIONxSEASON p=* but borderline
Anova(results, type = "III")
TukeyHSD(results)## IslaD differs from both Faro and Creek p=*** / Spring 16 differs from Summer 15 p=*** and Spring16 p=* but Stations in Spring 15 vs Spring16 don't show significant differences

results= aov(Taxa_SB_Yoldia$Yoldia_small ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SB_Yoldia)## STATION p=***
Anova(results, type = "III")
TukeyHSD(results)## Faro differs from both Isla D and Creek

results= aov(Taxa_SB_Yoldia$Yoldia_med ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SB_Yoldia)## STATION p=**, SEASON p=*
Anova(results, type = "III")
TukeyHSD(results)## Faro differs from Creek and Isla D / Spring 15 differs from Spring 16

results= aov(Taxa_SB_Yoldia$Yoldia_ss ~ factor_SB$STATION*factor_SB$SEASON, data=Taxa_SB_Yoldia)## STATION p=***
Anova(results, type = "III")
TukeyHSD(results)## Faro differs from Creek and Isla D 

class(Taxa_SB_Yoldia)
Taxa_SB_Yoldia <- as.data.frame(Taxa_SB_Yoldia)
class(Taxa_SB_Yoldia)

d_Taxa_SD_Yoldia_big=vegdist(Taxa_SD_Yoldia$Yoldia_big, method="euclidean")
adonis(d_Taxa_SD_Yoldia_big~ factor_SB$SEASON*factor_SB$STATION, permutations=999)
bd__SD_Yoldia_big=(betadisper(d_Taxa_SD_Yoldia_big,factor_SB$SEASON, type = "centroid"))
permutest(bd__SD_Yoldia_big)### Permdisp STATION p=**,  SEASON=**, STATIONSEASON p=***
boxplot(bd__SD_Yoldia_big, cex = 0.4, main="Dens Yoldia big (>2 cm)", las=2)


 ## Yoldia-Dens##
d_Taxa_SD_Yoldia_med=vegdist(Taxa_SD_Yoldia$Yoldia_med, method="euclidean")
adonis(d_Taxa_SD_Yoldia_med~ factor_SB$SEASON*factor_SB$STATION, permutations=999)
bd__SD_Yoldia_med=(betadisper(d_Taxa_SD_Yoldia_med,factor_SB$STATION, type = "centroid"))
permutest(bd__SD_Yoldia_med)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SD_Yoldia_med, cex = 0.4, main="Dens Yoldia medium (1.5-2 cm)", las=2)

d_Taxa_SD_Yoldia_sm=vegdist(Taxa_SD_Yoldia$Yoldia_small, method="euclidean")
adonis(d_Taxa_SD_Yoldia_sm~ factor_SB$SEASON*factor_SB$STATION, permutations=999)
bd__SD_Yoldia_sm=(betadisper(d_Taxa_SD_Yoldia_sm,factor_SB$STATION, type = "centroid"))
permutest(bd__SD_Yoldia_sm)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SD_Yoldia_sm, cex = 0.4, main="Dens Yoldia small (1-1.5 cm)", las=2)

d_Taxa_SB_Yoldia_ss=vegdist(Taxa_SB_Yoldia$Yoldia_ss, method="euclidean")
adonis(d_Taxa_SD_Yoldia_ss~ factor_SB$SEASON*factor_SB$STATION, permutations=999)
bd__SD_Yoldia_ss=(betadisper(d_Taxa_SD_Yoldia_ss,factor_SB$STATION, type = "centroid"))
permutest(bd__SD_Yoldia_ss)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SD_Yoldia_ss, cex = 0.4, main="Dens Yoldia recruits (<1 cm)", las=2)

### Yoldia Biom ##

d_Taxa_SB_Yoldia_big=vegdist(Taxa_SB_Yoldia$Yoldia_big, method="euclidean")
bd__SB_Yoldia_big=(betadisper(d_Taxa_SB_Yoldia_big,factor_SB$SEASON, type = "centroid"))
permutest(bd__SB_Yoldia_big)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SB_Yoldia_big, cex = 0.4, main="Biom Yoldia Big (>2 cm)", las=2)

d_Taxa_SB_Yoldia_med=vegdist(Taxa_SB_Yoldia$Yoldia_med, method="euclidean")
bd__SB_Yoldia_med=(betadisper(d_Taxa_SB_Yoldia_med,factor_SB$STATION, type = "centroid"))
permutest(bd__SB_Yoldia_med)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SB_Yoldia_med, cex = 0.4, main="Biom Yoldia medium (1.5-2 cm)", las=2)

d_Taxa_SB_Yoldia_sm=vegdist(Taxa_SB_Yoldia$Yoldia_small, method="euclidean")
bd__SB_Yoldia_sm=(betadisper(d_Taxa_SB_Yoldia_sm,factor_SB$STATION, type = "centroid"))
permutest(bd__SB_Yoldia_sm)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SB_Yoldia_sm, cex = 0.4, main="Biom Yoldia small (1-1.5 cm)", las=2)

d_Taxa_SB_Yoldia_ss=vegdist(Taxa_SB_Yoldia$Yoldia_ss, method="euclidean")
bd__SB_Yoldia_ss=(betadisper(d_Taxa_SB_Yoldia_ss,factor_SB$STATIONSEASON, type = "centroid"))
permutest(bd__SB_Yoldia_ss)### Permdisp STATION p=*,  SEASON=*, STATIONSEASON p=***
boxplot(bd__SB_Yoldia_ss, cex = 0.4, main="Biom Yoldia recruits (<1 cm)", las=2)

## Aequiyoldia ### looking at the species density more than its size classes

Aequi_Dens <- data.frame(MACRO_SEASON_DENS)
d_Aequi_SD <- vegdist(Aequi_Dens, method="euclidean")
phoc <- with(Aequi_Dens, betadisper(d_Aequi_SD,factor_SD$STATIONSEASON))
TukeyHSD(phoc) # Isla D Summer 15 vs Creek Summer 15 p=*/
adonis(d_Aequi_SD~ factor_SD$SEASON*factor_SD$STATION, permutations=999)
results= aov(Aequi_Dens$Aequiyoldia.eightsi ~ factor_SD$STATION*factor_SD$SEASON, data=Aequi_Dens)## NS
Anova(results, type = "III")
Aequi_Dens

###---------------------------

#####------ Laternula elliptica SEASON -----######

LAT <- data.frame(Laternula_SEASON)
LAT_factors <- data.frame(Laternula_SEASON_factors)
LAT_factors$STATION <- factor(LAT_factors$STATION, levels=c("FARO","ISLAD","CREEK"))
levels(LAT_factors$STATION)[levels(LAT_factors$STATION)=="ISLAD"] <- "ISLA D"
LAT_factors$SEASON <-  factor(LAT_factors$SEASON, levels=c("Spring15","Summer15","Spring16"))## define factor SEASON levels
LAT_factors$interaction  <- with(LAT_factors, interaction(LAT_factors$STATION,  LAT_factors$SEASON), drop = TRUE )#Define interaction factor

## I prepare to make a graph of the Densities of Laternula

library(ggplot2)       ##package with nice and easy graphics grammar
library(plyr)  ##package necessary to be able to use ddply - see further
library(car) ## for the stats
library(reshape2)

LAT_DensFactors <- data.frame(LAT_DensFactors)
LAT_DensFactors$STATION <- factor(LAT_DensFactors$STATION, levels=c("FARO","ISLAD","CREEK"))
LAT_DensFactors$SEASON <- factor(LAT_DensFactors$SEASON, levels=c("Summer15", "Spring15","Spring16"))
LAT.DENS.graph <- ddply(LAT_DensFactors, .(STATION, SEASON), summarise,                      # make summary with following info:    
                   N      = length(Densities),                            ##number of observations
                   mean = mean(Densities),                                 ##average
                   sd     = sd(Densities),                                 ##standard deviation
                   se     = sd(Densities) / sqrt(length(Densities)-1) )
LAT.DENS.graph$STATION <- factor(LAT.DENS.graph$STATION, levels=c("FARO","ISLAD","CREEK"))
station.names.D <- mapvalues(LAT.DENS.graph$STATION, from=c("FARO","ISLAD","CREEK"), to=c("Faro", "Isla D", "Creek"))
LAT.DENS.graph$SEASON <- factor (LAT.DENS.graph$SEASON, levels=c("Summer15", "Spring15","Spring16"))

Aequi <- data.frame (Yoldia_Summary)
library(dplyr)
Aequi$STATION <- factor(Aequi$STATION, levels=c("FARO", "ISLAD", "CREEK"))
Aequi$SEASON <- factor(Aequi$SEASON, levels=c("Summer15","Spring15","Spring16"))
Aequi.Dens <- ddply(Aequi, .(STATION, SEASON), summarise,
                    N      = length(A.eightsii_Dens),                            ##number of observations
                    mean = mean(A.eightsii_Dens),                                 ##average
                    sd     = sd(A.eightsii_Dens),                                 ##standard deviation
                    se     = sd(A.eightsii_Dens) / sqrt(length(A.eightsii_Dens)-1) )
Aequi.Biom <- ddply(Aequi, .(STATION, SEASON), summarise,
                    N      = length(A.eightsii_Biom),                            ##number of observations
                    mean = mean(A.eightsii_Biom),                                 ##average
                    sd     = sd(A.eightsii_Biom),                                 ##standard deviation
                    se     = sd(A.eightsii_Biom) / sqrt(length(A.eightsii_Biom)-1) )
Aequi.Dens_Size <- ddply(Taxa_SD_Yoldia, .(factor_SD$STATION, factor_SD$SEASON), summarise,
                    N      = length(Yoldia_ss),                            ##number of observations
                    mean = mean(Yoldia_ss),                                 ##average
                    sd     = sd(Yoldia_ss),                                 ##standard deviation
                    se     = sd(Yoldia_ss) / sqrt(length(Yoldia_ss)-1) )
names(Aequi.Dens_Size)[names(Aequi.Dens_Size)=="factor_SD$STATION"] <- "STATION"
names(Aequi.Dens_Size)[names(Aequi.Dens_Size)=="factor_SD$SEASON"] <- "SEASON"
Aequi.Dens_Size$Interaction <- with(Aequi.Dens_Size, interaction(Aequi.Dens_Size$STATION, Aequi.Dens_Size$SEASON), drop = TRUE)
Taxa_SD_Yoldia$Yoldia_ss
Taxa_SD_Yoldia
##
b <- ggplot(Aequi.Dens_Size, aes(x=station.names.D, y=mean, fill=SEASON)) +                           ##make plot and specify x and y axes, fill=Code to give a seperate color to each bar/factor level
  geom_bar(position=position_dodge(), stat = "identity", color="black")+ # bargraph, default is grey, but we prefer black and white (theme_bw)
  ggtitle(" Aequiyoldia eightsii recruits")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), colour = "black", size = 0.7, width = 0.3, position=position_dodge(.9))+
  scale_y_continuous(expand=c(0,0),limits=c(0,2500))+ ## here ypu put in the limits of the graph y axis --> if you put it in wrong, the plot will not show up
  ylab(bquote('Individuals '*m^-2*''))+ ## the special things in a title you need to put inside '* blablabla*' and then close the ''
  scale_fill_manual(values=c("#009E73","#000000", "#E69F00"))+
  theme(text=element_text(size=14,  family="Helvetica Neue"))+ ## to have helvetica neue
  theme(legend.position="bottom") +       ## puts legend at the bottom
  theme(plot.title = element_text(face="italic"))+ ##turns main title in italic
  theme(axis.title.x = element_blank())+ 
  
  theme(legend.title=element_blank())+
  theme(legend.background = element_rect())+
  theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"))
print(b)
##
theme(axis.text.x=element_text(angle=90,hjust = 1))+## To have x axis labels vertical

LAT.DENS.graph.B <- ddply(LAT_DensFactors, .(STATION, SEASON), summarise,                      # make summary with following info:    
                        N      = length(Total_Weight),                            ##number of observations
                        mean = mean(Total_Weight),                                 ##average
                        sd     = sd(Total_Weight),                                 ##standard deviation
                        se     = sd(Total_Weight) / sqrt(length(Total_Weight)-1) )

b <- ggplot(LAT.DENS.graph, aes(x=station.names.D, y=mean, fill=SEASON)) +                           ##make plot and specify x and y axes, fill=Code to give a seperate color to each bar/factor level
  geom_bar(position=position_dodge(), stat = "identity", color="black")+ # bargraph, default is grey, but we prefer black and white (theme_bw)
  ggtitle("Laternula elliptica")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), colour = "black", size = 0.7, width = 0.3, position=position_dodge(.9))+
  scale_y_continuous(expand=c(0,0),limits=c(0,300))+ ## here ypu put in the limits of the graph y axis --> if you put it in wrong, the plot will not show up
  ylab(bquote('Individuals '*m^-2*''))+ ## the special things in a title you need to put inside '* blablabla*' and then close the ''
  scale_fill_manual(values=c("#009E73","#000000", "#E69F00"))+
  theme(text=element_text(size=14,  family="Helvetica Neue"))+ ## to have helvetica neue
  theme(legend.position="bottom") +       ## puts legend at the bottom
  theme(plot.title = element_text(face="italic"))+ ##turns main title in italic
  theme(axis.title.x = element_blank())+ 
  theme(legend.title=element_blank())+
  theme(legend.background = element_rect())+
  theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"))
print(b)

ylab(bquote('g '*C[org]*' '*m^-2*''))+### for biomass

## get Helvetica font
library(extrafont) ## to have Helvetica neue installed as other fonts
font_import() ##
y ##
loadfonts(device = "win") ##

theme(text=element_text(size=16,  family="Helvetica Neue")) ## how to add the specific font above

library(ggplot2)
##plot the data
cbbPalette <- c("#000000","#E69FOO","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","CC79A7")

a <- ggplot(Aequi.Biom, aes(x=station.names.D, y=mean, fill=SEASON)) +                           ##make plot and specify x and y axes, fill=Code to give a seperate color to each bar/factor level
  geom_bar(position=position_dodge(), stat = "identity", color="black")+ # bargraph, default is grey, but we prefer black and white (theme_bw)
  ggtitle("Aequiyoldia eightsii")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), colour = "black", size = 0.7, width = 0.3, position=position_dodge(.9))+
  scale_y_continuous(expand=c(0,0),limits=c(0,160000))+ ## here ypu put in the limits of the graph y axis --> if you put it in wrong, the plot will not show up
  ylab(bquote('mg '*C[org]*' '*m^-2*''))+ ## the special things in a title you need to put inside '* blablabla*' and then close the ''
  scale_fill_manual(values=c("#009E73","#000000", "#E69F00"))+
  theme(text=element_text(size=14,  family="Helvetica Neue"))+ ## to have helvetica neue
  theme(legend.position="bottom") +       ## puts legend at the bottom
  theme(plot.title = element_text(face="italic"))+ ##turns main title in italic
  theme(axis.title.x = element_blank())+ 
  theme(legend.title=element_blank())+
  theme(legend.background = element_rect())+
  theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"))
print(a)


##perform PERMANOVA with interaction factor
adonis(LAT$Total_Weight~LAT_factors$STATION*LAT_factors$SEASON, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor STATION p=***, DEPTH p=***, STATION x DEPTH p=**
results_LAT=aov(LAT$Densities ~ LAT_factors$STATION*LAT_factors$SEASON, data=LAT)
summary(results_LAT)
TukeyHSD(results_LAT)##ISLA D:Spring16-ISLA D:Spring15  p=0.0000000/FARO:Spring16-FARO:Spring15 p=NS/
##CREEK:Summer15-CREEK:Spring15 p=0.0000000 /ISLA D:Spring16-ISLA D:Spring15 p=0.0000000/ISLA D:Summer15-FARO:Summer15  p=0.0000000/
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first. We use
#euclidean distance cause it is on one single variable Density
d_LAT=vegdist(LAT$Densities, method="euclidean")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_LAT=(betadisper(d_LAT,LAT_factors$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_LAT)
par(cex.axis=0.6)
boxplot(bd_LAT, cex=0.6,main="L.elliptica Biom SEASONAL analysis, PermDisp Inter p=***", names=c("FARO / SUM15","ISLA D / SUM15","CREEK / SUM15","FARO / SPR15","ISLA D / SPR15","CREEK / SPR15","FARO / SPR16","ISLA D / SPR16","CREEK / SPR16"), las=2)

## pairwise Test ##
phoc_LAT <- with(LAT, betadisper(d_LAT,LAT_factors$interaction))
TukeyHSD(phoc_LAT) ## BIOM ## ISLA D.Spring16-ISLA D.Spring15 p=0.0257617/ ISLA D.Spring16-FARO.Spring16  p= 0.0003379
## DENS ## CREEK.Spring16-ISLA D.Spring16  p=0.0148393 / ISLA D.Spring16-FARO.Spring16   p=0.0002173 /ISLA D.Spring16-ISLA D.Spring15 p=0.0076097

## LINEAR regression for Laternula and Aequiyoldia ##

cor(LAT$Densities,Aequi_Dens$Aequiyoldia.eightsi)
Aequi_Dens
factor_SD##33
LAT_factors## 135

###### Correlation Lat/Yoldia ##### I need to make to sets of data with same observations to ran the correlation on since
# I have many more Laternula obs than Yoldia

library(kimisc)
LAT_DensFactors$Interaction <- with(LAT_DensFactors, interaction(LAT_DensFactors$STATION,  LAT_DensFactors$SEASON), drop = TRUE )#Let's add an interaction factor
df1 <- sample.rows(subset(LAT_DensFactors, Interaction == "FARO.Spring16"), 3)##randomly select observations from Laternula to match the sample size of Yoldia
df2 <- sample.rows(subset(LAT_DensFactors, Interaction == "CREEK.Spring16"), 3)
df3 <- sample.rows(subset(LAT_DensFactors, Interaction == "ISLAD.Spring16"), 3)
df4 <- sample.rows(subset(LAT_DensFactors, Interaction == "ISLAD.Spring15"), 4)
df5 <- sample.rows(subset(LAT_DensFactors, Interaction == "CREEK.Spring15"), 4)
df6 <- sample.rows(subset(LAT_DensFactors, Interaction == "FARO.Spring15"), 4)
df7 <- sample.rows(subset(LAT_DensFactors, Interaction == "CREEK.Summer15"), 4)
df8 <- sample.rows(subset(LAT_DensFactors, Interaction == "FARO.Summer15"), 4)
df9 <- sample.rows(subset(LAT_DensFactors, Interaction == "ISLAD.Summer15"), 4)
LAT_cor <- Reduce(function(x, y) merge(x, y, all=TRUE), list(df1, df2, df3,df4,df5,df6,df7,df8,df9))### I merge the df1-df9 together
Yoldia_Summary$Interaction <- with(Yoldia_Summary, interaction(Yoldia_Summary$STATION,  Yoldia_Summary$SEASON), drop = TRUE )#I generate the commun column between the two df
write.csv(LAT_cor, file = "LAT_cor.csv")## I am having problems merging the df coz of different order so I do it manually in Numbers
write.csv(Yoldia_Summary, file = "Yoldia_Summary.csv")
Cor_YolLat <- data.frame(YolLat_Cor) 
write.csv(Cor_YolLat, file = "Cor_YolLat.csv")


## correlation by group STATION/SEASON/Interaction ---> Biomass seems to be more (negatively) 
## correlated than densities --> changes in individual biomass may explain 
func <- function(Cor_YolLat)
{
  return(data.frame(COR = cor(Cor_YolLat$L.elliptica_Biom, Cor_YolLat$A.eightsii_Biom, subset(Cor_YolLat, subset(SEASON=="Spring16")))))## using function cor but it does not give you the pvalue of the correlation
  
}

ddply(Cor_YolLat, .(SEASON=="Spring16", STATION), func)

## cor.test gives you also pvalue of the correlation -> I made manually subsets with single levels of SEASON and STATION  
cor.test(Cor_YolLat_2$L.elliptica_Biom,Cor_YolLat_2$A.eightsii_Biom, conf.level = 0.95, method="pearson")

Cor_YolLat_2 <- data.frame(Cor_YL_Spring16)
### Spring16
#Pearson's product-moment correlation

#data:  Cor_YolLat_2$L.elliptica_Biom and Cor_YolLat_2$A.eightsii_Biom
#t = -2.6291, df = 7, p-value = 0.03396
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.93247944 -0.07660924
#sample estimates:
#       cor 
#-0.7048687



corfun<-function(Cor_YolLat) {
  corr=(cor.test(Cor_YolLat$L.elliptica_Biom,Cor_YolLat$A.eightsii_Biom, conf.level = 0.95, method="pearson"))
}

ddply(Cor_YolLat, .(SEASON,STATION), summarise,z=corfun(SEASON=="Spring16", STATION=="FARO")$statistic,
      pval=corfun(SEASON=="Spring16", STATION=="FARO")$p.value,
      tau.est=corfun(SEASON=="Spring16", STATION=="FARO")$estimate,
      alt=corfun(SEASON=="Spring16", STATION=="FARO")$alternative
) 




library(vegan)



ggplot(Cor_YolLat, aes(x = L.elliptica_Biom, y = A.eightsii_Biom, color = SEASON, shape = STATION)) + geom_point(size=3) +
  labs(title="Biomass",
       x="Laternula elliptica", y = "Aequiyoldia eightsii") +
  theme(axis.title.x = element_text(face = "italic", size="14"))+
  theme(axis.title.y = element_text(face = "italic", size="14"))


#STATION        COR
#1    FARO  0.2442353
#2   ISLAD  0.1319616
#3   CREEK -0.1749662

## Correlation of Biomass
#STATION        COR
#1    FARO 0.27563946
#2   ISLAD 0.09806043
#3   CREEK 0.33658898
#    SEASON        COR
#1 Summer15 -0.3408261
#2 Spring15 -0.2581153
#3 Spring16 -0.5066449

library(vegan)
adonis(Cor_YolLat$A.eightsii_Biom~Cor_YolLat$STATION*Cor_YolLat$SEASON, permutation = 999, method='euclidean')
d_Yoldia_Biom=vegdist(Cor_YolLat$A.eightsii_Biom, method="euclidean")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd__DD=(betadisper(d_Yoldia_Biom,Cor_YolLat$Interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd__DD)
TukeyHSD(bd__DD)
###  Tukey multiple comparisons of means ## Yoldia DENSITIES
##95% family-wise confidence level

##Fit: aov(formula = distances ~ group, data = df)

##$group SEASON Spring15 differs from Summer15 and Spring16
#diff        lwr       upr     p adj
#Spring15-Summer15  287.9760   30.14115  545.8109 0.0245670
#Spring16-Summer15 -134.4275 -438.28879  169.4339 0.5470323
#Spring16-Spring15 -422.4035 -726.26483 -118.5422 0.0036325
# STATION CREEK differs from FARO and ISLAD
#             diff       lwr      upr     p adj
#ISLAD-FARO  -76.49754 -302.0903 149.0952 0.7008332
#CREEK-FARO  406.20234  180.6096 631.7951 0.0001145
#CREEK-ISLAD 482.69988  257.1071 708.2927 0.0000042

#within level Spring15 (All different)
#             diff       lwr      upr     p adj
#CREEK-FARO   690.583804   307.52482 1073.64279 0.0000033
#CREEK-ISLAD  789.595104   406.53612 1172.65409 0.0000001
#CREEK-FARO  690.583804   307.52482 1073.64279 0.0000033

#within level Summer15 (All NS)
#             diff       lwr      upr     p adj
#ISLAD-FARO  -235.310734  -618.36972  147.74825 0.5860737
#CREEK-FARO  -233.309793  -616.36878  149.74919 0.5973897
#CREEK-ISLAD    2.000942  -381.05805  385.05993 1.0000000

#Within level Spring16  (All NS)
#ISLAD-FARO   167.974070  -342.77125  678.71939 0.9810240
#CREEK-FARO  430.779403   -79.96591  941.52472 0.1716351
#CREEK-ISLAD  262.805333  -247.93998  773.55065 0.7879423

#  Within level FARO
#Spring15 vs Summer15 p=*
#Spring15 vs Spring16 p= NS
#Spring16 vs summer15 p=**

# Within level CREEK
#Spring15 vs Summer15 p=**
#Spring15 vs Spring16 p= NS
#Spring16 vs summer15 p= NS

# Within level ISLAD
#Spring15 vs Summer15 p=NS
#Spring15 vs Spring16 p= NS
#Spring16 vs summer15 p= NS

###  Tukey multiple comparisons of means ## Yoldia BIOMASS
##95% family-wise confidence level

##Fit: aov(formula = distances ~ group, data = df)

#$group SEASON 
#                   diff       lwr        upr     p adj
#Spring15-Summer15 -29107.38 -41225.80 -16988.960 0.0000003
#Spring16-Summer15 -39979.30 -54260.99 -25697.601 0.0000000
#Spring16-Spring15 -10871.92 -25153.61   3409.778 0.1717897

#$group STATION They all differ
#               diff       lwr        upr     p adj
#ISLAD-FARO  -27194.17 -42791.25 -11597.090 0.0001923
#CREEK-FARO  -16001.89 -31598.97   -404.808 0.0429522
#CREEK-ISLAD  11192.28  -4404.80  26789.364 0.2082898

#within level Spring15 (FARO vs CREEK p=* / FARO vs ISLAD p=**)
#             diff       lwr      upr     p adj
#ISLAD-FARO    -31939.8936 -56008.47  -7871.3147 0.0017202
#CREEK-ISLAD  3912.6480 -20155.93  27981.2269 0.9998674
#CREEK-FARO  -28027.2455 -52095.82  -3958.6667 0.0103304

#within level Summer15 (CREEK vs FARO p=**)
#             diff       lwr      upr     p adj
#ISLAD-FARO  -20468.8913 -44537.47   3599.6875 0.1634443
#CREEK-FARO  -28866.1052 -52934.68  -4797.5264 0.0071624
#CREEK-ISLAD    -8397.2139 -32465.79  15671.3649 0.9726605

#Within level Spring16  (All NS)
#ISLAD-FARO   NS
#CREEK-FARO NS
#CREEK-ISLAD NS

#  Within level FARO
#Spring15 vs Summer15 p=**
#Spring15 vs Spring16 p= NS
#Spring16 vs summer15 p=*

# Within level CREEK
#Spring15 vs Summer15 p=*
#Spring15 vs Spring16 p= NS
#Spring16 vs summer15 p= .

# Within level ISLAD
#Spring15 vs Summer15 p=**
#Spring15 vs Spring16 p= NS
#Spring16 vs summer15 p= **

##### DEPTH ANALYSIS MACROFAUNA #######

Taxa_DB <- data.frame(MACRO_DEPTH_BIOM)## load multivariate Biomass of Macrofauna from xcl imported dataset
Taxa_DD <- data.frame(MACRO_DEPTH_DENS)##load multivariate Densities of Macrofauna from xcl imported dataset
factor_DD <- data.frame(MACRO_DEPTH_DENS_factors)#load factors for Densities Seasonal
factor_DB <- data.frame(MACRO_DEPTH_BIOM_factors)


factor_DB$STATION <- factor(factor_DB$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(factor_DB$STATION)[levels(factor_DB$STATION)=="ISLAD"] <- "ISLA D"
factor_DB$DEPTH <-  factor(factor_DB$DEPTH, levels=c("9","15"))## define factor SEASON levels
factor_DB$STATIONDEPTH  <- with(factor_DB, interaction(factor_DB$STATION,  factor_DB$DEPTH), drop = TRUE )## DEfine interactions for PERMDISP later and fpr boxplot

factor_DB

factor_DD$STATION <- factor(factor_DD$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(factor_DD$STATION)[levels(factor_DD$STATION)=="ISLAD"] <- "ISLA D"
factor_DD$DEPTH <-  factor(factor_DD$DEPTH, levels=c("9","15"))## define factor SEASON levels
factor_DD$STATIONDEPTH  <- with(factor_DD, interaction(factor_DD$STATION,  factor_DD$DEPTH), drop = TRUE )## DEfine interactions for PERMDISP later and fpr boxplot

factor_DD

str(factor_DD$STATION)
factor_DD$STATION

##call the factors and check
factor_DD$STATION
factor_DB$STATION
factor_DB$DEPTH
factor_DD$DEPTH
factor_DD$STATIONDEPTH
factor_DB$STATIONDEPTH

## metaNMDS ##
## nMDS on the DENSITIES of Macrofauna DEPTH

metaMDS_DD <- metaMDS(Taxa_DD, distance="bray")## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 19)
# empty plot
plot(metaMDS_DD, type = 'n',main='Macrofauna Multivariate Densities DEPTH ANALYSIS')### Non-metric fit R^2=0.973, Linear fit, R^2=0.854
# add points
points(metaMDS_DD, col = cols[factor_DD$STATION], pch = shps[factor_DD$DEPTH])
# add legend
legend('topleft', col=cols, legend=levels(factor_DD$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor_DD$DEPTH), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_DD, factor_DD$STATIONDEPTH, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses


stressplot(metaMDS_DD) ##YES!!!

##let's transform the data fourth root
Taxa_DD_2 <- Taxa_DD^0.25
##perform PERMANOVA with interaction factor
adonis(Taxa_DD_2~factor_DD$STATION*factor_DD$DEPTH, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor STATION p=***, DEPTH p=***, STATION x DEPTH p=**

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Taxa_DD=vegdist(Taxa_DD_2, method="bray")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd__DD=(betadisper(d_Taxa_DD,factor_DD$STATIONDEPTH, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd__DD)## for factor STATION p=NS , DEPTH p=NS, STATION x DEPTH p=NS hence  OK! there are no differences in the dispersions we can trust the PERMANOVA 100% 


par(cex.axis=0.8)
boxplot(bd__DD, cex=0.8,main="MACRO Multivar Dens DEPTH analysis, PermDisp STxDEPTH p=NS")
permutest(bd__DD)

##let's do a pairwise comparison to see what groups differ significantly (we have a=groups=2 and n=4)
##number of unique permutations=(2*4)!/(2!(4!)^2))= 35 --> minimum significance level =1/35 = 0.028 --> it is ok
library(vegan)
library(RVAideMemoire)
pairwise.perm.manova(d_Taxa_DD, factor_DD$STATIONDEPTH, test = "Wilks", nperm = 9999, 
                     progress = TRUE, p.method = "bonferroni")###DEPTH is highly significant p=***, ISLAD differs from both FARO and CREEK p= **,  STATION effect with differences between Faro and IslaD and Faro and Creek, but 
##no clear pairwise result from STATIONDEPTH but from the nMDS we see that FAro is the one station with the least distance between centroids




simDD =simper(Taxa_DD_2, group = factor_DD$STATIONDEPTH, permutations= 9999)
simDD
summary(simDD)
SIM_DEPTH_MACRODENS <- as.list(summary(simDD))
lapply(SIM_DEPTH_MACRODENS, function(x) write.table( data.frame(x), 'SIM_DEPTH.csv'  , append= T, sep=',' ))
capture.output(summary(SIM_DEPTH_MACRODENS), file = "SIM_DEPTH.txt")



## nMDS on the BIOMASS of Macrofauna DEPTH

metaMDS_DB <- metaMDS(Taxa_DB, distance="bray", autotransform = T)## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 19)
# empty plot
plot(metaMDS_DB, type = 'n',main='Macrofauna Multivariate Biomass DEPTH ANALYSIS')### Non-metric fit R^2=0.974, Linear fit, R^2=0.866
# add points
points(metaMDS_DB, col = cols[factor_DB$STATION], pch = shps[factor_DB$DEPTH])
# add legend
legend('topleft', col=cols, legend=levels(factor_DB$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(factor_DB$DEPTH), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_DB, factor_DB$STATIONDEPTH, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses


stressplot(metaMDS_DB) ##YES!!!

##let's transform the data fourth root
Taxa_DB_2 <- Taxa_DB^0.25
##perform PERMANOVA with interaction factor
adonis(Taxa_DB_2~factor_DB$STATION*factor_DB$DEPTH, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor STATION p=***, DEPTH p=***, STATION x DEPTH p=**

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Taxa_DB=vegdist(Taxa_DB_2, method="bray")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd__DB=(betadisper(d_Taxa_DB,factor_DB$STATIONDEPTH, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd__DB)## for factor STATION p=**, DEPTH p=***, STATION x DEPTH p=*** hence  there are differences in the dispersions of the centroid we have to consider it part of the cause for the
##effect we see from the PERMANOVA 


par(cex.axis=0.8)
boxplot(bd__DB, cex=0.8,main="MACRO Multivar Biom DEPTH analysis, PermDisp STxDEPTH p=***")
permutest(bd__DB)

##let's do a pairwise comparison to see what groups differ significantly (we have a=groups=2 and n=4)
##number of unique permutations=(2*4)!/(2!(4!)^2))= 35 --> minimum significance level =1/35 = 0.028 --> it is ok
library(vegan)
library(RVAideMemoire)
pairwise.perm.manova(d_Taxa_DB, factor_DB$DEPTH, test = "Wilks", nperm = 9999, 
                     progress = TRUE, p.method = "bonferroni")###DEPTH is highly significant p=***, STATION effect with ISLAD differs from FARO,  
#no clear pairwise result from STATIONDEPTH but from the nMDS we see that FAro is the one station with the least distance between centroids




simDB =simper(Taxa_DB_2, group = factor_DB$STATIONDEPTH, permutations= 9999)
simDB
summary(simDB)
SIM_DEPTH_MACROBIOM <- as.list(summary(simDB))
lapply(SIM_DEPTH_MACROBIOM, function(x) write.table( data.frame(x), 'SIM_DEPTH_BIOM.csv'  , append= T, sep=',' ))
capture.output(summary(SIM_DEPTH_MACROBIOM), file = "SIM_DEPTH_BIOM.txt")


############## MEIOFAUNA ###############


###### SEASONAL MEIOFAUNA #########

Meio_Season_D <- data.frame(SEASON_Meio_Dens_Taxa)
Meio_Season_Factor <- data.frame(SEASON_Meio_Dens_Factors)
Meio_Season_Biom <-  data.frame(SEASON_MEIO_BIOM)
Meio_Season_Biom_Factor <- data.frame(SEASON_BIOM_MEIO_factors)


str(SEASON_Meio_Dens_Taxa)

Meio_Season_Factor$STATION <- factor(Meio_Season_Factor$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
Meio_Season_Factor$SEASON <-  factor(Meio_Season_Factor$SEASON, levels=c("Summer15","Spring15","Winter15","Spring16"))## define factor SEASON levels
## DEfine interactions for PERMDISP later and fpr boxplot
Meio_Season_Factor$interaction <- with(factor, interaction(Meio_Season_Factor$STATION,  Meio_Season_Factor$SEASON), drop = TRUE)## DEfine interactions for PERMDISP later and fpr boxplot
levels(Meio_Season_Factor$STATION)[levels(Meio_Season_Factor$STATION)=="ISLAD"] <- "ISLA D"

Meio_Season_Biom_Factor$STATION <- factor(Meio_Season_Biom_Factor$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(Meio_Season_Biom_Factor$STATION)[levels(Meio_Season_Biom_Factor$STATION)=="ISLAD"] <- "ISLA D"
Meio_Season_Biom_Factor$SEASON <-  factor(Meio_Season_Biom_Factor$SEASON, levels=c("Summer15","Spring15","Winter15","Spring16"))## define factor SEASON levels
## DEfine interactions for PERMDISP later and fpr boxplot
Meio_Season_Biom_Factor$interaction <- with(Meio_Season_Biom_Factor, interaction(Meio_Season_Biom_Factor$STATION,  Meio_Season_Biom_Factor$SEASON), drop = TRUE)## DEfine interactions for PERMDISP later and fpr boxplot



str(Meio_Season_Factor$interaction)
Meio_Season_Factor$STATION##call the factors to see they have been defined properly
Meio_Season_Factor$SEASON
Meio_Season_Biom_Factor

Meio_Season_Biom_Factor

library(vegan)
## nMDS on the Biomass of Macrofauna Seasonal
Meio_Season_D2 <- Meio_Season_D^0.25
Meio_Season_D2
metaMDS_SeaD <- metaMDS(Meio_Season_D2, distance="bray", autotransform = FALSE)## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 25, 18)
# empty plot
plot(metaMDS_SeaD, type = 'n',main='Meiofauna Densities SEASONAL analysis')
# add points
points(metaMDS_SeaD, col = cols[Meio_Season_Factor$STATION], pch = shps[Meio_Season_Factor$SEASON])
# add legend
legend('topleft', col=cols, legend=levels(Meio_Season_Factor$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(Meio_Season_Factor$SEASON), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_SeaD, Meio_Season_Factor$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

stressplot(metaMDS_SeaD) ##YES!!!

### Now some PERMANOVA on the densities of MEIOFAUNA 
##adonis relies on a long-understood phenomenon that allows one to partition 
## sums of squared deviations from a centroid in two different ways, he most widely recognized method, used, e.g., for ANOVA and MANOVA,
##is to first identify the relevant centroids and then to calculate the squared deviations from these points

##let's transform the data fourth root
Meio_Season_D2 <- Meio_Season_D^0.25
##perform PERMANOVA with interaction factor
adonis(Meio_Season_D2~Meio_Season_Factor$STATION*Meio_Season_Factor$SEASON, permutation = 9999, method='bray')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Season_D2=vegdist(Meio_Season_D2, method="bray")


##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioSeaD2=(betadisper(d_Meio_Season_D2,Meio_Season_Factor$STATION, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaD2)## for factor STATION p=NS  hence  OK! there are no differences in the dispersions we can trust the PERMANOVA 


##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioSeaD2=(betadisper(d_Meio_Season_D2,Meio_Season_Factor$SEASON, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaD2)## for factor SEASON p=NS  OK! there are no differences in the dispersions that can cause the significan PERMANOVA p
##Subsequently, we can conduct a PERMDISP on this distance matrix.


bd_MeioSeaD2=(betadisper(d_Meio_Season_D2,Meio_Season_Factor$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaD2)## for factor STATION*SEASON p=0.03  therefore there are differences in the dispersions 

## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioSeaD2, cex=0.6,las = 2 ,names = c("FARO / SUM15","ISLA D / SUM15","CREEK / SUM15","FARO / WIN15", "ISLA D / WIN15","FARO / SPR15","ISLA D / SPR15","CREEK / SPR15","FARO / SPR16","ISLA D / SPR16","CREEK / SPR16"), main="MEIO Multivar Dens SEASONAL analysis, PermDisp STATIONxSEASON p=*")


permutest(bd_MeioSeaD2)## for factor STATIONYEAR p=** hence there are differences in the dispersion of variances  at each STATION between YEARS and we should be careful considering the interaction term




metaMDS_SeaB <- metaMDS(Meio_Season_Biom, distance="bray", autotransform = T)## I prefer this type

cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 25, 18)
# empty plot
plot(metaMDS_SeaB, type = 'n',main='Meiofauna Biomass SEASONAL analysis')
# add points
points(metaMDS_SeaB, col = cols[Meio_Season_Biom_Factor$STATION], pch = shps[Meio_Season_Biom_Factor$SEASON])
# add legend
legend('topleft', col=cols, legend=levels(Meio_Season_Biom_Factor$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(Meio_Season_Biom_Factor$SEASON), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_SeaB, Meio_Season_Biom_Factor$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

stressplot(metaMDS_SeaB) #

d_Meio_Season_Biom=vegdist(Meio_Season_Biom, method="bray")

bd_MeioSeaBiom=(betadisper(d_Meio_Season_Biom,Meio_Season_Factor$interaction, type = "centroid"))
permutest(bd_MeioSeaBiom)## Permdisp was significant for both SEASON and Interaction factor p=**

simB =simper(Meio_Season_Biom, group = Meio_Season_Factor$interaction, permutations= 9999)
summary(simB)

##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioSeaBiom, cex=0.6,las = 2, 
        names = c("FARO / SUM15","ISLA D / SUM15","CREEK / SUM15","FARO / SPR15","ISLA D / SPR15","CREEK / SPR15","FARO / WIN15","ISLA D/ WIN15","FARO / SPR16","ISLA D / SPR16","CREEK / SPR16"),
        main="MEIO Multivar Biom SEASONAL analysis, PermDisp STATIONxSEASON p=*")
permutest(bd_MeioSeaBiom)## 

###anyway, since for the Meiofauna I only have 4 taxa I will run an anova on each taxon to check for effects in between groups

## Two-way ANOVA statistics on meiofauna Biomass Copepods
results_Meio_Sea_Biom = aov(Meio_Season_Biom$Copepoda ~ Meio_Season_Biom_Factor$interaction, data=Meio_Season_Biom)
aov(results_Meio_Sea_Biom, type = "III") # two way anova for unbalanced design with interaction on raw data selection (dataX$variable)
summary(results_Meio_Sea_Biom)
## ANOVA STATION*SEASON p=0.00737
TukeyHSD(results_Meio_Sea_Biom) ##Post hoc for two-way Anova
TK<-TukeyHSD(results_Meio_Sea_Biom)
## check assumptions
##2. Homogeneity of variances ##OK if p<0.05 the variance is not assumed to be equal, otherwise it is equal
plot(results_Meio_Sea_Biom, 1) ## the residuals versus fits plot is used to check the homogeneity of variances
library(car)
leveneTest(Meio_Season_Biom$Copepoda ~ Meio_Season_Biom_Factor$interaction, data = Meio_Season_Biom) ## p= *--> ?-NOT Ok cause the null Hp is rejected : Null Hp = variances are homogeneous
## 3. Normality
plot(results_Meio_Sea_Biom, 2) ## Normality plot of the residuals. The quantiles of the residuals are plotted against the quantiles of the normal distribution. A 45-degree reference line is also plotted.
# Extract the residuals
aov_residuals <- residuals(object = results_Meio_Sea_Biom)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals) # p= *** so signficant NOT OK, if p value is significant it rejects the null Hp of normal distribution
library(moments) ## for skewness and kurtosis
skewness(Meio_Season_Biom$Copepoda, na.rm = TRUE) ## = 2.44 which is rather high-- to analyse skewness (density of data left - or right +) removing missing NA data prior to modelling --> acceptable skewness range between -3 and +3
kurtosis(Meio_Season_Biom$Copepoda, na.rm = TRUE) ## = 8.31 TOO HIGH!
### I need to perform non-parametric tests 
kruskal.test(Meio_Season_Biom$Copepoda~Meio_Season_Biom_Factor$STATION, data = Meio_Season_Biom)
kruskal.test(Meio_Season_Biom$Copepoda~Meio_Season_Biom_Factor$SEASON, data = Meio_Season_Biom)
kruskal.test(Meio_Season_Biom$Copepoda~Meio_Season_Biom_Factor$interaction, data = Meio_Season_Biom)##  --> p=*** 
##Kruskal-Wallis chi-squared = 36.996, df = 10, p-value = 5.669e-05
## non parametric two-way anova Kruskall Wallis test on interaction factors STATION*SEASON
pairwise.wilcox.test(Meio_Season_Biom$Copepoda, Meio_Season_Biom_Factor$interaction, p.adj = "bonf", data=Meio_Season_Biom)
warnings()


## Adonis on UNOVARIATE data - Biomass Copepods
##perform PERMANOVA with interaction factor
adonis(Meio_Season_Biom$Copepoda~Meio_Season_Biom_Factor$STATION*Meio_Season_Biom_Factor$SEASON, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##Terms added sequentially (first to last)

##                                                                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##Meio_Season_Biom_Factor$STATION                                 2     84770   42385 25.1165 0.57507  0.001 ***
##Meio_Season_Biom_Factor$SEASON                                  3      1535     512  0.3031 0.01041  0.836    
##Meio_Season_Biom_Factor$STATION:Meio_Season_Biom_Factor$SEASON  5      7103    1421  0.8419 0.04819  0.536    
##Residuals                                                      32     54002    1688         0.36634           
##Total                                                          42    147410                 1.00000           
##---
 ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Season_Cope=vegdist(Meio_Season_Biom$Copepoda, method="euclidean")
####Subsequently, we can conduct a PERMDISP on this euclidean distance matrix.
bd_MeioSeaDCope=(betadisper(d_Meio_Season_Cope,Meio_Season_Biom_Factor$STATION, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaDCope)## for factor STATION p=***  hence there are  differences in the dispersions from centroids we have to consider this as a possible effect as well as the station effect trust the PERMANOVA 

par(cex.axis=0.6)
boxplot(bd_MeioSeaDCope, cex=0.6,las = 2, main="MEIO Copepoda Biom SEASONAL analysis, PermDisp STATION p=***")
permutest(bd_MeioSeaBiom)## 


#### Seasonal Biomass Nematoda ####
## adonis on univariate data - Biomass Nematoda
##perform PERMANOVA with interaction factor
adonis(Meio_Season_Biom$Nematoda~Meio_Season_Biom_Factor$STATION*Meio_Season_Biom_Factor$SEASON, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=**, SEASON p=NS, STATION x SEASON p=NS
##Terms added sequentially (first to last)

##                                                                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
##Meio_Season_Biom_Factor$STATION                                 2   4105256 2052628  8.7193 0.27195  0.005 **
##Meio_Season_Biom_Factor$SEASON                                  3   1871736  623912  2.6503 0.12399  0.057 . 
##Meio_Season_Biom_Factor$STATION:Meio_Season_Biom_Factor$SEASON  5   1585581  317116  1.3471 0.10503  0.254   
##Residuals                                                      32   7533218  235413         0.49903          
##Total                                                          42  15095791                 1.00000          
##---
  ##Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Season_Nema=vegdist(Meio_Season_Biom$Nematoda, method="euclidean")
####Subsequently, we can conduct a PERMDISP on this euclidean distance matrix.
bd_MeioSeaDNema=(betadisper(d_Meio_Season_Nema,Meio_Season_Biom_Factor$STATION, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaDNema)## for factor STATION p=***  hence there are  differences in the dispersions from centroids we have to consider this as a possible effect as well as the station effect trust the PERMANOVA 

par(cex.axis=0.6)
boxplot(bd_MeioSeaDNema, cex=0.6,las = 2, main="MEIO Nematoda Biom SEASONAL analysis, PermDisp STATION p=***")
permutest(bd_MeioSeaBiom)## 

## Two-way ANOVA statistics on meiofauna Biomass Nematodes
results_Meio_Sea_Biom = aov(Meio_Season_Biom$Nematoda ~ Meio_Season_Factor$STATION*Meio_Season_Factor$SEASON, data=Meio_Season_Biom)
aov(results_Meio_Sea_Biom, type = "III") # two way anova for unbalanced design with interaction on raw data selection (dataX$variable)
summary(results_Meio_Sea_Biom)
## ANOVA STATION p=0.0009***
TukeyHSD(results_Meio_Sea_Biom) ##Post hoc for two-way Anova
TK<-TukeyHSD(results_Meio_Sea_Biom)
## check assumptions
##2. Homogeneity of variances ##OK if p<0.05 the variance is not assumed to be equal, otherwise it is equal
plot(results_Meio_Sea_Biom, 1) ## the residuals versus fits plot is used to check the homogeneity of variances
library(car)
leveneTest(Meio_Season_Biom$Nematoda ~ Meio_Season_Factor$STATION*Meio_Season_Factor$SEASON, data = Meio_Season_Biom) ## p= NS-->  Ok cause the null Hp is rejected : Null Hp = variances are homogeneous
## 3. Normality
plot(results_Meio_Sea_Biom, 2) ## Normality plot of the residuals. The quantiles of the residuals are plotted against the quantiles of the normal distribution. A 45-degree reference line is also plotted.
# Extract the residuals
aov_residuals <- residuals(object = results_Meio_Sea_Biom)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals) # p= *** signficant NOT OK, if p value is significant it rejects the null Hp of normal distribution
library(moments) ## for skewness and kurtosis
skewness(Meio_Season_Biom$Nematoda, na.rm = TRUE) ## = 2.11 which is rather high-- to analyse skewness (density of data left - or right +) removing missing NA data prior to modelling --> acceptable skewness range between -3 and +3
kurtosis(Meio_Season_Biom$Nematoda, na.rm = TRUE) ## = 8.64 TOO HIGH!
### I need to perform non-parametric tests 
kruskal.test(Meio_Season_Biom$Nematoda~Meio_Season_Factor$STATION, data = Meio_Season_Biom)
kruskal.test(Meio_Season_Biom$Nematoda~Meio_Season_Factor$SEASON, data = Meio_Season_Biom)
kruskal.test(Meio_Season_Biom$Nematoda~Meio_Season_Factor$interaction, data = Meio_Season_Biom)##  --> p=*** 
##Kruskal-Wallis chi-squared = 24.437, df = 10, p-value = 0.006521
## non parametric two-way anova Kruskall Wallis test on interaction factors STATION*SEASON
pairwise.wilcox.test(Meio_Season_Biom$Nematoda, Meio_Season_Factor$interaction, p.adj = "bonf", data=Meio_Season_Biom)
warnings()

library(vegan)
Meio_Season_Biom$Nematoda2 <- Meio_Season_Biom$Nematoda^0.25
##perform PERMANOVA with interaction factor for univariate Biomass Nematoda
adonis(Meio_Season_Biom$Nematoda2~Meio_Season_Factor$STATION*Meio_Season_Factor$SEASON, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*

##I can't perform a PERMDISP since I can't make a Resemblance Matrix with only one species...I will do a pair-wise 
library(RVAideMemoire)
pairwise.perm.manova(Meio_Season_Biom,Meio_Season_Factor$interaction, test = c("Wilks"), nperm = 999, 
                     progress = TRUE, p.method = "fdr")

#### Seasonal Biomass Polychaeta ####
## adonis on univariate data - Biomass Polychaetes
##perform PERMANOVA with interaction factor
adonis(Meio_Season_Biom$Polychaeta~Meio_Season_Biom_Factor$STATION*Meio_Season_Biom_Factor$SEASON, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=**, SEASON p=NS, STATION x SEASON p=NS
##Terms added sequentially (first to last)

##                                                              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##Meio_Season_Biom_Factor$STATION                                 2    925144  462572  7.7575 0.27772  0.001 ***
##Meio_Season_Biom_Factor$SEASON                                  3    128716   42905  0.7195 0.03864  0.611    
##Meio_Season_Biom_Factor$STATION:Meio_Season_Biom_Factor$SEASON  5    369201   73840  1.2383 0.11083  0.306    
##Residuals                                                      32   1908126   59629         0.57281           
##Total                                                          42   3331187                 1.00000           
##---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Season_Poly=vegdist(Meio_Season_Biom$Polychaeta, method="euclidean")
####Subsequently, we can conduct a PERMDISP on this euclidean distance matrix.
bd_MeioSeaDPoly=(betadisper(d_Meio_Season_Poly,Meio_Season_Biom_Factor$STATION, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaDPoly)## for factor STATION p=***  hence there are  differences in the dispersions from centroids we have to consider this as a possible effect as well as the station effect trust the PERMANOVA 

par(cex.axis=0.6)
boxplot(bd_MeioSeaDPoly, cex=0.6,las = 2, main="MEIO Polychaeta Biom SEASONAL analysis, PermDisp STATION p=***")
permutest(bd_MeioSeaBiom)## 

#### Seasonal Biomass Cumacea ####
## adonis on univariate data - Biomass Polychaetes
##perform PERMANOVA with interaction factor
adonis(Meio_Season_Biom$Cumacea~Meio_Season_Biom_Factor$STATION*Meio_Season_Biom_Factor$SEASON, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=**, SEASON p=NS, STATION x SEASON p=NS
##Terms added sequentially (first to last)

##                                                              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##Meio_Season_Biom_Factor$STATION                                 2    925144  462572  7.7575 0.27772  0.001 ***
##Meio_Season_Biom_Factor$SEASON                                  3    128716   42905  0.7195 0.03864  0.611    
##Meio_Season_Biom_Factor$STATION:Meio_Season_Biom_Factor$SEASON  5    369201   73840  1.2383 0.11083  0.306    
##Residuals                                                      32   1908126   59629         0.57281           
##Total                                                          42   3331187                 1.00000           
##---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Season_Cuma=vegdist(Meio_Season_Biom$Cumacea, method="euclidean")
####Subsequently, we can conduct a PERMDISP on this euclidean distance matrix.
bd_MeioSeaDCuma=(betadisper(d_Meio_Season_Cuma,Meio_Season_Biom_Factor$STATIONSEASON, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioSeaDCuma)## for factor STATION p=***  hence there are  differences in the dispersions from centroids we have to consider this as a possible effect as well as the station effect trust the PERMANOVA 

par(cex.axis=0.6)
boxplot(bd_MeioSeaDCuma, cex=0.6,las = 2, main="MEIO Cumacea Biom SEASONAL analysis, PermDisp STATION p=***")
permutest(bd_MeioSeaBiom)## 


##### DEPTH ANALYSIS MEIOFAUNA #####

Meio_Depth_D <- data.frame(DEPTH_Meio_Taxa_Dens)
Meio_Depth_Factor <- data.frame(DEPTH_Meio_Factors_Dens)
Meio_Depth_Factor_B <- data.frame(Depth_Meio_Biom_Factors)
Meio_Depth_B <-  data.frame(DEPTH_Meio_Biom)

Meio_Depth_D
Meio_Depth_Factor

str(Meio_Depth_Factor)

Meio_Depth_Factor$STATION <- factor(Meio_Depth_Factor$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(Meio_Depth_Factor$STATION)[levels(Meio_Depth_Factor$STATION)=="ISLAD"] <- "ISLA D"
Meio_Depth_Factor$DEPTH <-  factor(Meio_Depth_Factor$DEPTH, levels=c("9","15"))## define factor DEPTH levels
## DEfine interactions for PERMDISP later and fpr boxplot
Meio_Depth_Factor$interaction <- with(factor, interaction(Meio_Depth_Factor$STATION,  Meio_Depth_Factor$DEPTH), drop = TRUE)## DEfine interactions for PERMDISP later and fpr boxplot



str(Meio_Depth_Factor$interaction)
Meio_Depth_Factor$interaction
Meio_Depth_Factor$STATION##call the factors to see they have been defined properly
Meio_Depth_Factor$DEPTH

library(vegan)

##metaNMDS

metaMDS_DepthD <- metaMDS(Meio_Depth_D, distance="bray", autotransform = TRUE)## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# empty plot
plot(metaMDS_DepthD, type = 'n',main='Meiofauna Densities DEPTH analysis')
# add points
points(metaMDS_DepthD, col = cols[Meio_Depth_Factor$STATION], pch = shps[Meio_Depth_Factor$DEPTH])## non metric fit R^2=0.963, linear fit R^2=0.776
# add legend
legend('topleft', col=cols, legend=levels(Meio_Depth_Factor$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(Meio_Depth_Factor$DEPTH), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_DepthD, Meio_Depth_Factor$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

stressplot(metaMDS_DepthD) ##YES!!!

##DEPTH MEIO BIOMASS

Meio_Depth_Factor_B$STATION <- factor(Meio_Depth_Factor_B$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(Meio_Depth_Factor_B$STATION)[levels(Meio_Depth_Factor_B$STATION)=="ISLAD"] <- "ISLA D"
Meio_Depth_Factor_B$DEPTH <-  factor(Meio_Depth_Factor_B$DEPTH, levels=c("9","15"))## define factor DEPTH levels
## DEfine interactions for PERMDISP later and fpr boxplot
Meio_Depth_Factor_B$interaction <- with(factor, interaction(Meio_Depth_Factor_B$STATION,  Meio_Depth_Factor_B$DEPTH), drop = TRUE)## DEfine interactions for PERMDISP later and fpr boxplot



metaMDS_DepthB <- metaMDS(Meio_Depth_B, distance="bray", autotransform = TRUE)## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# empty plot
plot(metaMDS_DepthB, type = 'n',main='Meiofauna Biomass DEPTH analysis')
# add points
points(metaMDS_DepthB, col = cols[Meio_Depth_Factor_B$STATION], pch = shps[Meio_Depth_Factor_B$DEPTH])
# add legend
legend('topleft', col=cols, legend=levels(Meio_Depth_Factor_B$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(Meio_Depth_Factor_B$DEPTH), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_DepthB, Meio_Depth_Factor_B$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

stressplot(metaMDS_DepthB) ##YES!!!


##DEPTH univariate PERMANOVA Copepoda Biomass##
Meio_Depth_B2Cope <- Meio_Depth_B$Copepoda^0.25
##perform PERMANOVA with interaction factor
adonis(Meio_Depth_B2Cope~Meio_Depth_Factor_B$STATION*Meio_Depth_Factor_B$DEPTH, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=NS / DEPTH p=NS / STATION*DEPTH p=*** SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Depth_BCope=vegdist(Meio_Depth_B2Cope, method="euclidean")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioDepthBCope=(betadisper(d_Meio_Depth_BCope,Meio_Depth_Factor_B$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioDepthBCope)## for factor interaction p=NS



##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioDepthBCope, cex=0.6,las = 2, 
        main="Meiofauna Copepoda Biomass DEPTH, PermDisp STATIONxSEASON p=NS")

##DEPTH univariate PERMANOVA Nematoda Biomass##
Meio_Depth_B2Nema <- Meio_Depth_B$Nematoda^0.25
##perform PERMANOVA with interaction factor
adonis(Meio_Depth_B2Nema~Meio_Depth_Factor_B$STATION*Meio_Depth_Factor_B$DEPTH, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=** / DEPTH p=NS / STATION*DEPTH p= ** SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Depth_BNema=vegdist(Meio_Depth_B2Nema, method="euclidean")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioDepthBNema=(betadisper(d_Meio_Depth_BNema,Meio_Depth_Factor_B$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioDepthBNema)## for factor interaction p=NS



##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioDepthBNema, cex=0.6,las = 2, 
        main="Meiofauna Nematoda Biomass DEPTH, PermDisp STATIONxSEASON p=NS")

##DEPTH univariate PERMANOVA Polychaeta Biomass##
Meio_Depth_B2Poly <- Meio_Depth_B$Polychaeta^0.25
##perform PERMANOVA with interaction factor
adonis(Meio_Depth_B2Poly~Meio_Depth_Factor_B$STATION*Meio_Depth_Factor_B$DEPTH, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=*** / DEPTH p=*** / STATION*DEPTH p= *** SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Depth_BPoly=vegdist(Meio_Depth_B2Poly, method="euclidean")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioDepthBPoly=(betadisper(d_Meio_Depth_BPoly,Meio_Depth_Factor_B$DEPTH, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioDepthBPoly)## for factor interaction p=** -> differences in the dispersion from groups centroids is also causing the significance
# in the Permanova univariate test



##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioDepthBPoly, cex=0.6,las = 2, 
        main="Meiofauna Polychaeta Biomass DEPTH, PermDisp STATIONxSEASON p=***")

##DEPTH univariate PERMANOVA Cumacea Biomass##
Meio_Depth_B2Cuma <- Meio_Depth_B$Cumacea^0.25
##perform PERMANOVA with interaction factor
adonis(Meio_Depth_B$Cumacea~Meio_Depth_Factor_B$STATION*Meio_Depth_Factor_B$DEPTH, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=* / DEPTH p=* / STATION*DEPTH p= NS SO we consider STATION and DEPTH separately
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Depth_BCuma=vegdist(Meio_Depth_B$Cumacea, method="euclidean")
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioDepthBCuma=(betadisper(d_Meio_Depth_BCuma,Meio_Depth_Factor_B$DEPTH, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioDepthBCuma)## for factor interaction p=** -> differences in the dispersion from groups centroids is also causing the significance
# in the Permanova univariate test



##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioDepthBCuma, cex=0.6,las = 2, 
        main="Meiofauna Cumacea Biomass DEPTH, PermDisp STATION p=NS")

## Pairwise univariate Meiofauna Biomass DEPTH analysis ##




## DEPTH PERMANOVA Densities
Meio_Depth_D2 <- Meio_Depth_D^0.25
##perform PERMANOVA with interaction factor
adonis(Meio_Depth_D2~Meio_Depth_Factor$STATION*Meio_Depth_Factor$DEPTH, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=*** / DEPTH p=*** / STATION*DEPTH p= *** SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Meio_Depth_D=vegdist(Meio_Depth_D2, method="bray")## distance matric on 4rth root transformed data
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_MeioDepthD2=(betadisper(d_Meio_Depth_D,Meio_Depth_Factor$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_MeioDepthD2)##  Permdisp for factor interaction = NS hence  OK! there are no differences in the dispersions we can trust the PERMANOVA 



d_Meio_Depth_D=vegdist(Meio_Depth_D2, method="bray")

bd_MeioDepthD=(betadisper(d_Meio_Depth_D,Meio_Depth_Factor$interaction, type = "centroid"))
permutest(bd_MeioDepthD)## Permdisp was significant STATION and DEPTH but not for the Interaction factor



##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_MeioDepthD, cex=0.6,las = 2, 
        main="Meiofauna Multivariate Densities DEPTH, PermDisp STATIONxSEASON p=NS")
permutest(bd_MeioDepthD)## 


###### Long-term analysis Meiofauna ######
Long_Meio_D <- data.frame(LONG_Meio_Dens)
Long_Meio_D_factors <- data.frame(LONG_Meio_Dens_factors)
Long_Meio_B <- data.frame(LONG_Meio_Biom)
Long_Meio_B_factors <- data.frame(LONG_Meio_Biom_Factors)


str(Long_Meio_D_factors)
Long_Meio_D_factors$STATION
Long_Meio_D_factors$YEAR
Long_Meio_D_factors$interaction


Long_Meio_D_factors$STATION <- factor(Long_Meio_D_factors$STATION, levels=c("Faro","IslaD","Creek"))## define factor STATION levels
Long_Meio_D_factors$YEAR <-  factor(Long_Meio_D_factors$YEAR, levels=c("2010","2015"))## define factor DEPTH levels
## DEfine interactions for PERMDISP later and fpr boxplot
Long_Meio_D_factors$interaction <- with(factor, interaction(Long_Meio_D_factors$STATION,  Long_Meio_D_factors$YEAR), levels= c("FARO.2010","ISLA D.2010","CREEK.2010", "FARO.2015", "ISLA D.2015", "CREEK.2015"), drop = TRUE)## DEfine interactions for PERMDISP later and fpr boxplot
levels(Long_Meio_D_factors$STATION)[levels(Long_Meio_D_factors$STATION)=="IslaD"] <- "ISLA D"
levels(Long_Meio_D_factors$STATION)[levels(Long_Meio_D_factors$STATION)=="Faro"] <- "FARO"
levels(Long_Meio_D_factors$STATION)[levels(Long_Meio_D_factors$STATION)=="Creek"] <- "CREEK"
Long_Meio_D_factors$interaction <- with(Long_Meio_D_factors, interaction(Long_Meio_D_factors$STATION,Long_Meio_D_factors$YEAR), drop = TRUE)#
##metaNMDS Densities Long-term
library(vegan)

metaMDS_LongD <- metaMDS(Long_Meio_D, distance="bray", autotransform = TRUE)## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# empty plot
plot(metaMDS_LongD, type = 'n',main='Meiofauna Densities FIVE YEAR TREND')
# add points
points(metaMDS_LongD, col = cols[Long_Meio_D_factors$STATION], pch = shps[Long_Meio_D_factors$YEAR])
# add legend
legend('topleft', col=cols, legend=levels(Long_Meio_D_factors$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(Long_Meio_D_factors$YEAR), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_LongD, Long_Meio_D_factors$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ploygon" allows to color fill the ellipses

stressplot(metaMDS_LongD) ##YES!!!

Long_Meio_D2 <- Long_Meio_D^0.25
##perform PERMANOVA with interaction factor
adonis(Long_Meio_D2~Long_Meio_D_factors$STATION*Long_Meio_D_factors$YEAR, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=*** / DEPTH p=*** / STATION*DEPTH p= *** SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Long_Meio_D2=vegdist(Long_Meio_D2, method="bray")## distance matric on 4rth root transformed data
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_Long_Meio_D2=(betadisper(d_Long_Meio_D2,Long_Meio_D_factors$YEAR, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_Long_Meio_D2)##  Permdisp for factor interaction = NS hence  OK! there are no differences in the dispersions we can trust the PERMANOVA 



d_Long_Meio_D2=vegdist(Long_Meio_D2, method="bray")

bd_Long_Meio_D2=(betadisper(Long_Meio_D2,Long_Meio_D_factors$interaction, type = "centroid"))
permutest(d_Long_Meio_D2)## Permdisp was significant STATION and DEPTH but not for the Interaction factor



##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.6)
boxplot(bd_Long_Meio_D2, cex=0.6,las = 2, 
        main="Meiofauna Multivar Dens FIVE YEAR TREND, PermDisp STxYEAR p=NS")
permutest(bd_Long_Meio_D2)## 

Long_Meio_B_factors$STATION <- factor(Long_Meio_B_factors$STATION, levels=c("Faro","IslaD","Creek"))## define factor STATION levels
levels(Long_Meio_B_factors$STATION)[levels(Long_Meio_B_factors$STATION)=="IslaD"] <- "ISLA D"
levels(Long_Meio_B_factors$STATION)[levels(Long_Meio_B_factors$STATION)=="Faro"] <- "FARO"
levels(Long_Meio_B_factors$STATION)[levels(Long_Meio_B_factors$STATION)=="Creek"] <- "CREEK"
Long_Meio_B_factors$YEAR <-  factor(Long_Meio_B_factors$YEAR, levels=c("2010","2015"))## define factor YEAR levels


Long_Meio_B_factors$interaction <- with(Long_Meio_B_factors, interaction(Long_Meio_B_factors$STATION,Long_Meio_B_factors$YEAR), drop = TRUE)### DEfine interactions for PERMDISP later and fpr boxplot
Long_Meio_B_factors
metaMDS_LongB <- metaMDS(Long_Meio_B, distance="bray", autotransform = TRUE)## I prefer this type

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21)
# empty plot
plot(metaMDS_LongB, type = 'n',main='Meiofauna Biomass FIVE YEAR TREND')
# add points
points(metaMDS_LongB, col = cols[Long_Meio_B_factors$STATION], pch = shps[Long_Meio_B_factors$YEAR])
# add legend
legend('topleft', col=cols, legend=levels(Long_Meio_B_factors$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(Long_Meio_B_factors$YEAR), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(metaMDS_LongB, Long_Meio_B_factors$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ploygon" allows to color fill the ellipses

stressplot(metaMDS_LongB) ##YES!!!




##perform PERMANOVA LONG TERM  with interaction factor
adonis(Long_Meio_B~Long_Meio_B_factors$STATION*Long_Meio_B_factors$YEAR, permutation = 999, method='bray')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=* / YEAR p=NS / STATION*YEAR p= * SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Long_Meio_B=vegdist(Long_Meio_B, method="bray")## distance matric on 4rth root transformed data
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_Long_Meio_B=(betadisper(d_Long_Meio_B,Long_Meio_B_factors$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_Long_Meio_B)##  PERMDISP STATIONxYEAR p=NS! OK 

##an F-test (ANOVA, if the assumptions are met), or a permutation test.
par(cex.axis=0.7)
boxplot(bd_Long_Meio_B, cex=0.7,las = 2, 
        main="Meiofauna Multivar Biom FIVE YEAR TREND, PermDisp STxYEAR p=NS")
permutest(bd_Long_Meio_B)## 

##perform PERMANOVA LONG TERM  with interaction factor Nematodes

adonis(Long_Meio_B$Copepoda~Long_Meio_B_factors$STATION*Long_Meio_B_factors$YEAR, permutation = 999, method='euclidean')## the PERMANOVA results give me a significant interaction factor SATION p=***, SEASON p=*, STATION x SEASON p=*
##STATION p=* / YEAR p=NS / STATION*YEAR p= * SO we consider interaction factor
##The PERMDISP uses the distance matrix (as also used by the PERMANOVA), but here we have to make it ourselves first
d_Long_Meio_BNema=vegdist(Long_Meio_B$Copepoda, method="euclidean")## distance matric on 4rth root transformed data
##Subsequently, we can conduct a PERMDISP on this distance matrix.
bd_Long_Meio_BNema=(betadisper(d_Long_Meio_BNema,Long_Meio_B_factors$interaction, type = "centroid"))
## To check if the dispersion of variances differ between the groups, we can use a visual method, 
##an F-test (ANOVA, if the assumptions are met), or a permutation test.
permutest(bd_Long_Meio_BNema)## 


########### DistLm Analysis on balanced designed files (I adapted the #of Rows between ENV and BIO datasets) #######

##### DistLm SEASON MEIO ####


Dist_ENV_SEA <- as.data.frame(DistLm_ENV_SEASON)
Dist_BIOM_SEA <- as.data.frame(BIOM_MEIO_SEASON_DistLm)
Dist_DENS_SEA <-  as.data.frame(SEASON_Meio_Dens_Taxa_DistLm)
Dist_ENV_SEA
Dist_BIOM_SEA
Dist_DENS_SEA

library(vegan) 
## We would like to quantify how much of the variation in the multivariate data cloud is explained by the environmental variables. 
## We therefore create a model that includes all our environmental variables.
SEA_DENS_model = adonis(vegdist(Dist_DENS_SEA, "bray") ~ Dist_ENV_SEA$RATIO + Dist_ENV_SEA$MGS + Dist_ENV_SEA$TN+ Dist_ENV_SEA$TOC)
## we call the model
SEA_DENS_model

##Call:
#adonis(formula = vegdist(Dist_DENS_SEA, "bray") ~ Dist_ENV_SEA$RATIO +      Dist_ENV_SEA$MGS + Dist_ENV_SEA$TN + Dist_ENV_SEA$TOC) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Dist_ENV_SEA$RATIO  1   0.34381 0.34381  4.5315 0.12679  0.025 *
  #Dist_ENV_SEA$MGS    1   0.19271 0.19271  2.5400 0.07107  0.081 .
#Dist_ENV_SEA$TN     1   0.00735 0.00735  0.0969 0.00271  0.949  
#Dist_ENV_SEA$TOC    1   0.04333 0.04333  0.5711 0.01598  0.544  
#Residuals          28   2.12436 0.07587         0.78345         
#Total              32   2.71156                 1.00000         
#---
  #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
### What we see is th&t only the RATIO seems to significantly explain the variation in the multivariate densited of the SEASONAL analysis 
#although it explains only 12% of the variation whereas most of the variation 78% is in the residuals


##dbRDA###
##Subsequently, we can VISUALIZE these results. We use a dbRDA for this, which is a constrained analysis that also uses the 
##environmental variables as input variables in the analysis. Because our response data are biological (count data), 
##we use a Bray-Curtis distance matrix. We only use the significant environmental variables from the previous analysis.

DistLm_ENV_SEASON$STATION <- factor(DistLm_ENV_SEASON$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(DistLm_ENV_SEASON$STATION)[levels(DistLm_ENV_SEASON$STATION)=="ISLAD"] <- "ISLA D"
DistLm_ENV_SEASON$SEASON <-  factor(DistLm_ENV_SEASON$SEASON, levels=c("Summer15","Spring15","Winter15","Spring16"))## define factor SEASON levels
## DEfine interactions for PERMDISP later and fpr boxplot
DistLm_ENV_SEASON$interaction <- with(DistLm_ENV_SEASON, interaction(DistLm_ENV_SEASON$STATION,  DistLm_ENV_SEASON$SEASON), drop = TRUE)## Define interactions for bdRDA
DistLm_ENV_SEASON$STATION
##capscale for distance-based redundancy analysis (db-RDA), based on metric multidimensional scaling,
# a.k.a. principal coordinates analysis (PCoA).
dbRDA_SEA_MEIO = capscale(vegdist(Dist_DENS_SEA, "bray") ~ Dist_ENV_SEA$RATIO + Dist_ENV_SEA$MGS) ## we only use the significant variables from the DistLm analysis
plot(dbRDA_SEA_MEIO)
summary(dbRDA_SEA_MEIO)

#Accumulated constrained eigenvalues
#Importance of components:
#                       CAP1     CAP2
#Eigenvalue            0.5308 0.008259
#Proportion Explained  0.9847 0.015323
#Cumulative Proportion 0.9847 1.000000


# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 17, 25)
# empty plot
plot(dbRDA_SEA_MEIO,type="points", xlab = "CAP1 [98%]", ylab = "CAP2 [1%]", main='MEIO DENS dbRDA SEASONAL analysis')
# add points
points(dbRDA_SEA_MEIO, col = cols[DistLm_ENV_SEASON$STATION], pch = shps[DistLm_ENV_SEASON$SEASON])## non metric fit R^2=0.963, linear fit R^2=0.776
# add legend
legend('topleft', col=cols, legend=levels(DistLm_ENV_SEASON$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(DistLm_ENV_SEASON$SEASON), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(dbRDA_SEA_MEIO, DistLm_ENV_SEASON$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

## Do the same for Biomass ## 
## We would like to quantify how much of the variation in the multivariate data cloud is explained by the environmental variables. 
## We therefore create a model that includes all our environmental variables.
SEA_BIOM_model = adonis(vegdist(Dist_BIOM_SEA, "bray") ~ Dist_ENV_SEA$RATIO + Dist_ENV_SEA$MGS + Dist_ENV_SEA$TN+ Dist_ENV_SEA$TOC)
## we call the model
SEA_BIOM_model

##Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Dist_ENV_SEA$RATIO  1    0.0816 0.08164  0.5776 0.01748  0.655  
#Dist_ENV_SEA$MGS    1    0.4609 0.46093  3.2610 0.09871  0.017 *
#  Dist_ENV_SEA$TN     1    0.0740 0.07400  0.5235 0.01585  0.682  
#Dist_ENV_SEA$TOC    1    0.0954 0.09540  0.6749 0.02043  0.585  
#Residuals          28    3.9577 0.14135         0.84753         
#Total              32    4.6696                 1.00000         
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
###only MGS explains the Biomass trends in the Meio dataset but only explains very little not even 10%, and 84% is unexplianed vairance (residuals)
dbRDA_BIOM_MEIO = capscale(vegdist(Dist_BIOM_SEA, "bray") ~ Dist_ENV_SEA$RATIO + Dist_ENV_SEA$MGS) 
plot(dbRDA_BIOM_MEIO)
summary(dbRDA_BIOM_MEIO)

#Accumulated constrained eigenvalues
#Importance of components:
 #                      CAP1    CAP2
#Eigenvalue            0.4857 0.08021
#Proportion Explained  0.8583 0.14173
#Cumulative Proportion 0.8583 1.00000


# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(15, 21, 17, 25)
# empty plot
plot(dbRDA_BIOM_MEIO,type="points",xlab = "CAP1 [85%]", ylab = "CAP2 [14%]",main='MEIO BIOM dbRDA SEASONAL analysis, MGS p=*')
# add points
points(dbRDA_BIOM_MEIO, col = cols[DistLm_ENV_SEASON$STATION], pch = shps[DistLm_ENV_SEASON$SEASON])## non metric fit R^2=0.963, linear fit R^2=0.776
# add legend
legend('topleft', col=cols, legend=levels(DistLm_ENV_SEASON$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(DistLm_ENV_SEASON$SEASON), pch = shps, cex = 1)
# draw dispersion ellipses around data points with standard error
ordiellipse(dbRDA_BIOM_MEIO, DistLm_ENV_SEASON$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses


###### DistLm DEPTH Analysis MEIOFAUNA  ######

ENV_DEP <- as.data.frame(DistLm_ENV_Depth)
Dist_BIOM_DEP <- as.data.frame(DEPTH_Meio_Biom_DistLm)
Dist_DENS_DEP <-  as.data.frame(DEPTH_Meio_Taxa_Dens_DistLm)
ENV_DEP
Dist_BIOM_DEP
Dist_DENS_DEP

library(vegan) 
## We would like to quantify how much of the variation in the multivariate data cloud is explained by the environmental variables. 
## We therefore create a model that includes all our environmental variables.
DEP_DENS_model = adonis(vegdist(Dist_DENS_DEP, "bray") ~ ENV_DEP$RATIO + ENV_DEP$MGS + ENV_DEP$TN+ ENV_DEP$TOC)
## we call the model
DEP_DENS_model

###Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Dist_ENV_DEP$RATIO  1   0.23057 0.23057  7.8414 0.20765  0.009 **
#Dist_ENV_DEP$MGS    1   0.06145 0.06145  2.0898 0.05534  0.139   
#Dist_ENV_DEP$TN     1   0.11654 0.11654  3.9634 0.10495  0.050 * 
#Dist_ENV_DEP$TOC    1   0.31956 0.31956 10.8679 0.28780  0.002 **
#Residuals          13   0.38226 0.02940         0.34426          
#Total              17   1.11038                 1.00000          
#---
  #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
### So it seems that the data is now explained by RATIO (20% of variation), TN (10% of variation) and TOC (28% of variation) whch means that
## summed up the three variables explain 68% of the total variation in meiofauna communities in the DEPTH anlaysis

##dbRDA###
##Subsequently, we can visualize these results. We use a dbRDA for this, which is a constrained analysis that also uses the 
##environmental variables as input variables in the analysis. Because our response data are biological (count data), 
##we use a Bray-Curtis distance matrix. We only use the significant environmental variables from the previous analysis.

ENV_DEP$STATION <- factor(ENV_DEP$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(ENV_DEP$STATION)[levels(ENV_DEP$STATION)=="ISLAD"] <- "ISLA D"
ENV_DEP$DEPTH <-  factor(ENV_DEP$DEPTH, levels=c("9","15"))## define factor SEASON levels
## DEfine interactions for PERMDISP later and fpr boxplot
ENV_DEP$interaction <- with(ENV_DEP, interaction(ENV_DEP$STATION,  ENV_DEP$DEPTH), drop = TRUE)## Define interactions for bdRDA
ENV_DEP$STATION
##capscale for distance-based redundancy analysis (db-RDA), based on metric multidimensional scaling,
# a.k.a. principal coordinates analysis (PCoA).
dbRDA_DEP_MEIO = capscale(vegdist(Dist_DENS_DEP, "bray") ~ ENV_DEP$RATIO + ENV_DEP$TN +  ENV_DEP$TOC) 
plot(dbRDA_DEP_MEIO)
summary(dbRDA_DEP_MEIO)
##Accumulated constrained eigenvalues
##Importance of components:
#                      CAP1    CAP2     CAP3
#Eigenvalue            0.494 0.03196 0.004666
#Proportion Explained  0.931 0.06024 0.008794
#Cumulative Proportion 0.931 0.99121 1.000000


ENV_DEP$DEPTH

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(17, 21)
# empty plot
plot(dbRDA_DEP_MEIO,xlab = "CAP1 [93%]", ylab = "CAP2 [6%]", type="points",main='MEIO DENS dbRDA DEPTH analysis RATIO p=**/ TN p=*/ TOC p=**')##how to add percentage to CAP from the eigenvalues analysis

# add points
points(dbRDA_DEP_MEIO, col = cols[ENV_DEP$STATION], pch = shps[ENV_DEP$DEPTH])## non metric fit R^2=0.963, linear fit R^2=0.776
# add legend
legend('topleft', col=cols, legend=levels(ENV_DEP$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(ENV_DEP$DEPTH), pch = shps, cex = 1)

# draw dispersion ellipses around data points with standard error
ordiellipse(dbRDA_DEP_MEIO, ENV_DEP$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses

## BIOMASS##
## We would like to quantify how much of the variation in the multivariate data cloud is explained by the environmental variables. 
## We therefore create a model that includes all our environmental variables.
DEP_BIOM_model = adonis(vegdist(Dist_BIOM_DEP, "bray") ~ ENV_DEP$RATIO + ENV_DEP$MGS + ENV_DEP$TN+ ENV_DEP$TOC)
## we call the model
DEP_BIOM_model

###Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#ENV_DEP$RATIO  1   0.20291 0.20291  3.0676 0.11051  0.031 * 
#ENV_DEP$MGS    1   0.13592 0.13592  2.0548 0.07403  0.103   
#ENV_DEP$TN     1   0.45465 0.45465  6.8733 0.24762  0.002 **
#ENV_DEP$TOC    1   0.18271 0.18271  2.7622 0.09951  0.051 . 
#Residuals     13   0.85991 0.06615         0.46833          
#Total         17   1.83609                 1.00000  of permutations: 999


#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
### So it seems that the data is now explained by RATIO (11% of variation), TN (24% of variation) and almost significant TOC (9% of variation) whch means that
## summed up the three variables explain only  35% of the total variation in meiofauna communities in the DEPTH anlaysis and most of it remaind unexplained

##dbRDA###
##Subsequently, we can visualize these results. We use a dbRDA for this, which is a constrained analysis that also uses the 
##environmental variables as input variables in the analysis. Because our response data are biological (count data), 
##we use a Bray-Curtis distance matrix. We only use the significant environmental variables from the previous analysis.
##Function: capscale() for distance-based redundancy analysis (db-RDA), based on metric multidimensional scaling,
# a.k.a. principal coordinates analysis (PCoA).first, a distance matrix is calculated using the distance measure of choice. 
#Next, a principle coordinates analysis (PCoA) is done on the matrix. 
#Finally, the eigenvalues obtained in the PCoA are plugged into an RDA. 

ENV_DEP$STATION <- factor(ENV_DEP$STATION, levels=c("FARO","ISLAD","CREEK"))## define factor STATION levels
levels(ENV_DEP$STATION)[levels(ENV_DEP$STATION)=="ISLAD"] <- "ISLA D"
ENV_DEP$DEPTH <-  factor(ENV_DEP$DEPTH, levels=c("9","15"))## define factor SEASON levels
## DEfine interactions for PERMDISP later and fpr boxplot
ENV_DEP$interaction <- with(ENV_DEP, interaction(ENV_DEP$STATION,  ENV_DEP$DEPTH), drop = TRUE)## Define interactions for bdRDA
ENV_DEP$STATION
##capscale for distance-based redundancy analysis (db-RDA), based on metric multidimensional scaling,
# a.k.a. principal coordinates analysis (PCoA).first, a distance matrix is calculated using the distance measure of choice. 
#Next, a principle coordinates analysis (PCoA) is done on the matrix. 
#Finally, the eigenvalues obtained in the PCoA are plugged into an RDA. 
dbRDA_DEP_MEIO_BIOM = capscale(vegdist(Dist_BIOM_DEP, "bray") ~ ENV_DEP$RATIO + ENV_DEP$TN +  ENV_DEP$TOC) 
plot(dbRDA_DEP_MEIO_BIOM)
summary(dbRDA_DEP_MEIO_BIOM)
##Accumulated constrained eigenvalues
##Importance of components:
  #                     CAP1   CAP2   CAP3
#Eigenvalue            0.4994 0.2140 0.1091
#Proportion Explained  0.6072 0.2602 0.1326
#Cumulative Proportion 0.6072 0.8674 1.0000


ENV_DEP$DEPTH

# set colors and shapes
cols = c('darkgreen', 'black', 'red')
shps = c(17, 21)
# empty plot
plot(dbRDA_DEP_MEIO_BIOM,xlab = "CAP1 [60%]", ylab = "CAP2 [26%]", type="points",cex= 1, main='MEIO BIOM dbRDA DEPTH analysis RATIO p=*/ TN p=**/ TOC p=.')##how to add percentage to CAP from the eigenvalues analysis

# add points
points(dbRDA_DEP_MEIO_BIOM, col = cols[ENV_DEP$STATION], pch = shps[ENV_DEP$DEPTH])## non metric fit R^2=0.963, linear fit R^2=0.776
# add legend
legend('topleft', col=cols, legend=levels(ENV_DEP$STATION), pch = 16, cex = 1)
legend('bottomleft', legend=levels(ENV_DEP$DEPTH), pch = shps, cex = 1)

# draw dispersion ellipses around data points with standard error
ordiellipse(dbRDA_DEP_MEIO_BIOM, ENV_DEP$interaction, display = "sites", kind = "se", conf=0.95,draw="polygon", border=c('darkgreen', 'black', 'red'), col=c('darkgreen', 'black', 'red') ,label = F)##draw="ployon" allows to color fill the ellipses



### DistLm Long Term Meiofauna ##
ENV_Long <- data.frame(DistLm_Meio_DENS_ENV)
Dist_Dens_MeioL <- data.frame(DistLm_Meio_DENS)
Dist_Biom_MeioL <- data.frame(LONG_Meio_Biom)
ENV_Long_B <-  data.frame(DistLm_Meio_BIOM_ENV)

  library(vegan) 
## We would like to quantify how much of the variation in the multivariate data cloud is explained by the environmental variables. 
## We therefore create a model that includes all our environmental variables.
Long_DENSMeio_model = adonis(vegdist(Dist_Dens_MeioL, "bray") ~ ENV_Long$RATIO + ENV_Long$MGS + ENV_Long$TN+ ENV_Long$TOC)
## we call the model
Long_DENSMeio_model  
    
### Permutation: free
##Number of permutations: 999

##Terms added sequentially (first to last)

##              Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
##ENV_Long$RATIO  1   0.26426 0.264264  3.9492 0.20262  0.056 .
##ENV_Long$MGS    1   0.01968 0.019681  0.2941 0.01509  0.749  
##ENV_Long$TN     1   0.13195 0.131948  1.9718 0.10117  0.160  
##ENV_Long$TOC    1   0.01841 0.018405  0.2750 0.01411  0.769  
##Residuals      13   0.86992 0.066917         0.66700         
##Total          17   1.30421                  1.00000         
##--- Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## We would like to quantify how much of the variation in the multivariate data cloud is explained by the environmental variables. 
## We therefore create a model that includes all our environmental variables.
Long_BIOMMeio_model = adonis(vegdist(Dist_Biom_MeioL, "bray") ~ ENV_Long_B$RATIO + ENV_Long_B$MGS + ENV_Long_B$TN+ ENV_Long_B$TOC)
## we call the model
Long_BIOMMeio_model  
## Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#ENV_Long_B$RATIO  1   0.09184 0.09184  0.8739 0.04608  0.404  
#ENV_Long_B$MGS    1   0.20902 0.20902  1.9888 0.10486  0.135  
#ENV_Long_B$TN     1   0.35413 0.35413  3.3695 0.17766  0.047 *  ## very borderline
#ENV_Long_B$TOC    1   0.07713 0.07713  0.7339 0.03869  0.518  
#Residuals        12   1.26116 0.10510         0.63271         
#Total            16   1.99328                 1.00000         
#---Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
