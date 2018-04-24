##############################################################################
# Urban stormwater run-off promotes compression of saltmarshes by freshwater #
#                           plants and mangrove forests                      #
#                                                                            #
#               Ina Geedicke, Jens Oldeland, Michelle R. Leishman            #
#                                                                            #
#                        Data analysis code, STOTEN 2018                     #
##############################################################################


### SEDIMENT ANALYSES -------------------------------------------------------

# 1. Load packages ---------------------------------------------------------

library(ggplot2)  # plotting
library(Rmisc)    # summary statistics
library(coin)      # Wilcoxon Rank Sum Test
library(multcompView) # letters out of significant differences for plotting
library(magrittr) # for piping
library(data.table) # for setDT
require(PMCMR) # Kruskal-Wallis Test


# 2. Load data -------------------------------------------------------------

data <- read.csv("IGCh12016_EnvData.csv")


# 3. Environmental pollutants statistics, table 1 ----------------------------

### Reshape data: group by Type and Urbanisation
Out <- subset(data, Type == 'Outlet')
Con <- subset(data, Type == 'Control')
Urb <- subset(data, Urbanisation == 'Industrial')
Sub <- subset(data, Urbanisation == 'Residential')

### Summary statistics on Type and Urbanisation
sink("Data_Suburban_summary_mean min max.txt")
summary(Out)
sink()

sink("Data_Control_summary_mean min max.txt")
summary(Con)
sink()

sink("Data_Urban_summary_mean min max.txt")
summary(Urb)
sink()

sink("Data_Suburban_summary_mean min max.txt")
summary(Sub)
sink()

### Calculate standard deviation on summary stattistics
sapply(Out, sd)
sapply(Con, sd)
sapply(Urb, sd)
sapply(Sub, sd)

### Compare summary statisics: Mann-Whitney-U/Wilcoxon Test

# Function
Wilcox <- function(data, resp, nom){
  wilcox_test(formula= as.formula(paste(resp, "~", nom)), data= data) 
}

# loop to automate output per environmental pollutant
vecW <- names(data[,22:51]) # all parameters I want to test
vecnom <- names(data[,5:6]) # against Urbanisation and Type
resW <- list() # save the result in a list
resp <- vector() # create vector for p-value output
resz <- vector() #create vector for z value output

for(i in vecW){
  print(Wilcox(data, i, "Urbanisation"))
  resW[[i]] <- Wilcox(data, i, "Urbanisation")
  resp[i] <- pvalue(resW[[i]])
  resz[i] <- statistic(resW[[i]])
}

# save output in CSV file
write.csv(resz, 'Wilcox_z_Urban.csv')
write.csv(resp, 'Wilcox_p_Urban.csv')


for(i in vecW){
  print(Wilcox(data, i, "Type"))
  resW[[i]] <- Wilcox(data, i, "Type")
  resp[i] <- pvalue(resW[[i]])
  resz[i] <- statistic(resW[[i]])
}

# save output in CSV file
write.csv(resz, 'Wilcox_z_Type.csv')
write.csv(resp, 'Wilcox_p_Type.csv')


# 4. Supplementary, Figure A1 ------------------------------------------------

### summarize statistics (mean, sd, ci) for a parameter and plot them for all sites
### as a barplot with error bars (Table 1, Supplementary Data Fig A2)

# Function:
SitePlot<-function(data,x){
  # Calculate summary statistics per study site
  Sumstat <- summarySEwithin(data, measurevar=x, withinvars="Site",
                             idvar="Plot.Number", na.rm=FALSE, conf.interval=.95)
  # assign urbanisation type (residential/industrial) to each study site
  Sumstat$Urbanisation <-c("Industrial","Residential","Industrial","Residential",
                           "Industrial", "Industrial","Industrial", "Residential") 
  # plot pollutant element per site: barchart and standard deviation (Figure A2)
  p <- ggplot(Sumstat, aes(x=as.factor(Site), y=Sumstat[,3], fill= Urbanisation)) +
    geom_bar(position=position_dodge(), stat="identity", colour='black') +
    geom_errorbar(aes(ymin=Sumstat[,3]-sd, ymax=Sumstat[,3]+sd), width=.2,
                  position=position_dodge(.2)) +
    ylab(paste(x, "[mg/kg]"))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black")) +
    scale_fill_manual(values=c("white", "lightgrey")) +
    theme(axis.title.x=element_blank(), legend.title = element_blank())
  
  # export barcharts as jpg
  png(paste("plot_", x, ".jpg", sep = ""), width=1100, height=500, res=120)
  print(p)
  dev.off()
  # print statistics
  print(Sumstat)
}

# Loop
results <- vector('list')
vec <- names(data[,22:51])

for(i in vec){
  results[[i]] <- SitePlot(data,i)
}

saveRDS(results, "SummaryStatEnv-Site.rds")

### Do sites differ significantly per environmental pollutant?

# Function
SignifCodes <- function(data,x){
  # Kruskal-Wallis Test
  kruskal.test(formula = as.formula(paste(x, "~ Site", sep = "")), data= data)
  # Post-hoc Test
  tuk <- posthoc.kruskal.nemenyi.test(formula = as.formula(paste(x, "~ Site", sep = "")),
                                      data= data, dist="Chisquare")
  out.p <- get.pvalues(tuk)
  # apply letter code to statistical differences
  out.mcV <- multcompLetters(out.p, threshold=0.05) 
  print(out.mcV)
}

# Loop to automate function
result_letters <- vector('list')
vec <- names(data[,22:51])

for(i in vec){
  result_letters[[i]] <- SignifCodes(data,i)
}

### Compile results in table
result_letters_table <- do.call(rbind, result_letters) %>% as.data.frame(row.names=NULL)
setDT(result_letters_table, keep.rownames = TRUE)[]


### IMPACT OF NUTRIENTS & POLLUTANTS ON VEGETATION WITHIN WETLANDS  --------

# 1. Load Packages -----------------------------------------------------------
library(vegan) # for ordination analysis
library(dplyr) # To structure the data
library(partykit) # Toolkit for recursive partitioning, e.g. conditional trees

# 2. Load data ---------------------------------------------------------------
# CCA:
Species <- read.csv("20161103Data_Species_CCA.csv", row.names = "Site.ID") 
Env <- read.csv("20161103Data_Environment_CCA.csv", row.names = "Site.ID")
# row.names = Site.ID very important, because it associates every row to the site ID

# Conditional regression tree:
Data <- read.csv("IGCh12016_EnvData.csv") # Environmental pollutants data , 
                                          #including vegetation categories 

# 3. Canonical Correlation Analysis (CCA) ------------------------------------

### Transform community data: Cover % data square root -------------------------
Sp.sqrt <- sqrt(Species)

## - for abundance data take hellinger transformation sp.hell <- decostand (species, 'hell')
## - for standardizing of data (i.e. soil paramter) <- (scale(Env[6:35])

### Drop correlated enviromental parameters -----------------------------

## I repeatedly called vif.cca + sort(mycca) to identify the parameter with 
## the highest VIF value.

mycca <- cca(Sp.sqrt ~ Urbanisation + 
               Type + 
               Distance + 
               Total.N + 
               # Total.P + 
               pH + 
               salinity + 
               # TDR + 
               # S +
               # Cl +
               K + 
               Ca + 
               #Ti +
               # V +
               # Cr +
               # Mn +
               # Fe +
               # Ni +
               Cu + 
               Zn + 
               As+ 
               # Rb+
               # Sr + 
               Y + 
               Zr + 
               #Nb + 
               Mo + 
               Cd + 
               Ba + 
               # Ta + 
               #Pb + 
               Th +
               U, data = Env )

mycca
Vif.sorting <- vif.cca(mycca) 
sort(Vif.sorting)

### Ordistep function to identify the relevant parameters ----------------------

mod <-ordistep(mycca, scope = formula(mycca),perm.max=999,steps=100, direction="both")

###### these parameters remain:
#     - Zn            1 394.12 1.6333  0.065 . 
#     - Type          1 394.14 1.6495  0.050 * 
#     - Distance      1 394.22 1.7213  0.045 * 
#     - Total.N       1 394.66 2.1247  0.005 **
#     - pH            1 394.82 2.2758  0.005 **
#     - Ca            1 395.83 3.2186  0.005 **
#     - Urbanisation  1 396.26 3.6225  0.005 **
#     ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### Final CCA with remaining parameters -------------------------------------

mycca.redux <- cca(Sp.sqrt ~ Zn + Distance + Type + Total.N + pH + Urbanisation + Ca, data = Env)
mycca.redux <- cca(Sp.sqrt ~ Distance + Type + Total.N + pH + Urbanisation + Ca, data = Env)

### Plotting of CCA with 'vegan' --------------------------------------------

par(oma=c(3,3,3,3)) # make a larger background around the figure
colvec<- c("black", "darkgrey") # assign colour vector
shapevec <- c(21, 24) # assign a shape vector

plot(mycca.redux, type = "n", 
     xlim = c(-4.5,4), ylim = c(-1.5,1),
     scaling = 1,
     ylab = "Axis 2", xlab = "Axis 1", 
     cex.lab=1.6, cex.axis=1.2, bty = "o")

with(Env, points(mycca.redux, display = "sites",  # plot points for each plot: diff. shape and color
                 col = "black", bg = colvec[Urbanisation],
                 scaling = 1, pch = shapevec[Type], cex = 1.7))

fit <- envfit(mycca.redux ~ Distance + Total.N + pH + Ca, data = Env) # plots singular paramters with arrows
plot(fit, cex=1.5, axis=TRUE, col = "black" )

with(Env, legend(2.6,2, legend = levels(Urbanisation),bty = "n", # plot legend
                 col = colvec, pch = 21, pt.bg = colvec, cex = 1.3))
with(Env, legend(2.6,1, legend = levels(Type),bty = "n",  
                 col = "black", pch = shapevec, cex = 1.3))

# for printing as jpeg used ratio: 1200 x 695

dev.off()


###############################################################################
# Comments ----------------------------------------------------------------

## CCA including TDR and transformed environmental parameters were tried 
## as well. VIF selection and pattern stays the same, only TDR is too strong
## and other effects are too week. TDR is left out, because we already had 
## doubts on the measurment method. Transforming of environmental data (EnvS) 
## resulted the same pattern, only upside down. For simplified interpretation 
## we decided to leave environmental data untransformed.

###############################################################################

# 4. Conditional Regression Tree --------------------------------------------

### Clean data ----------------------------------------------------------------

# Data sets for Cover data of Saltmarshes, Mangroves and Exotics 
# with all environmental parameters

dataSR.SM <-  Data[,c(3, 5:7, 10, 22:51)]
dataSR.M <-  Data[,c(3, 5:7, 11, 22:51)]
dataSR.E <-  Data[,c(3, 5:7, 12, 22:51)]

dataC.SM <-  Data[,c(3, 5:7, 16, 22:51)]
dataC.M <-  Data[,c(3, 5:7, 17, 22:51)]
dataC.E <-  Data[,c(3, 5:7, 18, 22:51)]


### Conditional Tree using important variables from CCA  -------------------------------

# Function
CondTree <- function(data, i){
  Ctree <- ctree(formula = as.formula(paste(i, "~ Type + Distance + Urbanisation + Zn + 
                                            Total.N + pH + Ca", sep = "")), data = data)
  p <- plot(Ctree, main=paste("Conditional Inference Tree for", i, 
                              "Parameters after VIF"), gp = gpar(fontsize = 22))
  print(p)
}

# Conditional trees for cover data of vegetation types
C.SM.vif<- CondTree(dataC.SM, "C.Native.Saltmarsh")
C.M.vif<- CondTree(dataC.M, "C.Mangrove")
C.E.vif<- CondTree(dataC.E, "C.Exotics")


