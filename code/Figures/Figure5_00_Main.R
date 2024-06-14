library(ggplot2)
library(ggbeeswarm)
library(lemon)

#set paths
wd <- "./"

setwd(paste0(primus, wd.path))

sc.path <- "code/Figures/"
data.path <- "data/Figure5/"
fig.path <- "Figures/Figure_5/"

#name of the summary table file
sumtab.name <- "KnockIn_Summary.txt"

#Set variables for samples
smpl.is.names <- c("DoT:3.1", "DoT:6.2", "DoT:6.3", "DoT:7.3", "IFT20")

#Create table with group names
group.tab <- data.frame(Gr = c("F", "M", "K", "A"),
                        Virus = c("FeLV", "MoMLV", "KoRV", "ASLV"))

#Load Summary table
print("loading summary table")
sum.df <- read.table(paste0(data.path, sumtab.name),
                     header = TRUE, sep = "\t")
#set levels
sum.df$Virus <- factor(sum.df$Virus, levels = group.tab$Virus)
sum.df$TargetSite <- factor(sum.df$TargetSite, levels = smpl.is.names)

print("Summary table loaded:")
head(sum.df)
print("...")
tail(sum.df)


#Create plot
#Main text and supplementary figures characterizing knock-in clones
source(paste0(sc.path, "Figure5_KI_plots.R"))

#Distance of target sites to nearest genomic segments
source(paste0(sc.path, "Figure5_Supp_Target-Segment_distance"))

#Target copy number plot
source(paste0(sc.path, "Figure5_Supp_Target_CN.R"))

#Bulk copy number of ASLV-nucleofected cells
source(paste0(sc.path, "Figure5_Supp_KI_ASLV_CN_bulk.R"))



