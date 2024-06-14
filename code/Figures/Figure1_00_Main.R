#Create plots for MLV bulk FACS measurements

#Plots created
#p.bulk.sp - scatter plots of cytometric data (GFP vs FCS)
#p.bulk - %GFP / time of all variants / stock volumes
#p.bulk.rel - p.bulk realtive change to expression at 3 dpi

library(ggplot2)
library(ggbeeswarm)
library(lemon)

#path to working data directory
wd.path <- "./"
setwd(wd.path)

data.path <- "data/Figure1/"
sc.path <- "code/Figures/"
fig.path <- "Figures/Figure_1/"


#Set variables
#maximum number of cells to be plotted
max.cells <- 10000

pos.value <- "GFP+"

#Set variables for plotting
x.values <- "FCS.A"
x.values <- "GFP..488B..A"

x.name <- "FCS.A"
y.name <- "GFP"

vir.lvls <- "LdG"
var.lvls <- c("wt", "W390A", "CBX")
lvls.pop <- c( "GFP+", "GFP-")

#Set color palettes
gr.cols <- c("gray90", "gray50", "gray10")
cb.cols.inv <- c("#004D40", "#1E88E5", "#D81B60")

cols <- c("darkgreen", "grey60")
names(cols) <- lvls.pop

ax.tex.size <- 3
ax.tit.size <- 4

#data files:
#export.file <- "MLV_3dpi_fcs.txt"
#bulk.file <- "MLV_bulk.txt"
#sort.file <- "MLV_sorted.txt"
#clone.file <- "MLV_clones.txt"


#1) FlowCyt Dot plots
source(paste0(sc.path, "Figure1_B1_load_data.R"))
source(paste0(sc.path, "Figure1_B2_fcs.R"))

#2) Bulk expression in time (absolute + relative)
source(paste0(sc.path, "Figure1_C_relbulkExpr.R"))

#3) Expression of the bulk-sorted population
source(paste0(sc.path, "Figure1_D_sort.R"))

#Expression of the clonal populations
source(paste0(sc.path, "Figure1_E_clones.R"))
