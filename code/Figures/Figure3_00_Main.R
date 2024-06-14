library(ggplot2)
library(ggbeeswarm)
library(lemon)

#path to working data directory
wd.path <- "./"
setwd(wd.path)

data.path <- "data/Figure3/
sc.path <- "code/Figures/"
fig.path <- "Figures/Figure_3/"

in.cat <- c("wt", "W390A", "CBX")

cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

#panel B GFP intensities
source(paste0(sc.path, "Figure3_B_FACS_3dpi.R"))

#panel C and D: ddPCR (copy number) and active copies
source(paste0(sc.path, "Figure3_C-D_FCS_vs_ddPCR_Main.R"))
source(paste0(sc.path, "Figure3_C-D_FCS_vs_ddPCR_Plot.R"))

#Supplementary figures
#data need to be loaded with Figure3_C-D_FCS_vs_ddPCR_Main.R
source(paste0(sc.path, "Figure3_Supp_Plot.R"))


