library(ggplot2)
library(ggbeeswarm)
library(lemon)

#path to working data directory
wd.path <- "./"
setwd(wd.path)

data.path <- "data/Figure2/"
sc.path <- "code/Figures/"
fig.path <- "Figures/Figure_2/"

in.cat <- c("wt", "WA", "CBX")

cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

#Figure size
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297

#ChromHMM distances
#Impact + distances
source(paste0(sc.path, "Figure2_A-B_ChromHMM.R"))

#Integration into different chromatin compartments: LADs + Hi-C subsomparments
source(paste0(sc.path, "Figure2_C-D_SubC.R"))

#Random targeting into different chromatin compartments: LADs + Hi-C subsomparments
source(paste0(sc.path, "Figure2_E_SubC_random.R"))


