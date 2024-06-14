library(ggplot2)
library(ggbeeswarm)
library(lemon)

#path to working data directory
wd <- "./"
setwd(wd)

#set Figure paths
data.path <- "data/Figure4/"
sc.path <- "code/Figures/"
fig.path <- "Figures/Figure_4/"

#Set general Figure variables
vir.lvls <- c("SFFV","MoMLV", "FeLV", "SNV", "KoRV", "CrERV")

#panel B: GFP intensities + intensity stats
source(paste0(sc.path, "Figure4_B1_FACS_dotplot.R"))
source(paste0(sc.path, "Figure4_B2_FACS_intensity_stat.R"))

#panel C: %GFP+ cells in time
source(paste0(sc.path, "Figure4_C_PercInTime.R"))

#panel D: ddPCR (copy number) 
source(paste0(sc.path, "Figure4_D_vCN_Barplot.R"))

#panel E: active copies (FCS vs CN)
source(paste0(sc.path, "Figure4_E_FCS_vs_ddPCR.R"))

#Supplementary figures
#SNV: expression intensity and active copies over time
source(paste0(sc.path, "Figure4_Supp_SNV_FACS_dotplot.R"))
source(paste0(sc.path, "Figure4_Supp_SNV_FCS_vs_ddPCR_time.R"))


