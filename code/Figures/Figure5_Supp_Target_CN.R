#Create supplementary figure
#Copy Number of target allele

#NormCN value was calculated as follows
#norm.smpl <- "K562-RPP30"
#norm.val <- mean(vcn$CopiesPer20uLWell[vcn$Sample == norm.smpl])
#vcn$NormCN <- vcn$CopiesPer20uLWell / norm.val


#Load data
dd.file <- "K562_target_CN_all.txt"

vcn <- read.table(paste0(data.path, dd.file),
                  sep = "\t",  dec = ".", header = TRUE)

#Set target site levels
target.lvls <- c("RPP30", "GAPDH", "DoT:3.1", "DoT:6.2", "DoT:6.3", "DoT:7.3", "IFT20")

#Colors for plots
#Set color palettes
gr.cols <- c("gray90", "gray50", "gray10", "gray0")
cb.cols.inv <- c("#004D40", "#1E88E5", "#D81B60", "#FFC107")

#Plot
data.to.plot <- vcn[vcn$Target %in% target.lvls &
                      grepl("K562", vcn$Sample),]

supp.ta.p <- ggplot(data = data.to.plot, aes(x = factor(Target, levels = target.lvls),
                                             y = log2(NormCN),
                                             color = Target)) +
  color_palette("viridis") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_beeswarm(shape = 21, size = .75, stroke = .5, fill = NA, cex = 3) +
  coord_cartesian(xlim = c(0.5, 7.5),
                  ylim = c(-2, 2.1),
                  expand = FALSE) +
  xlab("Target") +
  ylab("log2 relative copies") +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.text.x = element_text(color = "black", face = "bold", size = ax.tex.size, angle = 45, hjust = 1.1, vjust = 1.2),
        axis.title = element_text(color = "black", face = "bold", size = ax.tit.size),
  )

supp.ta.p

ggsave(paste0(fig.path, "Figure5_S8_Target_CN.png"), supp.ta.p,
       width = 40, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
