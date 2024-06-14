#Plot the GFP copy number per cell in bulk population
#Knock-in of ASLV vector

vcn.2 <- read.table(paste0(data.path, "ASLV_KI_CN_bulk.txt"),
                    sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)

neg.ctrl <- vcn.2[vcn.2$Sample == "K562",]

#copy number per cell
data.to.plot <- vcn.2[!is.na(vcn.2$Variant),
                      c("Sample", "Virus", "Variant", "Population", "Droplets", "Genomes", "perCell")]

#Format Factors
data.to.plot$Virus <- factor(data.to.plot$Virus,
                             levels = vir.lvls)
data.to.plot$Variant <- factor(data.to.plot$Variant,
                               levels = var.lvls)
data.to.plot$Population <- factor(data.to.plot$Population, levels = c("GFP", "nS"))

#Reorder data by
data.to.plot <- data.to.plot[order(data.to.plot$Virus, data.to.plot$Population),]
data.to.plot <- data.to.plot[order(data.to.plot$Virus, data.to.plot$Variant),]

data.to.plot$Sample <- factor(data.to.plot$Sample,
                              levels = unique(data.to.plot$Sample))

#Set colors according to populations
cols <- c("darkgreen", "gray")
names(cols) <- c("GFP", "nS")

#Create plot
p.vcn <- ggplot(data = data.to.plot, aes(x = Sample, y = perCell, group = Virus, fill = Population)) +
  scale_fill_manual(values = alpha(cols, 0.5)) +
  geom_hline(yintercept = seq(10, 100, 10), linetype = 2, color = "gray") +
  geom_hline(yintercept = c(0), linetype = 1, color = "black") +
  geom_hline(yintercept = mean(neg.ctrl$perCell), linetype = "dotted", color = "black") +
  geom_bar(color = "black", position = "dodge", stat = "summary", fun = "mean") +
  geom_beeswarm(pch = 21, fill = NA, size = 2, cex = 2) +
  coord_cartesian(ylim = c(0, 100),
                  expand = TRUE) +
  scale_y_continuous(name = "Copies per 100 Genomes",
                     breaks = seq(0, 100, 10)
                     ) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black")
  ) +
  facet_grid(~ Variant, scales = "free_x", space = "free")

p.vcn

#Save figure
ggsave(paste0(fig.path, "ddPCR_VCN.png"), p.vcn,
       width = 120, height = 80, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

