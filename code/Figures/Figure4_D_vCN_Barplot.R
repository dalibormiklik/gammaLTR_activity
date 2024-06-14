vcn.2 <- read.table(paste0(fig.data, "vCN_all.txt"),
                    sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)

#copy number per cell
data.to.plot <- vcn.2[!is.na(vcn.2$Virus),
                      c("Sample", "Virus", "Volume", "dpi", "Droplets", "Genomes", "perCell", "DropletQ")]

#Format Factors
#data.to.plot <- data.to.plot[order(match(data.to.plot$Variant, var.lvls)),]

data.to.plot$Virus <- factor(data.to.plot$Virus,
                             levels = vir.lvls)
data.to.plot$Volume <- factor(data.to.plot$Volume,
                               levels = var.lvls)
#Reorder data by 
data.to.plot <- data.to.plot[!is.na(data.to.plot$Virus) &
                               order(data.to.plot$Virus, data.to.plot$Volume),]

data.to.plot$Sample <- factor(data.to.plot$Sample,
                              levels = unique(data.to.plot$Sample))

#Set plot colors
cb.cols <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "#B9390C", "#C22DF3")
names(cb.cols) <- vir.lvls


ax.tex.size <- 4
ax.tit.size <- 4

p.gfp <- ggplot(data = data.to.plot, aes(x = Virus, y = GFP, fill = Virus)) +
  scale_fill_manual(values = alpha(cb.cols.inv, 0.5)) +
  scale_color_manual(values = cols.dropq) +
  geom_hline(yintercept = c(0), linetype = 1, color = "black") +
  geom_hline(yintercept = c(5, 10, 15), linetype = 3, color = "gray", size = .15) + 
  geom_bar(color = "black", position = "dodge", stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = DropletQ), pch = 21, fill = NA, size = 2) +
  coord_cartesian(ylim = c(0,20),
                  expand = TRUE) +
  scale_y_continuous(name = "Droplet Count (x10^3)",
                     ) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 1),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black")
  ) +
  facet_grid(~ Virus, scales = "free_x", space = "free")

p.cells <- ggplot(data = data.to.plot, aes(x = Volume, y = Genomes / 100, fill = Virus)) +
  #scale_fill_manual(values = alpha(cb.cols.inv, 0.5)) +
  scale_color_manual(values = cols.dropq) +
  geom_hline(yintercept = c(0), linetype = 1, color = "black") +
  geom_bar(color = "black", position = "dodge", stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = DropletQ), pch = 21, fill = NA, size = 2) +
  #coord_cartesian(ylim = c(0,5),
                  #xlim = c(0.25, 3.75),
  #                expand = TRUE) +
  scale_y_continuous(name = "Diploid Genomes (x10^2)",
                     #breaks = c(0, 0.25 , 0.5, 0.75, 1)
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 1),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black")
  ) +
  facet_grid(~ Virus, scales = "free_x", space = "free")


p.vcn <- ggplot(data = data.to.plot, aes(x = Virus, y = perCell, group = Virus, fill = Virus)) +
  scale_fill_manual(values = alpha(cb.cols, 0.5)) +
  scale_color_manual(values = cols.dropq) +
  geom_hline(yintercept = seq(5, 100, 5), linetype = 3, color = "gray", size = .15) +
  geom_hline(yintercept = c(0), linetype = 1, color = "black") +
  geom_bar(data = data.to.plot[data.to.plot$DropletQ == "OK" & data.to.plot$perCell > 0,],
           color = "black", position = "dodge", stat = "summary", fun = "mean", size = 0.25, width = .8) +
  geom_beeswarm(aes(color = DropletQ), cex = 2, pch = 21, fill = NA, size = 0.75, stroke = 0.25) +
  scale_y_continuous(name = "GFP copies / 200 RPP30 copies",
                     breaks = seq(0,100, 2)
                     ) +
  coord_cartesian(ylim = c(0, 10.5),
                  xlim = c(0.35, 6.65),
                  expand = FALSE) +
  guides(fill = "none",
         color = "none") +
  xlab("LTR U3") +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(color = "black", face = "bold", angle = 45, vjust = 1, hjust = 1, size = ax.tex.size + 2),
    axis.text.y = element_text(color = "black", size = ax.tex.size),
    axis.title.y = element_text(color = "black", face = "bold", size = ax.tit.size),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black")
  )

p.vcn

#Save figure
#A4 = 210 / 297 mm
a4w <- 180
a4h <- 267
ggsave(paste0(fig.save, "Figure4_D_vCN.png"), p.vcn,
       width = 27, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
