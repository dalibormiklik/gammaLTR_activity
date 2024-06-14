
#Load data
clone.file <- "MLV_clones.txt"
mlv.clone.df <- read.table(paste0(data.path, clone.file),
                           sep = "\t", dec = ".", header = TRUE)

#Modify data.frame
mlv.clone.df$Variant <- factor(mlv.clone.df$Varian, levels = var.lvls)
mlv.clone.count <- do.call(rbind,
                           lapply(levels(mlv.clone.df$Variant),
                                  function(x) {
                                    cl.count <- nrow(mlv.clone.df[mlv.clone.df$Variant == x,])
                                    data.frame(Variant = x,
                                               n = cl.count)
                                  }))
#Create Plots
pt.stroke <- 0.1

#Set data.to.plot for % of positives plot
data.to.plot <- mlv.clone.df[mlv.clone.df$Vector != "CTRL",]
p.clone.perc <- ggplot(data = data.to.plot, aes(x = Variant, y = GFP_Perc,
                                                color = Variant, fill = Variant)) +
  scale_fill_manual(values = alpha(select.pal, .25)) +
  scale_color_manual(values = select.pal) +
  scale_y_continuous(name = "% GFP+ cells", breaks = seq(0, 100, 25)) +
  geom_hline(yintercept = 90, linetype = 2, color = "gray") +
  geom_quasirandom(shape = 21, size = 1, color = alpha("black", .75), stroke = pt.stroke, groupOnX = TRUE) +
  coord_cartesian(ylim = c(0, 110), xlim = c(0.5, 3.5), expand = FALSE) +
  guides(fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = ax.tit.size))

#p.clone.perc

ggsave(paste0(fig.path, "Figure1_E_clones_perc.png"), p.clone.perc,
       width = 35, height = 25, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


#Set data.to.plot for MFI plot
data.to.plot <- mlv.clone.df[mlv.clone.df$Vector != "CTRL",]
p.clone.mfi <- ggplot(data = data.to.plot, aes(x = Variant, y = GFP_MFI,
                                               color = Variant, fill = Variant)) +
  scale_fill_manual(values = alpha(select.pal, .25)) +
  scale_color_manual(values = select.pal) +
  scale_y_log10(name = "MFI", labels = label_log()) +
  geom_quasirandom(shape = 21, size = 1, stroke = pt.stroke, color = alpha("black", .75), groupOnX = TRUE) +
  geom_violin(color = "black", linewidth = .3, draw_quantiles = 0.5) +
  coord_cartesian(ylim = c(70, 1.75*10^5),
                  xlim = c(0.5, 3.5),
                  expand = FALSE) +
  guides(fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = ax.tit.size),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size + 1),
        legend.key.size = unit(3, "pt"))

#p.clone.mfi

ggsave(paste0(fig.path, "Figure1_F_clones_MFI.png"), p.clone.mfi,
       width = 35, height = 25, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
