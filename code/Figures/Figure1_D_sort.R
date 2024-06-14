#Set point stroke
pt.stroke <- 0.1

#Load data
sort.file <- "MLV_sorted.txt"
mlv.sort.df <- read.table(paste0(data.path, sort.file),
                          sep = "\t", dec = ".", header = TRUE)
mlv.sort.df$Variant <- factor(mlv.sort.df$Varian, levels = var.lvls)

#Create Plots
data.to.plot <- mlv.sort.df[mlv.sort.df$dpi > 3,]
p.sort <- ggplot(data = data.to.plot, aes(x = dpi, y = GFP_Perc,
                                          color = Variant, fill = Variant)) +
  scale_fill_manual(values = alpha(select.pal, .9)) +
  scale_color_manual(values = alpha(select.pal, .5)) +
  scale_x_continuous(breaks = c(14, 23)) +
  scale_y_continuous(name = "% GFP+ cells") +
  geom_line(size = .5) +
  geom_point(shape = 21, size = 1, color = "black", stroke = pt.stroke) +
  coord_cartesian(xlim = c(12, 25),
                  ylim = c(0, 100),
                  expand = FALSE) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.title = element_text(color = "black", face = "bold", size = ax.tit.size),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size + 1),
        legend.key.size = unit(3, "pt"))


#p.sort

ggsave(paste0(fig.path, "Figure1_D_FACS_3dpi.png"), p.sort,
       width = 40, height = 25, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
