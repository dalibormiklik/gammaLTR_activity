#Select collor palette
select.pal <- cb.cols.inv


#Bulk expression
bulk.file <- "MLV_bulk.txt"
mlv.bulk.df <- read.table(paste0(data.path, bulk.file),
                          sep = "\t", dec = ".", header = TRUE)
mlv.bulk.df <- mlv.bulk.df[mlv.bulk.df$Vector != "CTRL",]

mlv.bulk.df$Variant <- factor(mlv.bulk.df$Varian, levels = var.lvls)
mlv.bulk.df$Stock_Volume <- factor(mlv.bulk.df$Stock_Volume, levels = sort(unique(mlv.bulk.df$Stock_Volume)))

#Create Sample column unique for each measured sample in time
mlv.bulk.df$Sample <- as.factor(paste(mlv.bulk.df$Vector,
                                      mlv.bulk.df$Variant,
                                      mlv.bulk.df$Stock_Volume,
                                      sep = "_"))

#Add column with relative %GFP to %GFP at 3 dpi
gfp.3dpi <- sapply(levels(mlv.bulk.df$Sample),
                   function(x) {
                     mlv.bulk.df$GFP_Perc[mlv.bulk.df$Sample == x & mlv.bulk.df$dpi == 3]
                   })
mlv.bulk.df$GFP_rPerc <- sapply(1:nrow(mlv.bulk.df),
                                function(wr) {
                                  smpl.wr <- mlv.bulk.df[wr,]
                                  smpl <- smpl.wr$Sample
                                  gfp.3 <- gfp.3dpi[names(gfp.3dpi) == smpl]
                                  gfp.x <- smpl.wr$GFP_Perc
                                  round(gfp.x / gfp.3, digits = 3)
                                })


#1) GFP expression in time
data.to.plot <- mlv.bulk.df
p.bulk <- ggplot(data = data.to.plot, aes(x = dpi, y = GFP_Perc,
                                          group = Sample)) +
  scale_fill_manual(values = alpha(cols, .5)) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_x_continuous(breaks = c(3, 10, 16, 30)) +
  scale_y_continuous(name = "% GFP+ cells") +
  geom_line(aes(color = Variant), size = 1) +
  geom_point(aes(fill = Variant, shape = Stock_Volume), size = 1, color = "black") +
  coord_cartesian(xlim = c(2, 32), ylim = c(0, 40), expand = FALSE) +
  guides(fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.title = element_text(color = "black", face = "bold", size = ax.tit.size),
        legend.justification = "top",
        legend.title = element_text(colour = "black", face = "bold", size = ax.tit.size),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(4, "pt")
  )

#p.bulk

ggsave(paste0(fig.path, "Figure1_S1_B_FACS_3-30dpi.png"), p.bulk,
       width = 50, height = 30, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


#3 fold change of GFP expression in time
data.to.plot <- mlv.bulk.df

p.bulk.rel <- ggplot(data = data.to.plot, aes(x = dpi, y = GFP_rPerc, color = Variant, fill = Variant, group = Variant)) +
  scale_fill_manual(values = alpha(select.pal, 1)) +
  scale_color_manual(values = select.pal) +
  scale_x_continuous(breaks = c(3, 10, 16, 30)) +
  scale_y_continuous(name = "% GFP / % GFP at 3 dpi", trans = "log2",
                     breaks = c(0.5, 1, 2),
                     minor_breaks = c(seq(0.5, 0.9, 0.1), seq(1.1, 1.9, 0.1))) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray") +
  geom_line(aes(group = Stock_Volume), size = .25) +
  geom_point(pch = 21, size = .8, stroke = 0.05) +
  stat_summary(fun = mean, geom="line", size = 0.3, color = "black") +
  stat_summary(fun = mean, geom="point", size = 1, shape = 21, color = "black") +
  coord_cartesian(xlim = c(0, 32),
                  ylim = c(0.5, 2),
                  expand = FALSE) +
  guides(color = "none",
         fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black" ,face = "bold", size = ax.tex.size),
        axis.title = element_text(color = "black",face = "bold", size = ax.tit.size),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = 5),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size + 1),
        legend.key.size = unit(3, "pt")
        ) +
  facet_wrap(~Variant, nrow = 1)

#p.bulk.rel

ggsave(paste0(fig.path, "Figure1_C_FACS_3-30dpi.png"), p.bulk.rel,
       width = 50, height = 30, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
