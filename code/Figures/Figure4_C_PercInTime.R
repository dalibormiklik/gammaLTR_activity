#Create plots of AS.gamma GFP expression in time

vir.lvls <- c("SFFV","MoMLV", "FeLV", "SNV", "KoRV", "CrERV")
virs <- c("AS-SF","AS-Mo", "AS-Fe", "AS-SN", "AS-Ko", "AS-Cr")

data.df <- read.table(paste0(data.path, "GFP_time.txt"),
    sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE
    )

#Define time points for plot
u.dpi <- sort(unique(data.df$dpi))
max.dpi <- max(u.dpi)
plot.dpi <- u.dpi[u.dpi <= max.dpi]

#Set color palette
cb.cols <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "#C22DF3", "#B9390C")
names(cb.cols) <- vir.lvls

#Size of axis font
ax.tex.size <- 3
ax.tit.size <- 4


#plot
vir.lvls2 <- vir.lvls
data.to.plot <- data.df[data.df$virus %in% vir.lvls2 &
                          data.df$dpi %in% plot.dpi,]


p.bulk <- ggplot(data = data.to.plot, aes(x = dpi, y = gfp_perc,
                                          group = ID)) +
  scale_fill_manual(values = alpha(cb.cols, .5)) +
  scale_color_manual(values = alpha(cb.cols, .5)) +
  scale_x_continuous(breaks = c(0, plot.dpi)) +
  scale_y_continuous(name = "% GFP+ cells", trans = "log10", breaks = c(1:5, 10, 20, 30, 50)) +
  geom_line(aes(color = virus), linewidth = 0.5) +
  geom_point(aes(fill = virus), shape = 21, size = 1, color = "black") +
  coord_cartesian(xlim = c(0, (max.dpi+max.dpi*0.1)),
                  ylim = c(0.4, 50),
                  expand = FALSE) +
  guides(fill = "none",
         color = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.title = element_text(color = "black", face = "bold", size = ax.tit.size),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = ax.tit.size+1),
        legend.justification = "top",
        legend.title = element_text(colour = "black", face = "bold", size = ax.tit.size),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(4, "pt")
  ) +
  facet_grid(~ factor(virus, levels = vir.lvls2))

#p.bulk

ggsave(paste0(fig.path,"Figure4_Supp01_AS-gamma_time.png"), p.bulk,
       width = 120, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


#Fold change in expression
vir.lvls2 <- vir.lvls

data.to.plot <- data.df[data.df$virus %in% vir.lvls2 &
                          data.df$dpi %in% plot.dpi,]

p.bulk.rel <- ggplot(data = data.to.plot, aes(x = as.numeric(dpi), y = GFP_rPerc, fill = virus, color = virus, group = virus)) +
  scale_fill_manual(values = alpha(cb.cols, .9)) +
  scale_color_manual(values = alpha(cb.cols, .5)) +
  scale_x_continuous(name = "dpi",
                     breaks = c(0, plot.dpi)) +
  scale_y_continuous(name = "% GFP+ / % GFP+ at 3 dpi", trans = "log2",
                     breaks = c(seq(0.1, 0.7, 0.1), seq(0.8, 2, 0.2)),
                     minor_breaks = c(seq(0.5, 0.9, 0.1), seq(1.1, 1.9, 0.1))) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray") +
  geom_line(aes(group = volume), size = .25) +
  geom_point(pch = 21, size = .8, stroke = 0.05) +
  stat_summary(fun = mean, geom="line", size = 0.3, color = "black") +
  stat_summary(fun = mean, geom="point", size = 1, shape = 21, color = "black") +
  coord_cartesian(xlim = c(0, (max.dpi+max.dpi*0.1)),
                  ylim = c(0.28, 2),
                  expand = FALSE) +
  guides(fill = "none",
         color = "none") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", face = "bold", size = ax.tex.size),
        axis.title = element_text(color = "black", face = "bold", size = ax.tit.size),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = ax.tit.size + 2),
        legend.justification = "top",
        legend.title = element_text(colour = "black", face = "bold", size = ax.tit.size),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(4, "pt")) +
        facet_grid(~ factor(virus, levels = vir.lvls2))

#p.bulk.rel

ggsave(paste0(fig.path,"Figure4_C_AS-gamma_time_fold.png"), p.bulk.rel,
       width = 120, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
