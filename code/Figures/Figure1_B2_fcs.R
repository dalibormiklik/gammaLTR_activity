#Create plots for Figure 1B + associated supplementary figures
#FACS dot plots

#biexp scale transform function
#Biexp transformation function
biexp_trans <- function(lim = 100, decade.size = lim){
  trans <- function(x){
    ifelse(x <= lim,
           x,
           lim + decade.size * (suppressWarnings(log(x, 10)) -
                                  log(lim, 10)))
  }
  inv <- function(x) {
    ifelse(x <= lim,
           x,
           10^(((x-lim)/decade.size) + log(lim,10)))
  }
  breaks <- function(x) {
    if (all(x <= lim)) {
      pretty_breaks()(x)
    } else if (all(x > lim)) {
      log_breaks(10)(x)
    } else {
      #unique(c(pretty_breaks()(c(x[1],lim)),
      #         log_breaks(10)(c(lim, x[2]))))
      #10^pretty(log10(100:max(x)))
      10^sort(unique(round(log10(100:max(x)))))
    }
  }
  labels <- label_log(base = 10, digits = 1)
  trans_new(paste0("biexp-",format(lim)), trans, inv, breaks)
} 

#-
#Main text figure: Fig 1 B
#-
sel.gid <- c("wt_50", "W390A_50", "CBX_50")
data.to.plot <- do.call(rbind,
                        lapply(sel.gid,
                               function(x) {
                                 df <- sel.export.df[sel.export.df$GID == x,]
                                 df[sample.int(n = nrow(df), size = min.cells),]
                               }))

data.to.plot$Volume <- factor(data.to.plot$Volume, levels = unique(data.to.plot$Volume))

pos.export.df <- data.to.plot[data.to.plot$Population == pos.value,]


sum.df <- export.sum.df[export.sum.df$GID %in% sel.gid,]

p01.b <-  ggplot(data = data.to.plot, aes(x = FSC.A, y = GFP..488B..A, colour = Population)) +
  scale_fill_manual(values = alpha(cols, c(.5, .5))) +
  scale_colour_manual(values = alpha(cols, c(.5, .5))) +
  scale_x_continuous(name = x.name, trans = log10_trans(), labels = label_log(base = 10, digits = 1)) +
  scale_y_continuous(name = y.name, trans = biexp_trans(lim = 700),
                     breaks = c(300, 500, 1000), labels = c(expression("3"%*%"10"^"2"), expression("5"%*%"10"^"2"), expression("1"%*%"10"^"3"))) +
  geom_point(shape = 1, size = 0.001) +
  geom_boxplot(data = pos.export.df,
               colour = alpha("black", .8), outlier.shape = NA, fill = NA, na.rm = TRUE, inherit.aes = TRUE) +
  geom_text(data = sum.df, aes(label = round(100 * Ratio), 1), inherit.aes = TRUE,
            x = 2.3, y = 8.1*10^2, colour = "black", size = 1.5, fontface = "bold") +
  coord_cartesian(ylim = c(250, 1*10^3),
                  #xlim = c(1, 10^5)
  ) +
  guides(color = "none"
         #colour = guide_legend(override.aes = list(size=1))
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black" ,face = "plain", size = ax.tex.size + 1),
        axis.text.x = element_blank(),
        axis.title = element_text(color = "black" , face = "bold", size = ax.tit.size),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = ax.tit.size + 1),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size - 1),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(3, "pt")) +
  facet_grid(~ Volume + factor(Variant, levels = var.lvls))


p01.b

ggsave(paste0(fig.path, "Figure1_B_FACS_3dpi.png"), p01.b,
       width = 40, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


#Supplementary figure with all "titration" measurements
data.to.plot <- do.call(rbind,
                        lapply(unique(sel.export.df$GID),
                               function(x) {
                                 df <- sel.export.df[sel.export.df$GID == x,]
                                 df[sample.int(n = nrow(df), size = min.cells),]
                               }))

pos.export.df <- sel.export.df[sel.export.df$Population == pos.value,]
sum.df <- export.sum.df

p.01.S1A <- ggplot(data = data.to.plot, aes(x = FSC.A, y = GFP..488B..A, colour = Population)) +
  scale_fill_manual(values = alpha(cols, c(.5, .5))) +
  scale_colour_manual(values = alpha(cols, c(.5, .5))) +
  scale_x_continuous(name = x.name, trans = log10_trans(), labels = label_log(base = 10, digits = 1)) +
  scale_y_continuous(name = y.name, trans = biexp_trans(lim = 700),
                     breaks = c(300, 500, 1000), labels = c(expression("3"%*%"10"^"2"), expression("5"%*%"10"^"2"), expression("1"%*%"10"^"3"))) +
  geom_point(shape = 1, size = 0.001) +
  geom_boxplot(data = pos.export.df,
               colour = alpha("black", .4), outlier.shape = NA, fill = NA, na.rm = TRUE, inherit.aes = TRUE) +
  geom_text(data = sum.df, aes(label = round(100 * Ratio), 1), inherit.aes = TRUE,
            x = 2.3, y = 8.1*10^2, colour = "black", size = 1.5, fontface = "bold") +
  coord_cartesian(ylim = c(250, 1*10^3),
                  #xlim = c(1, 10^5)
  ) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black" ,face = "plain", size = ax.tex.size),
        axis.text.x = element_blank(),
        axis.title = element_text(color = "black" , face = "bold", size = ax.tit.size),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = ax.tit.size),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size - 1),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(3, "pt")) +
  facet_grid(~ Volume + factor(Variant, levels = var.lvls))


p.01.S1A
ggsave(paste0(fig.path, "Figure1_S1_A_FACS_3dpi.png"), p.01.S1A,
       width = 120, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

