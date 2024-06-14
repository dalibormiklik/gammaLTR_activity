#Supplementary figures to Figure3
#01
#Create plot with % of GFP+ cells
data.to.plot <- fcs.tab[fcs.tab$Virus %in% vir.lvls,]

data.to.plot <- data.to.plot[-grep("_[25]00_",data.to.plot$Sample),]

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = vir.lvls)
data.to.plot$Variant <- factor(data.to.plot$Variant, levels = var.lvls)
data.to.plot$Experiment <- factor(data.to.plot$Experiment, levels = unique(data.to.plot$Experiment))

p.perc.time <- ggplot(data = data.to.plot, aes(x = dpi, y = GFP_Perc, color = Variant, fill = Variant, shape = Variant, group = Variant)) +
  scale_fill_manual(values = alpha(cb.cols.inv, 0.5)) +
  scale_color_manual(values = alpha(cb.cols.inv, 0.75)) +
  scale_shape_manual(values = c(21:23)) +
  scale_x_continuous(breaks = sort(unique(data.to.plot$dpi))) +
  geom_hline(yintercept = 0) +
  geom_line(size = 1.25) +
  geom_point(fill = alpha("gray90", 0.5), color = "black", size = 1.5, stroke = 1) +
  ylab("GFP+ [%]") +
  theme_classic() +
  coord_cartesian(expand = TRUE) +
  theme(axis.text = element_text(color = "black", size = 6, face = "bold"),
        axis.line.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        strip.text.y = element_text(angle = 0)
        ) +
  facet_wrap(Virus ~ Experiment, scales = "free_y")

#p.perc.time

#Save figure
ggsave(paste0(fig.path, "Figure3_Supp_GFP_perc.png"), p.perc.time,
       width = 216, height = 160, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#02
#Plot GFP intensity in time
in.cat <- c("wt", "W390A", "CBX")

cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat


relfi.df <- do.call(rbind,
                 lapply(unique(fcs.tab$Experiment),
                        function(expmnt) {
                          e.df <- fcs.tab[fcs.tab$Experiment == expmnt,]
                          do.call(rbind,
                                  lapply(unique(e.df$dpi),
                                         function(d) {
                                           d.df <- e.df[e.df$dpi == d,]
                                           v.df <- d.df[d.df$Virus != "Mock" & !is.na(d.df$Virus),]
                                           n.ctrl <- mean(d.df$GFP_p50[d.df$Virus == "Mock" & !is.na(d.df$Virus)])
                                           v.df$GFP_p50 <- v.df$GFP_p50 - n.ctrl
                                           v.df$GFP_p25 <- v.df$GFP_p25 - n.ctrl
                                           v.df$GFP_p75 <- v.df$GFP_p75 - n.ctrl
                                           v.df
                                         }))
                        }))

p.mfi.t <- ggplot(data = relfi.df, aes(x = dpi, y = log10(GFP_p50), fill = Variant)) +
  scale_y_continuous(name = bquote(log[10]~MFI)) +
  scale_fill_manual(values = alpha(cols, 0.5)) +
  scale_shape_manual(values = c(21, 22 ,23, 24)) +
  geom_hline(yintercept = seq(0,4,.5), color = "gray", size = .25) +
  geom_smooth(aes(color = Variant), color = "black", size = .15, se = FALSE) +
  geom_point(aes(shape = Experiment), size = 1, stroke = .15) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0,17), expand = FALSE) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme_classic() +
  theme(
    axis.title.y = element_text(color = "black", face = "bold", size = 4),
    axis.title.x = element_text(color = "black", face = "bold", size = 4),
    axis.text.x = element_text(color = "black", face = "plain", size = 4, angle = 0),
    axis.ticks.x = element_line(size = 0.25),
    axis.text.y = element_text(color = "black", face = "plain", size = 4, angle = 0),
    axis.line.x = element_line(size = 0.25),
    axis.line.y = element_line(size = 0.25),
    axis.ticks.y =  element_line(size = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black", size = 5, angle = 0),
    panel.spacing.x = unit(2, "pt"),
    panel.spacing.y = unit(10, "pt"),
    legend.justification = "top",
    legend.key.size = unit(4, "pt"),
    legend.spacing = unit(0, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4, face = "bold")
  ) +
  facet_grid(factor(Variant, levels = in.cat) ~ factor(Virus, levels = vir.lvls))

ggsave(paste0(fig.path, "Figure3_0_MFI-time.png"), p.mfi.t,
       width = 120, height = 90, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#03
#Supplementary Figure with ddPCR GFP/cell values
p.dd.all <- ggplot(vcn.2, aes(x = Experiment, y = perCell, fill = factor(Variant, levels = var.lvls), group = factor(Variant, levels = var.lvls))) +
  scale_fill_manual(values = alpha(cols, c(.5))) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  geom_bar(position = position_dodge2(preserve = "single"), stat = "summary", fun = "mean", color = "black") +
  geom_point(shape = 21, position = position_dodge(.9)) +
  ylab("GFP copies / 2(RPP30)") +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black" ,face = "plain", size = ax.tex.size + 1),
        axis.text.x = element_text(color = "black" ,face = "plain", size = ax.tex.size + 1, angle = 45),
        axis.title = element_text(color = "black" , face = "bold", size = ax.tit.size + 1),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = ax.tit.size + 1),
        legend.justification = "top",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(3, "pt")
  ) +
  facet_wrap(~ factor(Virus, levels = c(vir.lvls, c("KM_4E4", "K562"))), scales = "free")

ggsave(paste0(fig.path, "Figure3_Supp_ddPCR.png"), p.dd.all,
       width = 120, height = 90, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
