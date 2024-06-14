#Create plot with ddPCR and FACS data
#copy number per cell
select.virus <- c("MoMLV", "FeLV", "SNV", "KoRV", "CrERV")
select.virus <- c("MoMLV", "FeLV", "SNV", "KoRV")


data.to.plot <- fcs.dd.df[fcs.dd.df$Virus %in% select.virus,]

#Format Factors
data.to.plot$Virus <- factor(data.to.plot$Virus,
                             levels = vir.lvls)
data.to.plot$Variant <- factor(data.to.plot$Variant,
                               levels = var.lvls)
data.to.plot$Experiment <- factor(data.to.plot$Experiment,
                               levels = dd.exp.name)

#Reorder data by 
data.to.plot <- data.to.plot[order(data.to.plot$Virus, data.to.plot$Variant),]
data.to.plot$vCN_per_pGFP <- data.to.plot$vCN_per_pGFP# + 1
data.to.plot$Sample <- factor(data.to.plot$Sample,
                              levels = unique(data.to.plot$Sample))

#Create Barplot  CN / GFP%

stat.df <- do.call(rbind,
                   lapply(unique(data.to.plot$Sample),
                          function(x) {
                            smpl.df <- data.to.plot[data.to.plot$Sample == x,]
                            sngl.line <- smpl.df[1, -c(4, 5, 6)]
                            sngl.line$Mean <- mean(smpl.df$vCN_per_pGFP)
                            sngl.line$SD <- sd(smpl.df$vCN_per_pGFP)
                            sngl.line
                          }))


p.dd.fcs <- ggplot(data = data.to.plot, aes(x = Sample, y = vCN_per_pGFP, fill = Variant, group = Virus)) +
  scale_fill_manual(values = alpha(cols, 0.5)) +
  scale_shape_manual(values = c(21, 22 ,23)) +
  geom_hline(yintercept = c(0.1), linetype = 1, color = "black") +
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), size = .25, linetype = 2, color = "gray") +
  geom_errorbar(data = stat.df, aes(ymin=Mean, ymax=Mean+SD, x = Sample, y = Mean, group = Variant),
                size = .15, width = .4, position=position_dodge(.9), inherit.aes = FALSE) +
  geom_bar(data = stat.df, aes(y = Mean), color = "black", size = .15, position = "dodge", stat = "summary", fun = "mean", width = 0.8) +
  geom_point(aes(shape = Experiment), size = .75, stroke = .15, position = position_dodge(width = .9)) +
  coord_cartesian(ylim = c(0.1, 30),
                  xlim = c(0.25, 3.75),
                  expand = FALSE) +
  scale_y_continuous(name = "GFP CN / GFP+ cell", trans = pseudo_log_trans(base = 10),
                     breaks = c(1, 2, 3, seq(4, 10, 2), seq(15, 30, 5))
  ) +
  theme_classic() +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme(
    axis.title.y = element_text(color = "black", face = "bold", size = 4),
    axis.title.x = element_text(color = "black", face = "bold", size = 4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "plain", size = 4, angle = 0),
    axis.line.x = element_blank(),
    axis.line.y = element_line(size = 0.25),
    axis.ticks.y =  element_line(size = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black", size = 5, angle = 0),
    panel.spacing.x = unit(2, "pt"),
    legend.justification = "top",
    legend.key.size = unit(4, "pt"),
    legend.spacing = unit(0, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4, face = "bold")
  ) +
  facet_wrap(~ Virus, nrow = 1, scales = "free_x")

#p.dd.fcs

ggsave(paste0(fig.path, "Figure3_C_FCS-CN.png"), p.dd.fcs,
       width = 70, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Create Barplot for active copies (GFP% / CN)
#exclude CrERV from active copies plot
select.virus <- c("MoMLV", "FeLV", "SNV", "KoRV")

data.to.plot <- fcs.dd.df[fcs.dd.df$Virus %in% select.virus,]

#Format Factors
data.to.plot$Virus <- factor(data.to.plot$Virus,
                             levels = vir.lvls)
data.to.plot$Variant <- factor(data.to.plot$Variant,
                               levels = var.lvls)
data.to.plot$Experiment <- factor(data.to.plot$Experiment,
                               levels = dd.exp.name)

#Reorder data by 
data.to.plot <- data.to.plot[order(data.to.plot$Virus, data.to.plot$Variant),]
data.to.plot$vCN_per_pGFP <- data.to.plot$vCN_per_pGFP# + 1
data.to.plot$Sample <- factor(data.to.plot$Sample,
                              levels = unique(data.to.plot$Sample))

stat.df <- do.call(rbind,
                   lapply(unique(data.to.plot$Sample),
                          function(x) {
                            smpl.df <- data.to.plot[data.to.plot$Sample == x,]
                            sngl.line <- smpl.df[1, -c(4, 5, 6)]
                            sngl.line$Mean <- mean(smpl.df$ActiveCop_ratio)
                            sngl.line$SD <- sd(smpl.df$ActiveCop_ratio)
                            sngl.line
                          }))

p.act.cps <- ggplot(data = data.to.plot, aes(x = Sample, y = ActiveCop_ratio, fill = Variant, group = Virus)) +
  scale_fill_manual(values = alpha(cols, 0.5)) +
  scale_shape_manual(values = c(21, 22 ,23)) +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  geom_hline(yintercept = seq(0.1, 1, 0.1), size = .25, linetype = 2, color = "gray") +
  geom_errorbar(data = stat.df, aes(ymin=Mean, ymax=Mean+SD, x = Sample, y = Mean, group = Variant),
                size = .15, width = .4, position=position_dodge(.9), inherit.aes = FALSE) +
  geom_bar(data = stat.df, aes(y = Mean), color = "black", size = .15, position = "dodge", stat = "summary", fun = "mean", width = 0.8) +
  geom_point(aes(shape = Experiment), size = .75, stroke = .15, position = position_dodge(width = .9)) +
  coord_cartesian(ylim = c(0, 1.05),
                  xlim = c(0.25, 3.75),
                  expand = FALSE) +
  scale_y_continuous(name = "Genome CN / %GFP+ cells",
                     breaks = seq(0.2, 1, 0.2)
  ) +
  theme_classic() +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme(
    axis.title.y = element_text(color = "black", face = "bold", size = 4),
    axis.title.x = element_text(color = "black", face = "bold", size = 4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "plain", size = 4, angle = 0),
    axis.line.x = element_blank(),
    axis.line.y = element_line(size = 0.25),
    axis.ticks.y =  element_line(size = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black", size = 5, angle = 0),
    panel.spacing.x = unit(2, "pt"),
    legend.justification = "top",
    legend.key.size = unit(4, "pt"),
    legend.spacing = unit(0, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4, face = "bold")
  ) +
  facet_wrap(~ Virus, nrow = 1, scales = "free_x")

#p.act.cps

ggsave(paste0(fig.path, "Figure3_D_ActCop.png"), p.act.cps,
       width = 60, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

