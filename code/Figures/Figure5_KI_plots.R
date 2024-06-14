#Figure 5 panels

#set treshold for single copy number evaluation
scn.val <- 1.5

#Create subgroup data.frames
pcr.conf <- sum.df[sum.df$PCR == "++" & sum.df$Seq == "++" & !is.na(sum.df$PCR),]
pcr.unconf <- sum.df[sum.df$PCR != "++" & !is.na(sum.df$PCR),]

#Set variables for graphics
p.size <- 2
cols <- c("#09AB03", "#D81B60", "#1E88E5", "#FFC107")
bg.cols <- c(rep(c("gray80", NA), 2))

names(cols) <- group.tab$Virus

#Create plot(s)

#Figure5_C: % positive
data.to.plot <- sum.df
p.perc <- ggplot(data = data.to.plot, aes(x = Virus, y = mCherry_Freq, fill = Virus)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = alpha(cols, 0.5)) +
  geom_violin(trim = TRUE, scale = "width", adjust = 0.5, draw_quantiles = 0.5, show.legend = FALSE) +
  geom_quasirandom(shape = 21, size = p.size, stroke = .2) +
  coord_cartesian(ylim = c(0, 100)) +
  ylab("% mcherry+") +
  guides(color = "none") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "black")
  ) +
  facet_grid(~ TargetSite)

#p.perc

ggsave(paste0(fig.path, "Figure5_C_FACSperc.png"), p.perc,
       width = 200, height = 45, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Supplementary Figure: % GFP in confirmed clones
#Select only PCR-confirmed clones
## % positive
data.to.plot <- pcr.conf
p.perc.pcr <- ggplot(data = data.to.plot, aes(x = Virus, y = mCherry_Freq, fill = Virus)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = alpha(cols, 0.75)) +
  geom_quasirandom(data = pcr.unconf,
                   fill = alpha("grey", 0.1), shape = 21, size = p.size, stroke = .2) +
  geom_quasirandom(shape = 21, size = p.size, stroke = .2) +
  coord_cartesian(ylim = c(0, 100)) +
  ylab("% mcherry+") +
  guides(color = "none") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "black")
  ) +
  facet_grid(~ TargetSite)

p.perc.pcr

ggsave(paste0(fig.path, "Figure5_Supp_FACS-PCR.png"), p.perc.pcr,
       width = a4w, height = a4h * 1/6, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Figure_5_D: CN by TargetSite
data.to.plot <- pcr.conf

#Count of clones in each virus-site category
clone.count <- do.call(rbind,
                       lapply(smpl.is.names,
                              function(s) {
                                s.df <- data.to.plot[data.to.plot$TargetSite == s,]
                                data.frame(TargetSite = s,
                                           #mCherry_Freq = 100,
                                           clone_count = nrow(s.df)
                                )
                              }))
clone.count$TargetSite <- factor(clone.count$TargetSite, levels = smpl.is.names)


max.cn <- 14
y.breaks <- c(0.5, 1:max.cn)
y.breaks.select <- c(0.5, 1:4, c(6, 10, 14))
y.labs <- sapply(y.breaks, function(x) {ifelse(x %in% y.breaks.select, x, "")})

p.cn.ts <- ggplot(data = data.to.plot, aes(x = "X", y = CN, group = TargetSite)) +
  geom_hline(yintercept = y.breaks.select, linetype = 2, linewidth = .25, colour = "gray") +
  scale_fill_manual(values = alpha(cols, 0.75)) +
  scale_y_continuous(breaks = c(0.5, 1:max.cn), labels = y.labs, trans = "log2") +
  geom_text(data = clone.count, aes(label=clone_count, y = (max.cn + 2)), size = 3) +
  geom_violin(trim = TRUE, scale = "width", adjust = 0.5, draw_quantiles = 0.5,
              fill = alpha("gray", .1), colour = alpha("black", .4), show.legend = FALSE) +
  geom_beeswarm(aes(fill = Virus), shape = 21, size = 1.4, stroke = .2, cex = 5) +
  ylab("insert copies per cell genome") +
  xlab("target site") +
  guides(color = "none") +
  coord_cartesian(ylim = c(0.5, max.cn+4),
                  xlim = c(0.5, 1.5),
                  expand = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black")
  ) + facet_grid(~ TargetSite)

p.cn.ts

ggsave(paste0(fig.path ,"Figure5_D_FACS-PCR-CN-byTarget.png"), p.cn.ts,
       width = 110, height = 50, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Figure_5_E:
##Plot single knock-ins
pcr.conf2 <- pcr.conf[!is.na(pcr.conf$CN),]
pcr.conf2$CN_group <- sapply(pcr.conf2$CN,
                             function(cn) {
                               if(cn >= 2.5) {
                                 "CN_2+"
                               } else {
                                 if(cn >= 1.5) {
                                   "CN_2"
                                 } else {
                                   "CN_1"
                                 }
                               }
                             })
pcr.conf2$KI_group <- sapply(1:nrow(pcr.conf2),
                             function(w.r) {
                               r.smpl <- pcr.conf2[w.r,]
                               smpl.name <- r.smpl$TargetSite
                               smpl.cn <- r.smpl$CN_group
                               smpl.ta <- r.smpl$target_allele
                               
                               if(smpl.name == "IS11") {
                                 if(smpl.cn == "CN_1") {
                                   "single-allele KI"
                                 } else {
                                   if(smpl.cn == "CN_2") {
                                     if(smpl.ta %in% c("0","1")) {
                                       "double-allele KI"
                                     } else {
                                       "multi"
                                     }
                                   } else {
                                     "multi"
                                   }
                                 }
                               } else {
                                 if(smpl.cn == "CN_1") {
                                   "single-allele KI"
                                 } else {
                                   if(smpl.cn == "CN_2") {
                                     if(smpl.ta == "0") {
                                       "double-allele KI"
                                     } else {
                                       "multi"
                                     }
                                   } else {
                                     "multi"
                                   }
                                 }
                               }
                             })

data.to.plot <- pcr.conf2[pcr.conf2$KI_group == "single-allele KI",]

#data.frame for the backgroud rectangles (columns)
l <- nrow(group.tab)
bg.rects <- do.call(rbind,
                    lapply(smpl.is.names,
                           function(s) {
                             data.frame(TargetSite = factor(rep(s, l), levels = smpl.is.names),
                                        Virus = factor(group.tab$Virus, levels = group.tab$Virus),
                                        mCherry_Freq = rep(100, l),
                                        fill = bg.cols)
                           }))
#Count of clones in each virus-site category
clone.count <- do.call(rbind,
                       lapply(smpl.is.names,
                              function(s) {
                                s.df <- data.to.plot[data.to.plot$TargetSite == s,]
                                do.call(rbind,
                                        lapply(group.tab$Virus,
                                               function(v) {
                                                 data.frame(TargetSite = s,
                                                            Virus = v,
                                                            mCherry_Freq = 100,
                                                            clone_count = nrow(s.df[s.df$Virus == v,]))
                                                 
                                               }))
                              }))
clone.count$TargetSite <- factor(clone.count$TargetSite, levels = smpl.is.names)
clone.count$Virus <- factor(clone.count$Virus, levels = group.tab$Virus)

#Create plot
p.kig <- ggplot(data = data.to.plot, aes(x = Virus, y = mCherry_Freq, group = TargetSite)) +
  scale_fill_manual(values = alpha(cols, 0.75)) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  geom_col(data = bg.rects, color = NA, fill = alpha(bg.rects$fill, .75)) +
  geom_text(data = clone.count, aes(label=clone_count, y = mCherry_Freq + 6), size = 3) +
  geom_hline(yintercept = c(0, 100), linetype = 1, colour = alpha("black", 0.9), linewidth = .25) +
  geom_beeswarm(aes(fill = Virus), cex = 3, shape = 21, size = 1.7, stroke = .2) +
  ylab("% mcherry+") +
  xlab("target site") +
  guides(fill = "none",
         color = "none") +
  coord_cartesian(ylim = c(-3, 109) ,expand = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 5, face = "bold", angle = 45, hjust = 1.2, vjust = 1.4),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 8)
  ) + facet_grid(~ TargetSite)

p.kig

ggsave(paste0(fig.path, "Figure5_E_singleKI.png"), p.kig,
       width = 90, height = 54, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

##Supplementary Figures:

#Select only PCR-confirmed clones
## % positive
data.to.plot <- pcr.conf
p.perc.pcr <- ggplot(data = data.to.plot, aes(x = Virus, y = mCherry_Freq, fill = Virus)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = alpha(cols, 0.75)) +
  geom_quasirandom(data = pcr.unconf,
                   fill = alpha("grey", 0.1), shape = 21, size = p.size, stroke = .2) +
  geom_quasirandom(shape = 21, size = p.size, stroke = .2) +
  coord_cartesian(ylim = c(0, 100)) +
  ylab("% mcherry+") +
  guides(color = "none") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "black")
  ) +
  facet_grid(~ TargetSite)

p.perc.pcr

ggsave(paste0(fig.path ,"Figure5_Supp_FACS-PCR.png"), p.perc.pcr,
       width = 180, height = 45, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


##clones with CN == 2
data.to.plot <- pcr.conf2[pcr.conf2$KI_group == "double-allele KI",]

#data.frame for the backgroud rectangles (columns)
l <- nrow(group.tab)
bg.rects <- do.call(rbind,
                    lapply(smpl.is.names,
                           function(s) {
                             data.frame(TargetSite = factor(rep(s, l), levels = smpl.is.names),
                                        Virus = factor(group.tab$Virus, levels = group.tab$Virus),
                                        mCherry_Freq = rep(100, l),
                                        fill = bg.cols)
                           }))
#Count of clones in each virus-site category
clone.count <- do.call(rbind,
                       lapply(smpl.is.names,
                              function(s) {
                                s.df <- data.to.plot[data.to.plot$TargetSite == s,]
                                do.call(rbind,
                                        lapply(group.tab$Virus,
                                               function(v) {
                                                 data.frame(TargetSite = s,
                                                            Virus = v,
                                                            mCherry_Freq = 100,
                                                            clone_count = nrow(s.df[s.df$Virus == v,]))
                                                 
                                               }))
                              }))
clone.count$TargetSite <- factor(clone.count$TargetSite, levels = smpl.is.names)
clone.count$Virus <- factor(clone.count$Virus, levels = group.tab$Virus)

#Create plot
p.kig.d <- ggplot(data = data.to.plot, aes(x = Virus, y = mCherry_Freq, group = TargetSite)) +
  scale_fill_manual(values = alpha(cols, 0.75)) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  geom_col(data = bg.rects, color = NA, fill = alpha(bg.rects$fill, .75)) +
  geom_text(data = clone.count, aes(label=clone_count, y = mCherry_Freq + 6), size = 3) +
  geom_hline(yintercept = c(0, 100), linetype = 1, colour = alpha("black", 0.9), linewidth = .25) +
  geom_beeswarm(aes(fill = Virus), cex = 3, shape = 21, size = 1.7, stroke = .2) +
  ylab("% mcherry+") +
  xlab("target site") +
  guides(fill = "none",
         color = "none") +
  coord_cartesian(ylim = c(-3, 109) ,expand = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 5, face = "bold", angle = 45, hjust = 1.2, vjust = 1.4),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 8)
  ) + facet_grid(~ TargetSite)

#p.kig.d

ggsave(paste0(fig.path, "Figure5_Supp_doubleKI.png"), p.kig.d,
       width = 90, height = 55, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

