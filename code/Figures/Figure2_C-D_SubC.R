#Integration into different chromatin compartments
#Integration into LADs
lad.df <- read.table(paste0(data.path, "distance_LAD.txt"),
                     sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                     col.names = c("chrom", "start", "end", "IN", "Sel", "strand", "Fchrom", "Fstart", "Fend", "distance"))

lad.count.df <- do.call(rbind,
                        lapply(in.cat,
                               function(i) {
                                 i.df <- lad.df[lad.df$IN == i & lad.df$Sel == "GFP",]
                                 tot <- nrow(i.df)
                                 in.feat <- nrow(i.df[i.df$distance == 0,])
                                 data.frame(IN = i,
                                            Total = tot,
                                            inFeature = in.feat
                                 )
                               }))
lad.count.df$SubC = "LAD"

#Integration into chromatin subcompartments
subc.df <- read.table(paste0(data.path, "distance_SubC.txt"),
                      sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                      col.names = c("chrom", "start", "end", "IN", "Sel", "strand", "Fchrom", "Fstart", "Fend", "Fname", "Fscore", "Fstrand", "distance"))
subc.df$ISID <- paste(subc.df$chrom, subc.df$start, subc.df$end, subc.df$IN, subc.df$Sel, sep = "_")

in.cat <- c("wt", "W390A", "CBX")
cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

subc.names <- c("A1", "A2", "B1", "B2", "B3")
comp.lvls <- c(subc.names, "LAD")

subc.count.df <- do.call(rbind,
                         lapply(in.cat,
                                function(i) {
                                  i.df <- subc.df[subc.df$IN == i & subc.df$Sel == "GFP",]
                                  tot <- length(unique(i.df$ISID))
                                  
                                  do.call(rbind,
                                          lapply(subc.names,
                                                 function(s) {
                                                   s.df <- i.df[i.df$Fname == s & i.df$distance == 0,]
                                                   in.feat <- length(unique(s.df$ISID))
                                                   data.frame(IN = i,
                                                              Total = tot,
                                                              SubC = s,
                                                              inFeature = in.feat)
                                                   
                                                 }))
                                }))

comp.count.df <- rbind(subc.count.df, lad.count.df)
comp.count.df$Perc <- 100 * comp.count.df$inFeature / comp.count.df$Total
enrich.base <- sapply(comp.lvls, function(x) {comp.count.df$Perc[comp.count.df$IN == "wt" & comp.count.df$SubC == x]})
comp.count.df$foldChange <- sapply(1:nrow(data.to.plot),
                                   function(w.r) {
                                     r.df <- data.to.plot[w.r,]
                                     s <- r.df$SubC
                                     p <- r.df$Perc
                                     round(p / enrich.base[names(enrich.base) == s], 2)
                                   })

data.to.plot <- comp.count.df
data.to.plot$IN <- factor(data.to.plot$IN, levels = in.cat)
data.to.plot$SubC <- factor(data.to.plot$SubC, levels = comp.lvls)

#Set limits for y axis
cat.limits <- c(80, 25, 15, 4, 7, 10)

limit.df <- do.call(rbind,
                    lapply(in.cat,
                           function(smpl) {
                             do.call(rbind,
                                     lapply(1:length(comp.lvls),
                                            function(w.s) {
                                              subc <- comp.lvls[w.s]
                                              if(length(cat.limits) == length(comp.lvls)) {
                                                m <- cat.limits[w.s]
                                              } else {
                                                m <- ceiling(max(data.to.plot$Perc[data.to.plot$SubC == subc],
                                                                 na.rm = TRUE))
                                              }
                                              data.frame(IN = factor(smpl, levels = in.cat),
                                                         SubC = factor(subc, levels = comp.lvls),
                                                         Perc = m
                                              )
                                            })
                             )
                           })
                       )



p.subc.0 <- ggplot(data.to.plot, aes(x = IN, y = Perc, fill = IN, group = IN)) +
  scale_fill_manual(values = alpha(cols, 0.2)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_blank(data = limit.df) +
  geom_col(colour = alpha("black", .7), linewidth = .25) +
  labs(fill = "integrase") +
  ylab("% IS in compartment") +
  xlab("integrase") +
  guides(#fill = "none",
    color = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 6),
        axis.line.y = element_line(linewidth = .3),
        axis.ticks.y = element_line(linewidth = .3),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 6),
        legend.title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt")
  )

p.subc.i <- p.subc.0 + facet_rep_wrap(~ SubC, scales = "free_y", repeat.tick.labels = TRUE)

p.subc.ii <- p.subc.0 + facet_rep_wrap(~ SubC, scales = "free_y", nrow = 1, repeat.tick.labels = TRUE)

ggsave(paste0(fig.path, "Figure2_C_SubC.png"), p.subc.ii,
       width = 96, height = 28, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Plot fold change in targeting frequency
data.to.plot <- comp.count.df[comp.count.df$IN != "wt",]
data.to.plot$IN <- factor(data.to.plot$IN,levels =  rev(in.cat[in.cat != "wt"]))
data.to.plot$SubC <- factor(data.to.plot$SubC, levels = comp.lvls)

p.fe.bar <- ggplot(data.to.plot, aes(x = foldChange, y = factor(SubC, levels = rev(comp.lvls)), fill = IN)) +
  scale_fill_manual(values = alpha(cols, 0.2)) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(trans = "log2", position = "top",
                     breaks = c(0.5, 1:4), labels = c(0.5, 1:4)) +
  geom_vline(xintercept = 1, colour = "black") +
  geom_vline(xintercept = c(1/4:2, 2:4), colour = alpha("gray", 0.5), linewidth = 0.2, linetype = 2) +
  geom_col(position = "dodge", colour = alpha("black", .7), width = .75,  linewidth = .25) +
  labs(fill = "integrase") +
  xlab("fold frequency change") +
  guides(color = "none") +
  coord_cartesian(xlim = c(1/4, 4)) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black", size = 4, face = "bold", angle = 0, hjust = 0, vjust = 0.5),
        axis.text.x = element_text(color = "black", size = 5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 6),
        axis.line.x = element_line(linewidth = .3),
        axis.ticks.x = element_line(linewidth = .3),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x.bottom = element_text(face = "bold", colour = "black", size = 6)
  )

#p.fe.bar

ggsave(paste0(fig.path, "Figure2_D_SubC_fold_bar.png"), p.fe.bar,
       width = 40, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
