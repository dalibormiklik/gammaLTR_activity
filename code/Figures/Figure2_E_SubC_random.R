#Integration into different chromatin compartments
#Integration into LADs
lad.ran.df <- read.table(paste0(data.path, "distance_ran_LAD.txt"),
                     sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                     col.names = c("chrom", "start", "end", "IN", "Sel", "strand", "Fchrom", "Fstart", "Fend", "distance"))

lad.ran.count.df <- do.call(rbind,
                        lapply(in.cat,
                               function(i) {
                                 i.df <- lad.ran.df[lad.ran.df$IN == i & lad.ran.df$Sel == "GFP",]
                                 tot <- nrow(i.df)
                                 in.feat <- nrow(i.df[i.df$distance == 0,])
                                 data.frame(IN = i,
                                            Total = tot,
                                            inFeature = in.feat
                                 )
                               }))
lad.ran.count.df$SubC = "LAD"

#Integration into chromatin subcompartments
subc.ran.df <- read.table(paste0(data.path, "distance_ran_SubC.txt"),
                      sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                      col.names = c("chrom", "start", "end", "IN", "Sel", "strand", "Fchrom", "Fstart", "Fend", "Fname", "Fscore", "Fstrand", "distance"))
subc.ran.df$ISID <- paste(subc.ran.df$chrom, subc.ran.df$start, subc.ran.df$end, subc.ran.df$IN, subc.ran.df$Sel, sep = "_")

in.cat <- c("wt", "W390A", "CBX")
cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

subc.names <- c("A1", "A2", "B1", "B2", "B3")
comp.lvls <- c(subc.names, "LAD")

subc.ran.count.df <- do.call(rbind,
                         lapply(in.cat,
                                function(i) {
                                  i.df <- subc.ran.df[subc.ran.df$IN == i & subc.ran.df$Sel == "GFP",]
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

comp.ran.count.df <- rbind(subc.ran.count.df, lad.ran.count.df)
comp.ran.count.df$Perc <- 100 * comp.ran.count.df$inFeature / comp.ran.count.df$Total
enrich.base <- sapply(comp.lvls, function(x) {comp.ran.count.df$Perc[comp.ran.count.df$IN == "wt" & comp.ran.count.df$SubC == x]})
comp.ran.count.df$foldChange <- sapply(1:nrow(comp.ran.count.df),
                                   function(w.r) {
                                     r.df <- comp.ran.count.df[w.r,]
                                     s <- r.df$SubC
                                     p <- r.df$Perc
                                     round(p / enrich.base[names(enrich.base) == s], 2)
                                   })

data.to.plot <- comp.ran.count.df
data.to.plot$IN <- factor(data.to.plot$IN, levels = in.cat)
data.to.plot$SubC <- factor(data.to.plot$SubC, levels = comp.lvls)

p.subc.ran <- ggplot(data.to.plot, aes(x = SubC, y = Perc, fill = IN)) +
  scale_fill_manual(values = alpha(cols, 0.2)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", width = .75,
           fill = alpha("gray", 0.5), colour = alpha("black", .7), size = 0.3) +
  geom_beeswarm(shape = 21, size = .75, colour = alpha("black", 0.9), stroke = .1) +
  labs(fill = "integrase") +
  ylab("% shuffled sites in\ncompartment") +
  xlab("compartment") +
  guides(fill = "none",
    color = "none") +
  coord_cartesian(ylim =c(-0.5, 52), xlim = c(0.4, 6.6), expand = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", face = "bold", angle = 45, size = 4, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 4),
        axis.line.y = element_line(size = .3),
        axis.ticks.y = element_line(size = .3),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 6),
        legend.title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt")
  )

#p.subc.ran

ggsave(paste0(fig.path, "Figure2_E_Ran_dist_SubC.png"), p.subc.ran,
       width = 22, height = 28, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
