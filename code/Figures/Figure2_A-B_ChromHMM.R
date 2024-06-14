library(ImpactEffectsize)

#load ChromHMM distances
#Set variables
in.cat <- c("wt", "W390A", "CBX")

col.chrfam <- c("#D81B60", "#FFC107", "#CC79A7", "#09AB03", "#56B4E9", "#999999")
names(col.chrfam) <- c("Promoter", "Enhancer", "Poised", "Transcription", "Insulator", "Repressed")

#Load data
chrom.hmm <- read.table(paste0(data.path,"distance_ChromHMM.txt"),
                        sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                        col.names = c("chrom", "start", "end", "IN", "Sel", "strand", "Fchrom", "Fstart", "Fend", "Fname", "Fscore", "Fstrand", "distance"))

if(any(chrom.hmm$IN == "WA")) {
  w <- which(chrom.hmm$IN == "WA")
  chrom.hmm$IN[w] <- "W390A"
}

chrom.hmm$IN <- factor(chrom.hmm$IN, levels = in.cat)

#Set which chromHMM state correspods to which chromHMM family
chrom.states.df <- read.table(paste0(data.path, "chromHMM_segment_family.txt"),
                              sep = "\t", quote = "", header = TRUE)
chrom.states <- chrom.states.df$Segment
names(chrom.states) <- chrom.states.df$Family


#Create data.frame with distance medians
chrom.dist.stat <- do.call(rbind,
                          lapply(unique(chrom.hmm$Sel),
                                 function(s) {
                                   s.df <- chrom.hmm[chrom.hmm$Sel == s,]
                                   do.call(rbind,
                                           lapply(unique(s.df$IN),
                                                  function(i) {
                                                    i.df <- s.df[s.df$IN == i,]
                                                    do.call(rbind,
                                                            lapply(unique(i.df$Fname),
                                                                   function(f) {
                                                                     f.df <- i.df[i.df$Fname == f,]
                                                                     #Create dataframe with stats
                                                                     stat.df <- as.data.frame(t(unclass(summary(f.df$distance))))
                                                                     #Create dataframe with variable names
                                                                     var.df <- data.frame(IN = i,
                                                                                          Sel = s,
                                                                                          Fname = f)
                                                                     #Create final output dataframe
                                                                     cbind(var.df, stat.df)
                                                                   }))
                                                  }))
                                 }
                          ))

#Calculate Impact
in1 <- "wt"
in2 <- "CBX"

#Calculate impact
if(in1 != in2) {
  im.chrom <- do.call(rbind,
                      lapply(unique(chrom.hmm$Fname),
                             function(feat) {
                               
                               feat.chrom.hmm <- chrom.hmm[chrom.hmm$Fname == feat &
                                                             chrom.hmm$Sel == "GFP",]
                               
                               in.feat.chrom.hmm <- rbind(feat.chrom.hmm[feat.chrom.hmm$IN == in1,],
                                                          feat.chrom.hmm[feat.chrom.hmm$IN == in2,])
                               in.feat.chrom.hmm$IN <- factor(in.feat.chrom.hmm$IN, levels = c(in1, in2))
                               
                               feat.fam <- names(chrom.states[chrom.states == feat])
                               
                               
                               im <- Impact(in.feat.chrom.hmm$distance, Cl = in.feat.chrom.hmm$IN)
                               
                               if((im$Impact == 0)) {
                                 im$MorphDiff <- 0
                                 im$CTDiff <- 0
                               }
                               
                               im.df <- as.data.frame(t(unlist(im)))
                               
                               cbind(data.frame(Feature = feat,
                                                FeatFam = feat.fam,
                                                IN1 = in1,
                                                IN2 = in2),
                                     im.df
                                     )
                               
                             }))
}

#Plot pairwise impact

#Main Figure2 A:
#Impact vs CTdiff
data.to.plot <- im.chrom
p.impact <- ggplot(data.to.plot, aes(x = Impact, y = CTDiff, fill = FeatFam)) +
  scale_fill_manual(values = alpha(col.chrfam, 0.5)) +
  scale_y_continuous(breaks = seq(0, .5, 0.1)) +
  geom_hline(yintercept = -log10(0.01), colour = alpha("black", .1), linewidth = 0.2, linetype = 2) +
  geom_vline(xintercept = c(0,-0.2, 0.2), colour = c("black", "gray", "gray"), linewidth = 0.2, linetype = c(1, 2, 2)) +
  geom_point(shape = 21, size = 1.2, colour = "black", stroke = .2) +
  geom_text(data = data.to.plot[abs(as.numeric(data.to.plot$Impact)) > 0.5,],
            aes(label = Feature), fontface = 2, nudge_y = 0.006, nudge_x = -0.05, hjust = "right", size = 1.25) +
  ggtitle(paste0(in1, " to ", in2)) +
  guides(color = "none") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-0, .3), expand = FALSE) +
  theme_classic() +
  theme(plot.title = element_text(color="black", size = 4, face="bold", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_text(color = "black", size = 4, face = "bold"),
        axis.title.y = element_text(color = "black", size = 4, face = "bold"),
        axis.line = element_line(linewidth = .3),
        axis.ticks = element_line(linewidth = .3),
        legend.title = element_blank(),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt")
  )

p.impact

ggsave(paste0(fig.path, "Figure2_A_Impact_", in1, "-", in2, ".png"), p.impact,
       width = 55, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


#Supplemental Impact plots
#Impact distribution
data.to.plot <- im.chrom
p.impact.distr <- ggplot(data.to.plot, aes(x = IN1, y = abs(Impact))) +
  scale_fill_manual(values = alpha(col.chrfam, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), name = "Impact (absolute value)") +
  geom_violin(trim = TRUE, scale = "width", adjust = 0.5,
              fill = alpha("gray", .1), colour = alpha("black", .4), show.legend = FALSE) +
  geom_beeswarm(aes(fill = factor(FeatFam, levels = names(col.chrfam))),
                shape = 21, size = 1.2, colour = "black", stroke = .2, cex = 4) +
  ggtitle(paste0(in1, " to ", in2)) +
  guides(color = "none") +
  coord_cartesian(xlim = c(.53, 1.47),
                  ylim = c(0, .75),
                  expand = FALSE
                  ) +
  theme_classic() +
  theme(plot.title = element_text(color="black", size = 4, face="bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 4, face = "bold"),
        axis.line = element_line(linewidth = .3),
        axis.ticks = element_line(linewidth = .3),
        legend.title = element_blank(),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt")
  )

p.impact.distr

ggsave(paste0(fig.path, "Figure2_Supp_A_Impact_", in1, "-", in2, "_distribution.png"), p.impact.distr,
       width = 55, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Supplementary
#Central vs shape
data.to.plot <- im.chrom
p.ct.morph <- ggplot(data.to.plot, aes(x = CTDiff, y = MorphDiff)) +
  scale_fill_manual(values = alpha(col.chrfam, 0.5)) +
  geom_point(aes(fill = factor(FeatFam, levels = names(col.chrfam))),
                shape = 21, size = 1.2, colour = "black", stroke = .2) +
  xlab("Central tendency difference") +
  ylab("Morphology difference") +
  ggtitle(paste0(in1, " to ", in2)) +
  guides(color = "none") +
  coord_cartesian(xlim = c(0, .4),
                  ylim = c(0, .8),
                  expand = FALSE
  ) +
  theme_classic() +
  theme(plot.title = element_text(color="black", size = 4, face="bold", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_text(color = "black", size = 4, face = "bold"),
        axis.title.y = element_text(color = "black", size = 4, face = "bold"),
        axis.line = element_line(linewidth = .3),
        axis.ticks = element_line(linewidth = .3),
        legend.title = element_blank(),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt")
  )

p.ct.morph

ggsave(paste0(fig.path, "Figure2_Supp_B_CT-Morph_", in1, "-", in2, ".png"), p.ct.morph,
       width = 55, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#
#Plot distance distributions

cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

#Select plotted features
select.feat <- c("Tss", "Enh")
data.to.plot <- chrom.hmm[chrom.hmm$Sel == "GFP" & chrom.hmm$Fname %in% select.feat,]

p.chst.0 <- ggplot(data.to.plot, aes(x = IN, y = distance / 1000, fill = IN, group = IN)) +
  scale_fill_manual(values = alpha(cols, 0.2)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_quasirandom(shape = 21, size = 0.5, colour = alpha("black", 0.25), stroke = .1) +
  geom_boxplot(size = .1, outlier.shape = NA) +
  ylab("Distance (kb)") +
  xlab("integrase") +
  guides(fill = "none",
         color = "none") +
  coord_cartesian(ylim = c(-1, 101), xlim = c(0.5, 3.5), expand = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 5, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 6),
        axis.line.y = element_line(size = .3),
        axis.ticks.y = element_line(size = .3),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 6)
  )
p.chst <- p.chst.0 +
  facet_rep_grid(~ factor(Fname, levels = select.feat))

#p.chst

ggsave(paste0(fig.path, "Figure2_B_chromHMM_distance.png"), p.chst,
       width = 45, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Plot all feature distances
cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

select.feat <- chrom.states
data.to.plot <- chrom.hmm[chrom.hmm$Sel == "GFP" & chrom.hmm$Fname %in% select.feat,]

p.chst.all <- p.chst.0 +
  facet_wrap(~ factor(Fname, levels = select.feat), ncol = 5)

ggsave(paste0(fig.path, "Figure2_Supplement_chromHMM.png"), p.chst.all,
       width = 100, height = 120, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
                                            
