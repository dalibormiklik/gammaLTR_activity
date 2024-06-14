#Load data
chrom.hmm <- read.table(paste0(data.path, "distance_ChromHMM_hg38.txt"),
                        sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                        col.names = c("chrom", "start", "end", "name", "seq", "strand", "Fchrom", "Fstart", "Fend", "Fname", "Fscore", "Fstrand", "FthickS", "FthickE", "Frgb", "distance"))

chrom.states.df <- read.table(paste0(data.path, "chromHMM_segment_family.txt"),
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

target.names <- c("DoT:3.1", "DoT:6.2", "DoT:6.3", "DoT:7.3", "IFT20")

chrom.hmm$name <- factor(chrom.hmm$name, levels = target.names)

col.chrfam <- c("#D81B60", "#FFC107", "#CC79A7", "#09AB03", "#56B4E9", "#999999")
names(col.chrfam) <- c("Promoter", "Enhancer", "Poised", "Transcription", "Insulator", "Repressed")

#Set which chromHMM state correspods to which chromHMM family
chrom.states <- chrom.states.df$Segment
names(chrom.states) <- chrom.states.df$Family

chrom.hmm$FeatureFam <- factor(sapply(chrom.hmm$Fname,
                                      function(x) {
                                        names(chrom.states[chrom.states == x])
                                        }),
                               levels = names(col.chrfam))

data.to.plot <- chrom.hmm

p.chst.0 <- ggplot(data.to.plot, aes(x = name, y = log10((distance + 1)), fill = FeatureFam)) +
  scale_fill_manual(values = alpha(col.chrfam, .6)) +
  scale_y_continuous(breaks = 0:7, labels = c(0, "", 2, "", 4, "", 6, "")) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .3) +
  geom_beeswarm(shape = 21, size = 1, cex = 1.5, colour = alpha("black", 1), stroke = .2) +
  geom_text(data = data.to.plot[data.to.plot$distance == 0,],
            aes(label = Fname), size = 1.2, angle = 45, hjust = 0, nudge_x = 0, nudge_y = 0.2) +
  ylab(bquote(log[10]~distance)) +
  xlab("target site") +
  guides(color = "none") +
  coord_cartesian(ylim = c(-.2, 7.5), xlim = c(0.5, 5.5), expand = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 5, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 6),
        axis.line.y = element_line(size = .3),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = .3),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 6),
        legend.text = element_text(size = 5, face = "bold", colour = "black"),
        legend.key.size = unit(5, "pt")
  )

#p.chst.0

ggsave(paste0(fig.path, "Figure5_S7B_KItarget_chromHMM.png"), p.chst.0,
       width = 60, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
