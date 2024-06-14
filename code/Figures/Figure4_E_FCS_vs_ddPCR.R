#Connect FACS nd ddPCR data

#General variables
dd.file <- "vCN_all.txt"
fcs.file <- "GFP_time.txt"

#Selected sample names
vir.lvls2 <- c("MoMLV", "SFFV", "FeLV", "SNV", "CrERV", "KoRV")

#Set objects w/ info to experiments
dd.dpi <- c(14)
dd.exp.name <- c("E1")

#Colors for plots
#Set color palettes
cb.cols <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "#B9390C", "#C22DF3")
names(cb.cols) <- vir.lvls

#Load ddPCR data
vcn.2 <- read.table(paste0(data.path, dd.file),
                    sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)

#Time points of measurements
fcs.tab <- read.table(paste0(data.path, fcs.file),
                      sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)

#Join data from FCS and ddPCR
fcs.sel.tab <- fcs.tab[fcs.tab$virus %in% vir.lvls2 &
                       fcs.tab$dpi %in% dd.dpi &
                       fcs.tab$volume %in% unique(vcn.2$Volume), c("virus", "volume", "dpi", "gfp_perc")]
fcs.sel.tab$Sample <- with(fcs.sel.tab, paste0(virus, "_", volume, "_", dpi))

dd.sel.tab <- vcn.2[vcn.2$Virus %in% vir.lvls2 &
                    vcn.2$perCell > 0,]
dd.sel.tab$Sample <- with(dd.sel.tab, paste0(Virus, "_", Volume, "_", dpi))

smpl.names <- unique(dd.sel.tab$Sample)

fcs.dd.df <- do.call(rbind,
                     lapply(smpl.names,
                            function(s.nam) {
                              
                              s.df <- fcs.sel.tab[fcs.sel.tab$Sample == s.nam,]
                              
                              s.vcn <- dd.sel.tab$perCell[dd.sel.tab$Sample == s.nam]
                              
                              df <- do.call(rbind,
                                            lapply(s.vcn,
                                                   function(n) {
                                                     s.df$vcn <- n
                                                     s.df
                                                   }))
                              #active copies ratio (acr)
                              # = %GFP / VCN  
                              df$acr <- with(df, gfp_perc / vcn)

                              df
                              
                            }))
head(fcs.dd.df)

data.to.plot <- fcs.dd.df[fcs.dd.df$acr != Inf,]
data.to.plot$virus <- factor(data.to.plot$virus, levels = vir.lvls2)

p.acr <- ggplot(data = data.to.plot, aes(x = virus, y = acr, group = virus, fill = virus)) +
  scale_fill_manual(values = alpha(cb.cols, 0.5)) +
  geom_hline(yintercept = seq(0, 1, 0.1), linetype = 3, color = "gray", size = .15) +
  geom_bar(color = "black", position = "dodge", stat = "summary", fun = "mean", size = 0.25, width = .8) +
  geom_beeswarm(cex = 2, pch = 21, fill = NA, size = 0.75, stroke = 0.25) +
  coord_cartesian(ylim = c(0, 1.05),
                  xlim = c(0.35, 6.65),
                  expand = FALSE) +
  scale_y_continuous(name = "Active Genome Ratio",
                     breaks = c(seq(0, 1, 0.1), seq(1.2, 2, 0.2))
  ) +
  xlab("LTR U3") +
  guides(fill = "none",
         color = "none") +
  theme_classic() +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(color = "black", face = "bold", angle = 45, vjust = 1, hjust = 1, size = ax.tex.size + 2),
    axis.text.y = element_text(color = "black", size = ax.tex.size),
    axis.title.y = element_text(color = "black", face = "bold", size = ax.tit.size),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = 2, color = "black")
  )


#p.acr

#Save figure
ggsave(paste0(fig.path, "Figure4_E_Active_Copies_Ratio.png"), p.acr,
       width = 27, height = 35, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
