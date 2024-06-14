#Create plots for AS.gamma bulk FACS measurements

#Plots created
#p.bulk.sp - scatter plots of cytometric data (GFP vs FCS)
#p.bulk - %GFP / time of all variants / stock volumes
#p.bulk.rel - p.bulk realtive change to expression at 3 dpi

#Intensity plot
#Extract data for each measured cell from FACS
export.dir <- data.path
export.file <- "Cell_Export.txt"

#load tables

#Set variables
#maximum number of cells to be plotted
max.cells <- 10000

pos.value <- "GFP+"
lvls.pop <- c( "GFP+", "GFP-")

cols <- c("darkgreen", "grey60")
names(cols) <- lvls.pop

#Set variables for plotting
x.values <- "FSC-A"
y.values <- "B.530_30.A"

x.name <- "FCS.A"
y.name <- "GFP"


ax.tex.size <- 3
ax.tit.size <- 4

#Run 'join_FCS_data.R' on raw FlowJo-exported files to create export.file
#if(!file.exists(gfpint.file)) {
#  source(paste0(fig.scr, "Figure4_join_FCS_data.R"))
#}

gfpint.df <- read.table(paste0(data.path, gfpint.file),
                        sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE,
                        )
gfpint.df$GID <- gfpint.df$Virus

sel.gfpint.df <- gfpint.df[gfpint.df$Virus %in% vir.lvls,]
pos.gfpint.df <- sel.gfpint.df[sel.gfpint.df$Population == pos.value,]

#Calculate global number of cells to be plotted
#equals to min number of cells per sample | to limit set in max.cells object
min.cells <- floor(min(sapply(unique(sel.gfpint.df$GID),
                              function(x) {
                                nrow(sel.gfpint.df[sel.gfpint.df$GID == x,])
                              }))/100) * 100
if(min.cells > max.cells) {min.cells <- max.cells}

data.to.plot <- do.call(rbind,
                        lapply(unique(sel.gfpint.df$GID),
                               function(x) {
                                 df <- sel.gfpint.df[sel.gfpint.df$GID == x,]
                                 df[sample.int(n = nrow(df), size = min.cells),]
                               }))

#Calculate "summary data" - % of GFP positive cells
gfpint.sum.df <- do.call(rbind,
                         lapply(unique(gfpint.df$GID)[unique(gfpint.df$GID) %in% unique(data.to.plot$GID)],
                                function(v) {
                                  v.df <- gfpint.df[gfpint.df$GID == v,]
                                  tot <- nrow(v.df)
                                  vrs <- vir.lvls[which(vir.lvls == unique(v.df$Virus))]
                                  vrs.var <- unique(v.df$Variant)
                                  vrs.vol <- unique(v.df$Volume)
                                  v.pos <- nrow(v.df[v.df$Population == pos.value,])
                                  df <- data.frame(GID = v,
                                             Virus = vrs,
                                             Cell_Count = tot,
                                             Positive = v.pos,
                                             Ratio = round(v.pos / tot, 3))
                                  df
                                }))
gfpint.sum.df$Virus <- factor(gfpint.sum.df$Virus, levels = vir.lvls)

meds <- sapply(gfpint.sum.df$Virus,
               function(v) {
                 v.df <- pos.gfpint.df[pos.gfpint.df$Virus == v,]
                 median(v.df[,which(colnames(v.df) == y.values)], na.rm = TRUE)
               })
names(meds) <- gfpint.sum.df$Virus
vir.lvls2 <- names(meds)[order(meds, decreasing = TRUE)]

#biexp scale transform function
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
      10^sort(unique(round(log10(100:max(x)))))
    }
  }
  labels <- label_log(base = 10, digits = 1)
  trans_new(paste0("biexp-",format(lim)), trans, inv, breaks)
} 


y.breaks = 10^pretty(log10(data.to.plot$B.530_30.A))
y.labels = formatC(y.breaks, format = "e", digits = 2)

p.f04.b <- ggplot(data = data.to.plot, aes(x = FSC.A, y = B.530_30.A, colour = Population, group = Virus)) +
  scale_fill_manual(values = alpha(cols, c(.5, .5))) +
  scale_colour_manual(values = alpha(cols, c(.5, .5))) +
  scale_x_continuous(name = x.name, trans = log10_trans(), labels = label_log(base = 10, digits = 1)) +
  scale_y_continuous(name = y.name, trans = biexp_trans(lim = 300), labels = label_log(base = 10, digits = 1)) +
  geom_point(shape = 1, size = 0.001) +
  geom_boxplot(data = pos.gfpint.df,
               colour = alpha("black", .8), outlier.shape = NA, fill = NA, na.rm = TRUE, inherit.aes = TRUE) +
  geom_text(data = gfpint.sum.df, aes(label = round(100 * Ratio), 1), inherit.aes = TRUE,
            x = 4.8, y = 920, colour = "black", size = 2.5, fontface = "bold") +
  coord_cartesian(ylim = c(1, 3*10^4)
  ) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black" ,face = "plain", size = ax.tex.size + 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(color = "black" , face = "bold", size = ax.tit.size),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = ax.tit.size + 2),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size - 1),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(3, "pt")) +
  facet_grid(~ factor(Virus, levels = vir.lvls2))

#p.f04.b

ggsave(paste0(fig.path, "/Figure4_B_FACS-intens_3dpi.png"), p.f04.b,
       width = 120, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

