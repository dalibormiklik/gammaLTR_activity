#Create scatter plot representing the FACS measurement

#maximum number of cells to be plotted
max.cells <- 10000

pos.value <- "GFP+"

#Set variables for plotting
x.values <- "FCS.A"
x.values <- "GFP"

x.name <- "FCS.A"
y.name <- "GFP"

lvls.pop <- c( "GFP+", "GFP-")

cols <- c("darkgreen", "grey60")
names(cols) <- lvls.pop

#Set variables
vir <- c("MoM", "FL", "SN", "Ko", "Cr")
vir.lvls <- c("MLV", "FeLV", "SNV", "KoRV", "CrERV")
var.lvls <- c("wt", "W390A", "CBX")

#load tables

#Extract data for each measured cell from FACS
export.dir <- data.path
export.file <- "Cell_Export.txt"

#Run 'join_FCS_data.R' on raw FlowJo-exported files to create export.file
#if(!file.exists(paste0(export.dir, export.file))) {
#  source(paste0(export.dir, "/join_FCS_data.R"))
#}

export.df <- read.table(paste0(export.dir, export.file),
                        sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE,
                        col.names = c("FSC.A", "FSC.H", "SSC.A", "SSC.H", "GFP", "Hoechst", "mCherry", "Time", "Virus", "Variant", "dpi", "Population"))

#set data.to.plot
sel.var <- c("wt", "W390A", "CBX")
sel.dpi <- 3
export.df$GID <- paste0(export.df$Virus, "_", export.df$Variant, "_", export.df$dpi)

sel.export.df <- export.df[!is.na(export.df$Variant) &
                             !is.na(export.df$Variant) &
                             export.df$Variant %in% sel.var &
                             export.df$Virus %in% vir &
                             export.df$dpi == sel.dpi,]
sel.export.df$Virus <- factor(sapply(sel.export.df$Virus,
                                    function(v) {
                                      vir.lvls[which(vir == v)]
                                    }),
                             levels = vir.lvls)

pos.export.df <- sel.export.df[sel.export.df$Population == pos.value,]
min.max.fsca <- c(min(sel.export.df$FSC.A), max(sel.export.df$FSC.A))
pos.export.df$FSC.A <- sample(min.max.fsca, size = nrow(pos.export.df), replace = TRUE)

#Calculate global number of cells to be plotted
#equals to min number of cells per sample | to limit set in max.cells object
min.cells <- floor(min(sapply(unique(sel.export.df$GID),
                              function(x) {
                                nrow(sel.export.df[sel.export.df$GID == x,])
                              }))/100) * 100
if(min.cells > max.cells) {min.cells <- max.cells}

#Set data.to.plot
data.to.plot <- do.call(rbind,
                        lapply(unique(sel.export.df$GID),
                               function(x) {
                                 df <- sel.export.df[sel.export.df$GID == x,]
                                 df[sample.int(n = nrow(df), size = min.cells),]
                               }))
                               
#Calculate "summary data" - % of GFP positive cells
export.sum.df <- do.call(rbind,
                         lapply(unique(sel.export.df$GID)[unique(sel.export.df$GID) %in% unique(data.to.plot$GID)],
                                function(v) {
                                  v.df <- sel.export.df[sel.export.df$GID == v,]
                                  tot <- nrow(v.df)
                                  vrs <- unique(v.df$Virus)
                                  vrs.var <- unique(v.df$Variant)
                                  v.pos <- nrow(v.df[v.df$Population == pos.value,])
                                  data.frame(GID = v,
                                             Virus = vrs,
                                             Variant = vrs.var,
                                             Cell_Count = tot,
                                             Positive = v.pos,
                                             Ratio = round(v.pos / tot, 3))
                                }))
export.sum.df$Virus <- factor(export.sum.df$Virus, levels = vir.lvls)

ax.tex.size <- 5
ax.tit.size <- 6

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

#Create plot
p.bulk.sp <- ggplot(data = data.to.plot, aes(x = FSC.A, y = GFP, colour = Population)) +
  scale_fill_manual(values = alpha(cols, c(.5, .5))) +
  scale_colour_manual(values = alpha(cols, c(.5, .5))) +
  scale_x_continuous(name = x.name, trans = log10_trans(), labels = label_log(base = 10, digits = 1)) +
  scale_y_continuous(name = y.name, trans = biexp_trans(lim = 700), labels = label_log(base = 10, digits = 1)) +
  geom_point(shape = 1, size = 0.001) +
  geom_boxplot(data = pos.export.df,
               colour = alpha("black", .8), outlier.shape = NA, fill = NA, na.rm = TRUE, inherit.aes = TRUE) +
  geom_text(data = export.sum.df, aes(label = 100 * Ratio), inherit.aes = TRUE,
            x = 5, y = 2200, colour = "black", size = 1.5, fontface = "bold") +
  coord_cartesian(ylim = c(1, 8*10^4),
                  #xlim = c(4.7, 5.3)
                  ) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black" ,face = "plain", size = ax.tex.size),
        #axis.text.x = element_text(color = "black" ,face = "plain", size = ax.tex.size, angle = 45),
        axis.text.x = element_blank(),
        axis.title = element_text(color = "black" , face = "bold", size = ax.tit.size),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "bold", size = 5),
        legend.justification = "top",
        legend.title = element_text(colour = "black", size = ax.tit.size - 1),
        legend.text = element_text(colour = "black", face = "bold", size = ax.tex.size),
        legend.key.size = unit(3, "pt")) +
  facet_grid(~ Virus + factor(Variant, levels = var.lvls))

#p.bulk.sp

ggsave(paste0(fig.path, "Figure3_B_FACS_3dpi.png"), p.bulk.sp,
       width = 150, height = 40, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

