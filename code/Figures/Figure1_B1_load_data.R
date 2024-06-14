#Load fcs data for Figure 1 A (fcs dot plot)
#Extract data for each measured cell from FACS
export.file <- "MLV_3dpi_fcs.txt"

#Create new "export.file" if one does not exist
if(!file.exists(paste0(data.path, export.file))) {
  source(paste0(export.dir, "/join_FCS_data.R"))
}

export.df <- read.table(paste0(data.path, export.file),
                        sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE,
                        col.names = c("FSC.A", "SSC.A", "GFP..488B..A", "mCherry..561D..A", "Cell_ID", "Variant", "Volume", "Population"))
export.df$GID <- paste0(export.df$Variant, "_", export.df$Volume)
export.df$Virus <- vir.lvls

#Select data from the joined fcs data
sel.export.df <- export.df


#Calculate global number of cells to be plotted
#equals to min number of cells per sample | to limit set in max.cells object
min.cells <- floor(min(sapply(unique(sel.export.df$GID),
                              function(x) {
                                nrow(sel.export.df[sel.export.df$GID == x,])
                              }))/100) * 100
if(min.cells > max.cells) {min.cells <- max.cells}

#Calculate "summary data" - % of GFP positive cells
export.sum.df <- do.call(rbind,
                         lapply(unique(export.df$GID)[unique(export.df$GID) %in% unique(data.to.plot$GID)],
                                function(v) {
                                  v.df <- export.df[export.df$GID == v,]
                                  tot <- nrow(v.df)
                                  vrs <- vir.lvls[which(vir.lvls == unique(v.df$Virus))]
                                  vrs.var <- unique(v.df$Variant)
                                  vrs.vol <- unique(v.df$Volume)
                                  v.pos <- nrow(v.df[v.df$Population == pos.value,])
                                  data.frame(GID = v,
                                             Virus = vrs,
                                             Variant = vrs.var,
                                             Volume = vrs.vol,
                                             Cell_Count = tot,
                                             Positive = v.pos,
                                             Ratio = round(v.pos / tot, 3))
                                }))
export.sum.df$Virus <- factor(export.sum.df$Virus, levels = vir.lvls)
