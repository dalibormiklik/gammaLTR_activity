#Connect FACS nd ddPCR data

#General variables
#Set working directory and experiment namewd
fcs.file <- "KgLTR_FCS_allE.txt"
dd.file <- "vCN_all.txt"

#Sample names
vir.lvls.short <- c("MoM", "FL", "SN", "Ko", "Cr")
vir.lvls <- c("MoMLV", "FeLV", "SNV", "KoRV", "CrERV")
var.lvls <- c("wt", "W390A", "CBX")

#Set objects w/ info to experiments
dd.dpi <- c(16,16,14)
dd.exp.name <- c("E1", "E2", "E4")

#Colors for plots
#Set color palettes
gr.cols <- c("gray90", "gray50", "gray10")
cb.cols.inv <- c("#004D40", "#1E88E5", "#D81B60")

select.pal <- cb.cols.inv

cols <- c("#09AB03", "#FFC107", "#D81B60")
names(cols) <- in.cat

#Load FCS data
#Time points of measurements
fcs.tab <- read.table(fcs.file,
                      sep = "\t", dec = ".", stringsAsFactors = FALSE, header = TRUE)

fcs.tab.dpi <- do.call(rbind,
                       lapply(1:length(dd.exp.name),
                              function(w.exp) {
                                fcs.tab[fcs.tab$Experiment == dd.exp.name[w.exp] & fcs.tab$dpi == dd.dpi[w.exp],]
                              }))

#Load ddPCR data
vcn.2 <- read.table(dd.file,
            quote = "", sep = "\t", header = TRUE)

#Join data from FCS and ddPCR
fcs.sel.tab <- fcs.tab.dpi[, c("Virus", "Variant", "GFP_Perc", "Experiment")]
fcs.sel.tab$Sample2 <- paste0(fcs.sel.tab$Virus, "-", fcs.sel.tab$Variant)

smpl.names <- unique(vcn.2$Sample[vcn.2$Virus %in% vir.lvls])

fcs.dd.df <- do.call(rbind,
                     lapply(smpl.names,
                            function(s.nam) {
                              dd.s.nam <- vcn.2[vcn.2$Sample == s.nam,]
                              do.call(rbind,
                                      lapply(unique(dd.s.nam$Experiment),
                                              function(dd.exp) {
                                                
                                                dd.s.nam.e <- dd.s.nam[dd.s.nam$Experiment == dd.exp,]
                                                
                                                vcn <- 100 * mean(dd.s.nam.e$perCell, na.rm = TRUE)
                                                fcs.s.nam <- fcs.sel.tab[fcs.sel.tab$Sample2 == s.nam & fcs.sel.tab$Experiment == dd.exp,]
                                                gfp <- fcs.s.nam$GFP_Perc
                                                #Create data frame for s.nam (Sample NAMe)
                                                df <- data.frame(Sample = s.nam,
                                                                 Virus = fcs.s.nam$Virus,
                                                                 Variant = fcs.s.nam$Variant,
                                                                 vCN = vcn,
                                                                 pGFP = gfp,
                                                                 vCN_per_pGFP = vcn / gfp,
                                                                 ActiveCop_ratio = gfp / vcn,
                                                                 Experiment = dd.exp)
                                                df
                                                
                                        
                                        
                                      }))
                            }))
#head(fcs.dd.df)
