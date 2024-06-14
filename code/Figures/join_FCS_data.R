#Create single file from csv files exported from the Flow Jo

#Each sample (.fcs file) is exported into two .csv files according to gates:
#one .csv for positive (here GFP+), one .csv for negative (GFP-) population
#Here, the data from each sample are joined
#All samples (csv files) in 'export.dir' are joined into single .txt file

#variables like 'export.file' name are defined in the 'FigureX_00_Main.R' code

#Path to directory with exported fcs files
export.dir <- "./"

#Set variables describing samples
vars <- c("wt", "W390A", "CBX")

vir.lvls <- c("MoMLV", "FeLV", "SNV", "KoRV", "CrERV")
virs <- c("MoM", "FL", "SN", "Ko", "Cr")

c.names.sep <- "_"
c.name.pos <- c(3, 4, 5)
c.name.names <- c("virus", "variant", "dpi")

#Set population names
pos.name <- "GFP+"
neg.name <- "GFP-"

#Load names of files in 'export.dir'
f.names <- list.files(export.dir)

f.names.neg <- f.names[grep(neg.name, f.names)]

export.df <- do.call(rbind,
                     lapply(f.names.neg,
                            function(f.name) {
                              
                              p.name <- gsub(neg.name, pos.name, f.name, fixed = TRUE)
                              
                              #Load data from parent and daughter populations
                              neg.df <- read.csv(paste0(export.dir, "/", f.name))
                              
                              pos.df <- read.csv(paste0(export.dir, "/", p.name))
                              
                              neg.df$Population <- neg.name
                              pos.df$Population <- pos.name
                              
                              df <- rbind(neg.df, pos.df)
                              
                              sub.names <- unlist(strsplit(f.name, "_"))
                              
                              df$Cell <- sub.names[2]
                              print(paste0(sub.names[3], "_", sub.names[4]))
                              if(sub.names[3] %in% virs) {
                                df$Virus <- vir.lvls[which(virs == sub.names[3])]
                              } else {
                                df$Virus <- sub.names[3]
                              }
                              df$Variant <- sub.names[4]
                              df$dpi <- sub.names[5]
                              
                              df
                              
                            }))

write.table(export.df,
            paste0(data.path, export.file),
            sep = "\t", quote = FALSE, row.names = FALSE)
