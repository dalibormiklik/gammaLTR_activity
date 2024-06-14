#Path to the data
export.dir <- gfpint.dir


#Set variables describing samples
vars <- c("wt")

vir.lvls <- c("SFFV","MoMLV", "FeLV", "SNV", "KoRV", "CrERV")
virs <- c("AS-SF","AS-Mo", "AS-Fe", "AS-SN", "AS-Ko", "AS-Cr")

c.names.sep <- "_"
c.name.pos <- c(4)
c.name.names <- c("virus")

#Set population names
parent.name <- "alive"
daughter.name <- "GFP+"
pos.name <- "GFP+"
neg.name <- "GFP-"

#Load file names
f.names <- list.files(export.dir)

f.names.parent <- f.names[grep(parent.name, f.names)]
f.names.neg <- f.names[grep(neg.name, f.names)]

export.df <- do.call(rbind,
                     lapply(f.names.neg,
                            function(f.name) {
                              
                              p.name <- gsub(neg.name, pos.name, f.name, fixed = TRUE)
                              
                              #Load data from parent and daughter populations
                              neg.df <- read.csv(paste0(export.dir, "/", f.name))
                              
                              pos.df <- read.csv(paste0(export.dir, "/", p.name))
                              
                              if(nrow(neg.df) > 0) {neg.df$Population <- neg.name}
                              if(nrow(pos.df) > 0) {pos.df$Population <- pos.name}
                              
                              df <- rbind(neg.df, pos.df)
                              
                              sub.names <- unlist(strsplit(f.name, "_"))
                              
                              df$Cell <- sub.names[2]
                              print(paste0(sub.names[3], "_", sub.names[4]))
                              if(sub.names[3] %in% virs) {
                                df$Virus <- vir.lvls[which(virs == sub.names[3])]
                              } else {
                                df$Virus <- sub.names[3]
                              }
                              df$well <- sub.names[6]
                              df$dpi <- sub.names[4]
                              
                              df
                              
                            }))



write.table(export.df,
            paste0(export.dir, "/", gfpint.file),
            sep = "\t", quote = FALSE, row.names = FALSE)
