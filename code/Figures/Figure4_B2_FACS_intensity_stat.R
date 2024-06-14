MarkSignif <- function(log.pval) {
  if(log.pval > -log10(0.001)) {
    "***"
  } else {
    if(log.pval > -log10(0.01)) {
      "**"
    } else {
      if(log.pval > -log10(0.05)) {
        "*"
      } else {
        "NS"
      }
    }
  }
}
SignifDF <- function(input.df, features) {
  
  out.df <- do.call(rbind,
                    lapply(features,
                           function(f) {
                             
                             in.df <- input.df[input.df$Fname == f,]
                             
                             meds <- sapply(unique(in.df$GID), function(x) {median(in.df$GFP[in.df$GID == x], na.rm = TRUE)})
                             
                             k <- kruskal.test(GFP ~ GID, data = in.df)
                             kp <- -log10(k$p.value)
                             kp.df <-  data.frame(Test = "Kruskal",
                                                  Feature = f,
                                                  Var1 = "all",
                                                  Var2 = "all",
                                                  logP = kp,
                                                  sgnf = MarkSignif(kp),
                                                  foldChange = "_",
                                                  foldChange2 = "_")
                             
                             if(k$p.value < 0.05) {
                               kp <- -log10(k$p.value)
                               pwt <- pairwise.wilcox.test(in.df$GFP, in.df$GID, p.adjust.method = "BH")
                               pwt.m <- pwt$p.value
                               nr <- nrow(pwt.m)
                               cnames <- colnames(pwt.m)
                               rnames <- rownames(pwt.m)
                               pwt.df <- do.call(rbind,
                                                 lapply(1:ncol(pwt.m),
                                                        function(w.c) {
                                                          lp <- -log10(pwt.m[,w.c])
                                                          df <- data.frame(Test = "pWilcox",
                                                                           Feature = f,
                                                                           Var1 = rep(cnames[w.c], nr),
                                                                           Var2 = rnames,
                                                                           logP = -log10(pwt.m[,w.c]))
                                                          df <- df[!is.na(df$logP),]
                                                          df$sgnf <- sapply(df$logP, MarkSignif)
                                                          #fold change in median of a distance
                                                          df$foldChange <- sapply(1:nrow(df),
                                                                                  function(w.r) {
                                                                                    r.df <- df[w.r,]
                                                                                    v1 <- r.df$Var1
                                                                                    v2 <- r.df$Var2
                                                                                    round(meds[names(meds) == v2] / meds[names(meds) == v1], 4)
                                                                                  })
                                                          df$foldChange2 <- sapply(df$foldChange,
                                                                                   function(x) {
                                                                                     if(x < 1) {
                                                                                       round(1/x, 1)
                                                                                       } else {
                                                                                         round(x, 1)
                                                                                       }
                                                                                     })
                                                          df
                                                        }))
                               
                               rbind(kp.df, pwt.df)
                               
                             } else {
                               #output only output of Kruskal test
                               kp.df
                               
                             }
                             
                             
                           }))
  
  return(out.df)
  
}
FoldChDF <- function(input.df, col.group, col.values) {
  
  gid.col <- which(colnames(input.df) == col.group)
  val.col <- which(colnames(input.df) == col.values)
  
  gids <- unique(input.df[,gid.col])
  meds <- sapply(gids, function(x) {median(input.df[input.df[,gid.col] == x, val.col], na.rm = TRUE)})
  
  
  out.df <- do.call(rbind,
                    lapply(meds,
                           function(g) {
                             sapply(meds, function(m) {
                               if(m > g) {
                                 round(m / g, 1)
                               } else {
                                 round(g / m, 1)
                               }
                             }) 
                           })
                    )
  

  return(out.df)
  
}


features <- "GFP+"
v.df <- pos.gfpint.df[,c(y.values, "Virus", "Population", "GID")]
colnames(v.df) <- c("GFP", "Virus", "Population", "GID")
v.df$Fname <- v.df$Population

#calculate fold-change
f.df <- FoldChDF(v.df, "GID", "GFP")
#save table
write.table(f.df,
            paste0(data.path, "Figure4_SupTab_AS-gamma_intensity.txt"),
            quote = FALSE, sep = "\t")

#calculate p-values
p.df <- SignifDF(v.df, features)

p.mat.df <- pairwise.wilcox.test(v.df$GFP, v.df$GID, p.adjust.method = "BH")

#save table
write.table(round(-log10(p.mat.df$p.value),1),
            paste0(data.path, "Figure4_SupTab_AS-gamma_intensity_pair-pVal.txt"),
            quote = FALSE, sep = "\t")

