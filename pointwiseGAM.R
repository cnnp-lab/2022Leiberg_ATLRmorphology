
library(mgcv)
library(data.table)

# Onset side, measure, hemi
#args <- commandArgs(TRUE)
args <- c("R", "T", "r") # run for onset side R/L, measure T/At/Ae, hemi r/l

dataset <- fread(paste("not_corrected/", args[2], "_", args[3], "h_", args[1], "TLE_full.csv", sep=""), header = F)
meta <- read.csv(paste("not_corrected/", args[1], "TLE_fullMeta.csv", sep=""))

meta$isFemale = as.factor(meta$isFemale)
meta$SP = as.factor(meta$SP)

corrected <- matrix(, nrow(dataset), ncol(dataset))

cols <- names(dataset)

# GAM correction in each point
for (kk in 1:163842){
  
  response = dataset[[cols[kk]]]
  
  tmp = data.frame(response = response[meta$IsControl==1], isFemale = meta[meta$IsControl==1,"isFemale"],
                   Age = meta[meta$IsControl==1,"Age"], SP = meta[meta$IsControl==1,"SP"])
  
  # run GAM if there are at least 5 controls of each SP
  if (sum(!is.nan(tmp[tmp$SP==1 & tmp$isFemale==1,]$response))>4 & sum(!is.nan(tmp[tmp$SP==1 & tmp$isFemale==0,]$response))>4 &
      sum(!is.nan(tmp[tmp$SP==2 & tmp$isFemale==1,]$response))>4 & sum(!is.nan(tmp[tmp$SP==2 & tmp$isFemale==0,]$response))>4 &
      sum(!is.nan(tmp[tmp$SP==3 & tmp$isFemale==1,]$response))>4 & sum(!is.nan(tmp[tmp$SP==3 & tmp$isFemale==0,]$response))>4){
    
    m = gam(response ~ s(Age, by=isFemale) + s(SP, bs = "re"), data=tmp, method = "REML")
    
    corrected[,kk] <- response - predict(m, data.frame(isFemale = meta$isFemale, Age = meta$Age, SP = meta$SP))
  }
}

fwrite(corrected, paste("age_sex_corrected/", args[2], "_", args[3], "h_", args[1], "TLE_full_corrected.csv", sep=""), sep=",",  col.names=FALSE, row.names = F)

