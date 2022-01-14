library(missMDA)

# load data
contdata<-read.csv("processed_continuous_data.csv")
catdata<-read.csv("processed_categorical_data.csv")

# Merge data
data <- merge(contdata, catdata, by="prim_key")

# Impute data
res_filled <- imputeFAMD(data, ncp=5)
res<-res_filled$completeObs
row.names(res) <- res$prim_key
res <- res[,-1]

# Save data
write.csv(res,"processed_filled.csv",row.names=TRUE)
