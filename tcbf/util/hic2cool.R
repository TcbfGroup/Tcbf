library(data.table)
data <- fread("A2_40000.matrix",header=F)
bed <-  fread("A2_40000.bed",header = F)

bed <- bed[,.(m1 = min(V4),m2 = max(V4)),by = .(V1)]
bed$index <-  1:dim(bed)[1]
for (line in  1:dim(bed)[1]){
  start = as.integer(bed[index==line][,2])
  end <- as.integer(bed[index==line][,3])
  result <- data[V1>=start&V1<=end&V2>=start&V2<=end]
  fwrite(result,paste("A2/",line,"_",line,".txt",sep = ""),sep = "\t")
}