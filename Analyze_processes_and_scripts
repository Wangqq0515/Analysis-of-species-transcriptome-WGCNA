1. 首先整合拟南芥各个数据的FPKM表达文件到一个文件里
WGCNA代码：
```console
goodSample <- c("FPKM.H10m.1","FPKM.H10m.2","FPKM.H10m.3","FPKM.H1h.1","FPKM.H1h.2","FPKM.H1h.3","FPKM.H30m.1","FPKM.H30m.2","FPKM.H30m.3","FPKM.H4h.1","FPKM.H4h.2","FPKM.H4h.3","FPKM.S10m.1","FPKM.S10m.2","FPKM.S10m.3","FPKM.S1h.1","FPKM.S1h.2","FPKM.S1h.3","FPKM.S30m.1","FPKM.S30m.2","FPKM.S30m.3","FPKM.S4h.1","FPKM.S4h.2","FPKM.S4h.3")
cot_goodSample <- total_FPKM[,goodSample]
cot_goodSample_FPKM <- total_FPKM[,goodSample]
cot_goodSample_phone <- cot_pheno[goodSample,]
datExpr <- cot_goodSample_FPKM
```

# 筛选方法：mad>1且top5000
```console
WGCNA_matrix = t(datExpr[order(apply(datExpr,1,mad),decreasing = T)[1:5000],])
datExpr_filted <- WGCNA_matrix
```

检查离群样品
如果你的是TCGA上下载的数据，几百个临床样品，这一步是必须的，可能会有离群值，需要把它们去掉。
```console
library(WGCNA)
gsg = goodSamplesGenes(datExpr_filted, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr_filted), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
save(datExpr_filted, file = "datExpr_filted_hclust.RData")
```
