## Fig.S1E ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_scRNA1_GeneralCelltype_20240902.Rdata") # Control&PD
QC_cluster_vlnplot1 = VlnPlot(MOU_scRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "Sample", ncol = 2, pt.size = 0) 
ggsave(QC_cluster_vlnplot1,file="./MOU_scRNA1_Vlnplot_QC1_Sample.jpg",width=10,height=10,dpi=600)
QC_cluster_vlnplot2 = VlnPlot(MOU_scRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "Sample", ncol = 2, pt.size = 0.01) 
ggsave(QC_cluster_vlnplot2,file="./MOU_scRNA1_Vlnplot_QC2_Sample.jpg",width=10,height=10,dpi=600)


## Fig.S1F-G ######################
# Mouse
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_scRNA1_GeneralCelltype_20240902.Rdata") # Control&PD
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
#绘制细胞种类百分比 Excel表及条形图（各样本）
project_percent <- with(MOU_scRNA1@meta.data, prop.table(table(general_celltype, Sample), margin = 2) * 100)
project_percent <- as.data.frame(project_percent)
project_percent$general_celltype <- as.character(project_percent$general_celltype)
project_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project_percent$general_celltype <- as.factor(c(project_percent$general_celltype))
project_percent$general_celltype <- fct_relevel(project_percent$general_celltype, c("Fibroblast","Epithelial","Endothelial","Immune", "Myocyte","Schwann_Neuron"))
project_percent$Sample <- fct_relevel(project_percent$Sample, c("Control-1","Control-2","Control-3","PD-1","PD-2","PD-3","PD-4","PD-5","PD-6"))
project_percent$general_celltype
project_percent$Sample
project_percent <- project_percent[project_percent$general_celltype != "x", ] #删除x行
project_percent
write.table(project_percent, "./MOU_general_celltype_Bar2_percentage_samples.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project_percent, aes(x = Sample, y = Freq, fill = general_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./MOU_general_celltype_Bar2_percentage_samples.jpg",width = 8, height=4.5, dpi=600)

# Human
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step2_HU_scRNA_GeneralCelltype_20240620.Rdata")
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
#绘制细胞种类百分比 Excel表及条形图（各样本）
project_percent <- with(HU_scRNA@meta.data, prop.table(table(general_celltype, Sample), margin = 2) * 100)
project_percent <- as.data.frame(project_percent)
project_percent$general_celltype <- as.character(project_percent$general_celltype)
project_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project_percent$general_celltype <- as.factor(c(project_percent$general_celltype))
project_percent$general_celltype <- fct_relevel(project_percent$general_celltype, c("Fibroblast","Epithelial","Endothelial","Immune", "Myocyte","Schwann_Neuron"))
project_percent$Sample <- fct_relevel(project_percent$Sample, c("Control-1","Control-2","Control-3","PD-1","PD-2","PD-3","PD-4","PD-5","PD-6"))
project_percent$general_celltype
project_percent$Sample
project_percent <- project_percent[project_percent$general_celltype != "x", ] #删除x行
project_percent
write.table(project_percent, "./HU_general_celltype_Bar2_percentage_samples.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project_percent, aes(x = Sample, y = Freq, fill = general_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./HU_general_celltype_Bar2_percentage_samples.jpg",width = 8, height=4.5, dpi=600)


## Fig.S2A ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_Immune1_Celltype_20240907.Rdata") # Control&PD
#绘制Marker基因表达的Dot图
Immune_celltype_markers =c("Cd3d","Cd3e","Trac","Trbc2",  #T & αβT ("Cd2","Cd3","Trac","Trbc","Cd28")
                           "Cd4",  #αβT_CD4
                           "Cd8b1",  #αβT_CD8
                           "Trdc","Tcrg-C1",  #γδT ("Trdc","Tcrg")
                           "Ms4a1","Bank1","Mef2c","Jchain","Ighg1","Ighg3",  #B/Plasma (all)
                           "Cd86","H2-Aa","H2-Eb1","Cd209a", #DC ("Cd40","Cd80","Cd86","Cd209a","H2-Aa","H2-Eb1")
                           "C1qa","C1qb","Ctsb","Adgre1",  #Macrophage(Adgre1=F4/80,Itgam=Cd11b)
                           "S100a8","S100a9","G0s2","Trem1", #Neutrophil
                           "Tpsab1","Tpsb2","Cpa3","Kit")  #Mast (all)
p1 <-  DotPlot(object =MOU_Immune1, features = Immune_celltype_markers, group.by = "Immune_celltype") + 
  scale_color_gradient2(low = "white", high = "black") +  #点的颜色 
  labs(x = " ", y = " ") + coord_flip() +
  theme(plot.background = element_rect(fill = "white"), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(p1, file = "./Dotplot_MOU_Immune1_Immune_celltype_MarkerGenes.jpg",width = 8, height=13, dpi=600)


## Fig.S2B ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step4_HU_Immune_celltype_20240907.Rdata") 
#绘制Marker基因表达的Dot图
HU_Immune_celltype_markers =c("CD3D","CD3E","TRAC","TRBC2",  #T & αβT ("CD2","CD3","TRAC","TRBC")
                              "CD4",  #αβT_CD4C
                              "CD8B", #αβT_CD8
                              "TRDC","TRGC1",  #γδT ("TRDC","TCRG","IL17";多数γδT为CD4和CD8双阴性，少数γδT表达CD8)
                              "MS4A1","BANK1","MEF2C","JCHAIN","IGHG1","IGHG3", #B/PLASMA (ALL)
                              "CD86","HLA-DRA","HLA-DRB1",  #DC ("CD40","CD80","CD86","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1")
                              "C1QA","C1QB","CTSB","CD163",  #MACROPHAGE
                              "S100A8","S100A9","G0S2","TREM1", #NEUTROPHIL
                              "TPSAB1","TPSB2","CPA3","KIT")  #MAST (ALL)
p1 <-  DotPlot(object =HU_Immune, features = HU_Immune_celltype_markers, group.by = "Immune_celltype") + 
  scale_color_gradient2(low = "white", high = "black") +  #点的颜色 
  labs(x = " ", y = " ") + coord_flip() +
  theme(plot.background = element_rect(fill = "white"), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(p1, file = "./Dotplot_HU_Immune_Immune_celltype_MarkerGenes.jpg",width = 8, height=13, dpi=600)



## Fig.S2F-G ######################
setwd(dir = "/public/home/yinqi/zhaoyi/periodontitis/Review_20250303/Review文章")
getwd()
load("/public/home/yinqi/zhaoyi/periodontitis/MOU_HU_Fig_20240902_now/Step8_Neutro_HU_MOU_celltype_20241117.Rdata")

# 鼠Neutro_celltype
Idents(Neutro) <- "origin"
MOU_Neutro=subset(Neutro,idents=c("MOU-Ctrl","MOU-PD"))

colors <- c("#50013b", "#64004a","#8F006D","#a20063","#AF4A92", "#D2A6C7","#E8D3F3")
group_percent <- with(MOU_Neutro@meta.data, prop.table(table(Neutro_celltype, project), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$Neutro_celltype <- as.character(group_percent$Neutro_celltype)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$Neutro_celltype <- as.factor(c(group_percent$Neutro_celltype))
group_percent$Neutro_celltype <- fct_relevel(group_percent$Neutro_celltype, c("G1","G2","G3","G4","G5a","G5b","G5c"))
group_percent$project <- fct_relevel(group_percent$project, c("Control","PD"))
group_percent$Neutro_celltype
group_percent$project
group_percent <- group_percent[group_percent$Neutro_celltype != "x", ] #删除x行
group_percent
Bar_celltype <- ggplot(group_percent, aes(x = project, y = Freq, fill = Neutro_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_MOU_Neutro_celltypes_percentage.jpg",width = 4, height=4.5, dpi=600)

# 人Neutro_celltype
Idents(Neutro) <- "origin"
HU_Neutro=subset(Neutro,idents=c("HU-Ctrl","HU-PD"))

colors <- c("#50013b", "#64004a","#8F006D","#a20063","#AF4A92", "#D2A6C7","#E8D3F3")
group_percent <- with(HU_Neutro@meta.data, prop.table(table(Neutro_celltype, project), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$Neutro_celltype <- as.character(group_percent$Neutro_celltype)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$Neutro_celltype <- as.factor(c(group_percent$Neutro_celltype))
group_percent$Neutro_celltype <- fct_relevel(group_percent$Neutro_celltype, c("G1","G2","G3","G4","G5a","G5b","G5c"))
group_percent$project <- fct_relevel(group_percent$project, c("Control","PD"))
group_percent$Neutro_celltype
group_percent$project
group_percent <- group_percent[group_percent$Neutro_celltype != "x", ] #删除x行
group_percent
Bar_celltype <- ggplot(group_percent, aes(x = project, y = Freq, fill = Neutro_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_HU_Neutro_celltypes_percentage.jpg",width = 4, height=4.5, dpi=600)

# origin
Idents(Neutro) <- "origin"
Neutro2=subset(Neutro,idents=c("MOU-Ctrl","MOU-PD", "HU-Ctrl","HU-PD"))

origin_data <- data.frame(Group = Neutro2@meta.data$origin)
origin_data_summary <- origin_data %>% group_by(Group) %>% summarise(Freq = n()) %>% mutate(Freq = Freq / sum(Freq))
origin_data_summary$ALL <- "ALL"
origin_data_summary$Group <- fct_relevel(origin_data_summary$Group, c("MOU-Ctrl","MOU-PD","HU-Ctrl","HU-PD"))
ggplot(origin_data_summary, aes(x = ALL, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Proportion", fill = " ") + 
  scale_fill_manual(values = c("HU-Ctrl" = "#B7B2D0", "HU-PD" = "#7DA6C6",
                               "MOU-Ctrl" = "#EAAA60", "MOU-PD" = "#E68B81")) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
# 保存图像
ggsave("Bar_Neutro_origin_percentage.jpg", width = 2.5, height = 4.5, dpi = 600)


## Fig.S2H-I ######################
load("./Step8_Neutro_HU_MOU_celltype_20241117.Rdata")
Idents(Neutro) <- "origin"
MOU_Neutro=subset(Neutro,idents=c("MOU-Ctrl","MOU-PD"))

## monocle
# 提取表达矩阵、细胞信息和特征信息
expr_matrix <- as(as.matrix(MOU_Neutro@assays$RNA@data), "sparseMatrix")
cell_metadata <- MOU_Neutro@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expr_matrix))
rownames(gene_annotation) <- rownames(expr_matrix)
# 创建Monocle对象
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
mono <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
# 进行数据标准化和预处理
mono <- estimateSizeFactors(mono)
mono <- estimateDispersions(mono) #用时较长
# 过滤掉低表达基因和低质量细胞
mono <- mono[rowMeans(exprs(mono)) > 0.1, ]
mono <- mono[, colSums(exprs(mono) > 0) > 200]
# 选择用于拟时序分析的高变基因
disp_table <- dispersionTable(mono)
vari_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mono <- setOrderingFilter(mono, vari_clustering_genes$gene_id)
# 降维和拟时序分析
mono <- reduceDimension(mono, max_components = 2, method = 'DDRTree') #用时很长
mono <- orderCells(mono) #用时很长
save(mono, file="./Monocle2_MOU_Neutro_20241118.Rdata")
head(mono@phenoData@data)

# MOU 细胞密度图
load("/public/home/yinqi/zhaoyi/periodontitis/MOU_HU_Fig_20240902_now/Monocle2_MOU_Neutro_20241118.Rdata")
bar_data <- mono@phenoData@data[, c("Pseudotime", "origin", "celltype_col")]
bar_data_order <- bar_data[order(bar_data$Pseudotime), ] #数据框按Pseudotime列由小到大排序

colors1 <- c("#6EC3C9","#F19373")
gg <- ggplot(bar_data_order, aes(Pseudotime, colour = origin, fill=origin)) +
  geom_density(bw=0.5,size=1,alpha = 0.5) +
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除主要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(gg,file="./Density_plot_MOU_Neutro_origin.jpg",width=8,height=4,dpi=600)


## Human
load("./Step8_Neutro_HU_MOU_celltype_20241117.Rdata")
Idents(Neutro) <- "origin"
HU_Neutro=subset(Neutro,idents=c("HU-Ctrl","HU-PD"))

## monocle
# 提取表达矩阵、细胞信息和特征信息
expr_matrix <- as(as.matrix(HU_Neutro@assays$RNA@data), "sparseMatrix")
cell_metadata <- HU_Neutro@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expr_matrix))
rownames(gene_annotation) <- rownames(expr_matrix)
# 创建Monocle对象
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
mono <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
# 进行数据标准化和预处理
mono <- estimateSizeFactors(mono)
mono <- estimateDispersions(mono) #用时较长
# 过滤掉低表达基因和低质量细胞
mono <- mono[rowMeans(exprs(mono)) > 0.1, ]
mono <- mono[, colSums(exprs(mono) > 0) > 200]
# 选择用于拟时序分析的高变基因
disp_table <- dispersionTable(mono)
vari_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mono <- setOrderingFilter(mono, vari_clustering_genes$gene_id)
# 降维和拟时序分析
mono <- reduceDimension(mono, max_components = 2, method = 'DDRTree') #用时很长
mono <- orderCells(mono) #用时很长
save(mono, file="./Monocle2_HU_Neutro_20241118.Rdata")
head(mono@phenoData@data)

# HU 细胞密度图
load("/public/home/yinqi/zhaoyi/periodontitis/MOU_HU_Fig_20240902_now/Monocle2_HU_Neutro_20241118.Rdata")
bar_data <- mono@phenoData@data[, c("Pseudotime", "origin", "celltype_col")]
bar_data_order <- bar_data[order(bar_data$Pseudotime), ] #数据框按Pseudotime列由小到大排序

colors1 <- c("#6EC3C9","#F19373")
gg <- ggplot(bar_data_order, aes(Pseudotime, colour = origin, fill=origin)) +
  geom_density(bw=0.5,size=1,alpha = 0.5) +
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除主要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(gg,file="./Density_plot_HU_Neutro_origin.jpg",width=8,height=4,dpi=600)



## Fig.S3E ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_scRNA3_GeneralCelltype_20240907.Rdata") # Lig$LigPg
Idents(MOU_scRNA3) <- "Sample"
MOU_scRNA3=subset(MOU_scRNA3,idents=c("LigPg-1","LigPg-2","LigPg-3"))
QC_cluster_vlnplot1 = VlnPlot(MOU_scRNA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "Sample", ncol = 4, pt.size = 0) 
ggsave(QC_cluster_vlnplot1,file="./MOU_scRNA3_Vlnplot_QC1_Sample.jpg",width=10,height=5,dpi=600)


