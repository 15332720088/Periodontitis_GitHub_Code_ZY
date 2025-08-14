# 用于绘图的Seurat文件
load("/public/home/yinqi/zhaoyi/periodontitis/MOU_HU_Fig_20240902_now/Step8_Macro_HU_MOU_celltype_20241101.Rdata")

load("/public/home/yinqi/zhaoyi/periodontitis/Review_20250303/Review文章/Macro_celltypes_20250321.Rdata") 
setwd(dir = "/public/home/yinqi/zhaoyi/periodontitis/Review_20250303/Review文章")
getwd()


## Fig.4A-B ######################
## 绘制UMAP图
colors <- c("#38044B", "#511985","#79378B",
            "#50013c", "#78004d", "#a20063","#AF4A92", "#C57CAC", "#D2A6C7")
DefaultAssay(Macro) <- "integrated" 
Idents(Macro) <- "Macro_types"
UMAP_plot2 = DimPlot(Macro, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./Macro_UMAP_celltypes_20250321.jpg",width=8,height=7,dpi=600)

colors <- c("#8EA0CC", "#377BAB", "#CCB4D7","#956A88","transparent")
DefaultAssay(Macro) <- "integrated" 
Idents(Macro) <- "origin"
UMAP_plot2 = DimPlot(Macro, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./Macro_UMAP_origin_20250321.jpg",width=8,height=7,dpi=600)


## 绘制细胞比例条形图
# 鼠Macro_types
Idents(Macro) <- "origin"
MOU_Macro=subset(Macro,idents=c("MOU_Ctrl","MOU_Lig"))

colors <- c("#38044B", "#511985","#79378B",
            "#50013c", "#78004d", "#a20063","#AF4A92", "#C57CAC", "#D2A6C7")
group_percent <- with(MOU_Macro@meta.data, prop.table(table(Macro_types, project), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$Macro_types <- as.character(group_percent$Macro_types)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$Macro_types <- as.factor(c(group_percent$Macro_types))
group_percent$Macro_types <- fct_relevel(group_percent$Macro_types, c("Classical Monocytes", "Intermediate Monocytes","Non-classical Monocytes",
                                                                      "M1 Macro 1","M1 Macro 2",
                                                                      "M2 Macro 1","M2 Macro 2","M2 Macro 3","M2 Macro 4"))
group_percent$project <- fct_relevel(group_percent$project, c("Control","PD"))
group_percent$Macro_types
group_percent$project
group_percent <- group_percent[group_percent$Macro_types != "x", ] #删除x行
group_percent
Bar_celltype <- ggplot(group_percent, aes(x = project, y = Freq, fill = Macro_types)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_MOU_Macro_celltypes_percentage.jpg",width = 4, height=4.5, dpi=600)


# 人Macro_types
Idents(Macro) <- "origin"
HU_Macro=subset(Macro,idents=c("HU_Ctrl","HU_PD"))

colors <- c("#38044B", "#511985","#79378B",
            "#50013c", "#78004d", "#a20063","#AF4A92", "#C57CAC", "#D2A6C7")
group_percent <- with(HU_Macro@meta.data, prop.table(table(Macro_types, project), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$Macro_types <- as.character(group_percent$Macro_types)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$Macro_types <- as.factor(c(group_percent$Macro_types))
group_percent$Macro_types <- fct_relevel(group_percent$Macro_types, c("Classical Monocytes", "Intermediate Monocytes","Non-classical Monocytes",
                                                                      "M1 Macro 1","M1 Macro 2",
                                                                      "M2 Macro 1","M2 Macro 2","M2 Macro 3","M2 Macro 4"))
group_percent$project <- fct_relevel(group_percent$project, c("Control","PD"))
group_percent$Macro_types
group_percent$project
group_percent <- group_percent[group_percent$Macro_types != "x", ] #删除x行
group_percent
Bar_celltype <- ggplot(group_percent, aes(x = project, y = Freq, fill = Macro_types)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_HU_Macro_celltypes_percentage.jpg",width = 4, height=4.5, dpi=600)


# origin
Idents(Macro) <- "origin"
Macro2=subset(Macro,idents=c("MOU_Ctrl","MOU_Lig", "HU_Ctrl","HU_PD"))

origin_data <- data.frame(Group = Macro2@meta.data$origin)
origin_data_summary <- origin_data %>% group_by(Group) %>% summarise(Freq = n()) %>% mutate(Freq = Freq / sum(Freq))
origin_data_summary$ALL <- "ALL"
origin_data_summary$Group <- fct_relevel(origin_data_summary$Group, c("MOU_Ctrl","MOU_Lig","HU_Ctrl","HU_PD"))
ggplot(origin_data_summary, aes(x = ALL, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Proportion", fill = " ") + 
  scale_fill_manual(values = c("HU_Ctrl" = "#B7B2D0", "HU_PD" = "#7DA6C6",
                               "MOU_Ctrl" = "#EAAA60", "MOU_Lig" = "#E68B81")) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
# 保存图像
ggsave("Bar_Macro_origin_percentage.jpg", width = 2.5, height = 4.5, dpi = 600)

table(Macro2$origin)
# HU_Ctrl  HU_PD   MOU_Ctrl   MOU_Lig  MOU_LigPg 
# 199      1070    59         1719      0




## Fig.4C-F ######################
library(monocle)
# 提取表达矩阵、细胞信息和特征信息
expr_matrix <- as(as.matrix(Macro@assays$integrated@data), "sparseMatrix")
cell_metadata <- Macro@meta.data
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
genes <- c("CD14","CCR2", # Classical Monocytes 1
           "CCR5","HLA-DQB2","NRROS","NR4A1", # Intermediate Monocytes
           "CX3CR1","ITGAX","IL1RN","SERPINB9","NR4A3", # Non-classical Monocytes
           "CD80","CD86", # M1
           "DENND4A","IRAK2","CCRL2", # M1 Macro 1
           "TNF","ALCAM","CLEC4D","ARHGAP10", # M1 Macro 2
           "CD163","MRC1", # M2
           "RETNLB","FOLR3","LYVE1","CCL24", # M2 Macro 1
           "RDH11","CREG1","NCOA4","DUSP6") # M2 Macro 4
mono <- setOrderingFilter(mono, genes)
# 降维和拟时序分析
mono <- reduceDimension(mono, max_components = 2, method = 'DDRTree') #用时很长
mono <- orderCells(mono) #用时很长
head(mono@phenoData@data)
save(mono, file="./Monocle2_Macro_20250322.Rdata")

####
load("./Monocle2_Macro_20250322.Rdata")

## 绘制分支图
trajectory_plot <- plot_cell_trajectory(mono, color_by = "Pseudotime", cell_size = 1, show_branch_points = T) + scale_colour_viridis(option = "inferno")
ggsave(trajectory_plot,file="./Trajectory_plot_Macro_Pseudotime.jpg",width=5,height=5,dpi=600)

trajectory_plot <- plot_cell_trajectory(mono, color_by = "State", cell_size = 1, show_branch_points = T)
ggsave(trajectory_plot,file="./Trajectory_plot_Macro_State.jpg",width=5,height=5,dpi=600)

colors <- c("#38044B", "#511985","#79378B",
            "#50013c", "#78004d", "#a20063","#AF4A92", "#C57CAC", "#D2A6C7")
trajectory_plot <- plot_cell_trajectory(mono, color_by = "Macro_types", cell_size = 1, show_branch_points = T) +
  scale_colour_manual(values = colors)
ggsave(trajectory_plot,file="./Trajectory_plot_Macro_celltypes.jpg",width=5,height=5,dpi=600)


## plot_cell_trajectory()函数无法更改点的覆盖顺序
## 用ggplot分层绘图

# 人
df <- as.data.frame(pData(mono))  # 提取细胞元数据
df$Dim1 <- reducedDimS(mono)[1, ]  # 第一维
df$Dim2 <- reducedDimS(mono)[2, ]  # 第二维
df1 <- df[df$origin == "HU_PD", ]   # 低层（粉色点）
df2 <- df[df$origin == "HU_Ctrl", ] # 高层（蓝色点）

colors <- c("transparent", "transparent", "transparent","transparent","transparent")
trajectory_plot <- plot_cell_trajectory(mono, show_branch_points = T) +
  scale_colour_manual(values = colors) +  # 原图像仅保留trajectory分支，点都为透明
  geom_point(data = df1, aes(x = Dim1, y = Dim2, color = origin), 
             color = "#F19373", size = 1, alpha = 1) +  # 先绘制df1
  geom_point(data = df2, aes(x = Dim1, y = Dim2, color = origin), 
             color = "#6EC3C9", size = 1, alpha = 1)  # 再绘制df2
ggsave(trajectory_plot,file="./Trajectory2_plot_HU_Macro_origin.jpg",width=5,height=5,dpi=600)

# 鼠
df <- as.data.frame(pData(mono))  # 提取细胞元数据
df$Dim1 <- reducedDimS(mono)[1, ]  # 第一维
df$Dim2 <- reducedDimS(mono)[2, ]  # 第二维
df1 <- df[df$origin == "MOU_Lig", ]   # 低层（粉色点）
df2 <- df[df$origin == "MOU_Ctrl", ] # 高层（蓝色点）

colors <- c("transparent", "transparent", "transparent","transparent","transparent")
trajectory_plot <- plot_cell_trajectory(mono, show_branch_points = T) +
  scale_colour_manual(values = colors) +  # 原图像仅保留trajectory分支，点都为透明
  geom_point(data = df1, aes(x = Dim1, y = Dim2, color = origin), 
             color = "#F19373", size = 0.8, alpha = 1) +  # 先绘制df1
  geom_point(data = df2, aes(x = Dim1, y = Dim2, color = origin), 
             color = "#6EC3C9", size = 1, alpha = 1)  # 再绘制df2
ggsave(trajectory_plot,file="./Trajectory2_plot_MOU_Macro_origin.jpg",width=5,height=5,dpi=600)



## 绘制细胞密度图（Macro_types）
bar_data <- mono@phenoData@data[, c("Pseudotime", "origin", "Macro_types")]
bar_data_order <- bar_data[order(bar_data$Pseudotime), ] #数据框按Pseudotime列由小到大排序

table(Macro$Macro_types)

colors1 <- c("#38044B", "#511985","#79378B",
             "#50013c", "#78004d", "#a20063","#AF4A92", "#C57CAC", "#D2A6C7")
gg <- ggplot(bar_data_order, aes(Pseudotime, colour = Macro_types, fill=Macro_types)) +
  geom_density(bw=0.07, size=1, alpha = 0.5) +  # 调整bw值可以将挤在一起的平滑波峰分散开
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  scale_y_continuous(trans = scales::log1p_trans()) +  # 让高的峰更高，小的峰更低
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", colour = "white"))
ggsave(gg, file="./Density2_plot_Macro_celltypes.jpg", width=10, height=4, dpi=600)



## 绘制细胞密度图（人origin）
bar_data <- mono@phenoData@data[, c("Pseudotime", "origin", "Macro_types")]
bar_data_order <- bar_data[order(bar_data$Pseudotime), ] #数据框按Pseudotime列由小到大排序

# 增加高密度细胞区域权重，使波峰更易清晰
bar_data_order$weight <- 1  
bar_data_order$weight[bar_data_order$Pseudotime > 0 & bar_data_order$Pseudotime < 0.5] <- 2.5
bar_data_order$weight[bar_data_order$Pseudotime > 2 & bar_data_order$Pseudotime < 2.3] <- 2.5

bar_data_order <- bar_data_order %>% filter(!origin %in% c("MOU_Ctrl","MOU_Lig","MOU_LigPg"))
colors1 <- c("#6EC3C9","#F19373", "transparent","transparent","transparent")
gg <- ggplot(bar_data_order, aes(Pseudotime, colour = origin, fill=origin, weight=weight)) +
  geom_density(bw=0.15,size=1,alpha = 0.5) +  # 调整bw值可以将挤在一起的平滑波峰分散开
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(gg,file="./Density2_plot_HU_Macro_origin.jpg",width=10,height=4,dpi=600)


## 绘制细胞密度图（鼠origin）
bar_data <- mono@phenoData@data[, c("Pseudotime", "origin", "Macro_types")]
bar_data_order <- bar_data[order(bar_data$Pseudotime), ] #数据框按Pseudotime列由小到大排序

# 增加高密度细胞区域权重，使波峰更易清晰
bar_data_order$weight <- 1  
bar_data_order$weight[bar_data_order$Pseudotime > 0 & bar_data_order$Pseudotime < 0.5] <- 2.5
bar_data_order$weight[bar_data_order$Pseudotime > 2 & bar_data_order$Pseudotime < 2.3] <- 2.5

bar_data_order <- bar_data_order %>% filter(origin != "MOU_LigPg")
colors1 <- c("transparent", "transparent", "#6EC3C9","#F19373","transparent")
gg <- ggplot(bar_data_order, aes(Pseudotime, colour = origin, fill=origin, weight=weight)) +
  geom_density(bw=0.15,size=1,alpha = 0.5) +  # 调整bw值可以将挤在一起的平滑波峰分散开
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(gg,file="./Density2_plot_MOU_Macro_origin.jpg",width=10,height=4,dpi=600)




## Fig.4G-H ######################
##恢复过滤基因
load("./Macro_celltypes_20250321.Rdata")
load("./Monocle2_Macro_20250322.Rdata")

# 提取Macro中的表达矩阵，重新构建exprs(mono)
new_expr_matrix <- as(as.matrix(Macro@assays$RNA@data), "sparseMatrix")
all_genes <- union(rownames(new_expr_matrix), rownames(mono))
new_expr_matrix_full <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(mono), sparse = TRUE)
rownames(new_expr_matrix_full) <- all_genes
colnames(new_expr_matrix_full) <- colnames(mono)
new_expr_matrix_full[rownames(new_expr_matrix), ] <- new_expr_matrix

dim(new_expr_matrix_full)
[1] 16542  3560
dim(exprs(mono))
[1] 5417 3560

# 重新构建和new_expr_matrix_full基因名对应的基因信息slot：featureData(mono)
new_gene_annotation <- data.frame(gene_short_name = rownames(new_expr_matrix_full))
rownames(new_gene_annotation) <- rownames(new_expr_matrix_full)
new_fd <- new("AnnotatedDataFrame", data = new_gene_annotation)

## 整合为mono对象
mono_new <- newCellDataSet(
  new_expr_matrix_full,
  phenoData = phenoData(mono),  # 细胞信息保持不变
  featureData = new_fd,  # 更新基因信息
  lowerDetectionLimit = mono@lowerDetectionLimit,
  expressionFamily = negbinomial.size()
)
mono_new <- estimateSizeFactors(mono_new)
mono_new <- estimateDispersions(mono_new) 

dim(exprs(mono_new))
[1] 16542  3560

save(mono_new, file="./Monocle2_Macro_mono_new_heatmap_20250329.Rdata")


## 绘制热图
load("./Monocle2_Macro_mono_new_heatmap_20250329.Rdata")
# HU
gene_to_cluster <- c("MT-ND6","HUS1","RACK1","ATP5F1E","SELENOK",
                     "ATP5MC2","ATP5MC3","NDUFS7","TENT5A","SPARC",
                     "HMOX1","CD86","CXCL9","ATP5F1","PLA2G2D")
trajectory_plot <- plot_pseudotime_heatmap(mono_new[gene_to_cluster,],
                                           num_clusters = 3, 
                                           show_rownames = TRUE,
                                           return_heatmap = TRUE,
                                           hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
ggsave(trajectory_plot,file="./Heatmap_HU_Macro_monocle_genes.jpg",width=5,height=3,dpi=600)

# MOU
gene_to_cluster <- c("SELL","VCAN","EREG","FCN1","ITGAX",
                     "MAFB","C1QB","C1QA","C1QC","CTSD",
                     "CXCL12","PLA2G2D","CD163L1","CTSL","MS4A7")
trajectory_plot <- plot_pseudotime_heatmap(mono_new[gene_to_cluster,],
                                           num_clusters = 3, 
                                           show_rownames = TRUE,
                                           return_heatmap = TRUE,
                                           hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
ggsave(trajectory_plot,file="./Heatmap_MOU_Macro_monocle_genes.jpg",width=5,height=3,dpi=600)


