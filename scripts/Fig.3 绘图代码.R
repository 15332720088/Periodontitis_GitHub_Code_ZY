# 用于绘图的Seurat文件
load("/public/home/yinqi/zhaoyi/periodontitis/Review_20250303/Review文章/B_cell_Picture_20250326.Rdata")

setwd(dir = "/public/home/yinqi/zhaoyi/periodontitis/Review_20250303/Review文章")
getwd()

## Fig.3A ##################################################
colors <- c("#103667", "#1B4F93","#426EB4","#635BA2", "#8273B0", "#A095C4","#7388d0" ,
            "#00676B", "#008489", "#00B2BF","#6EC3C9","#99D1D3")
Idents(B_cell) <- "B_celltypes"
UMAP_plot2 = DimPlot(B_cell, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./UMAP4_B_cell_celltypes.jpg",width=8,height=7,dpi=600)


## Fig.3B ##################################################
colors <- c("#B7B2D0", "#7DA6C6", "#EAAA60", "#E68B81","transparent")
Idents(B_cell) <- "origin"
UMAP_plot2 = DimPlot(B_cell, reduction = "umap", label=F, cols = colors, raster=FALSE, pt.size = 2) 
ggsave(UMAP_plot2,file="./UMAP4_B_cell_origin.jpg",width=8,height=7,dpi=600)


## 绘制细胞比例条形图 ##################################################
# 鼠B_celltypes
Idents(B_cell) <- "origin"
MOU_B_cell=subset(B_cell,idents=c("MOU-Ctrl","MOU-Lig"))

colors <- c("#103667", "#1B4F93","#426EB4","#635BA2", "#8273B0", "#A095C4","#7388d0" ,
            "#00676B", "#008489", "#00B2BF","#6EC3C9","#99D1D3")
group_percent <- with(MOU_B_cell@meta.data, prop.table(table(B_celltypes, project), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$B_celltypes <- as.character(group_percent$B_celltypes)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$B_celltypes <- as.factor(c(group_percent$B_celltypes))
group_percent$B_celltypes <- fct_relevel(group_percent$B_celltypes, c("Naive B","Activated B","Transitional B","Conventional Memory B",
                                                                      "CD11c+CD27- Memory B","Atypical Memory B","Breg",
                                                                      "Plasma celltype 1","Plasma celltype 2","Plasma celltype 3",
                                                                      "Plasma celltype 4","Plasma celltype 5"))
group_percent$project <- fct_relevel(group_percent$project, c("Control","PD"))
group_percent$B_celltypes
group_percent$project
group_percent <- group_percent[group_percent$B_celltypes != "x", ] #删除x行
group_percent
Bar_celltype <- ggplot(group_percent, aes(x = project, y = Freq, fill = B_celltypes)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_MOU_B_cell_celltypes_percentage.jpg",width = 4, height=4.5, dpi=600)


# 人B_celltypes
Idents(B_cell) <- "origin"
HU_B_cell=subset(B_cell,idents=c("HU-Ctrl","HU-PD"))

colors <- c("#103667", "#1B4F93","#426EB4","#635BA2", "#8273B0", "#A095C4","#7388d0" ,
            "#00676B", "#008489", "#00B2BF","#6EC3C9","#99D1D3")
group_percent <- with(HU_B_cell@meta.data, prop.table(table(B_celltypes, project), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$B_celltypes <- as.character(group_percent$B_celltypes)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$B_celltypes <- as.factor(c(group_percent$B_celltypes))
group_percent$B_celltypes <- fct_relevel(group_percent$B_celltypes, c("Naive B","Activated B","Transitional B","Conventional Memory B",
                                                                      "CD11c+CD27- Memory B","Atypical Memory B","Breg",
                                                                      "Plasma celltype 1","Plasma celltype 2","Plasma celltype 3",
                                                                      "Plasma celltype 4","Plasma celltype 5"))
group_percent$project <- fct_relevel(group_percent$project, c("Control","PD"))
group_percent$B_celltypes
group_percent$project
group_percent <- group_percent[group_percent$B_celltypes != "x", ] #删除x行
group_percent
Bar_celltype <- ggplot(group_percent, aes(x = project, y = Freq, fill = B_celltypes)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_HU_B_cell_celltypes_percentage.jpg",width = 4, height=4.5, dpi=600)


# origin
Idents(B_cell) <- "origin"
B_cell2=subset(B_cell,idents=c("HU-Ctrl","HU-PD", "MOU-Ctrl","MOU-Lig"))

origin_data <- data.frame(Group = B_cell2@meta.data$origin)
origin_data_summary <- origin_data %>% group_by(Group) %>% summarise(Freq = n()) %>% mutate(Freq = Freq / sum(Freq))
origin_data_summary$ALL <- "ALL"
origin_data_summary$Group <- fct_relevel(origin_data_summary$Group, c("MOU-Ctrl","MOU-Lig","HU-Ctrl","HU-PD"))
ggplot(origin_data_summary, aes(x = ALL, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Proportion", fill = " ") + 
  scale_fill_manual(values = c("HU-Ctrl" = "#B7B2D0", "HU-PD" = "#7DA6C6",
                               "MOU-Ctrl" = "#EAAA60", "MOU-Lig" = "#E68B81")) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
# 保存图像
ggsave("Bar_B_cell_origin_percentage.jpg", width = 2.5, height = 4.5, dpi = 600)


table(B_cell2$origin)


## Fig.3C-E ######################################################
load("./B_cell_celltypes_20250314.Rdata") 

library(monocle)
# 提取表达矩阵、细胞信息和特征信息
expr_matrix <- as(as.matrix(B_cell@assays$RNA@data), "sparseMatrix")
cell_metadata <- B_cell@meta.data
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
save(mono, file="./Monocle2_B_cell_20250318.Rdata")
head(mono@phenoData@data)

load("./Monocle2_B_cell_20250318.Rdata")
trajectory_plot <- plot_cell_trajectory(mono, color_by = "Pseudotime", cell_size = 1, show_branch_points = T) + scale_colour_viridis(option = "inferno")
ggsave(trajectory_plot,file="./Trajectory_plot_B_cell_Pseudotime.jpg",width=5,height=5,dpi=600)

colors1 <- c("#103667", "#1B4F93","#426EB4","#635BA2", "#8273B0", "#A095C4","#7388d0" ,
             "#00676B", "#008489", "#00A6AD","#00B2BF","#6EC3C9","#99D1D3","#99D3F0")
trajectory_plot <- plot_cell_trajectory(mono, color_by = "B_celltypes", cell_size = 1, show_branch_points = T)+  
  scale_colour_manual(values = colors1)
ggsave(trajectory_plot,file="./Trajectory_plot_B_cell_celltype.jpg",width=5,height=5,dpi=600)

colors2 <- c("#B7B2D0", "#7DA6C6", "#EAAA60", "#E68B81","transparent")
trajectory_plot <- plot_cell_trajectory(mono, color_by = "origin", cell_size = 1, show_branch_points = T) +
  scale_colour_manual(values = colors2)
ggsave(trajectory_plot,file="./Trajectory_plot_B_cell_origin.jpg",width=5,height=5,dpi=600)


## Fig.3F ######################################################
## 细胞密度图
bar_data <- mono@phenoData@data[, c("Pseudotime", "origin", "B_celltypes")]
bar_data_order <- bar_data[order(bar_data$Pseudotime), ] #数据框按Pseudotime列由小到大排序

table(B_cell$B_celltypes)

colors1 <- c("#103667", "#1B4F93","#426EB4","#635BA2", "#8273B0", "#A095C4","#7388d0" ,
             "#00676B", "#008489", "#00A6AD","#00B2BF","#6EC3C9","#99D1D3","#99D3F0")
gg <- ggplot(bar_data_order, aes(Pseudotime, colour = B_celltypes, fill=B_celltypes)) +
  geom_density(bw=0.5,size=1,alpha = 0.5) +
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(gg,file="./Density3_plot_B_cell_celltypes.jpg",width=8,height=4,dpi=600)


bar_data_order <- bar_data_order %>% filter(origin != "MOU-LigPg")
colors1 <- c("#B7B2D0", "#7DA6C6", "#EAAA60", "#E68B81","transparent")
trajectory_plot <- plot_cell_trajectory(mono, color_by = "origin", cell_size = 1, show_branch_points = T) +
  scale_colour_manual(values = colors2)
ggsave(trajectory_plot,file="./Trajectory_plot_B_cell_origin.jpg",width=5,height=5,dpi=600)

gg <- ggplot(bar_data_order, aes(Pseudotime, colour = origin, fill=origin)) +
  geom_density(bw=0.5,size=1,alpha = 0.5) +
  scale_fill_manual(name = "", values = colors1) +
  scale_color_manual(name= "", values = colors1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(gg,file="./Density3_plot_B_cell_origin.jpg",width=8,height=4,dpi=600)


## BEAM算法绘制与分支高度相关的基因热图
load("/public/home/yinqi/zhaoyi/periodontitis/MOU_HU_Fig_20240902_now/Monocle2_B_cell_20241031.Rdata")
beam_results <- BEAM(mono, branch_point = 1)
save(beam_results, file="./B_cell_beam_results.Rdata")
gene_to_cluster = beam_results %>% arrange(qval) %>% head(100) %>% pull(gene_short_name) # 选出用于绘图的top80基因
trajectory_plot <- plot_pseudotime_heatmap(mono[gene_to_cluster,],
                                           num_clusters = 3, 
                                           show_rownames = TRUE,
                                           return_heatmap = TRUE,
                                           hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
ggsave(trajectory_plot,file="./Heatmap3_B_cell_genes.jpg",width=8,height=5,dpi=600)
## 调整目的基因排列顺序
gene_to_cluster = c("LYN","CD24","MS4A1","SELL","CD83","CXCR4","CD74",
                    "HLA-DRB5","GLTSCR2","FTL","CYBA","PDIA6","TMEM59","MZB1",
                    "IGKC","IGHG1","IGHG2","JCHAIN")
trajectory_plot <- plot_pseudotime_heatmap(mono[gene_to_cluster,], cluster_rows = FALSE,
                                           num_clusters = 3, 
                                           show_rownames = TRUE,
                                           return_heatmap = TRUE,
                                           hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
ggsave(trajectory_plot,file="./Heatmap3_B_cell_genes.jpg",width=8,height=3,dpi=600)



## Fig.3G-H ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step8_B_cell_HU_MOU_celltypes_20241031.Rdata")
colors <- c("#B7B2D0", "#7DA6C6","#EAAA60", "#E68B81", "#84C3B7")
VlnPlot <- VlnPlot(B_cell, features = "BTG1", group.by = "origin", pt.size = 0, combine = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black",alpha = 0.5) + #箱式图
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        legend.position = "right", legend.title = element_blank())
ggsave(VlnPlot, file = "./VlnPlot_B_cell_TFs_BTG1.jpg",width = 5, height=5, dpi=600)

colors <- c("#B7B2D0", "#7DA6C6", "#EAAA60", "#E68B81","#84C3B7")
VlnPlot <- VlnPlot(B_cell, features = "XBP1", group.by = "origin", pt.size = 0, combine = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black",alpha = 0.5) + #箱式图
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        legend.position = "right", legend.title = element_blank())
ggsave(VlnPlot, file = "./VlnPlot_B_cell_TFs_XBP1.jpg",width = 5, height=5, dpi=600)


# 计算差异基因表达是否有统计学意义
## HU
Idents(B_cell) <- "origin"
HU_B_cell=subset(B_cell,idents=c("HU-Ctrl","HU-PD"))

expression <- HU_B_cell@assays$RNA@data["BTG1", ]
group_info <- HU_B_cell@meta.data$origin
gene_data <- data.frame(Expression = expression, Group = group_info)
result_HU_BTG1 <- wilcox.test(Expression ~ Group, data = gene_data)
result_HU_BTG1

expression <- HU_B_cell@assays$RNA@data["XBP1", ]
group_info <- HU_B_cell@meta.data$origin
gene_data <- data.frame(Expression = expression, Group = group_info)
result_HU_XBP1 <- wilcox.test(Expression ~ Group, data = gene_data)
result_HU_XBP1


## MOU
Idents(B_cell) <- "origin"
MOU_B_cell=subset(B_cell,idents=c("MOU-Ctrl","MOU-Lig"))

expression <- MOU_B_cell@assays$RNA@data["BTG1", ]
group_info <- MOU_B_cell@meta.data$origin
gene_data <- data.frame(Expression = expression, Group = group_info)
result_MOU_BTG1 <- wilcox.test(Expression ~ Group, data = gene_data)
result_MOU_BTG1

expression <- MOU_B_cell@assays$RNA@data["XBP1", ]
group_info <- MOU_B_cell@meta.data$origin
gene_data <- data.frame(Expression = expression, Group = group_info)
result_MOU_XBP1 <- wilcox.test(Expression ~ Group, data = gene_data)
result_MOU_XBP1







