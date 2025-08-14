# 用于绘图的Seurat文件
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step4_HU_Immune_celltype_20240907.Rdata") 
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_Immune1_Celltype_20240907.Rdata") # Control&PD


## Fig.2A ######################
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
Idents(MOU_Immune1) <- "Immune_celltype"
UMAP_plot2 = DimPlot(MOU_Immune1, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./MOU_Immune1_UMAP.plot2_clusters.jpg",width=5,height=5,dpi=600)


## Fig.2B ######################
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
#绘制细胞种类百分比 Excel表及条形图(各组)
project_percent <- with(MOU_Immune1@meta.data, prop.table(table(Immune_celltype, project), margin = 2) * 100)
project_percent <- as.data.frame(project_percent)
project_percent$Immune_celltype <- as.character(project_percent$Immune_celltype)
project_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project_percent$Immune_celltype <- as.factor(c(project_percent$Immune_celltype))
project_percent$Immune_celltype <- fct_relevel(project_percent$Immune_celltype, c("αβT_CD4","αβT_CD8","αβT_DP","αβT_DN","γδT","B/Plasma","DC","Macrophage","Neutrophil","Mast"))
project_percent$project <- fct_relevel(project_percent$project, c("Control","PD"))
project_percent$Immune_celltype
project_percent$project
project_percent <- project_percent[project_percent$Immune_celltype != "x", ] #删除x行
project_percent
write.table(project_percent, "./MOU_Immune1_celltype_Bar2_percentage_projects.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project_percent, aes(x = project, y = Freq, fill = Immune_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./MOU_Immune1_celltype_Bar2_percentage_projects.jpg",width = 3.5, height=4.5, dpi=600)


## Fig.2C ######################
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
Idents(HU_Immune) <- "Immune_celltype"
UMAP_plot2 = DimPlot(HU_Immune, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./HU_Immune_UMAP.plot2_clusters.jpg",width=5,height=5,dpi=600)


## Fig.2D ######################
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
#绘制细胞种类百分比 Excel表及条形图(各组)
project_percent <- with(HU_Immune@meta.data, prop.table(table(Immune_celltype, project), margin = 2) * 100)
project_percent <- as.data.frame(project_percent)
project_percent$Immune_celltype <- as.character(project_percent$Immune_celltype)
project_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project_percent$Immune_celltype <- as.factor(c(project_percent$Immune_celltype))
project_percent$Immune_celltype <- fct_relevel(project_percent$Immune_celltype, c("αβT_CD4","αβT_CD8","αβT_DP","αβT_DN","γδT","B/Plasma","DC","Macrophage","Neutrophil","Mast"))
project_percent$project <- fct_relevel(project_percent$project, c("Control","PD"))
project_percent$Immune_celltype
project_percent$project
project_percent <- project_percent[project_percent$Immune_celltype != "x", ] #删除x行
project_percent
write.table(project_percent, "./HU_Immune_celltype_Bar2_percentage_projects.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project_percent, aes(x = project, y = Freq, fill = Immune_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./HU_Immune_celltype_Bar2_percentage_projects.jpg",width = 3.5, height=4.5, dpi=600)




