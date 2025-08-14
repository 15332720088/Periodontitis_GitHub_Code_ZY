# 用于绘图的Seurat文件
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_scRNA3_general_celltype_20250104.Rdata") # Lig$LigPg
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_Immune2_Celltype_20250105.Rdata") # Control&Lig$LigPg
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_Immune3_Immune_celltype_20250104.Rdata") # Lig$LigPg


## Fig.5A ######################
MOU_scRNA3$general_celltype <- as.factor(MOU_scRNA3$general_celltype)
MOU_scRNA3$general_celltype <- fct_relevel(MOU_scRNA3$general_celltype, c("Fibroblast","Epithelial","Endothelial","Immune","Myocyte","Schwann_Neuron"))
levels(MOU_scRNA3$general_celltype)

colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
Idents(MOU_scRNA3) <- "general_celltype"
UMAP_plot2 = DimPlot(MOU_scRNA3, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./UMAP_MOU_scRNA3_general_celltype.jpg",width=5,height=3,dpi=600)


## Fig.5B ######################
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
#绘制细胞种类百分比 Excel表及条形图(各组)
group_percent <- with(MOU_scRNA3@meta.data, prop.table(table(general_celltype, project2), margin = 2) * 100)
group_percent <- as.data.frame(group_percent)
group_percent$general_celltype <- as.character(group_percent$general_celltype)
group_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
group_percent$general_celltype <- as.factor(c(group_percent$general_celltype))
group_percent$general_celltype <- fct_relevel(group_percent$general_celltype, c("Fibroblast","Epithelial","Endothelial","Immune", "Myocyte","Schwann_Neuron"))
group_percent$Group <- fct_relevel(group_percent$Group, c("Lig","LigPg"))
group_percent$general_celltype
group_percent$Group
group_percent <- group_percent[group_percent$general_celltype != "x", ] #删除x行
group_percent
write.table(group_percent, "./Bar_MOU_scRNA3_general_celltype_Propotion_project2.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(group_percent, aes(x = project2, y = Freq, fill = general_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./Bar_MOU_scRNA3_general_celltype_Propotion_project2.jpg",width = 3.5, height=4.5, dpi=600)


## Fig.5C ######################
MOU_Immune3$Immune_celltype <- as.factor(MOU_Immune3$Immune_celltype)
MOU_Immune3$Immune_celltype <- fct_relevel(MOU_Immune3$Immune_celltype, c("αβT_CD4","αβT_CD8","αβT_DP","αβT_DN","γδT","B/Plasma","DC","Macrophage","Neutrophil","Mast"))
levels(MOU_Immune3$Immune_celltype)

colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
Idents(MOU_Immune3) <- "Immune_celltype"
UMAP_plot2 = DimPlot(MOU_Immune3, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(UMAP_plot2,file="./UMAP_MOU_Immune3_Immune_celltype.jpg",width=5,height=3,dpi=600)


## Fig.5D ######################
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
#绘制细胞种类百分比 Excel表及条形图(各组)
project2_percent <- with(MOU_Immune3@meta.data, prop.table(table(Immune_celltype, project2), margin = 2) * 100)
project2_percent <- as.data.frame(project2_percent)
project2_percent$Immune_celltype <- as.character(project2_percent$Immune_celltype)
project2_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project2_percent$Immune_celltype <- as.factor(c(project2_percent$Immune_celltype))
project2_percent$Immune_celltype <- fct_relevel(project2_percent$Immune_celltype, c("αβT_CD4","αβT_CD8","αβT_DP","αβT_DN","γδT","B/Plasma","DC","Macrophage","Neutrophil","Mast"))
project2_percent$project2 <- fct_relevel(project2_percent$project2, c("Control","PD"))
project2_percent$Immune_celltype
project2_percent$project2
project2_percent <- project2_percent[project2_percent$Immune_celltype != "x", ] #删除x行
project2_percent
write.table(project2_percent, "./MOU_Immune3_celltype_Bar2_percentage_project2s.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project2_percent, aes(x = project2, y = Freq, fill = Immune_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./MOU_Immune3_celltype_Bar2_percentage_project2s.jpg",width = 3.5, height=4.5, dpi=600)



## Fig.5H ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_Immune2_Celltype_20250105.Rdata") # Control&Lig$LigPg
##总的Immune
## Lig&Ctrl
# (小鼠)找出两组间差异基因（只比较指定两组）
Idents(MOU_Immune2) <- "project2"
DEGs1 = FindMarkers(MOU_Immune2, ident.1 = "Lig", ident.2 = c("Control"), logfc.threshold = 0, only.pos = FALSE) #（only.pos = TRUE，只返回阳性结果）
save(DEGs1, file="./MOU_Immune2_Lig&Control_DEGs1.Rdata")
write.table(DEGs1, "./MOU_Immune2_Lig&Control_DEGs1.xls", col.names = T, row.names = F, sep = "\t")

#命名基因列名，计算表达分布Difference
DEGs1 <- DEGs1 %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2)
#根据阈值筛选差异显著的上下调基因，与差异不显著的基因
DEGs1$Sig_group = "no-significant"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC > 0.25))] = "up-regulated"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC < -0.25))] = "down-regulated"
table(DEGs1$Sig_group)
DEGs_up_Lig <- DEGs1[DEGs1$Sig_group == "up-regulated", ]
table(DEGs_up_Lig$Sig_group)
save(DEGs_up_Lig, file="./MOU_Immune2_Lig&Control_DEGs_up_Lig.Rdata")
DEGs_down_Lig <- DEGs1[DEGs1$Sig_group == "down-regulated", ]
table(DEGs_down_Lig$Sig_group)
save(DEGs_down_Lig, file="./MOU_Immune2_Lig&Control_DEGs_down_Lig.Rdata")


## LigPg&Ctrl
# (小鼠)找出两组间差异基因（只比较指定两组）
Idents(MOU_Immune2) <- "project2"
DEGs1 = FindMarkers(MOU_Immune2, ident.1 = "LigPg", ident.2 = c("Control"), logfc.threshold = 0, only.pos = FALSE) #（only.pos = TRUE，只返回阳性结果）
save(DEGs1, file="./MOU_Immune2_LigPg&Control_DEGs1.Rdata")
write.table(DEGs1, "./MOU_Immune2_LigPg&Control_DEGs1.xls", col.names = T, row.names = F, sep = "\t")

#命名基因列名，计算表达分布Difference
DEGs1 <- DEGs1 %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2)
#根据阈值筛选差异显著的上下调基因，与差异不显著的基因
DEGs1$Sig_group = "no-significant"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC > 0.25))] = "up-regulated"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC < -0.25))] = "down-regulated"
table(DEGs1$Sig_group)
DEGs_up_LigPg <- DEGs1[DEGs1$Sig_group == "up-regulated", ]
table(DEGs_up_LigPg$Sig_group)
save(DEGs_up_LigPg, file="./MOU_Immune2_LigPg&Control_DEGs_up_LigPg.Rdata")
DEGs_down_LigPg <- DEGs1[DEGs1$Sig_group == "down-regulated", ]
table(DEGs_down_LigPg$Sig_group)
save(DEGs_down_LigPg, file="./MOU_Immune2_LigPg&Control_DEGs_down_LigPg.Rdata")


load("./MOU_Immune2_Lig&Control_DEGs_up_Lig.Rdata")
load("./MOU_Immune2_LigPg&Control_DEGs_up_LigPg.Rdata")
#（两交集）
venn_list <- list(Lig = DEGs_up_Lig$gene, LigPg = DEGs_up_LigPg$gene)
venn.diagram(venn_list, filename = 'Venn_Lig&LigPg_Immune_DEGs_up.tiff', imagetype = 'tiff', 
             fill = c("#103667","#8B0016"), # 圆圈内部充填颜色
             category.names = c("Lig", "LigPg"),
             alpha = 0.8, # 圆圈内部充填颜色的透明度
             label.col = "white", # 圈内基因数目标签的颜色
             cex = 2.2, # 圈内基因数目标签的字体大小
             col = "transparent", # 圆圈线的颜色
             lwd = 1, # 圆圈线的宽度
             cat.col = "black", # 圈外样本名标签的颜色
             cat.cex = 2.2, # 圈外样本名标签的字体大小
             cat.dist = c(0.04, 0.04), # 圈外样本名标签距圆圈的距离（几个标签几个距离）
             cat.pos = c(-225, 225), # 圈外样本名标位于圆圈的角度
             main = "Immune Up",  # 标题              
             main.col = "black",  # 标题颜色              
             main.cex = 2.2,  # 标题字体大小              
             ext.text = F, # 取消引导线注释              
             rotation.degree = 180, # 调整旋转角度              
             scaled = TRUE) # TRUE按比例绘制圆圈大小；FALSE圆圈大小一致

#导出韦恩图交集基因
inter <- get.venn.partitions(venn_list)
save(inter, file="./Venn_Lig&LigPg_Immune_DEGs_up_interGenes.Rdata")
str(inter)
# int  379 318 35

load("./MOU_Immune2_Lig&Control_DEGs_down_Lig.Rdata")
load("./MOU_Immune2_LigPg&Control_DEGs_down_LigPg.Rdata")
#（两交集）
venn_list <- list(Lig = DEGs_down_Lig$gene, LigPg = DEGs_down_LigPg$gene)
venn.diagram(venn_list, filename = 'Venn_Lig&LigPg_Immune_DEGs_down.tiff', imagetype = 'tiff', 
             fill = c("#103667","#8B0016"), # 圆圈内部充填颜色
             category.names = c("Lig", "LigPg"),
             alpha = 0.8, # 圆圈内部充填颜色的透明度
             label.col = "white", # 圈内基因数目标签的颜色
             cex = 2.2, # 圈内基因数目标签的字体大小
             col = "transparent", # 圆圈线的颜色
             lwd = 1, # 圆圈线的宽度
             cat.col = "black", # 圈外样本名标签的颜色
             cat.cex = 2.2, # 圈外样本名标签的字体大小
             cat.dist = c(0.04, 0.04), # 圈外样本名标签距圆圈的距离（几个标签几个距离）
             cat.pos = c(-225, 225), # 圈外样本名标位于圆圈的角度
             main = "Immune Down",  # 标题              
             main.col = "black",  # 标题颜色              
             main.cex = 2.2,  # 标题字体大小              
             ext.text = T, # F取消引导线注释              
             rotation.degree = 180, # 调整旋转角度              
             scaled = TRUE) # TRUE按比例绘制圆圈大小；FALSE圆圈大小一致
#导出韦恩图交集基因
inter <- get.venn.partitions(venn_list)
save(inter, file="./Venn_Lig&LigPg_Immune_DEGs_down_interGenes.Rdata")
str(inter)
# int  677 428 25

## Fig.5I ######################
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_Immune3_Immune_celltype_20250104.Rdata") # Lig$LigPg
###鼠Lig&Pg亚群差异基因对比
##找出各免疫亚群差异基因
immune_celltypes <- levels(MOU_Immune3$Immune_celltype)
celltype_subsets <- list() #建一个空列表存储每个细胞类型的子集对象
DEGs_subsets <- list() #建一个空列表存储每个亚群差异基因的子集对象
Idents(MOU_Immune3) <- "Immune_celltype"
for (celltype in immune_celltypes) {
  subset = subset(MOU_Immune3,idents = c(as.character(celltype)))
  celltype_subsets[[celltype]] <- subset
  Idents(subset) <- "project2"
  # 检查LigPg和Lig组中的细胞数目
  cells_LigPg <- WhichCells(subset, idents = "LigPg")
  cells_Lig <- WhichCells(subset, idents = "Lig")
  # 仅在LigPg和Lig组中的细胞数目都不少于3个时才执行FindMarkers操作
  if (length(cells_LigPg) >= 3 && length(cells_Lig) >= 3) {
    DEGs <- FindMarkers(subset, ident.1 = "LigPg", ident.2 = c("Lig"), logfc.threshold = 0, only.pos = FALSE)
    DEGs_subsets[[celltype]] <- DEGs
  } else {
    message(paste("Skipping cell type", celltype, "due to insufficient cell numbers in LigPg or Lig group"))
  }}
#去除列表中各数据框名字中的空格
new_names <- gsub(" ", "", names(DEGs_subsets))
names(DEGs_subsets) <- new_names
#显示并保存结果
lapply(DEGs_subsets, function(df) head(df, 5))
save(DEGs_subsets, file=("./MOU_Immune3_subset_LigPg&Lig_DEGs.Rdata"))

##调整及标注列表
DEGs_subsets2 <- list() #建一个空列表用于存储结果
for (name in names(DEGs_subsets)) { #对列表中的每个数据框执行操作
  df <- DEGs_subsets[[name]]
  df$group = as.character(name) #df$group中充填亚群名
  df <- df %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2) %>% filter(abs(Difference) > 0.1)
  df$P_group <- ""
  df$P_group[which((df$p_val_adj >= 0.05))] = "no significant"
  df$P_group[which((df$p_val_adj < 0.05 & df$avg_log2FC > 0))] = "up-regulated"
  df$P_group[which((df$p_val_adj < 0.05 & df$avg_log2FC < 0))] = "down-regulated"
  DEGs_subsets2[[name]] <- df
}
lapply(DEGs_subsets2, function(df) head(df, 5))
lapply(DEGs_subsets2, function(df) table(df$P_group))
# ！ 截图P_group细胞比例
save(DEGs_subsets2, file=("./MOU_Immune3_subset2_LigPg&Lig_DEGs.Rdata"))


##筛选图片中标注基因名的top基因
TopGene <- list() #建一个空列表用于存储结果
for (name in names(DEGs_subsets2)) { #对列表中的每个数据框执行操作
  df <- DEGs_subsets2[[name]]
  top5_log2FC_up = df %>% filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% top_n(n = 8, wt = abs(avg_log2FC))
  top5_log2FC_down = df %>% filter(p_val_adj < 0.05 & avg_log2FC < 0) %>% top_n(n = 8, wt = abs(avg_log2FC))
  df <- bind_rows(top5_log2FC_up, top5_log2FC_down)
  TopGene[[name]] <- df
}
TopGene <- bind_rows(TopGene) #合并列表中的所有数据框
TopGene
table(TopGene$group)
TopGene$gene
#调整TopGene$group顺序
TopGene$group <- as.factor(TopGene$group)
TopGene$group <- fct_relevel(TopGene$group, c("αβT_CD4","αβT_CD8","αβT_DP","αβT_DN","γδT","B/Plasma","DC","Macrophage","Neutrophil","Mast"))
levels(TopGene$group)
save(TopGene, file=("./MOU_Immune3_TopGene_LigPg&Lig_DEGs.Rdata"))

##背景柱状图数据准备
DEGs_subsets3 <- bind_rows(DEGs_subsets2) #合并列表中的所有数据框
DEGs_subsets3$group <- as.factor(DEGs_subsets3$group)
DEGs_subsets3$group <- fct_relevel(DEGs_subsets3$group, c("αβT_CD4","αβT_CD8","αβT_DP","αβT_DN","γδT","B/Plasma","DC","Macrophage","Neutrophil","Mast"))
levels(DEGs_subsets3$group)
dbar <- 
  DEGs_subsets3 %>% 
  group_by(group) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  select(group, avg_log2FC_min, avg_log2FC_max, P_group_min) %>% 
  rename(P_group = P_group_min)
dbar


## 有文字标注
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
ggplot()+
  geom_col(data = dbar,  # 绘制负向背景柱状图
           mapping = aes(x = group,y = avg_log2FC_min-0.3), fill = "#BFCAE6",
           alpha = 0.7, width = 0.9) + # alpha透明度，width柱子的宽度
  geom_col(data = dbar, # 绘制正向背景柱状图
           mapping = aes(x = group,y = avg_log2FC_max+0.3), fill = "#F5A89A",
           alpha = 0.7, width = 0.9) +
  geom_jitter(data = DEGs_subsets3, # 绘制所有数据点（geom_jitter对点的位置进行小幅度随机偏移以减少点的重叠）
              aes(x = group, y = avg_log2FC, color = P_group), 
              size = 0.15, width =0.38) + # 点的大小、偏移幅度
  scale_color_manual(values = c("up-regulated" = "#EB7153", "down-regulated" = "#205AA7", "no significant" = "#B7B7B7")) + #根据分类设置散点的颜色
  geom_jitter(data = TopGene, # 绘制TopGene数据点
              aes(x = group, y = avg_log2FC), color = "transparent", # 所有TopGene点的颜色
              size = 0.2, width =0.3) + 
  geom_text_repel(data = TopGene, # 绘制TopGene的标签
                  aes(x = group, y = avg_log2FC, label = gene),
                  size = 2, color = 'black', # 标签大小、颜色
                  segment.size = 0.2,  # 指示线粗细
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50), force = 0.6) + # default，最大重叠数; force, 排斥力，越高越推开标签防止重叠
  geom_tile(data = DEGs_subsets3, # 绘制中心分组名称条框
            aes(x = group, y = 0, fill = group), 
            color = "transparent", height=0.2, alpha = 0.6, show.legend = F) + #边框颜色、条框高度、透明度、不显示图例
  scale_fill_manual(values = colors) +  #设置各组名称条框内部充填颜色
  labs(x="General celltypes", y="Average log2FC") +
  geom_text(data=DEGs_subsets3, # 绘制中心分组名称文本
            aes(x=group,  y=0,  label=group),
            size = 3, color ="black") + # 文字大小及颜色
  theme_minimal() + 
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 13))
ggsave('MOU_Immune3_subset_LigPg&Lig_DEGs_lable.jpg', width = 10, height = 6, bg = 'white')


## 没有文字标注
colors <- c("#211551", "#2c1e69", "#392884","#625ba1","#8272b0", "#008489", "#C71585", "#D94D9C", "#EC87B3","#D2A6C7")
ggplot()+
  geom_col(data = dbar,  # 绘制负向背景柱状图
           mapping = aes(x = group,y = avg_log2FC_min-0.3), fill = "#BFCAE6",
           alpha = 0.7, width = 0.9) + # alpha透明度，width柱子的宽度
  geom_col(data = dbar, # 绘制正向背景柱状图
           mapping = aes(x = group,y = avg_log2FC_max+0.3), fill = "#F5A89A",
           alpha = 0.7, width = 0.9) +
  geom_jitter(data = DEGs_subsets3, # 绘制所有数据点（geom_jitter对点的位置进行小幅度随机偏移以减少点的重叠）
              aes(x = group, y = avg_log2FC, color = P_group), 
              size = 0.15, width =0.38) + # 点的大小、偏移幅度
  scale_color_manual(values = c("up-regulated" = "#EB7153", "down-regulated" = "#205AA7", "no significant" = "#B7B7B7")) + #根据分类设置散点的颜色
  geom_jitter(data = TopGene, # 绘制TopGene数据点
              aes(x = group, y = avg_log2FC), color = "transparent", # 所有TopGene点的颜色
              size = 0.2, width =0.3) + 
  geom_text_repel(data = TopGene, # 绘制TopGene的标签
                  aes(x = group, y = avg_log2FC, label = gene),
                  size = 2, color = 'transparent', # 标签大小、颜色
                  segment.size = 0.2,  # 指示线粗细
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50), force = 0.6) + # default，最大重叠数; force, 排斥力，越高越推开标签防止重叠
  geom_tile(data = DEGs_subsets3, # 绘制中心分组名称条框
            aes(x = group, y = 0, fill = group), 
            color = "transparent", height=0.25, alpha = 0.6, show.legend = F) + #边框颜色、条框高度、透明度、不显示图例
  scale_fill_manual(values = colors) +  #设置各组名称条框内部充填颜色
  labs(x="General celltypes", y="Average log2FC") +
  geom_text(data=DEGs_subsets3, # 绘制中心分组名称文本
            aes(x=group,  y=0,  label=group),
            size = 3, color ="transparent") + # 文字大小及颜色
  theme_minimal() + 
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 13))
ggsave('MOU_Immune3_subset_LigPg&Lig_DEGs_noLable.jpg', width = 10, height = 6, bg = 'white')








