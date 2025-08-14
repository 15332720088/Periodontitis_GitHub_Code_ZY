# 用于绘图的Seurat对象
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step2_HU_scRNA_GeneralCelltype_20240620.Rdata")
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_scRNA1_GeneralCelltype_20240902.Rdata") # Control&PD


## Fig.1A ######################
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
Idents(MOU_scRNA1) <- "general_celltype"
p <-  DimPlot(MOU_scRNA1, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(p, file="./UMAP_MOU_scRNA1_general_celltype.tiff",width=5,height=3,dpi=600)


## Fig.1B ######################
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
#绘制细胞种类百分比 Excel表及条形图(各组)
project_percent <- with(MOU_scRNA1@meta.data, prop.table(table(general_celltype, project), margin = 2) * 100)
project_percent <- as.data.frame(project_percent)
project_percent$general_celltype <- as.character(project_percent$general_celltype)
project_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project_percent$general_celltype <- as.factor(c(project_percent$general_celltype))
project_percent$general_celltype <- fct_relevel(project_percent$general_celltype, c("Fibroblast","Epithelial","Endothelial","Immune", "Myocyte","Schwann_Neuron"))
project_percent$project <- fct_relevel(project_percent$project, c("Control","PD"))
project_percent$general_celltype
project_percent$project
project_percent <- project_percent[project_percent$general_celltype != "x", ] #删除x行
project_percent
write.table(project_percent, "./MOU_general_celltype_Bar2_percentage_projects.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project_percent, aes(x = project, y = Freq, fill = general_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./MOU_general_celltype_Bar2_percentage_projects.jpg",width = 3.5, height=4.5, dpi=600)



## Fig.1C ######################
gene_names=c("Col1a2")
p2 <- FeaturePlot(MOU_scRNA1, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#008489")) #pt.size选择点的大小
ggsave(p2, file = "./UMAP1_MOU_scRNA1_general_celltype_Fibroblast_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("Krt5")
p2 <- FeaturePlot(MOU_scRNA1, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#ff8b00")) #pt.size选择点的大小
ggsave(p2, file = "./UMAP2_MOU_scRNA1_general_celltype_Epithelial_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("Pecam1")
p2 <- FeaturePlot(MOU_scRNA1, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#DF0029")) #pt.size选择点的大小
ggsave(p2, file = "./UMAP3_MOU_scRNA1_general_celltype_Endothelial_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("Ptprc")
p2 <- FeaturePlot(MOU_scRNA1, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#103667")) #pt.size选择点的大小
ggsave(p2, file = "./UMAP4_MOU_scRNA1_general_celltype_Immune_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("Acta2")
p2 <- FeaturePlot(MOU_scRNA1, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#511F90")) #pt.size选择点的大小
ggsave(p2, file = "./UMAP5_MOU_scRNA1_general_celltype_Myocyte_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("Plp1")
p2 <- FeaturePlot(MOU_scRNA1, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#BD6B09")) #pt.size选择点的大小
ggsave(p2, file = "./UMAP6_MOU_scRNA1_general_celltype_Schwann_Neuron_MarkerGenes.jpg", width = 5, height=5, dpi=600)


## Fig.1D ######################
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")
Idents(HU_scRNA) <- "general_celltype"
p <-  DimPlot(HU_scRNA, reduction = "umap", label=F, cols = colors, raster=FALSE) 
ggsave(p, file="./UMAP_HU_scRNA_general_celltype.tiff",width=5,height=3,dpi=600)


## Fig.1E ######################
colors <- c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2")

#绘制细胞种类百分比 Excel表及条形图(各组)
project_percent <- with(HU_scRNA@meta.data, prop.table(table(general_celltype, project), margin = 2) * 100)
project_percent <- as.data.frame(project_percent)
project_percent$general_celltype <- as.character(project_percent$general_celltype)
project_percent
library(forcats) #调整因子level的排列顺序→调整Bar图中细胞种类的排列顺序
project_percent$general_celltype <- as.factor(c(project_percent$general_celltype))
project_percent$general_celltype <- fct_relevel(project_percent$general_celltype, c("Fibroblast","Epithelial","Endothelial","Immune", "Myocyte","Others"))
project_percent$project <- fct_relevel(project_percent$project, c("Control","PD"))
project_percent$general_celltype
project_percent$project
project_percent <- project_percent[project_percent$general_celltype != "x", ] #删除x行
project_percent
write.table(project_percent, "./HU_general_celltype_Bar2_percentage_projects.xls", col.names = T, row.names = F, sep = "\t")
Bar_celltype <- ggplot(project_percent, aes(x = project, y = Freq, fill = general_celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0, 1, by = 0.5), labels = c("0", "50%", "100%")) + 
  labs(x = " ", y = "Propotion", fill = " ") + 
  scale_fill_manual(values = colors) + 
  theme(panel.background = element_rect(fill = "white", linewidth = 2), axis.text.x = element_text(angle = 30, hjust = 1), axis.ticks = element_blank()) 
ggsave(Bar_celltype, file = "./HU_general_celltype_Bar2_percentage_projects.jpg",width = 3.5, height=4.5, dpi=600)


## Fig.1F ######################
gene_names=c("COL1A2")
p2 <- FeaturePlot(HU_scRNA, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#008489"), raster=FALSE) #pt.size选择点的大小
ggsave(p2, file = "./UMAP1_HU_scRNA_general_celltype_Fibroblast_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("KRT5")
p2 <- FeaturePlot(HU_scRNA, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#ff8b00"), raster=FALSE) #pt.size选择点的大小
ggsave(p2, file = "./UMAP2_HU_scRNA_general_celltype_Epithelial_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("PECAM1")
p2 <- FeaturePlot(HU_scRNA, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#DF0029"), raster=FALSE) #pt.size选择点的大小
ggsave(p2, file = "./UMAP3_HU_scRNA_general_celltype_Endothelial_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("PTPRC")
p2 <- FeaturePlot(HU_scRNA, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#103667"), raster=FALSE) #pt.size选择点的大小
ggsave(p2, file = "./UMAP4_HU_scRNA_general_celltype_Immune_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("CD79A")
p2 <- FeaturePlot(HU_scRNA, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#103667"), raster=FALSE) #pt.size选择点的大小
ggsave(p2, file = "./UMAP4_HU_scRNA_general_celltype_Immune2_MarkerGenes.jpg", width = 5, height=5, dpi=600)
gene_names=c("ACTA2")
p2 <- FeaturePlot(HU_scRNA, features = gene_names, pt.size = 0.5, reduction = "umap", cols = c("#D7D7D7", "#BD6B09"), raster=FALSE) #pt.size选择点的大小
ggsave(p2, file = "./UMAP5_HU_scRNA_general_celltype_Myocyte_MarkerGenes.jpg", width = 5, height=5, dpi=600)


## Fig.1G ######################
# HU
load("./HU_scRNA_PD&Control_DEGs1.Rdata") #绘制火山图的DEGs
#命名基因列名，计算表达分布Difference
DEGs1 <- DEGs1 %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2)
#根据阈值筛选差异显著的上下调基因，与差异不显著的基因
DEGs1$Sig_group = "no-significant"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC > 0.25))] = "up-regulated"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC < -0.25))] = "down-regulated"
table(DEGs1$Sig_group)
DEGs_up <- DEGs1[DEGs1$Sig_group == "up-regulated", ]
table(DEGs_up$Sig_group)
save(DEGs_up, file="./HU_scRNA_PD&Control_DEGs_up.Rdata")
DEGs_down <- DEGs1[DEGs1$Sig_group == "down-regulated", ]
table(DEGs_down$Sig_group)
save(DEGs_down, file="./HU_scRNA_PD&Control_DEGs_down.Rdata")

# MOU
load("./MOU_scRNA1_PD&Control_DEGs1.Rdata") #绘制火山图的DEGs
#命名基因列名，计算表达分布Difference
DEGs1 <- DEGs1 %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2)
#根据阈值筛选差异显著的上下调基因，与差异不显著的基因
DEGs1$Sig_group = "no-significant"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC > 0.25))] = "up-regulated"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC < -0.25))] = "down-regulated"
table(DEGs1$Sig_group)
DEGs_up <- DEGs1[DEGs1$Sig_group == "up-regulated", ]
table(DEGs_up$Sig_group)
save(DEGs_up, file="./MOU_scRNA1_PD&Control_DEGs_up.Rdata")
DEGs_down <- DEGs1[DEGs1$Sig_group == "down-regulated", ]
table(DEGs_down$Sig_group)
save(DEGs_down, file="./MOU_scRNA1_PD&Control_DEGs_down.Rdata")

load("./MOU_scRNA1_PD&Control_DEGs_up.Rdata")
#将小鼠基因名转化为人类基因ENTREZID（导入的是基因名集合，导出的是列表）
library(biomaRt)
load("/public/home/yinqi/zhaoyi/Aging/RData/Ref/Ensembl_human_gene.Rdata") 
load("/public/home/yinqi/zhaoyi/Aging/RData/Ref/Ensembl_mouse_gene.Rdata") 
gene_symbol <- as.character(DEGs_up$gene)
tDEGs_up <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol", #指定了输入和筛选的属性
                   values = gene_symbol, mart = mouse, #输入的待转换基因集及物种（注意values输入的数据结构应为基因的集合，不是一个数据框或列表）
                   attributesL = c("ensembl_gene_id","hgnc_symbol","entrezgene_id"), martL = human,  #输出的属性及物种
                   uniqueRows = T) #TRUE为即使输入中可能存在重复的基因也仅返回一行；FALSE则返回重复的行
save(tDEGs_up, file="./MOU_transHU_scRNA1_PD&Control_tDEGs_up.Rdata")
#列名： MGI.symbol    Gene.stable.ID    HGNC.symbol     NCBI.gene..formerly.Entrezgene..ID
#判断转化中丢失基因名数目
length(DEGs_up$gene)
length(tDEGs_up$MGI.symbol)

gene_symbol <- as.character(DEGs_down$gene)
tDEGs_down <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol", #指定了输入和筛选的属性
                     values = gene_symbol, mart = mouse, #输入的待转换基因集及物种（注意values输入的数据结构应为基因的集合，不是一个数据框或列表）
                     attributesL = c("ensembl_gene_id","hgnc_symbol","entrezgene_id"), martL = human,  #输出的属性及物种
                     uniqueRows = T) #TRUE为即使输入中可能存在重复的基因也仅返回一行；FALSE则返回重复的行
save(tDEGs_down, file="./MOU_transHU_scRNA1_PD&Control_tDEGs_down.Rdata")
#列名： MGI.symbol    Gene.stable.ID    HGNC.symbol     NCBI.gene..formerly.Entrezgene..ID
#判断转化中丢失基因名数目
length(DEGs_down$gene)
length(tDEGs_down$MGI.symbol)

load("./MOU_transHU_scRNA1_PD&Control_tDEGs_up.Rdata")
load("./HU_scRNA_PD&Control_DEGs_up.Rdata")
#（两交集）
venn_list <- list(MOU = tDEGs_up$HGNC.symbol, HU = DEGs_up$gene)
venn.diagram(venn_list, filename = 'Venn_HU&MOU_scRNA_PD&Control_DEGs_up.tiff', imagetype = 'tiff', 
             fill = c("#103667","#8B0016"), # 圆圈内部充填颜色
             category.names = c("HU", "MOU"),
             alpha = 0.8, # 圆圈内部充填颜色的透明度
             label.col = "transparent", # 圈内基因数目标签的颜色
             cex = 1.8, # 圈内基因数目标签的字体大小
             col = "transparent", # 圆圈线的颜色
             lwd = 1, # 圆圈线的宽度
             cat.col = "black", # 圈外样本名标签的颜色
             cat.cex = 1.8, # 圈外样本名标签的字体大小
             cat.dist = c(0.07, 0.07), # 圈外样本名标签距圆圈的距离（几个标签几个距离）
             cat.pos = c(180, 180), # 圈外样本名标位于圆圈的角度
             main = "PD&Control Up",  # 标题              
             main.col = "black",  # 标题颜色              
             main.cex = 1.8,  # 标题字体大小              
             ext.text = F, # 取消引导线注释              
             rotation.degree = 180, # 调整旋转角度              
             scaled = TRUE) # TRUE按比例绘制圆圈大小；FALSE圆圈大小一致

#导出韦恩图交集基因
inter <- get.venn.partitions(venn_list)
save(inter, file="./Venn_HU&MOU_scRNA_PD&Control_DEGs_up_interGenes.Rdata")
str(inter)

load("./MOU_transHU_scRNA1_PD&Control_tDEGs_down.Rdata")
load("./HU_scRNA_PD&Control_DEGs_down.Rdata")
#（两交集）
venn_list <- list(HU = DEGs_down$gene, MOU = tDEGs_down$HGNC.symbol)
venn.diagram(venn_list, filename = 'Venn_HU&MOU_scRNA_PD&Control_DEGs_down.tiff', imagetype = 'tiff', 
             fill = c("#103667","#8B0016"), # 圆圈内部充填颜色
             category.names = c("HU", "MOU"),
             alpha = 0.8, # 圆圈内部充填颜色的透明度
             label.col = "transparent", # 圈内基因数目标签的颜色
             cex = 1.8, # 圈内基因数目标签的字体大小
             col = "transparent", # 圆圈线的颜色
             lwd = 1, # 圆圈线的宽度
             cat.col = "black", # 圈外样本名标签的颜色
             cat.cex = 1.8, # 圈外样本名标签的字体大小
             cat.dist = c(0.07, 0.07), # 圈外样本名标签距圆圈的距离（几个标签几个距离）
             cat.pos = c(180, 180), # 圈外样本名标位于圆圈的角度
             main = "PD&Control Down",  # 标题              
             main.col = "black",  # 标题颜色              
             main.cex = 1.8,  # 标题字体大小              
             ext.text = F, # 取消引导线注释              
             rotation.degree = 180, # 调整旋转角度              
             scaled = TRUE) # TRUE按比例绘制圆圈大小；FALSE圆圈大小一致
#导出韦恩图交集基因
inter <- get.venn.partitions(venn_list)
save(inter, file="./Venn2_HU&MOU_scRNA_PD&Control_DEGs_down_interGenes.Rdata")
str(inter)



## Fig.1H-I ######################
##小鼠
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step7_MOU_scRNA1_GeneralCelltype_20240902.Rdata")
# 找出两组间差异基因（只比较指定两组）
Idents(MOU_scRNA1) <- "project"
DEGs1 = FindMarkers(MOU_scRNA1, ident.1 = "PD", ident.2 = c("Control"), logfc.threshold = 0.1, only.pos = FALSE) #（only.pos = TRUE，只返回阳性结果）
save(DEGs1, file="./MOU_scRNA1_PD&Control_DEGs1.Rdata")
write.table(DEGs1, "./MOU_scRNA1_PD&Control_DEGs1.xls", col.names = T, row.names = F, sep = "\t")

load("./MOU_scRNA1_PD&Control_DEGs1.Rdata")

#命名基因列名，计算表达分布Difference
DEGs1 <- DEGs1 %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2)
#对p.adjust列进行-log10转换，并存储在新列log10_p.adjust中
DEGs1$log10_P_adj <- -log10(DEGs1$p_val_adj) 
#根据阈值筛选差异显著的上下调基因，与差异不显著的基因
DEGs1$Sig_group = "no-significant"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC > 0.1) & (DEGs1$Difference > 0.1))] = "up-regulated"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC < -0.1) & (DEGs1$Difference < -0.1))] = "down-regulated"
table(DEGs1$Sig_group)

#添加Label基因
#解释，例：取出DEGs1中up-regulated组log10_P_adj前8名，前8名可能会有很多重复相同的log10_P_adj → 再按avg_log2FC绝对值abs()由大到小排序，取出前8行
# top8_log10P_up = DEGs1 %>% filter(Sig_group == "up-regulated") %>% top_n(n = 8, wt = log10_P_adj) %>% arrange(desc(abs(avg_log2FC)))  %>% slice_head(n = 8)
# top8_log10P_down = DEGs1 %>% filter(Sig_group == "down-regulated") %>% top_n(n = 8, wt = log10_P_adj) %>% arrange(desc(abs(avg_log2FC)))  %>% slice_head(n = 8)
top8_log2FC_up = DEGs1 %>% filter(Sig_group == "up-regulated") %>% top_n(n = 5, wt = abs(avg_log2FC))
top8_log2FC_down = DEGs1 %>% filter(Sig_group == "down-regulated") %>% top_n(n = 5, wt = abs(avg_log2FC))
Lable_gene <- unique(as.character(c(top8_log2FC_up$gene,top8_log2FC_down$gene)))
# 删掉过于拥挤的基因标签
# Lable_gene <- Lable_gene[Lable_gene != "AY036118"]

volcano_plot <- ggplot(DEGs1, aes(x=Difference, y=avg_log2FC, color = Sig_group)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c( "#103667","grey","#8B0016") ) + 
  geom_label_repel(data=subset(DEGs1, gene %in% Lable_gene), 
                   aes(label=gene),  #添加label
                   color="transparent", #设置label中文字的颜色
                   size = 4, # 设置label中文字的大小
                   label.size = 0,
                   label.padding = 0.1,  # 设置标签框与内部文本之间的距离
                   segment.colour = "black", #设置label框的颜色
                   fill = "transparent", # 设置标签的背景色
                   segment.size = 0.25, #设置label框的大小
                   box.padding = 0.3, # 调整标签和其他数据框及点之间的距离
                   point.padding = 0.3, # 调整标签和对应数据点之间的距离
                   max.overlaps = 50) + #允许重叠的最大标签数量
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()+
  scale_x_continuous(limits = c(-1.5, 1.5)) +  # 设置x轴范围
  scale_y_continuous(limits = c(-3.5, 3.5))  # 设置y轴范围
ggsave(volcano_plot, file = "./Volcano_MOU_scRNA1_PD&Control_DEGs1.jpg",width = 8, height=4.5, dpi=600)

# 共用上述绘图代码
ggsave(volcano_plot, file = "./Volcano2_HU_scRNA_PD&Control_DEGs1.jpg",width = 9, height=4.5, dpi=600)


##人类
load("/public/home/yinqi/zhaoyi/periodontitis/Rdata/Step2_HU_scRNA_GeneralCelltype_20240620.Rdata")
# 找出两组间差异基因（只比较指定两组）
Idents(HU_scRNA) <- "project"
DEGs1 = FindMarkers(HU_scRNA, ident.1 = "PD", ident.2 = c("Control"), logfc.threshold = 0.1, only.pos = FALSE) #（only.pos = TRUE，只返回阳性结果）
save(DEGs1, file="./HU_scRNA_PD&Control_DEGs1.Rdata")
write.table(DEGs1, "./HU_scRNA_PD&Control_DEGs1.xls", col.names = T, row.names = F, sep = "\t")

load("./HU_scRNA_PD&Control_DEGs1.Rdata")

#命名基因列名，计算表达分布Difference
DEGs1 <- DEGs1 %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2)
#对p.adjust列进行-log10转换，并存储在新列log10_p.adjust中
DEGs1$log10_P_adj <- -log10(DEGs1$p_val_adj) 
#根据阈值筛选差异显著的上下调基因，与差异不显著的基因
DEGs1$Sig_group = "no-significant"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC > 0.1) & (DEGs1$Difference > 0.1))] = "up-regulated"
DEGs1$Sig_group[which((DEGs1$p_val_adj < 0.05) & (DEGs1$avg_log2FC < -0.1) & (DEGs1$Difference < -0.1))] = "down-regulated"
table(DEGs1$Sig_group)

#添加Label基因
#解释，例：取出DEGs1中up-regulated组log10_P_adj前5名，前5名可能会有很多重复相同的log10_P_adj → 再按avg_log2FC绝对值abs()由大到小排序，取出前5行
# top5_log10P_up = DEGs1 %>% filter(Sig_group == "up-regulated") %>% top_n(n = 5, wt = log10_P_adj) %>% arrange(desc(abs(avg_log2FC)))  %>% slice_head(n = 5)
# top5_log10P_down = DEGs1 %>% filter(Sig_group == "down-regulated") %>% top_n(n = 5, wt = log10_P_adj) %>% arrange(desc(abs(avg_log2FC)))  %>% slice_head(n = 5)
top5_log2FC_up = DEGs1 %>% filter(Sig_group == "up-regulated") %>% top_n(n = 6, wt = abs(avg_log2FC))
top5_log2FC_down = DEGs1 %>% filter(Sig_group == "down-regulated") %>% top_n(n = 5, wt = abs(avg_log2FC))
Lable_gene <- unique(as.character(c(top5_log2FC_up$gene,top5_log2FC_down$gene)))
# 删掉过于拥挤的基因标签
# Lable_gene <- Lable_gene[Lable_gene != "AY036115"]

# （人）绘图共用上述（小鼠）绘图代码



## Fig.1J ######################
##找出各免疫亚群差异基因
general_celltypes <- levels(MOU_scRNA1$general_celltype)
celltype_subsets <- list() #建一个空列表存储每个细胞类型的子集对象
DEGs_subsets <- list() #建一个空列表存储每个亚群差异基因的子集对象
Idents(MOU_scRNA1) <- "general_celltype"
for (celltype in general_celltypes) {
  subset = subset(MOU_scRNA1,idents = c(as.character(celltype)))
  celltype_subsets[[celltype]] <- subset
  Idents(subset) <- "project"
  DEGs <- FindMarkers(subset, ident.1 = "PD", ident.2 = c("Control"), logfc.threshold = 0, only.pos = FALSE)
  DEGs_subsets[[celltype]] <- DEGs
}
#去除列表中各数据框名字中的空格
new_names <- gsub(" ", "", names(DEGs_subsets))
names(DEGs_subsets) <- new_names
#显示并保存结果
lapply(DEGs_subsets, function(df) head(df, 5))
save(DEGs_subsets, file=("./MOU_scRNA1_subset_PD&Control_DEGs.Rdata"))

###有分组细胞数小于3时，报错用以下代码：
general_celltypes <- levels(MOU_Immune$general_celltype)
celltype_subsets <- list() #建一个空列表存储每个细胞类型的子集对象
DEGs_subsets <- list() #建一个空列表存储每个亚群差异基因的子集对象
Idents(MOU_Immune) <- "general_celltype"
for (celltype in general_celltypes) {
  subset = subset(MOU_Immune,idents = c(as.character(celltype)))
  celltype_subsets[[celltype]] <- subset
  Idents(subset) <- "project"
  # 检查PD和Control组中的细胞数目
  cells_PD <- WhichCells(subset, idents = "PD")
  cells_Control <- WhichCells(subset, idents = "Control")
  # 仅在PD和Control组中的细胞数目都不少于3个时才执行FindMarkers操作
  if (length(cells_PD) >= 3 && length(cells_Control) >= 3) {
    DEGs <- FindMarkers(subset, ident.1 = "PD", ident.2 = c("Control"), logfc.threshold = 0, only.pos = FALSE)
    DEGs_subsets[[celltype]] <- DEGs
  } else {
    message(paste("Skipping cell type", celltype, "due to insufficient cell numbers in PD or Control group"))
  }}
#去除列表中各数据框名字中的空格
new_names <- gsub(" ", "", names(DEGs_subsets))
names(DEGs_subsets) <- new_names
#显示并保存结果
lapply(DEGs_subsets, function(df) head(df, 5))
save(DEGs_subsets, file=("./MOU_scRNA1_subset_PD&Control_DEGs.Rdata"))


##对上述结果列表进行调整
DEGs_subsets2 <- list() #建一个空列表用于存储结果
for (name in names(DEGs_subsets)) { #对列表中的每个数据框执行操作
  df <- DEGs_subsets[[name]]
  df <- df %>% filter(avg_log2FC > 0) #!筛选上调基因
  df$group = as.character(name) #df$group中充填亚群名
  df <- df %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2) %>% filter(abs(Difference) > 0.2 & p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
  DEGs_subsets2[[name]] <- df
}
lapply(DEGs_subsets2, function(df) head(df, 5))
DEGs2 <- bind_rows(DEGs_subsets2) #合并列表中的所有数据框
head(DEGs2)
save(DEGs2, file=("./MOU_scRNA1_GeneralCelltype_PD&Control_up_DEGs2.Rdata"))

##对各亚群差异基因进行ENTREZID_ID转换
table(DEGs2$group)
general_celltype <- c("Endothelial","Epithelial","Fibroblast","Immune","Myocyte","Schwann_Neuron")
ENTREZID_ID <- list() #循环分别取出每个组的gene名，转换为ENTREZID形式，最后存在 “ENTREZID_ID” list中
for (name in general_celltype) { #对列表中的每个数据框执行操作
  sub <- subset(DEGs2, group==name, select = gene)
  gene_symbol <- as.character(sub$gene)
  gene_ENTREZID <-  bitr(gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  result <- as.character(gene_ENTREZID[,2])
  ENTREZID_ID[[name]] <- result
}
save(ENTREZID_ID, file=("./MOU_scRNA1_GeneralCelltype_PD&Control_up_ENTREZID_ID.Rdata"))

##进行GO分析
GO_DEGs <- compareCluster(ENTREZID_ID, fun="enrichGO", OrgDb="org.Mm.eg.db", ont= "ALL", pvalueCutoff=0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
# ！人类为org.Mm.eg.db; ont = "ALL",对所有三个分类进行富集分析，也可单独填写“BP/MF/CC”; readable = TRUE,将ID转换为基因的名字
save(GO_DEGs, file="./MOU_scRNA1_GeneralCelltype_PD&Control_up_GO_DEGs.Rdata")
write.table(GO_DEGs, "./MOU_scRNA1_GeneralCelltype_PD&Control_up_GO_DEGs.xls", col.names = T, row.names = F, sep = "\t")


##总的细胞差异基因对比及GO分析
load("./MOU_scRNA1_PD&Control_DEGs1.Rdata")
DEGs <- DEGs1 %>% filter(avg_log2FC > 0) #!筛选上调基因
DEGs <- DEGs  %>% mutate(Difference = pct.1 - pct.2) %>% filter(abs(Difference) > 0.2 & p_val_adj < 0.05 & avg_log2FC > 0.5)
head(DEGs)
gene_symbol <- as.character(rownames(DEGs))
gene_ENTREZID <-  bitr(gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ENTREZID <- as.character(gene_ENTREZID[,2])
GO_all_DEGs <- enrichGO(gene = ENTREZID, OrgDb = "org.Mm.eg.db", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
save(GO_all_DEGs, file="./MOU_scRNA1_PD&Control_up_GO_all_DEGs.Rdata")
write.table(GO_all_DEGs, "./MOU_scRNA1_PD&Control_up_GO_all_DEGs.xls", col.names = T, row.names = F, sep = "\t")
# ONTOLOGY; ID; Description; GeneRatio; BgRatio; pvalue; p.adjust; qvalue; geneID; Count


load("./MOU_scRNA1_GeneralCelltype_PD&Control_up_GO_DEGs.Rdata")
##制作数据框包括以下几列：总细胞的top GO条目(不够时加一些挑选各亚群的标志性top条目)、细胞亚群名、Count数、-log10(P.adj)
#筛选特定的GO条目行：总GO分析的前4条+各亚群GO分析的代表条目
GO_DEGs <- as.data.frame(GO_DEGs)
GO_clusters_1 <- GO_DEGs %>%
  filter(Cluster == "Fibroblast" & ID %in% c("GO:0046034", "GO:0006091", "GO:0050727", "GO:2001233", "GO:0050900", "GO:0009611", "GO:0050678", "GO:0042060"))
GO_clusters_2 <- GO_DEGs %>%
  filter(Cluster == "Epithelial" & ID %in% c("GO:0046034", "GO:0006091", "GO:0050727", "GO:2001233", "GO:0050900", "GO:0009611", "GO:0050678", "GO:0042060"))
GO_clusters_3 <- GO_DEGs %>%
  filter(Cluster == "Endothelial" & ID %in% c("GO:0046034", "GO:0006091", "GO:0050727", "GO:2001233", "GO:0050900", "GO:0009611", "GO:0050678", "GO:0042060"))
GO_clusters_4 <- GO_DEGs %>%
  filter(Cluster == "Immune" & ID %in% c("GO:0046034", "GO:0006091", "GO:0050727", "GO:2001233", "GO:0050900", "GO:0009611", "GO:0050678", "GO:0042060"))
GO_clusters_6 <- GO_DEGs %>%
  filter(Cluster == "Myocyte" & ID %in% c("GO:0046034", "GO:0006091", "GO:0050727", "GO:2001233", "GO:0050900", "GO:0009611", "GO:0050678", "GO:0042060"))
# 合并筛选后的GO条目行
GO_clusters_specific_terms <- bind_rows(GO_clusters_1, GO_clusters_2, GO_clusters_3, GO_clusters_4, GO_clusters_6)
#调整细胞亚群的名称排列顺序
GO_clusters_specific_terms$Cluster <- fct_relevel(GO_clusters_specific_terms$Cluster, c("Fibroblast","Epithelial","Endothelial","Immune","Schwann_Neuron","Myocyte"))
table(GO_clusters_specific_terms$Cluster)

#调整GO条目的level顺序
GO_clusters_specific_terms$Description <- as.factor(GO_clusters_specific_terms$Description)
GO_clusters_specific_terms$Description <- fct_relevel(GO_clusters_specific_terms$Description, rev(c(
  "ATP metabolic process", 
  "generation of precursor metabolites and energy",
  "regulation of inflammatory response",
  "regulation of apoptotic signaling pathway", 
  "leukocyte migration",
  "response to wounding",
  "regulation of epithelial cell proliferation",
  "wound healing")))
levels(GO_clusters_specific_terms$Description)

#添加-log10(P.adj)列
GO_clusters_specific_terms$log10_p_adjust <- -log10(GO_clusters_specific_terms$p.adjust)
#保存用于绘图的数据框文件
save(GO_clusters_specific_terms, file="./MOU_scRNA1_PD&Control_general_celltype_up_GO_specific_terms.Rdata")


##绘制GO分析结果Dot图
Dotplot <- ggplot(GO_clusters_specific_terms, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = log10_p_adjust)) +
  scale_color_gradient2(low = "#FCD9C4", mid = "#EB7153", high = "#C50023", #点的颜色
                        limits = c(2, 8), midpoint = 3, # 设置颜色渐变的上下限；中间值的位置
                        oob = scales::squish, # 将超出范围的值挤压到边界
                        breaks = seq(2, 8, by = 2)) +  # 设置颜色渐变的间隔  
  scale_size_continuous(range = c(3, 15)) +  # 设置点大小的间隔
  theme_minimal() +  # 使用简洁的主题
  labs(title = "", x = "", y = "",
       color = "-log10(p.adjust)", size = "Count") +
  theme(text = element_text(size = 15),  # 统一设置所有文本的字体大小
        axis.title = element_text(size = 20,color = "black",face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.justification = c(1,0),
        legend.title = element_text(size = 15),  # 设置图例标题字体大小
        legend.text = element_text(size = 12),  # 设置图例文本字体大小
        panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        # panel.grid.major = element_blank(),  # 移除主要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(Dotplot, file = "./MOU_scRNA1_PD&Control_general_celltype_up_GO_specific_terms.jpg",width = 15, height=10, dpi=600)



## Fig.1K ######################
##找出各免疫亚群差异基因
general_celltypes <- levels(HU_scRNA$general_celltype)
celltype_subsets <- list() #建一个空列表存储每个细胞类型的子集对象
DEGs_subsets <- list() #建一个空列表存储每个亚群差异基因的子集对象
Idents(HU_scRNA) <- "general_celltype"
for (celltype in general_celltypes) {
  subset = subset(HU_scRNA,idents = c(as.character(celltype)))
  celltype_subsets[[celltype]] <- subset
  Idents(subset) <- "project"
  DEGs <- FindMarkers(subset, ident.1 = "PD", ident.2 = c("Control"), logfc.threshold = 0, only.pos = FALSE)
  DEGs_subsets[[celltype]] <- DEGs
}
#去除列表中各数据框名字中的空格
new_names <- gsub(" ", "", names(DEGs_subsets))
names(DEGs_subsets) <- new_names
#显示并保存结果
lapply(DEGs_subsets, function(df) head(df, 5))
save(DEGs_subsets, file=("./HU_scRNA_subset_PD&Control_DEGs.Rdata"))

###有分组细胞数小于3时，报错用以下代码：
general_celltypes <- levels(MOU_Immune$general_celltype)
celltype_subsets <- list() #建一个空列表存储每个细胞类型的子集对象
DEGs_subsets <- list() #建一个空列表存储每个亚群差异基因的子集对象
Idents(MOU_Immune) <- "general_celltype"
for (celltype in general_celltypes) {
  subset = subset(MOU_Immune,idents = c(as.character(celltype)))
  celltype_subsets[[celltype]] <- subset
  Idents(subset) <- "project"
  # 检查PD和Control组中的细胞数目
  cells_PD <- WhichCells(subset, idents = "PD")
  cells_Control <- WhichCells(subset, idents = "Control")
  # 仅在PD和Control组中的细胞数目都不少于3个时才执行FindMarkers操作
  if (length(cells_PD) >= 3 && length(cells_Control) >= 3) {
    DEGs <- FindMarkers(subset, ident.1 = "PD", ident.2 = c("Control"), logfc.threshold = 0, only.pos = FALSE)
    DEGs_subsets[[celltype]] <- DEGs
  } else {
    message(paste("Skipping cell type", celltype, "due to insufficient cell numbers in PD or Control group"))
  }}
#去除列表中各数据框名字中的空格
new_names <- gsub(" ", "", names(DEGs_subsets))
names(DEGs_subsets) <- new_names
#显示并保存结果
lapply(DEGs_subsets, function(df) head(df, 5))
save(DEGs_subsets, file=("./HU_scRNA_subset_PD&Control_DEGs.Rdata"))


##对上述结果列表进行调整
DEGs_subsets2 <- list() #建一个空列表用于存储结果
for (name in names(DEGs_subsets)) { #对列表中的每个数据框执行操作
  df <- DEGs_subsets[[name]]
  df <- df %>% filter(avg_log2FC > 0) #!筛选上调基因
  df$group = as.character(name) #df$group中充填亚群名
  df <- df %>% rownames_to_column("gene") %>% mutate(Difference = pct.1 - pct.2) %>% filter( p_val_adj < 0.05 & avg_log2FC > 0.5)
  DEGs_subsets2[[name]] <- df
}
lapply(DEGs_subsets2, function(df) head(df, 5))
DEGs2 <- bind_rows(DEGs_subsets2) #合并列表中的所有数据框
head(DEGs2)
save(DEGs2, file=("./HU_scRNA_GeneralCelltype_PD&Control_up_DEGs2.Rdata"))

DEGs2 <- DEGs2[!DEGs2$group %in% c("Others"), ]

##对各亚群差异基因进行ENTREZID_ID转换
table(DEGs2$group)
general_celltype <- c("Endothelial","Epithelial","Fibroblast","Immune","Myocyte")
ENTREZID_ID <- list() #循环分别取出每个组的gene名，转换为ENTREZID形式，最后存在 “ENTREZID_ID” list中
for (name in general_celltype) { #对列表中的每个数据框执行操作
  sub <- subset(DEGs2, group==name, select = gene)
  gene_symbol <- as.character(sub$gene)
  gene_ENTREZID <-  bitr(gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  result <- as.character(gene_ENTREZID[,2])
  ENTREZID_ID[[name]] <- result
}
save(ENTREZID_ID, file=("./HU_scRNA_GeneralCelltype_PD&Control_up_ENTREZID_ID.Rdata"))

##进行GO分析
GO_DEGs <- compareCluster(ENTREZID_ID, fun="enrichGO", OrgDb="org.Hs.eg.db", ont= "ALL", pvalueCutoff=0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
# ！人类为org.Hs.eg.db; ont = "ALL",对所有三个分类进行富集分析，也可单独填写“BP/MF/CC”; readable = TRUE,将ID转换为基因的名字
save(GO_DEGs, file="./HU_scRNA_GeneralCelltype_PD&Control_up_GO_DEGs.Rdata")
write.table(GO_DEGs, "./HU_scRNA_GeneralCelltype_PD&Control_up_GO_DEGs.xls", col.names = T, row.names = F, sep = "\t")


##总的细胞差异基因对比及GO分析
Idents(HU_scRNA) <- "project"
DEGs = FindMarkers(HU_scRNA, ident.1 = "PD", ident.2 = c("Control"), min.pct = 0.25, logfc.threshold = 0.5, only.pos = TRUE) #（only.pos = TRUE，只返回阳性结果）
save(DEGs, file="./HU_scRNA_PD&Control_all_up_DEGs.Rdata")
write.table(DEGs, "./HU_scRNA_PD&Control_all_up_DEGs.xls", col.names = T, row.names = F, sep = "\t")
DEGs <- DEGs  %>% mutate(Difference = pct.1 - pct.2) %>% filter(abs(Difference) > 0.2 & p_val_adj < 0.05 & avg_log2FC > 0.5)
head(DEGs)
gene_symbol <- as.character(rownames(DEGs))
gene_ENTREZID <-  bitr(gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ENTREZID <- as.character(gene_ENTREZID[,2])
GO_all_DEGs <- enrichGO(gene = ENTREZID, OrgDb = "org.Mm.eg.db", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
save(GO_all_DEGs, file="./HU_scRNA_PD&Control_all_up_GO_all_DEGs.Rdata")
write.table(GO_all_DEGs, "./HU_scRNA_PD&Control_all_up_GO_all_DEGs.xls", col.names = T, row.names = F, sep = "\t")
# ONTOLOGY; ID; Description; GeneRatio; BgRatio; pvalue; p.adjust; qvalue; geneID; Count

load("./HU_scRNA_GeneralCelltype_PD&Control_up_GO_DEGs.Rdata")

##制作数据框包括以下几列：总细胞的top GO条目及各亚群的top条目、细胞亚群名、Count数、-log10(P.adj)
#筛选特定的GO条目行：总GO分析的前4条+各亚群GO分析的代表条目
GO_DEGs <- as.data.frame(GO_DEGs)
GO_clusters_1 <- GO_DEGs %>%
  filter(Cluster == "Fibroblast" & ID %in% c("GO:0050867", "GO:0002696", "GO:0002697", "GO:0045785", "GO:0042113", "GO:0042110", "GO:0002283", "GO:0002237"))
GO_clusters_2 <- GO_DEGs %>%
  filter(Cluster == "Epithelial" & ID %in% c("GO:0050867", "GO:0002696", "GO:0002697", "GO:0045785", "GO:0042113", "GO:0042110", "GO:0002283", "GO:0002237"))
GO_clusters_3 <- GO_DEGs %>%
  filter(Cluster == "Endothelial" & ID %in% c("GO:0050867", "GO:0002696", "GO:0002697", "GO:0045785", "GO:0042113", "GO:0042110", "GO:0002283", "GO:0002237"))
GO_clusters_4 <- GO_DEGs %>%
  filter(Cluster == "Immune" & ID %in% c("GO:0050867", "GO:0002696", "GO:0002697", "GO:0045785", "GO:0042113", "GO:0042110", "GO:0002283", "GO:0002237"))
GO_clusters_6 <- GO_DEGs %>%
  filter(Cluster == "Myocyte" & ID %in% c("GO:0050867", "GO:0002696", "GO:0002697", "GO:0045785", "GO:0042113", "GO:0042110", "GO:0002283", "GO:0002237"))
# 合并筛选后的GO条目行
GO_clusters_specific_terms <- bind_rows(GO_clusters_1, GO_clusters_2, GO_clusters_3, GO_clusters_4, GO_clusters_6)
#调整细胞亚群的名称排列顺序
GO_clusters_specific_terms$Cluster <- fct_relevel(GO_clusters_specific_terms$Cluster, c("Fibroblast","Epithelial","Endothelial","Immune","Schwann_Neuron","Myocyte"))
table(GO_clusters_specific_terms$Cluster)

#调整GO条目的level顺序
GO_clusters_specific_terms$Description <- as.factor(GO_clusters_specific_terms$Description)
GO_clusters_specific_terms$Description <- fct_relevel(GO_clusters_specific_terms$Description, rev(c(
  "positive regulation of cell activation", 
  "positive regulation of leukocyte activation",
  "regulation of immune effector process",
  "positive regulation of cell adhesion", 
  "B cell activation",
  "T cell activation",
  "neutrophil activation involved in immune response",
  "response to molecule of bacterial origin")))
levels(GO_clusters_specific_terms$Description)

#添加-log10(P.adj)列
GO_clusters_specific_terms$log10_p_adjust <- -log10(GO_clusters_specific_terms$p.adjust)
#保存用于绘图的数据框文件
save(GO_clusters_specific_terms, file="./HU_scRNA_PD&Control_general_celltype_up_GO_specific_terms.Rdata")


##绘制GO分析结果Dot图
Dotplot <- ggplot(GO_clusters_specific_terms, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = log10_p_adjust)) +
  scale_color_gradient2(low = "#FCD9C4", mid = "#EB7153", high = "#C50023", #点的颜色
                        limits = c(2, 8), midpoint = 3, # 设置颜色渐变的上下限；中间值的位置
                        oob = scales::squish, # 将超出范围的值挤压到边界
                        breaks = seq(2, 8, by = 2)) +  # 设置颜色渐变的间隔  
  scale_size_continuous(range = c(3, 15)) +  # 设置点大小的间隔
  theme_minimal() +  # 使用简洁的主题
  labs(title = "", x = "", y = "",
       color = "-log10(p.adjust)", size = "Count") +
  theme(text = element_text(size = 15),  # 统一设置所有文本的字体大小
        axis.title = element_text(size = 20,color = "black",face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.justification = c(1,0),
        legend.title = element_text(size = 15),  # 设置图例标题字体大小
        legend.text = element_text(size = 12),  # 设置图例文本字体大小
        panel.background = element_rect(fill = "white", colour = "white"),  # 设置背景为白色
        # panel.grid.major = element_blank(),  # 移除主要网格线
        plot.background = element_rect(fill = "white", colour = "white")) # 设置绘图区域背景为白色
ggsave(Dotplot, file = "./HU_scRNA_PD&Control_general_celltype_up_GO_specific_terms.jpg",width = 15, height=10, dpi=600)






