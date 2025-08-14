# Distinct Immune Landscapes between Murine and Human Periodontitis

**Abstract**  

**Aim:** Periodontitis (PD) affects over 1 billion people worldwide, emphasizing the need for suitable animal models to study its pathogenesis. While murine ligature models are widely used to study periodontitis (PD), their immunological fidelity to human PD remains unclear. This study aimed to (1) systematically compare single-cell immune landscapes between human PD and murine ligature models (± _Porphyromonas gingivalis_, Pg), and (2) evaluate whether Pg supplementation improves model translatability.

**Materials and methods:** In this study, we performed scRNA-seq on murine gingival tissues under three conditions: ligature-induced PD, ligature+Pg, and healthy controls. Cross-species integration analysis was conducted with human PD datasets.

**Results:** scRNA-seq revealed distinct immune profiles between murine and human periodontitis: (1) Murine models showed significantly lower B/plasma cell Δ(PD-Ctrl) proportions than human (murine 0.10% vs. human 22.39%, _P_<0.001) and attenuated activation states; (2) Conversely, the change in macrophages Δ(PD-Ctrl) is significantly greater in mice (murine 28.53% vs. human -6.603%, _P_<0.0001) and exhibited hyperactivation; (3) Pg co-administration failed to rescue human-like immune signatures (_P_>0.05).

**Conclusion:** Our findings highlight critical immune cellular and transcriptomic differences between the murine ligature models and human PD, particularly in B/plasma cells and macrophages, suggesting to approach with caution when using these models for studying immune interactions in PD.

**Keywords***：

Periodontitis; Single-Cell Gene Expression Analysis; Immunity; B-Lymphocytes; Macrophages, _Porphyromonas gingivalis_

---

## 📂 文件结构  
├── data/ # 存放原始或处理后的数据  
├── scripts/ # 分析代码（R）  
└── README.md # 本文件

---
## 🛠️ 安装与运行  
### **依赖环境**  
- R ≥ 4.0.0
- 安装及载入依赖包：  
```R
# 提前安装并载入依赖包 
library(Seurat)
library(SeuratWrappers)
library(DoubletFinder)
library(batchelor)
library(harmony)
library(stringr)
library(GO.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(SingleR)
library(celldex)
library(ggplot2) 
library(ggrepel)
library(patchwork)
library(pheatmap)
library(magrittr)
library(tidyr) 
library(dplyr)
library(tibble)
library(reshape2)
library(scales)
library(ggraph)
library(cowplot)
library(data.table)
library(viridis)
library(ComplexHeatmap)
library(knitr)
library(purrr)
library(forcats)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(gridExtra)
library (VennDiagram)
library(monocle)
library(tidyverse)

```



