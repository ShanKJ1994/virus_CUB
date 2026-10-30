
library(ComplexHeatmap)
library(tidyverse)
library(circlize)

ht_opt$message = FALSE

# ==========================================
# 0. 参数开关
# ==========================================
remove_Met_Trp <- T   # 保留 Met/Trp

# 对于行 dendrogram，"bottom" 表示把 Met/Trp 推到该 slice 的最后，也就是热图最下方
met_trp_position <- "bottom"

# ==========================================
# 1. 读取数据
# ==========================================



library(ComplexHeatmap)
library(tidyverse)
library(circlize)

ht_opt$message = FALSE

# ==========================================
# 0. 参数开关
# ==========================================
remove_Met_Trp <- T   # 保留 Met/Trp

# 对于行 dendrogram，"bottom" 表示把 Met/Trp 推到该 slice 的最后，也就是热图最下方
met_trp_position <- "bottom"

# ==========================================
# 1. 读取数据
# ==========================================
file_path <- "C:/Users/Ke-jia Shan/Desktop/Virus_CUB/Fig1/Results/RSCU.csv"

if (!file.exists(file_path)) {
  stop("找不到文件：", file_path)
}

raw_data <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ==========================================
# 2. 纵轴属性：codon 信息
# ==========================================
codon_cols <- colnames(raw_data)[3:ncol(raw_data)]

row_meta <- data.frame(Codon_Raw = codon_cols, stringsAsFactors = FALSE) %>%
  mutate(
    Codon = sub(".*:", "", Codon_Raw),
    Third_Base = toupper(substr(Codon, 3, 3)),
    Third_Base = ifelse(Third_Base == "T", "U", Third_Base),
    Amino_Acid = sub(":.*", "", Codon_Raw),
    
    # 这里标签写 A/T-ending，但由于上面 T 已转成 U，
    # 实际判断是 A/U-ending vs G/C-ending
    Ending_Group = factor(
      ifelse(Third_Base %in% c("A", "U"), "A/T-ending", "G/C-ending"),
      levels = c("A/T-ending", "G/C-ending")
    )
  )

if (remove_Met_Trp) {
  row_meta <- row_meta %>%
    filter(!Amino_Acid %in% c("Met", "Trp"))
}

# 先保证 A/T-ending 在上，G/C-ending 在下
row_meta <- row_meta %>%
  arrange(Ending_Group, Codon_Raw)

final_codon_cols <- row_meta$Codon_Raw

# Met 和 Trp 的行名
met_trp_rows <- row_meta$Codon_Raw[row_meta$Amino_Acid %in% c("Met", "Trp")]

message("Met/Trp rows: ", paste(met_trp_rows, collapse = ", "))

# ==========================================
# 3. 横轴顺序：按第二个 code
# ==========================================
desired_order <- c(
  "human", "bat", "Other mammals", "bird", "reptile", "fish",
  "Liliopsida_Plant", "eudicotyledons_Plant", "OtherPlants",
  "DNA", "RNA"
)

# 重要：先排序 raw_data，再生成矩阵和 col_meta
# 这样 Species 名称、数值矩阵、顶部注释不会错位
raw_data_ordered <- raw_data %>%
  mutate(Class_Factor = factor(Class, levels = desired_order)) %>%
  arrange(Class_Factor)

col_meta <- raw_data_ordered %>%
  select(Class, Species, Class_Factor)

# ==========================================
# 4. 矩阵
# ==========================================
plot_mat <- t(as.matrix(raw_data_ordered[, final_codon_cols]))

colnames(plot_mat) <- raw_data_ordered$Species
plot_mat <- plot_mat[row_meta$Codon_Raw, , drop = FALSE]

# ==========================================
# 4.1 聚类函数：保持默认聚类方式，只旋转 Met/Trp 到 G/C-ending 最下方
# ==========================================

rotate_targets_to_bottom <- function(dend, targets) {
  
  contains_target <- function(x) {
    any(labels(x) %in% targets)
  }
  
  swap_children <- function(d) {
    d2 <- d
    d2[[1]] <- d[[2]]
    d2[[2]] <- d[[1]]
    return(d2)
  }
  
  orient_node <- function(d) {
    if (is.leaf(d)) {
      return(d)
    }
    
    d[[1]] <- orient_node(d[[1]])
    d[[2]] <- orient_node(d[[2]])
    
    left_has_target <- contains_target(d[[1]])
    right_has_target <- contains_target(d[[2]])
    
    # 对行来说，labels(dend) 的后面通常对应热图更靠下的位置
    # 因此把包含 Met/Trp 的分支放到右侧，也就是 order 的末端
    if (left_has_target && !right_has_target) {
      d <- swap_children(d)
    }
    
    return(d)
  }
  
  orient_node(dend)
}

cluster_rows_like_default_with_mettrp_bottom <- function(m) {
  
  # 这一步对应 ComplexHeatmap 默认行聚类：
  # dist = euclidean, hclust method = complete
  hc <- hclust(dist(m), method = "complete")
  dend <- as.dendrogram(hc)
  
  # 模拟 ComplexHeatmap 默认的 dendrogram reorder
  dend <- reorder(dend, wts = rowMeans(m, na.rm = TRUE))
  
  # 当前 slice 内是否包含 Met/Trp
  current_targets <- intersect(labels(dend), met_trp_rows)
  
  # 只有 G/C-ending slice 里会包含 Met/Trp
  if (length(current_targets) > 0) {
    dend <- rotate_targets_to_bottom(
      dend = dend,
      targets = current_targets
    )
  }
  
  return(dend)
}

# ==========================================
# 5. 色阶：按第二个 code
# ==========================================
col_fun <- colorRamp2(
  c(0, 1, 4),
  c("#2166AC", "#f7f7f7", "#B2182B")
)

# ==========================================
# 6. 顶部注释颜色：按第二个 code
# ==========================================
class_colors <- c(
  Liliopsida_Plant = '#006400',
  eudicotyledons_Plant = '#00FF00',
  OtherPlants = '#b8ebb0',
  fish = '#B09C85FF',
  reptile = 'black',
  bird = '#FFC0CB',
  `Other mammals` = '#ADB6B6FF',
  human = 'red',
  bat = '#631879FF',
  RNA = '#57C3F3',
  DNA = '#0000FF'
)

top_anno <- HeatmapAnnotation(
  `Detailed Class` = col_meta$Class_Factor,
  col = list(`Detailed Class` = class_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "Class")
)

# ==========================================
# 7. 导出
# ==========================================
output_pdf <- "C:/Users/Ke-jia Shan/Desktop/RSCU_Heatmap_All_Codons_MetTrp_GC_Bottom.pdf"

pdf(output_pdf, width = 22, height = 10)

ht <- Heatmap(
  plot_mat,
  col = col_fun,
  name = "RSCU Value",
  top_annotation = top_anno,
  
  # --- 纵轴 ---
  row_split = row_meta$Ending_Group,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_gap = unit(6, "mm"),
  
  # 分组后，每个 slice 内部聚类；
  # 其中 G/C-ending slice 里额外旋转 Met/Trp 所在枝到最下方
  cluster_rows = cluster_rows_like_default_with_mettrp_bottom,
  
  # 防止 ComplexHeatmap 在我们旋转后再次自动 reorder
  row_dend_reorder = FALSE,
  
  row_names_gp = gpar(fontsize = 8.5),
  
  # --- 横轴：按第二个 code ---
  column_split = col_meta$Class_Factor,
  column_gap = unit(c(rep(1, 8), 4, 1), "mm"),
  column_title = NULL,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  
  use_raster = TRUE,
  raster_quality = 4
)

draw(ht, merge_legends = TRUE)

dev.off()

ht_opt$message = TRUE

message("绘图完成：", output_pdf)
message("密码子设置：全部61个密码子；Met/Trp保留在G/C-ending中，并旋转到G/C-ending最下方")


if (!file.exists(file_path)) {
  stop("找不到文件：", file_path)
}

raw_data <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ==========================================
# 2. 纵轴属性：codon 信息
# ==========================================
codon_cols <- colnames(raw_data)[3:ncol(raw_data)]

row_meta <- data.frame(Codon_Raw = codon_cols, stringsAsFactors = FALSE) %>%
  mutate(
    Codon = sub(".*:", "", Codon_Raw),
    Third_Base = toupper(substr(Codon, 3, 3)),
    Third_Base = ifelse(Third_Base == "T", "U", Third_Base),
    Amino_Acid = sub(":.*", "", Codon_Raw),
    
    # 这里标签写 A/T-ending，但由于上面 T 已转成 U，
    # 实际判断是 A/U-ending vs G/C-ending
    Ending_Group = factor(
      ifelse(Third_Base %in% c("A", "U"), "A/T-ending", "G/C-ending"),
      levels = c("A/T-ending", "G/C-ending")
    )
  )

if (remove_Met_Trp) {
  row_meta <- row_meta %>%
    filter(!Amino_Acid %in% c("Met", "Trp"))
}

# 先保证 A/T-ending 在上，G/C-ending 在下
row_meta <- row_meta %>%
  arrange(Ending_Group, Codon_Raw)

final_codon_cols <- row_meta$Codon_Raw

# Met 和 Trp 的行名
met_trp_rows <- row_meta$Codon_Raw[row_meta$Amino_Acid %in% c("Met", "Trp")]

message("Met/Trp rows: ", paste(met_trp_rows, collapse = ", "))

# ==========================================
# 3. 横轴顺序：按第二个 code
# ==========================================
desired_order <- c(
  "human", "bat", "Other mammals", "bird", "reptile", "fish",
  "Liliopsida_Plant", "eudicotyledons_Plant", "OtherPlants",
  "DNA", "RNA"
)

# 重要：先排序 raw_data，再生成矩阵和 col_meta
# 这样 Species 名称、数值矩阵、顶部注释不会错位
raw_data_ordered <- raw_data %>%
  mutate(Class_Factor = factor(Class, levels = desired_order)) %>%
  arrange(Class_Factor)

col_meta <- raw_data_ordered %>%
  select(Class, Species, Class_Factor)

# ==========================================
# 4. 矩阵
# ==========================================
plot_mat <- t(as.matrix(raw_data_ordered[, final_codon_cols]))

colnames(plot_mat) <- raw_data_ordered$Species
plot_mat <- plot_mat[row_meta$Codon_Raw, , drop = FALSE]

# ==========================================
# 4.1 聚类函数：保持默认聚类方式，只旋转 Met/Trp 到 G/C-ending 最下方
# ==========================================

rotate_targets_to_bottom <- function(dend, targets) {
  
  contains_target <- function(x) {
    any(labels(x) %in% targets)
  }
  
  swap_children <- function(d) {
    d2 <- d
    d2[[1]] <- d[[2]]
    d2[[2]] <- d[[1]]
    return(d2)
  }
  
  orient_node <- function(d) {
    if (is.leaf(d)) {
      return(d)
    }
    
    d[[1]] <- orient_node(d[[1]])
    d[[2]] <- orient_node(d[[2]])
    
    left_has_target <- contains_target(d[[1]])
    right_has_target <- contains_target(d[[2]])
    
    # 对行来说，labels(dend) 的后面通常对应热图更靠下的位置
    # 因此把包含 Met/Trp 的分支放到右侧，也就是 order 的末端
    if (left_has_target && !right_has_target) {
      d <- swap_children(d)
    }
    
    return(d)
  }
  
  orient_node(dend)
}

cluster_rows_like_default_with_mettrp_bottom <- function(m) {
  
  # 这一步对应 ComplexHeatmap 默认行聚类：
  # dist = euclidean, hclust method = complete
  hc <- hclust(dist(m), method = "complete")
  dend <- as.dendrogram(hc)
  
  # 模拟 ComplexHeatmap 默认的 dendrogram reorder
  dend <- reorder(dend, wts = rowMeans(m, na.rm = TRUE))
  
  # 当前 slice 内是否包含 Met/Trp
  current_targets <- intersect(labels(dend), met_trp_rows)
  
  # 只有 G/C-ending slice 里会包含 Met/Trp
  if (length(current_targets) > 0) {
    dend <- rotate_targets_to_bottom(
      dend = dend,
      targets = current_targets
    )
  }
  
  return(dend)
}

# ==========================================
# 5. 色阶：按第二个 code
# ==========================================
col_fun <- colorRamp2(
  c(0, 1, 4),
  c("#2166AC", "#f7f7f7", "#B2182B")
)

# ==========================================
# 6. 顶部注释颜色：按第二个 code
# ==========================================
class_colors <- c(
  Liliopsida_Plant = '#006400',
  eudicotyledons_Plant = '#00FF00',
  OtherPlants = '#b8ebb0',
  fish = '#B09C85FF',
  reptile = 'black',
  bird = '#FFC0CB',
  `Other mammals` = '#ADB6B6FF',
  human = 'red',
  bat = '#631879FF',
  RNA = '#57C3F3',
  DNA = '#0000FF'
)

top_anno <- HeatmapAnnotation(
  `Detailed Class` = col_meta$Class_Factor,
  col = list(`Detailed Class` = class_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "Class")
)

# ==========================================
# 7. 导出
# ==========================================
output_pdf <- "C:/Users/Ke-jia Shan/Desktop/RSCU_Heatmap_All_Codons_MetTrp_GC_Bottom.pdf"

pdf(output_pdf, width = 22, height = 10)

ht <- Heatmap(
  plot_mat,
  col = col_fun,
  name = "RSCU Value",
  top_annotation = top_anno,
  
  # --- 纵轴 ---
  row_split = row_meta$Ending_Group,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_gap = unit(6, "mm"),
  
  # 分组后，每个 slice 内部聚类；
  # 其中 G/C-ending slice 里额外旋转 Met/Trp 所在枝到最下方
  cluster_rows = cluster_rows_like_default_with_mettrp_bottom,
  
  # 防止 ComplexHeatmap 在我们旋转后再次自动 reorder
  row_dend_reorder = FALSE,
  
  row_names_gp = gpar(fontsize = 8.5),
  
  # --- 横轴：按第二个 code ---
  column_split = col_meta$Class_Factor,
  column_gap = unit(c(rep(1, 8), 4, 1), "mm"),
  column_title = NULL,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  
  use_raster = TRUE,
  raster_quality = 4
)

draw(ht, merge_legends = TRUE)

dev.off()

ht_opt$message = TRUE

message("绘图完成：", output_pdf)
message("密码子设置：59个密码子")








