library(ComplexHeatmap)
library(tidyverse)
library(circlize)

ht_opt$message = FALSE

# ==========================================
# 0. 参数开关
# ==========================================
remove_Met_Trp <- FALSE   # FALSE=保留Met/Trp；TRUE=去掉Met/Trp

# ==========================================
# 1. 读取数据
# ==========================================
file_path <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/RSCU_FigS1.csv"

if (!file.exists(file_path)) {
  stop("找不到文件：", file_path)
}

raw_data <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ==========================================
# 2. 纵轴属性
# ==========================================
codon_cols <- colnames(raw_data)[3:ncol(raw_data)]

row_meta <- data.frame(Codon_Raw = codon_cols, stringsAsFactors = FALSE) %>%
  mutate(
    Codon = sub(".*:", "", Codon_Raw),
    Third_Base = toupper(substr(Codon, 3, 3)),
    Third_Base = ifelse(Third_Base == "T", "U", Third_Base),
    Amino_Acid = sub(":.*", "", Codon_Raw),
    
    # 注意这里的顺序：
    # G/C-ending -> Met/Trp -> A/T-ending
    Ending_Group = factor(
      case_when(
        Amino_Acid %in% c("Met", "Trp") ~ "Met/Trp",
        Third_Base %in% c("A", "U") ~ "A/T-ending",
        TRUE ~ "G/C-ending"
      ),
      levels = c("G/C-ending", "Met/Trp", "A/T-ending")
    )
  )

if (remove_Met_Trp) {
  row_meta <- row_meta %>% filter(!Amino_Acid %in% c("Met", "Trp"))
  row_meta$Ending_Group <- droplevels(row_meta$Ending_Group)
}

row_meta <- row_meta %>% arrange(Ending_Group, Codon_Raw)
final_codon_cols <- row_meta$Codon_Raw

# ==========================================
# 3. 动物→植物→病毒
# ==========================================
desired_order <- c(
  "human", "bat", "Other mammals", "reptile", "bird", "fish",
  "eudicotyledons_Plant", "Liliopsida_Plant", "OtherPlants",
  "DNA", "RNA"
)

col_meta <- raw_data %>%
  select(Class, Species) %>%
  mutate(Class_Factor = factor(Class, levels = desired_order)) %>%
  arrange(Class_Factor)

# ==========================================
# 4. 矩阵
# ==========================================
plot_mat <- t(as.matrix(raw_data[, final_codon_cols]))
colnames(plot_mat) <- raw_data$Species

plot_mat <- plot_mat[row_meta$Codon_Raw, ]
plot_mat <- plot_mat[, col_meta$Species]

# ==========================================
# 4.1 构建自定义 row tree
# 目的：
# 1. 热图不分成3块，是一个整体
# 2. 左侧保留完整聚类树
# 3. 树结构强制为：
#    (G/C-ending + Met/Trp) + A/T-ending
# 4. 行顺序为：
#    G/C-ending -> Met/Trp -> A/T-ending
# ==========================================

make_leaf_dend <- function(label) {
  structure(
    list(),
    members = 1L,
    height = 0,
    label = label,
    leaf = TRUE,
    class = "dendrogram"
  )
}

make_group_dend <- function(mat, rows, method = "complete") {
  rows <- intersect(rows, rownames(mat))
  
  if (length(rows) == 0) {
    return(NULL)
  }
  
  if (length(rows) == 1) {
    return(make_leaf_dend(rows))
  }
  
  dend <- as.dendrogram(
    hclust(
      dist(mat[rows, , drop = FALSE]),
      method = method
    )
  )
  
  return(dend)
}

get_dend_height <- function(dend) {
  if (is.null(dend)) {
    return(0)
  } else {
    return(attr(dend, "height"))
  }
}

gc_rows <- row_meta$Codon_Raw[row_meta$Ending_Group == "G/C-ending"]
mettrp_rows <- row_meta$Codon_Raw[row_meta$Ending_Group == "Met/Trp"]
at_rows <- row_meta$Codon_Raw[row_meta$Ending_Group == "A/T-ending"]

gc_dend <- make_group_dend(plot_mat, gc_rows)
mettrp_dend <- make_group_dend(plot_mat, mettrp_rows)
at_dend <- make_group_dend(plot_mat, at_rows)

base_height <- max(
  get_dend_height(gc_dend),
  get_dend_height(mettrp_dend),
  get_dend_height(at_dend)
)

height_gap <- ifelse(base_height == 0, 1, base_height * 0.08)

if (remove_Met_Trp) {
  row_dend_use <- merge(
    gc_dend,
    at_dend,
    height = max(get_dend_height(gc_dend), get_dend_height(at_dend)) + height_gap
  )
} else {
  # 先让 G/C-ending 和 Met/Trp 聚成一支
  gc_mettrp_dend <- merge(
    gc_dend,
    mettrp_dend,
    height = max(get_dend_height(gc_dend), get_dend_height(mettrp_dend)) + height_gap
  )
  
  # 然后再和 A/T-ending 聚成完整的一棵树
  row_dend_use <- merge(
    gc_mettrp_dend,
    at_dend,
    height = max(get_dend_height(gc_mettrp_dend), get_dend_height(at_dend)) + height_gap
  )
}

# 按自定义 tree 的 leaf 顺序重排矩阵
row_order_use <- labels(row_dend_use)

plot_mat <- plot_mat[row_order_use, , drop = FALSE]
row_meta <- row_meta[match(row_order_use, row_meta$Codon_Raw), ]

# ==========================================
# 5. 色阶固定：0 → 1（白色）→ 4
# ==========================================
col_fun <- colorRamp2(c(0, 1, 4), c("#2166AC", "#f7f7f7", "#B2182B"))

# ==========================================
# 6. 顶部注释
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
if (remove_Met_Trp) {
  output_pdf <- "C:/Users/Ke-jia Shan/Desktop/RSCU_Heatmap_No_MetTrp.pdf"
  codon_msg <- "59个密码子（已去掉Met和Trp）"
} else {
  output_pdf <- "C:/Users/Ke-jia Shan/Desktop/RSCU_Heatmap_OneTree_GC_MetTrp_AT.pdf"
  codon_msg <- "全部61个密码子（一个整体；有完整row tree；顺序为G/C-ending -> Met/Trp -> A/T-ending）"
}

pdf(output_pdf, width = 22, height = 10)

ht <- Heatmap(
  plot_mat,
  col = col_fun,
  name = "RSCU Value",
  top_annotation = top_anno,
  
  # --- 纵轴 ---
  # 不使用 row_split，因此热图是一个整体
  # 但使用自定义 row dendrogram，因此左侧仍然有完整聚类树
  cluster_rows = row_dend_use,
  show_row_dend = TRUE,
  row_dend_width = unit(25, "mm"),
  row_names_gp = gpar(fontsize = 8.5),
  
  # --- 横轴 ---
  column_split = col_meta$Class_Factor,
  column_gap = unit(c(rep(1, 8), 1, 1), "mm"),
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
message("密码子设置：", codon_msg)






library(ComplexHeatmap)
library(tidyverse)
library(circlize)

ht_opt$message = FALSE

# ==========================================
# 0. 参数开关
# ==========================================
remove_Met_Trp <- FALSE   # 保留 Met/Trp

# 对于行 dendrogram，"bottom" 表示把 Met/Trp 推到该 slice 的最后，也就是热图最下方
met_trp_position <- "bottom"

# ==========================================
# 1. 读取数据
# ==========================================
file_path <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/RSCU_FigS1.csv"

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










library(ComplexHeatmap)
library(tidyverse)
library(circlize)

ht_opt$message = FALSE

# ==========================================
# 0. 参数开关
# ==========================================
remove_Met_Trp <- F   # FALSE=保留全部，TRUE=去掉Met/Trp

# Met/Trp 放到 G/C-ending 的哪一侧
# "right" 通常对应 G/C-ending 的下边缘，也就是靠近 A/T-ending 的那一侧
# 如果跑出来方向反了，改成 "left"
met_trp_side <- "right"

# ==========================================
# 1. 读取数据
# ==========================================
file_path <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/RSCU_FigS1.csv"

if (!file.exists(file_path)) {
  stop("找不到文件：", file_path)
}

raw_data <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
head(raw_data)

# raw_data[raw_data$Class %in% c("human","bat"),]$Class="Other mammals"

# ==========================================
# 2. 纵轴属性
# ==========================================
codon_cols <- colnames(raw_data)[3:ncol(raw_data)]

row_meta <- data.frame(Codon_Raw = codon_cols, stringsAsFactors = FALSE) %>%
  mutate(
    Codon = sub(".*:", "", Codon_Raw),
    Third_Base = toupper(substr(Codon, 3, 3)),
    Third_Base = ifelse(Third_Base == "T", "U", Third_Base),
    Amino_Acid = sub(":.*", "", Codon_Raw),
    Ending_Group = factor(ifelse(Third_Base %in% c("A", "U"), "A/T-ending", "G/C-ending"))
  )

if (remove_Met_Trp) {
  row_meta <- row_meta %>% filter(!Amino_Acid %in% c("Met", "Trp"))
}

row_meta <- row_meta %>% arrange(Ending_Group, Codon_Raw)
final_codon_cols <- row_meta$Codon_Raw

# Met 和 Trp 的行名
met_trp_rows <- row_meta$Codon_Raw[row_meta$Amino_Acid %in% c("Met", "Trp")]

# ==========================================
# 3. 动物→植物→病毒
# ==========================================
desired_order <- c(
  "human", "bat", "Other mammals","bird", "reptile", "fish", 
  "Liliopsida_Plant","eudicotyledons_Plant",  "OtherPlants",
  "DNA", "RNA"
)

col_meta <- raw_data %>%
  select(Class, Species) %>%
  mutate(Class_Factor = factor(Class, levels = desired_order)) %>%
  arrange(Class_Factor)

# ==========================================
# 4. 矩阵
# ==========================================
plot_mat <- t(as.matrix(raw_data[, final_codon_cols]))

# 为了尽量保持你原图一致，这里保留你的原始写法
colnames(plot_mat) <- col_meta$Species

plot_mat <- plot_mat[row_meta$Codon_Raw, ]
plot_mat <- plot_mat[, col_meta$Species]

# ==========================================
# 4.1 保持原聚类方式，只旋转 Met/Trp 所在枝
# ==========================================

rotate_targets_to_edge <- function(dend, targets, side = c("right", "left")) {
  side <- match.arg(side)
  
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
    
    if (side == "right") {
      if (left_has_target && !right_has_target) {
        d <- swap_children(d)
      }
    }
    
    if (side == "left") {
      if (!left_has_target && right_has_target) {
        d <- swap_children(d)
      }
    }
    
    return(d)
  }
  
  orient_node(dend)
}

cluster_rows_original_plus_mettrp_rotate <- function(m) {
  # 这一步等价于 ComplexHeatmap 默认 row clustering：
  # 欧氏距离 + complete linkage
  hc <- hclust(dist(m), method = "complete")
  dend <- as.dendrogram(hc)
  
  # 尽量模拟 ComplexHeatmap 默认的 dendrogram reorder
  # 这样图形会更接近你原始代码的结果
  dend <- reorder(dend, wts = rowMeans(m, na.rm = TRUE))
  
  # 只有 G/C-ending 这个 slice 里面会包含 Met/Trp
  current_targets <- intersect(labels(dend), met_trp_rows)
  
  if (length(current_targets) > 0) {
    dend <- rotate_targets_to_edge(
      dend = dend,
      targets = current_targets,
      side = met_trp_side
    )
  }
  
  return(dend)
}

# ==========================================
# 5. 色阶固定：0 → 1（白色）→ 4
# ==========================================
col_fun <- colorRamp2(c(0, 1, 4), c("#2166AC", "#f7f7f7", "#B2182B"))

# ==========================================
# 6. 顶部注释
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
if (remove_Met_Trp) {
  output_pdf <- "C:/Users/Ke-jia Shan/Desktop/RSCU_Heatmap_No_MetTrp.pdf"
  codon_msg <- "59个密码子（已去掉Met和Trp）"
} else {
  output_pdf <- "C:/Users/Ke-jia Shan/Desktop/RSCU_Heatmap_All_Codons_MetTrp_Outer.pdf"
  codon_msg <- "全部61个密码子；Met/Trp保留在G/C-ending中，仅旋转到该枝最外侧"
}

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
  
  # 原来是 cluster_rows = TRUE
  # 这里仍然使用默认聚类逻辑：
  # hclust(dist(...), method = "complete")
  # 只是额外对 Met/Trp 所在枝进行旋转
  cluster_rows = cluster_rows_original_plus_mettrp_rotate,
  
  # 防止 ComplexHeatmap 在我们旋转后再次自动重排 dendrogram
  row_dend_reorder = FALSE,
  
  row_names_gp = gpar(fontsize = 8.5),
  
  # --- 横轴 ---
  # 第9个间隙（OtherPlants和DNA之间）加大到4mm，其余0.3mm
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
message("密码子设置：", codon_msg)