rm(list = ls())
library(dplyr)
library(ArrayExpress)
library(RCurl)
# 1. 读取数据文件
patient_data <- read.csv("./data/information.csv")
head(patient_data)

# 2. 数据预处理
# 解析 File Names 列，提取每个样本的第一个文件名，并创建新列 `Selected_File`
patient_data <- patient_data %>%
  rowwise() %>%
  mutate(
    # 提取第一个样本名（从字符串中以 OR 分隔）
    Selected_File = strsplit(`X.search.`, " OR ")[[1]][1],
    # 提取EventFreeSurvival的数值部分（天数）
    EventFreeSurvival_days = as.numeric(sub(" day", "", `EventFreeSurvival`)),
    # 提取OverallSurvival的数值部分（天数）
    OverallSurvival_days = as.numeric(sub(" day", "", `OverallSurvival`))
  ) %>%
  ungroup()


head(patient_data %>% select(Selected_File, EventFreeSurvival_days, OverallSurvival_days, ClinicalInformation))

# 创建分组因子列
patient_data$Group <- NA

# # 正响应组：3年内没有事件发生
# patient_data$Group[patient_data$EventFreeSurvival >= 1095 & patient_data$ClinicalInformation == "no event"] <- "Positive"
# # 负响应组：3年内或3年后发生事件的样本
# patient_data$Group[patient_data$EventFreeSurvival < 1095 | 
#                      (patient_data$EventFreeSurvival >= 1095 & patient_data$ClinicalInformation %in% c("relapse of disease", "death from disease"))] <- "Negative"


patient_data$Group[patient_data$ClinicalInformation == "no event"] <- "Negative"
patient_data$Group[patient_data$ClinicalInformation != "no event"] <- "Positive"

table(patient_data$Group)

final_data <- patient_data %>% select(c(Selected_File, Group))
head(final_data)
write.csv(final_data, "./data/group_info.csv", row.names = TRUE)




# 推荐使用 “gBGSubSignal”（绿色通道背景减除信号）或 “rBGSubSignal”（红色通道背景减除信号）
# 作为表达量列。如果需要控制染料效应，也可以选择 “gDyeNormSignal” 或 “rDyeNormSignal”。
# 本文中使用了背景减除信号，因此使用 “gBGSubSignal” 或 “rBGSubSignal” 更为常见。
data_folder <- "./data/data" # 设置表达量数据文件夹路径
signal_column <- "Feature.Extraction.Software.gBGSubSignal"

# 提取所需样本文件名和分组信息
sample_files <- final_data$Selected_File
sample_groups <- final_data$Group

# 初始化表达矩阵和探针信息
expression_matrix <- NULL
probe_ids <- NULL

# 遍历每个样本文件并提取表达量
for (file in sample_files) {
  # 构建完整的文件路径
  file_path <- file.path(data_folder, file)
  
  # 读取样本表达量文件（假设为制表符分隔的文本文件）
  sample_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 提取探针ID和表达量列（假设探针ID列为 "Reporter identifier")
  if (is.null(probe_ids)) {
    probe_ids <- sample_data[["Reporter identifier"]]
  }
  
  # 提取指定信号列作为表达量
  expression_values <- sample_data[[signal_column]]
  
  # 将表达量合并到表达矩阵中
  expression_matrix <- cbind(expression_matrix, expression_values)
}

# 为表达矩阵添加行名（探针ID）和列名（样本文件名）
rownames(expression_matrix) <- probe_ids
colnames(expression_matrix) <- sample_files

# 保存表达矩阵为 CSV 文件
write.csv(expression_matrix, "./data/expression_matrix.csv", row.names = TRUE)

# 显示表达矩阵的前几行
head(expression_matrix)





rm(list = ls())
library(dplyr)
expression_matrix <- read.csv("./data/expression_matrix.csv", row.names = 1, stringsAsFactors = FALSE)
expression_samples <- colnames(expression_matrix)

final_data <- read.csv("./data/group_info.csv", row.names = 1, stringsAsFactors = FALSE)
print(setdiff(final_data$Selected_File, expression_samples))  # 如果为空，表示完全匹配



min_value <- min(expression_matrix, na.rm = TRUE)
if (min_value <= 0) {
  # 添加常数（比如 min_value 的绝对值加 1），确保所有值都大于零
  shift_constant <- abs(min_value) + 1
  expression_matrix_shifted <- expression_matrix + shift_constant
} else {
  expression_matrix_shifted <- expression_matrix
}

expression_matrix <- log2(expression_matrix_shifted)
write.csv(expression_matrix, "./data/log2_expression_matrix.csv", row.names = TRUE)



positive_samples <- final_data %>%
  filter(Group == "Positive") %>%
  pull(Selected_File)
positive_expression <- expression_matrix[, colnames(expression_matrix) %in% positive_samples]
write.csv(positive_expression, "./data/positive_expression_matrix.csv", row.names = TRUE)


negative_samples <- final_data %>%
  filter(Group == "Negative") %>%
  pull(Selected_File)
negative_expression <- expression_matrix[, colnames(expression_matrix) %in% negative_samples]
write.csv(negative_expression, "./data/negative_expression_matrix.csv", row.names = TRUE)